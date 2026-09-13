# ============================================================
# cell.py — COMMON DEFINITIONS
#
# Definitions only.
#
# No expensive calculations are performed on import.
# ============================================================

import mpmath as mp

from connes_cvs import (
    build_galerkin_matrix,
    compute_ground_state,
)
from connes_cvs.operator import (
    _compute_psi_pair,
    prime_powers_up_to,
    HAS_FLINT,
)

try:
    from flint import ctx as flint_ctx
except ImportError:
    flint_ctx = None

import hashlib
import json
import time
from pathlib import Path

# ============================================================
# DEFAULT NUMERICAL PARAMETERS
# ============================================================

DEFAULT_DPS = 80

# ============================================================
# GENERIC JSON RESULT CACHE
# ============================================================
#
# Persistent, content-addressed cache for expensive calculations.
#
# The cache is deliberately generic. It knows only:
#
#   namespace
#   parameters
#   results
#
# Parameters and results must be JSON-compatible.
#
# The cache key is SHA-256 over a canonical JSON representation of
# the namespace and parameters.
#
# Cache semantics:
#
#   lookup -> hit:
#       decode -> validate -> return
#
#   lookup -> miss:
#       generate -> write -> lookup again -> decode -> validate -> return
#
# The second path deliberately rejoins the first path after writing.
# Thus cache hit and cache miss have identical observable numerical
# semantics. The cache changes speed, not results.
# ============================================================

CACHE_SCHEMA_VERSION = 3

CELL_CACHE_DIR = (
    Path(__file__).resolve().parent / ".cell_cache"
)


def _cache_canonical_json(obj):
    """
    Return the canonical JSON representation used for hashing.

    Sorting keys and using fixed separators makes the representation
    independent of dictionary insertion order.
    """
    return json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    )


def cache_key(namespace, parameters):
    """
    Return:

        (sha256_hex_digest, canonical_parameter_record)

    for a cache entry.

    `parameters` must be JSON-compatible.
    """
    record = {
        "schema_version": CACHE_SCHEMA_VERSION,
        "namespace": str(namespace),
        "parameters": parameters,
    }

    canonical = _cache_canonical_json(record)

    digest = hashlib.sha256(
        canonical.encode("utf-8")
    ).hexdigest()

    return digest, record


def _cache_namespace_dir(namespace):
    """
    Return the cache directory for a namespace.
    """
    path = CELL_CACHE_DIR / str(namespace)
    path.mkdir(
        parents=True,
        exist_ok=True,
    )
    return path


def _cache_path(namespace, digest):
    """
    Return the path for a cache entry.
    """
    return (
        _cache_namespace_dir(namespace)
        / f"{digest}.json"
    )


def cache_save(
    namespace,
    parameters,
    results,
    *,
    timing=None,
):
    """
    Save a JSON-compatible result under its content-derived key.

    Existing entries are never silently replaced by a different
    calculation. If the path already exists, its identity must match
    the requested namespace and parameters.

    Returns:

        (digest, path)
    """
    digest, identity = cache_key(
        namespace,
        parameters,
    )

    path = _cache_path(
        namespace,
        digest,
    )

    # If an entry already exists, verify that it is the exact same
    # cache identity. Never mutate an existing entry.
    if path.exists():
        with path.open(
            "r",
            encoding="utf-8",
        ) as f:
            existing = json.load(f)

        if (
            existing.get("cache_schema_version")
            != CACHE_SCHEMA_VERSION
        ):
            raise ValueError(
                "existing cache entry has incompatible schema"
            )

        if existing.get("cache_key") != digest:
            raise ValueError(
                "existing cache entry has inconsistent key"
            )

        if existing.get("namespace") != str(namespace):
            raise ValueError(
                "existing cache entry has inconsistent namespace"
            )

        if existing.get("parameters") != parameters:
            raise ValueError(
                "existing cache entry has inconsistent parameters"
            )

        # Exact same immutable entry already exists.
        return digest, path

    payload = {
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "cache_key": digest,
        "namespace": str(namespace),
        "parameters": parameters,
        "results": results,
    }

    if timing is not None:
        payload["timing"] = timing

    # Atomic write. The temporary file is in the same directory so
    # filesystem rename semantics remain atomic. Use a unique filename
    # per PID and timestamp to avoid race conditions during concurrent runs.
    import os
    temporary = path.with_suffix(f".tmp.{os.getpid()}_{time.time_ns()}")

    with temporary.open(
        "w",
        encoding="utf-8",
    ) as f:
        json.dump(
            payload,
            f,
            indent=2,
            sort_keys=True,
            ensure_ascii=False,
        )
        f.write("\n")

    temporary.replace(path)

    return digest, path


def cache_load(
    namespace,
    parameters,
):
    """
    Load a cache entry.

    Returns:

        (results, metadata)

    Raises FileNotFoundError if the entry does not exist.
    Raises ValueError if the entry exists but its identity is
    inconsistent with the requested namespace/parameters.
    """
    digest, identity = cache_key(
        namespace,
        parameters,
    )

    path = _cache_path(
        namespace,
        digest,
    )

    if not path.exists():
        raise FileNotFoundError(path)

    try:
        with path.open(
            "r",
            encoding="utf-8",
        ) as f:
            payload = json.load(f)
    except (SystemError, Exception):
        raise FileNotFoundError(path)

    # Validate the cache record itself rather than trusting the
    # filename.
    if (
        payload.get("cache_schema_version")
        != CACHE_SCHEMA_VERSION
    ):
        raise ValueError(
            "cache schema version mismatch"
        )

    if payload.get("cache_key") != digest:
        raise ValueError(
            "cache key mismatch"
        )

    if payload.get("namespace") != str(namespace):
        raise ValueError(
            "cache namespace mismatch"
        )

    if payload.get("parameters") != parameters:
        raise ValueError(
            "cache parameters do not match request"
        )

    if "results" not in payload:
        raise ValueError(
            "cache entry contains no results"
        )

    return (
        payload["results"],
        {
            "cache_hit": True,
            "cache_key": digest,
            "cache_path": path,
            "cache_record": identity,
            "stored_timing": payload.get("timing"),
        },
    )


# ============================================================
# GALERKIN MATRIX CACHE WRAPPER
# ============================================================

GALERKIN_MATRIX_OPERATOR_VERSION = (
    "cell.py-galerkin-matrix-v1"
)


def _galerkin_matrix_parameters(
    c,
    N,
    T,
    dps,
    flint_bits=None,
):
    """
    Construct the complete identity of a Galerkin matrix calculation.

    `dps` is the generation/certification precision and forms part
    of the immutable cache identity.
    """
    if isinstance(c, float):
        raise TypeError(
            "Galerkin matrix c must not be a Python float; "
            "use an integer, decimal string, or mp.mpf."
        )

    c_mp = mp.mpf(c)

    if flint_bits is None:
        flint_bits = int(
            int(dps) * 3.5
        )

    return {
        "operator_version": (
            GALERKIN_MATRIX_OPERATOR_VERSION
        ),
        "c": mp.nstr(
            c_mp,
            max(50, int(dps) + 10),
        ),
        "N": int(N),
        "T": int(T),
        "dps": int(dps),
        "flint_bits": int(flint_bits),
    }


def _assemble_Q_from_psi(psi_vals, psi_deriv_vals, N):
    """
    Assemble the (2N+1) x (2N+1) Galerkin matrix Q from cached
    basis functional values psi(n) and psi'(n) for n in [0, N].

    Uses the exact parity identities:
        psi(-n)  = -psi(n)
        psi'(-n) =  psi'(n)
    and the difference quotient:
        Q[m, n] = psi'(n)                   if m == n
        Q[m, n] = (psi(m) - psi(n))/(m - n) if m != n
    for indices m, n in [-N, N].
    """
    dim = 2 * N + 1
    Q = mp.matrix(dim, dim)

    full_psi = [mp.mpf(0)] * dim
    full_psi_d = [mp.mpf(0)] * dim

    for n in range(N + 1):
        p = psi_vals[n]
        pd = psi_deriv_vals[n]
        full_psi[N + n] = p
        full_psi[N - n] = -p
        full_psi_d[N + n] = pd
        full_psi_d[N - n] = pd

    for i in range(dim):
        m = i - N
        p_m = full_psi[i]
        for j in range(i, dim):
            n = j - N
            if m == n:
                val = full_psi_d[j]
            else:
                val = (p_m - full_psi[j]) / (m - n)
            Q[i, j] = val
            Q[j, i] = val

    return Q


# ============================================================
# MATRIX ALGEBRA & INVARIANTS
# ============================================================

def frobenius_norm(A: mp.matrix) -> mp.mpf:
    """
    Compute the Frobenius norm ||A||_F = sqrt(sum_{i,j} |A_{ij}|^2) of an mpmath matrix.

    Uses explicit accumulation to ensure universal mpmath portability without
    relying on ambiguous string order parameters in mp.norm (which raise TypeError
    in mpmath when passed 'fro' or 'f').
    """
    s = mp.mpf(0)
    for i in range(A.rows):
        for j in range(A.cols):
            val = A[i, j]
            s += val * val
    return mp.sqrt(s)


def _galerkin_matrix_encode(
    psi_vals,
    psi_deriv_vals,
    N,
    dps,
    trace_val,
    frob_val,
):
    """
    Convert arbitrary-precision Galerkin basis values into a
    JSON-compatible fragment, retaining guard digits.
    """
    digits = int(dps) + 10

    return {
        "N": int(N),
        "psi_vals": [
            mp.nstr(psi_vals[n], digits)
            for n in range(N + 1)
        ],
        "psi_deriv_vals": [
            mp.nstr(psi_deriv_vals[n], digits)
            for n in range(N + 1)
        ],
        "trace": mp.nstr(trace_val, digits),
        "frobenius_norm": mp.nstr(frob_val, digits),
    }


def _galerkin_matrix_decode(
    results,
    N,
):
    """
    Reconstruct mpmath Galerkin matrix Q and basis values from cached JSON.

    Decoding occurs at caller's current mp.mp.dps.
    """
    cached_N = int(results["N"])
    if cached_N < N:
        raise ValueError(
            f"cached Galerkin data has N={cached_N}, expected at least {N}"
        )

    raw_psi = results["psi_vals"]
    raw_psi_d = results["psi_deriv_vals"]

    psi_vals = [mp.mpf(raw_psi[n]) for n in range(N + 1)]
    psi_deriv_vals = [mp.mpf(raw_psi_d[n]) for n in range(N + 1)]

    Q = _assemble_Q_from_psi(psi_vals, psi_deriv_vals, N)

    return Q, psi_vals, psi_deriv_vals


def _validate_galerkin_matrix_structure(
    Q,
    N,
    stored_trace=None,
):
    """
    Cheap intrinsic validation of the reconstructed Galerkin matrix Q.
    Checks dimensions, symmetry, finiteness, realness, and trace checksum.
    """
    expected_dim = 2 * N + 1

    if Q.rows != expected_dim or Q.cols != expected_dim:
        raise ValueError(
            f"Galerkin matrix dimensions ({Q.rows}, {Q.cols}) "
            f"do not match expected ({expected_dim}, {expected_dim})"
        )

    trace = mp.mpf(0)
    for i in range(expected_dim):
        diag_val = Q[i, i]
        if not mp.isfinite(diag_val):
            raise ValueError(f"diagonal entry Q[{i}, {i}] is not finite")
        trace += diag_val

        for j in range(i + 1, expected_dim):
            val_ij = Q[i, j]
            val_ji = Q[j, i]
            if not mp.isfinite(val_ij):
                raise ValueError(f"matrix entry Q[{i}, {j}] is not finite")
            if val_ij != val_ji:
                raise ValueError(f"matrix asymmetry detected at ({i}, {j})")

    trace_error = None
    if stored_trace is not None:
        expected_trace = mp.mpf(stored_trace)
        trace_error = abs(trace - expected_trace)
        tolerance = mp.sqrt(mp.eps) * max(1, abs(expected_trace))
        if trace_error > tolerance:
            raise ValueError(
                f"Galerkin matrix trace validation failed: "
                f"{mp.nstr(trace_error, 6)} > {mp.nstr(tolerance, 6)}"
            )

    return {
        "dimension": expected_dim,
        "trace": trace,
        "trace_error": trace_error,
    }


def _find_compatible_psi_cache(
    c_str,
    T,
    dps,
    flint_bits,
    target_N,
):
    """
    Search .cell_cache/galerkin_matrix for entries with matching
    (operator_version, c, T, dps, flint_bits).

    Returns:
        (best_entry_results, best_N) or (None, -1)
    """
    namespace = "galerkin_matrix"
    cache_dir = _cache_namespace_dir(namespace)

    best_results = None
    best_N = -1

    try:
        cache_paths = list(cache_dir.glob("*.json"))
    except (SystemError, Exception):
        cache_paths = []

    for path in cache_paths:
        try:
            with path.open("r", encoding="utf-8") as f:
                payload = json.load(f)
            params = payload.get("parameters", {})
            if (
                params.get("operator_version") == GALERKIN_MATRIX_OPERATOR_VERSION
                and params.get("c") == c_str
                and params.get("T") == int(T)
                and params.get("dps") == int(dps)
                and params.get("flint_bits") == int(flint_bits)
            ):
                entry_N = int(params.get("N", -1))
                if entry_N >= target_N:
                    return payload.get("results"), entry_N
                if entry_N > best_N:
                    best_N = entry_N
                    best_results = payload.get("results")
        except (SystemError, Exception):
            continue

    return best_results, best_N


def _generate_galerkin_matrix(
    c,
    N,
    T,
    dps,
    flint_bits,
):
    """
    Generate Galerkin matrix data at generation precision, reusing any
    available points from existing cached entries for the same parameters.
    """
    c_mp = mp.mpf(c)
    L = mp.log(c_mp)
    prime_data, _ = prime_powers_up_to(int(mp.floor(c_mp)))

    import connes_cvs.operator as _op
    if _op.HAS_FLINT and flint_ctx is not None:
        flint_ctx.prec = flint_bits
        flint_ctx.threads = 1

    c_str = mp.nstr(c_mp, max(50, int(dps) + 10))
    existing_results, avail_N = _find_compatible_psi_cache(
        c_str=c_str,
        T=T,
        dps=dps,
        flint_bits=flint_bits,
        target_N=N,
    )

    psi_vals = {}
    psi_deriv_vals = {}
    reused_count = 0

    if existing_results is not None and avail_N >= 0:
        raw_p = existing_results["psi_vals"]
        raw_pd = existing_results["psi_deriv_vals"]
        max_take = min(avail_N, N)
        for n in range(max_take + 1):
            psi_vals[n] = mp.mpf(raw_p[n])
            psi_deriv_vals[n] = mp.mpf(raw_pd[n])
        reused_count = max_take + 1

    start_idx = reused_count
    quad_start = time.perf_counter()

    for n_idx in range(start_idx, N + 1):
        psi, psi_d = _compute_psi_pair(n_idx, L, T, dps, prime_data)
        psi_vals[n_idx] = psi
        psi_deriv_vals[n_idx] = psi_d

    quad_elapsed = time.perf_counter() - quad_start

    assembly_start = time.perf_counter()
    Q = _assemble_Q_from_psi(psi_vals, psi_deriv_vals, N)
    assembly_elapsed = time.perf_counter() - assembly_start

    trace_val = mp.mpf(0)
    dim = 2 * N + 1
    for i in range(dim):
        trace_val += Q[i, i]
    frob_val = frobenius_norm(Q)

    encoded = _galerkin_matrix_encode(
        psi_vals=psi_vals,
        psi_deriv_vals=psi_deriv_vals,
        N=N,
        dps=dps,
        trace_val=trace_val,
        frob_val=frob_val,
    )

    return (
        encoded,
        {
            "quadrature_seconds": quad_elapsed,
            "assembly_seconds": assembly_elapsed,
            "reused_points": reused_count,
            "computed_points": (N + 1) - reused_count,
            "trace": trace_val,
            "frobenius_norm": frob_val,
        },
    )


def get_galerkin_matrix(
    c,
    N,
    T,
    dps,
    *,
    flint_bits=None,
    verbose=True,
):
    """
    Obtain the truncated Weil Galerkin matrix Q(c, N) through the persistent cache.

    Cache semantics
    ---------------
    The cache is always enabled.
    `dps` is the generation/certification precision and forms part
    of the immutable cache identity.
    The caller's current `mp.mp.dps` is the working precision.

    A cache hit:
        lookup -> decode -> validate -> return Q, metadata

    A cache miss:
        generate -> write -> fresh lookup -> decode -> validate -> return Q, metadata

    Sub-dimension reuse:
        If an existing cache entry for the same (c, T, dps) has dimension N' >= N,
        the matrix Q is assembled instantly from the existing basis points without
        re-evaluating any quadratures.
        If N' < N, existing points are reused and only the missing points
        n in [N' + 1, N] are evaluated.

    Returns
    -------
    Q, metadata
        Q : mpmath.matrix of dimension (2N+1) x (2N+1), symmetric and real.
        metadata : dict of cache and timing metrics.
    """
    working_dps = int(mp.mp.dps)
    if working_dps <= 0:
        raise ValueError("current mp.mp.dps must be positive")

    requested_generation_dps = int(dps)
    if requested_generation_dps <= 0:
        raise ValueError("dps must be positive")

    if flint_bits is None:
        flint_bits = int(requested_generation_dps * 3.5)

    namespace = "galerkin_matrix"

    generation_dps = max(requested_generation_dps, working_dps)
    if generation_dps != requested_generation_dps and verbose:
        print()
        print("GALERKIN MATRIX CACHE: GENERATION PRECISION PROMOTED")
        print(f"  requested dps = {requested_generation_dps}")
        print(f"  working dps   = {working_dps}")
        print(f"  generation dps = {generation_dps}")

    parameters = _galerkin_matrix_parameters(
        c=c,
        N=N,
        T=T,
        dps=generation_dps,
        flint_bits=flint_bits,
    )

    lookup_start = time.perf_counter()

    try:
        results, cache_meta = cache_load(namespace, parameters)
        lookup_elapsed = time.perf_counter() - lookup_start
        cache_hit = True
        generation_metadata = None

        if verbose:
            print()
            print("GALERKIN MATRIX CACHE: HIT")
            print(f"  key            = {cache_meta['cache_key']}")
            print(f"  N              = {N} (dim = {2*N+1})")
            print(f"  generation dps = {generation_dps}")
            print(f"  lookup         = {lookup_elapsed:.6f} s")

    except FileNotFoundError:
        lookup_elapsed = time.perf_counter() - lookup_start
        cache_hit = False

        if verbose:
            print()
            print("GALERKIN MATRIX CACHE: MISS")
            print(f"  N              = {N} (dim = {2*N+1})")
            print(f"  generation dps = {generation_dps}")
            print(f"  working dps    = {working_dps}")

        caller_dps = mp.mp.dps
        generation_start = time.perf_counter()
        mp.mp.dps = generation_dps

        try:
            generated_results, generation_metadata = _generate_galerkin_matrix(
                c=c,
                N=N,
                T=T,
                dps=generation_dps,
                flint_bits=flint_bits,
            )

            save_start = time.perf_counter()
            cache_save(
                namespace,
                parameters,
                generated_results,
                timing={
                    **{
                        k: v for k, v in generation_metadata.items()
                        if k.endswith("_seconds")
                    },
                    "generation_dps": generation_dps,
                    "reused_points": generation_metadata["reused_points"],
                    "computed_points": generation_metadata["computed_points"],
                },
            )
            save_elapsed = time.perf_counter() - save_start
            generation_elapsed = time.perf_counter() - generation_start

        finally:
            mp.mp.dps = caller_dps

        if verbose:
            print()
            print("GALERKIN MATRIX CACHE: GENERATED")
            print(f"  reused points  = {generation_metadata['reused_points']}")
            print(f"  computed pts   = {generation_metadata['computed_points']}")
            print(f"  quadrature     = {generation_metadata['quadrature_seconds']:.6f} s")
            print(f"  assembly       = {generation_metadata['assembly_seconds']:.6f} s")
            print(f"  save           = {save_elapsed:.6f} s")
            print(f"  generation     = {generation_elapsed:.6f} s")

        lookup_start = time.perf_counter()
        results, cache_meta = cache_load(namespace, parameters)
        second_lookup_elapsed = time.perf_counter() - lookup_start

        cache_meta = dict(cache_meta)
        cache_meta["initial_cache_miss_seconds"] = lookup_elapsed
        cache_meta["generation_seconds"] = generation_elapsed
        cache_meta["save_seconds"] = save_elapsed
        cache_meta["final_lookup_seconds"] = second_lookup_elapsed
        cache_meta["cache_hit"] = False
        cache_meta["generated"] = True
        cache_meta["reused_points"] = generation_metadata["reused_points"]
        cache_meta["computed_points"] = generation_metadata["computed_points"]

    decode_start = time.perf_counter()
    Q, psi_vals, psi_deriv_vals = _galerkin_matrix_decode(results, N)
    decode_elapsed = time.perf_counter() - decode_start

    validation_start = time.perf_counter()
    structural = _validate_galerkin_matrix_structure(
        Q,
        N,
        stored_trace=results.get("trace"),
    )
    validation_elapsed = time.perf_counter() - validation_start

    cache_meta = dict(cache_meta)
    cache_meta["generation_dps"] = generation_dps
    cache_meta["working_dps"] = working_dps
    cache_meta["decode_seconds"] = decode_elapsed
    cache_meta["validation_seconds"] = validation_elapsed
    cache_meta["structural_validation"] = structural

    if cache_hit:
        cache_meta["total_seconds"] = (
            lookup_elapsed + decode_elapsed + validation_elapsed
        )
    else:
        cache_meta["total_seconds"] = (
            cache_meta["initial_cache_miss_seconds"]
            + cache_meta["generation_seconds"]
            + cache_meta["save_seconds"]
            + cache_meta["final_lookup_seconds"]
            + decode_elapsed
            + validation_elapsed
        )

    if verbose:
        print()
        print("GALERKIN MATRIX CACHE: RETURN")
        print(f"  decode         = {decode_elapsed:.6f} s")
        print(f"  validation     = {validation_elapsed:.6f} s")
        print(f"  total          = {cache_meta['total_seconds']:.6f} s")

    return Q, cache_meta


get_cached_galerkin_matrix = get_galerkin_matrix


# ============================================================
# GROUND-STATE CACHE WRAPPER
# ============================================================

GROUND_STATE_OPERATOR_VERSION = (
    "cell.py-ground-state-v2"
)


def _ground_state_parameters(
    c,
    N,
    T,
    dps,
    flint_bits=None,
):
    """
    Construct the complete identity of a ground-state calculation.

    `dps` is the generation/certification precision and therefore
    forms part of the immutable cache identity.
    """
    if isinstance(c, float):
        raise TypeError(
            "Ground-state c must not be a Python float; "
            "use an integer, decimal string, or mp.mpf."
        )

    c_mp = mp.mpf(c)

    if flint_bits is None:
        flint_bits = int(
            int(dps) * 3.5
        )

    return {
        "operator_version": (
            GROUND_STATE_OPERATOR_VERSION
        ),
        "c": mp.nstr(
            c_mp,
            max(50, int(dps) + 10),
        ),
        "N": int(N),
        "T": int(T),
        "dps": int(dps),
        "flint_bits": int(flint_bits),
    }


def _ground_state_encode(
    lambda_min,
    v_canonical,
    dps,
):
    """
    Convert arbitrary-precision numerical results into a
    JSON-compatible fragment.

    The serialised result retains guard digits beyond generation dps.
    """
    digits = int(dps) + 10

    return {
        "lambda_min": mp.nstr(
            lambda_min,
            digits,
        ),
        "v_canonical": [
            mp.nstr(
                v_canonical[i, 0],
                digits,
            )
            for i in range(v_canonical.rows)
        ],
    }


def _ground_state_decode(
    results,
    N,
):
    """
    Reconstruct mpmath numerical objects from cached JSON.

    Decoding occurs at the caller's current mp.mp.dps. This is
    deliberate: the cached artefact may have been generated at a
    substantially higher precision than the current calculation
    requires.
    """
    lambda_min = mp.mpf(
        results["lambda_min"]
    )

    values = results["v_canonical"]

    expected = N + 1

    if len(values) != expected:
        raise ValueError(
            "cached ground-state vector has "
            f"length {len(values)}, expected {expected}"
        )

    v_canonical = mp.matrix(
        expected,
        1,
    )

    for i, value in enumerate(values):
        v_canonical[i, 0] = mp.mpf(value)

    return lambda_min, v_canonical


def _validate_ground_state_structure(
    lambda_min,
    v_canonical,
    N,
):
    """
    Cheap intrinsic validation performed on every retrieval.

    This validates the integrity and basic numerical structure of the
    cached result. It does NOT rebuild the Galerkin matrix.

    Mathematical correctness of the ground state is deliberately not
    part of ordinary cache retrieval. If that is ever in doubt, an
    explicit audit should be performed.
    """
    expected = N + 1

    if v_canonical.rows != expected:
        raise ValueError(
            "ground-state vector has incorrect row count"
        )

    if v_canonical.cols != 1:
        raise ValueError(
            "ground-state vector must be a column vector"
        )

    if not mp.isfinite(lambda_min):
        raise ValueError(
            "ground-state eigenvalue is not finite"
        )

    for i in range(expected):
        value = v_canonical[i, 0]

        if not mp.isfinite(value):
            raise ValueError(
                f"ground-state vector entry {i} "
                "is not finite"
            )

        if mp.im(value) != 0:
            raise ValueError(
                f"ground-state vector entry {i} "
                "is not real"
            )

    norm = mp.sqrt(
        mp.fdot(v_canonical, v_canonical)
    )

    norm_error = abs(norm - 1)

    tolerance = mp.sqrt(
        mp.eps
    )

    if norm_error > tolerance:
        raise ValueError(
            "ground-state norm validation failed: "
            f"{mp.nstr(norm_error, 10)} > "
            f"{mp.nstr(tolerance, 10)}"
        )

    return {
        "norm": norm,
        "norm_error": norm_error,
        "tolerance": tolerance,
    }


def _generate_ground_state(
    c,
    N,
    T,
    dps,
    flint_bits,
):
    """
    Generate a ground state at exactly the requested generation dps.

    The caller is responsible for establishing the generation
    precision context before calling this function.
    """
    Q_start = time.perf_counter()

    Q, _ = get_galerkin_matrix(
        c=c,
        N=N,
        T=T,
        dps=dps,
        flint_bits=flint_bits,
        verbose=False,
    )

    Q_build_elapsed = (
        time.perf_counter()
        - Q_start
    )

    eig_start = time.perf_counter()

    lambda_min, u_full = (
        compute_ground_state(Q)
    )

    v_canonical = full_to_canonical(u_full)

    eig_elapsed = (
        time.perf_counter()
        - eig_start
    )

    # Generation-time structural validation.
    validation_start = time.perf_counter()

    structural = (
        _validate_ground_state_structure(
            lambda_min,
            v_canonical,
            N,
        )
    )

    validation_elapsed = (
        time.perf_counter()
        - validation_start
    )

    results = _ground_state_encode(
        lambda_min,
        v_canonical,
        dps,
    )

    return (
        results,
        {
            "Q_build_seconds": Q_build_elapsed,
            "eigensolve_seconds": eig_elapsed,
            "validation_seconds": validation_elapsed,
            "structural_validation": structural,
        },
    )


def get_ground_state(
    c,
    N,
    T,
    dps,
    *,
    flint_bits=None,
    verbose=True,
):
    """
    Obtain the CvS ground state through the persistent cache.

    Cache semantics
    ---------------
    The cache is always enabled.

    `dps` is the generation/certification precision and forms part
    of the immutable cache identity.

    The caller's current `mp.mp.dps` is the working precision.

    A cache hit:
        lookup -> decode -> validate -> return

    A cache miss:
        generate -> write -> fresh lookup -> decode -> validate -> return

    The deliberate second lookup after generation means that cache
    hit and cache miss use exactly the same decode/validation/return
    path.

    If the requested generation precision is below the caller's
    working precision, that cache entry is not sufficient. A new
    calculation is generated at the caller's working precision and
    stored under its own immutable cache key.

    Existing cache entries are never mutated.

    Mathematical validation of the cached eigenpair against a freshly
    rebuilt Galerkin matrix is NOT performed here. The cache is
    assumed mathematically correct unless there is reason to suspect
    otherwise. Such an expensive check belongs in a separate explicit
    audit operation.

    Returns
    -------
    lambda_min, v_canonical, metadata
    """

    working_dps = int(mp.mp.dps)

    if working_dps <= 0:
        raise ValueError(
            "current mp.mp.dps must be positive"
        )

    requested_generation_dps = int(dps)

    if requested_generation_dps <= 0:
        raise ValueError(
            "dps must be positive"
        )

    if flint_bits is None:
        flint_bits = int(
            requested_generation_dps * 3.5
        )

    namespace = "ground_state"

    # --------------------------------------------------------
    # Resolve the generation precision.
    #
    # If the requested generation precision is insufficient for
    # the current working precision, promote the calculation to
    # working precision. The lower-precision cache entry is left
    # untouched.
    # --------------------------------------------------------

    generation_dps = max(
        requested_generation_dps,
        working_dps,
    )

    if generation_dps != requested_generation_dps:
        if verbose:
            print()
            print(
                "GROUND STATE CACHE: "
                "GENERATION PRECISION PROMOTED"
            )
            print(
                f"  requested dps = "
                f"{requested_generation_dps}"
            )
            print(
                f"  working dps   = "
                f"{working_dps}"
            )
            print(
                f"  generation dps = "
                f"{generation_dps}"
            )
            print(
                "  existing lower-precision "
                "cache entry will not be modified"
            )

    parameters = _ground_state_parameters(
        c=c,
        N=N,
        T=T,
        dps=generation_dps,
        flint_bits=flint_bits,
    )

    # --------------------------------------------------------
    # LOOKUP
    #
    # This is the only route to the returned result.
    # --------------------------------------------------------

    lookup_start = time.perf_counter()

    try:
        results, cache_meta = cache_load(
            namespace,
            parameters,
        )

        lookup_elapsed = (
            time.perf_counter()
            - lookup_start
        )

        cache_hit = True
        generation_metadata = None

        if verbose:
            print()
            print(
                "GROUND STATE CACHE: HIT"
            )
            print(
                f"  key            = "
                f"{cache_meta['cache_key']}"
            )
            print(
                f"  generation dps = "
                f"{generation_dps}"
            )
            print(
                f"  working dps    = "
                f"{working_dps}"
            )
            print(
                f"  lookup         = "
                f"{lookup_elapsed:.6f} s"
            )

    except FileNotFoundError:

        lookup_elapsed = (
            time.perf_counter()
            - lookup_start
        )

        cache_hit = False

        if verbose:
            print()
            print(
                "GROUND STATE CACHE: MISS"
            )
            print(
                f"  generation dps = "
                f"{generation_dps}"
            )
            print(
                f"  working dps    = "
                f"{working_dps}"
            )

        # ----------------------------------------------------
        # GENERATION
        #
        # Do the entire generation and generation-time
        # validation at generation precision.
        # ----------------------------------------------------

        caller_dps = mp.mp.dps

        generation_start = time.perf_counter()

        mp.mp.dps = generation_dps

        try:
            (
                generated_results,
                generation_metadata,
            ) = _generate_ground_state(
                c=c,
                N=N,
                T=T,
                dps=generation_dps,
                flint_bits=flint_bits,
            )

            save_start = time.perf_counter()

            cache_save(
                namespace,
                parameters,
                generated_results,
                timing={
                    **{
                        key: value
                        for key, value
                        in generation_metadata.items()
                        if key.endswith("_seconds")
                    },
                    "generation_dps": generation_dps,
                },
            )

            save_elapsed = (
                time.perf_counter()
                - save_start
            )

            generation_elapsed = (
                time.perf_counter()
                - generation_start
            )

        finally:
            # The generated result must NOT be returned directly.
            # The next operation is a completely fresh cache lookup,
            # and therefore the returned result will have exactly the
            # same semantics as a cache hit.
            mp.mp.dps = caller_dps

        if verbose:
            print()
            print(
                "GROUND STATE CACHE: GENERATED"
            )
            print(
                f"  Q build        = "
                f"{generation_metadata['Q_build_seconds']:.6f} s"
            )
            print(
                f"  eigensolve     = "
                f"{generation_metadata['eigensolve_seconds']:.6f} s"
            )
            print(
                f"  validation     = "
                f"{generation_metadata['validation_seconds']:.6f} s"
            )
            print(
                f"  save           = "
                f"{save_elapsed:.6f} s"
            )
            print(
                f"  generation     = "
                f"{generation_elapsed:.6f} s"
            )

        # ----------------------------------------------------
        # CRITICAL:
        #
        # Do NOT decode generated_results here.
        #
        # Start again from cache_load(), exactly as though this
        # had been a cache hit from the beginning.
        # ----------------------------------------------------

        lookup_start = time.perf_counter()

        results, cache_meta = cache_load(
            namespace,
            parameters,
        )

        second_lookup_elapsed = (
            time.perf_counter()
            - lookup_start
        )

        cache_meta = dict(cache_meta)
        cache_meta["initial_cache_miss_seconds"] = (
            lookup_elapsed
        )
        cache_meta["generation_seconds"] = (
            generation_elapsed
        )
        cache_meta["save_seconds"] = (
            save_elapsed
        )
        cache_meta["final_lookup_seconds"] = (
            second_lookup_elapsed
        )
        cache_meta["cache_hit"] = False
        cache_meta["generated"] = True

    # --------------------------------------------------------
    # DECODE
    #
    # Always performed after the final cache lookup.
    # Therefore hit and miss are identical from here onward.
    #
    # Decoding occurs at the caller's current mp.mp.dps.
    # --------------------------------------------------------

    decode_start = time.perf_counter()

    lambda_min, v_canonical = (
        _ground_state_decode(
            results,
            N,
        )
    )

    decode_elapsed = (
        time.perf_counter()
        - decode_start
    )

    # --------------------------------------------------------
    # STRUCTURAL VALIDATION
    #
    # Always performed at working precision.
    # --------------------------------------------------------

    validation_start = time.perf_counter()

    structural = (
        _validate_ground_state_structure(
            lambda_min,
            v_canonical,
            N,
        )
    )

    validation_elapsed = (
        time.perf_counter()
        - validation_start
    )

    # --------------------------------------------------------
    # FINAL METADATA
    # --------------------------------------------------------

    cache_meta = dict(cache_meta)

    cache_meta["generation_dps"] = generation_dps
    cache_meta["working_dps"] = working_dps
    cache_meta["decode_seconds"] = decode_elapsed
    cache_meta["validation_seconds"] = (
        validation_elapsed
    )
    cache_meta["structural_validation"] = structural

    if cache_hit:
        cache_meta["total_seconds"] = (
            lookup_elapsed
            + decode_elapsed
            + validation_elapsed
        )

    else:
        cache_meta["total_seconds"] = (
            cache_meta["initial_cache_miss_seconds"]
            + cache_meta["generation_seconds"]
            + cache_meta["save_seconds"]
            + cache_meta["final_lookup_seconds"]
            + decode_elapsed
            + validation_elapsed
        )

    if verbose:
        print()
        print(
            "GROUND STATE CACHE: RETURN"
        )
        print(
            f"  generation dps = "
            f"{generation_dps}"
        )
        print(
            f"  working dps    = "
            f"{working_dps}"
        )
        print(
            f"  decode         = "
            f"{decode_elapsed:.6f} s"
        )
        print(
            f"  validation     = "
            f"{validation_elapsed:.6f} s"
        )

        if not cache_hit:
            print(
                f"  final lookup   = "
                f"{cache_meta['final_lookup_seconds']:.6f} s"
            )

        print(
            f"  total          = "
            f"{cache_meta['total_seconds']:.6f} s"
        )

    return (
        lambda_min,
        v_canonical,
        cache_meta,
    )


# ============================================================
# CANONICAL FORENSIC GROUND-STATE CONFIGURATION
# ============================================================
#
# All diagnostic cells investigating the current Cell-5
# discrepancy should use these parameters unless they are
# explicitly performing a separate convergence experiment.
#
# `dps` is generation/certification precision.
# The caller's current mp.mp.dps is working precision.
# ============================================================

FORENSIC_GROUND_STATE = {
    "c": 13,
    "N": 8,
    "T": 400,
    "dps": 150,
}

SURVEY_GROUND_STATE = {
    "c": 13,
    "N": 8,
    "T": 400,
    "dps": 50,
}


# ============================================================
# BASIC GEOMETRIC / FOURIER PARAMETERS
# ============================================================

def compute_L(c):
    """
    L = log(c)
    """
    return mp.log(mp.mpf(c))


def compute_beta(L):
    """
    beta = L / (4*pi)
    """
    return L / (4 * mp.pi)


def compute_delta(L):
    """
    Delta = L / (2*pi)
    """
    return L / (2 * mp.pi)


# ============================================================
# PRIME POWERS
# ============================================================

def prime_power_terms(c):
    """
    Return

        [(q, Lambda(q)), ...]

    for all prime powers q <= c.

    Here Lambda(p^k) = log(p).
    """

    c_int = int(mp.floor(c))

    terms = []

    for p in range(2, c_int + 1):

        is_prime = True

        for d in range(2, int(mp.sqrt(p)) + 1):

            if p % d == 0:
                is_prime = False
                break

        if not is_prime:
            continue

        q = p

        while q <= c_int:

            terms.append(
                (
                    mp.mpf(q),
                    mp.log(p)
                )
            )

            q *= p

    return terms


# ============================================================
# CANONICAL <-> FULL SYMMETRIC COORDINATES
# ============================================================

def canonical_to_full(v):
    """
    Convert canonical real-even coordinates

        v = (v_0, v_1, ..., v_N)

    to the full symmetric coefficient vector

        u_{-N}, ..., u_0, ..., u_N

    with

        u_0 = v_0
        u_{+k} = u_{-k} = v_k / sqrt(2).
    """

    N = len(v) - 1
    u = mp.matrix(2 * N + 1, 1)

    for m in range(-N, N + 1):

        if m == 0:
            u[m + N] = v[0]

        else:
            u[m + N] = (
                v[abs(m)]
                / mp.sqrt(2)
            )

    return u


def full_to_canonical(u):
    """
    Convert a symmetric full-space vector to canonical
    real-even coordinates.
    """

    N = int((len(u) - 1) / 2)
    v = mp.matrix(N + 1, 1)

    v[0] = u[N]

    for k in range(1, N + 1):

        v[k] = (
            mp.sqrt(2)
            * u[N + k]
        )

    return v


# ============================================================
# CANONICAL BASIS PAIRS
# ============================================================

def canonical_pairs(k):
    """
    Return the full-space Fourier coefficient pairs belonging
    to canonical basis vector e_k.

    k = 0:
        [(0, 1)]

    k > 0:
        [(+k, 1/sqrt(2)), (-k, 1/sqrt(2))]
    """

    if k == 0:

        return [
            (0, mp.mpf(1))
        ]

    ck = 1 / mp.sqrt(2)

    return [
        ( k, ck),
        (-k, ck),
    ]


# ============================================================
# F BASIS RESPONSE
# ============================================================

def F_basis(k, tau, L):
    """
    Canonical basis response F_k(tau).
    """

    tau = mp.mpf(tau)

    exp_tL = mp.exp(
        -1j * tau * L
    )

    total = mp.mpc(0)

    for kk, ck in canonical_pairs(k):

        denom = (
            2 * mp.pi * kk / L
            - tau
        )

        if denom == 0:

            term = mp.mpc(L)

        elif abs(denom * L) < mp.sqrt(mp.eps):

            term = (
                mp.expm1(1j * denom * L)
                / (1j * denom)
            )

        else:

            term = (
                exp_tL - 1
            ) / (1j * denom)

        total += ck * term

    return mp.re(
        mp.exp(1j * tau * L / 2)
        * total
        / mp.sqrt(L)
    )


def sum_v_F(v, tau, L):
    """
    Coefficient-weighted sum of the canonical F basis responses:

        F_v(tau) = sum_k v_k F_k(tau).

    This is the construction of F_v from the canonical basis
    responses F_k using the entries of v as coefficients.
    """

    return sum(
        v[k] * F_basis(k, tau, L)
        for k in range(len(v))
    )


# ============================================================
# F' BASIS RESPONSE
# ============================================================

def Fprime_basis(k, tau, L):
    """
    Analytic derivative F'_k(tau).
    """

    tau = mp.mpf(tau)

    exp_tL = mp.exp(
        -1j * tau * L
    )

    H = mp.mpc(0)
    Hp = mp.mpc(0)

    for kk, ck in canonical_pairs(k):

        a = 2 * mp.pi * kk / L
        denom = a - tau

        g = (
            exp_tL - 1
        ) / (1j * denom)

        gp = (
            -L * exp_tL / denom
            - 1j * (exp_tL - 1)
            / denom**2
        )

        H += ck * g
        Hp += ck * gp

    return mp.re(
        mp.exp(1j * tau * L / 2)
        * (
            1j * L / 2 * H
            + Hp
        )
        / mp.sqrt(L)
    )


def sum_v_Fprime(v, tau, L):
    """
    Coefficient-weighted sum of the canonical F' basis responses:

        F'_v(tau) = sum_k v_k F'_k(tau).

    This is the construction of F'_v from the canonical basis
    responses F'_k using the entries of v as coefficients.
    """

    return sum(
        v[k]
        * Fprime_basis(k, tau, L)
        for k in range(len(v))
    )


# ============================================================
# G BASIS RESPONSE
# ============================================================

def G_basis_complex(k, z, L):
    """
    Complex analytic response G_k(z).
    """

    z = mp.mpc(z)

    exp_zL = mp.exp(
        -1j * z * L
    )

    total = mp.mpc(0)

    for kk, ck in canonical_pairs(k):

        a = 2 * mp.pi * kk / L
        denom = a - z

        if denom == 0:

            term = mp.mpc(L)

        elif abs(denom * L) < mp.sqrt(mp.eps):

            term = (
                mp.expm1(1j * denom * L)
                / (1j * denom)
            )

        else:

            term = (
                exp_zL - 1
            ) / (1j * denom)

        total += ck * term

    return (
        mp.exp(1j * z * L / 2)
        * total
        / mp.sqrt(L)
    )


def sum_v_G(v, z, L):
    """
    Coefficient-weighted sum of the complex G basis responses:

        G_v(z) = sum_k v_k G_k(z).

    This is the construction of G_v from the canonical basis
    responses G_k using the entries of v as coefficients.

    IMPORTANT:
        G_v(z) is a coefficient-weighted basis sum. It is therefore
        linear in v, but this function is not itself the quadratic
        functional represented by a Galerkin matrix quadratic form
        such as v* Q v.

        Quadratic quantities involving v must be constructed by the
        appropriate quadratic machinery.
    """

    return sum(
        v[k]
        * G_basis_complex(k, z, L)
        for k in range(len(v))
    )


# ============================================================
# POLE FUNCTIONAL
# ============================================================

def pole_basis(k, L):
    """
    P(e_k) for the canonical basis.
    """

    beta = compute_beta(L)

    if k == 0:

        return 1 / beta**2

    return (
        mp.sqrt(2)
        / (k**2 + beta**2)
    )


def pole_row(N, L):
    """
    Canonical coordinate row representing the pole functional P:

        [P(e_0), ..., P(e_N)].
    """

    return mp.matrix([
        pole_basis(k, L)
        for k in range(N + 1)
    ])


# ============================================================
# TRIGONOMETRIC POLYNOMIAL
# ============================================================

def T_canonical(v, t):
    """
    Trigonometric polynomial T_v(t) corresponding to canonical
    coefficient vector v.

    The dependence of T_v on v is through the canonical coefficient
    combination defining the trigonometric polynomial.
    """

    N = len(v) - 1

    total = mp.mpc(0)

    total += v[0]

    for k in range(1, N + 1):

        uk = v[k] / mp.sqrt(2)

        total += (
            uk * mp.exp(
                2j * mp.pi * k * t
            )
            +
            uk * mp.exp(
                -2j * mp.pi * k * t
            )
        )

    return total


# ============================================================
# VOLTERRA SINE-CHORD KERNEL
# ============================================================

def K_canonical(v, omega):
    """
    K_v(omega) =
        2 int_0^omega
            T_v(t) T_v(omega-t) dt

    for 0 <= omega <= 1.

    K_v is constructed from the product of two copies of T_v and
    is therefore quadratic in the coefficient vector v.
    """

    omega = mp.mpf(omega)

    if omega <= 0:
        return mp.mpf(0)

    if omega >= 1:
        raise ValueError(
            "K_canonical expects 0 <= omega <= 1"
        )

    integrand = lambda t: (
        T_canonical(v, t)
        * T_canonical(v, omega - t)
    )

    return 2 * mp.quad(
        integrand,
        [0, omega]
    )


# ============================================================
# ANALYTIC FOURIER REPRESENTATION OF K_v
# ============================================================
#
# The Archimedean calculation can be reduced analytically to
# a single r-integral by expressing the Fourier transform of
# the quadratic K_v kernel as a finite sum over the Fourier
# modes.
#
# For canonical coefficient vector v, let
#
#     u_{-N}, ..., u_0, ..., u_N
#
# be the corresponding full symmetric coefficients.
#
# The following functions implement the analytic ingredients
# of that finite Fourier representation.
#
# These expressions are algebraically equivalent to the
# numerical inner integration used in Cells 21/22.
#
# No mathematical change is intended by this implementation.
# ============================================================


def sinc(z):
    """
    Stable sinc function:

        sinc(z) = sin(z) / z

    with the removable value

        sinc(0) = 1.
    """

    z = mp.mpf(z)

    if z == 0:
        return mp.mpf(1)

    return mp.sin(z) / z


def one_minus_cos(x):
    """
    Stable evaluation of

        1 - cos(x)

    using

        1 - cos(x) = 2 sin^2(x/2).
    """

    x = mp.mpf(x)

    return (
        2
        * mp.sin(x / 2) ** 2
    )


def S_mode(m, r, L):
    """
    Fourier-mode sine/cosine integral:

        S_m(r)
          = integral_0^L
                sin(a_m y) cos(r y)
            dy

    where

        a_m = 2*pi*m/L.

    For m = 0 the value is exactly zero.

    The removable limit when a_m +/- r -> 0 is also handled
    explicitly.
    """

    m = int(m)
    r = mp.mpf(r)
    L = mp.mpf(L)

    a_m = (
        2
        * mp.pi
        * m
        / L
    )

    if m == 0:
        return mp.mpf(0)

    k_plus = a_m + r
    k_minus = a_m - r

    plus = (
        one_minus_cos(
            k_plus * L
        )
        / k_plus
        if k_plus != 0
        else mp.mpf(0)
    )

    minus = (
        one_minus_cos(
            k_minus * L
        )
        / k_minus
        if k_minus != 0
        else mp.mpf(0)
    )

    return (
        plus + minus
    ) / 2


def W_triangular(k, L):
    """
    Fourier transform of the triangular window:

        W(k)
          = integral_0^L
                (1 - y/L) cos(k y)
            dy

    with the stable representation

        W(k)
          = L/2 * sinc(k L/2)^2.

    In particular,

        W(0) = L/2.
    """

    k = mp.mpf(k)
    L = mp.mpf(L)

    return (
        L
        / 2
        * sinc(
            k * L / 2
        ) ** 2
    )


def C_mode(m, r, L):
    """
    Fourier-mode triangular-window integral:

        C_m(r)
          = integral_0^L
                (1 - y/L)
                cos(a_m y)
                cos(r y)
            dy

    where

        a_m = 2*pi*m/L.

    Using the product-to-sum identity,

        C_m(r)
          = 1/2 [
                W_triangular(a_m-r, L)
                + W_triangular(a_m+r, L)
            ].
    """

    m = int(m)
    r = mp.mpf(r)
    L = mp.mpf(L)

    a_m = (
        2
        * mp.pi
        * m
        / L
    )

    return (
        W_triangular(
            a_m - r,
            L,
        )
        + W_triangular(
            a_m + r,
            L,
        )
    ) / 2


def K_fourier(v, r, L):
    """
    Analytic Fourier-side representation of the quadratic K_v
    kernel.

    Parameters
    ----------
    v:
        Canonical real-even coefficient vector

            v = (v_0, ..., v_N).

    r:
        Real spectral variable.

    L:
        Logarithmic interval length, L = log(c).

    Returns
    -------
    mp.mpf
        The value J_v(r) of the Fourier-side representation of the
        quadratic K_v kernel entering the explicit Archimedean
        calculation.

    Construction
    ------------
    Starting from

        K_v(omega)
            = 2 * integral_0^omega
                T_v(t) T_v(omega - t) dt,

    the convolution is evaluated analytically in Fourier modes and
    then transformed to the r-variable. The resulting expression is
    evaluated directly in canonical coordinates:

        J_v(r)
          = 2 sum_{m=0}^N v_m^2 C_m(r)

            - (2 sqrt(2) v_0 / pi)
                sum_{m=1}^N v_m S_m(r) / m

            - (1 / pi)
                sum_{m=1}^N v_m^2 S_m(r) / m

            + (4 / pi)
                sum_{1 <= m < n <= N}
                    v_m v_n
                    * [m S_m(r) - n S_n(r)]
                    / (n^2 - m^2).

    Here

        C_m(r)
          = integral_0^L
                (1 - y/L)
                cos(a_m y) cos(r y) dy,

        S_m(r)
          = integral_0^L
                sin(a_m y) cos(r y) dy,

        a_m = 2 pi m / L.

    The formula exploits the real-even symmetry of the underlying
    full Fourier coefficients. The zero-mode, opposite-sign, and
    positive-mode pair contributions are combined analytically.

    For m < n, the triangular sum represents both signs of the
    Fourier modes after symmetry reduction. The factor 4/pi is the
    complete multiplicity factor for the resulting canonical
    positive-mode pair contribution; no additional factor is applied
    merely because the sum is restricted to m < n.

    Each S_m(r) is evaluated once and reused throughout the
    triangular sum.

    The detailed derivation from the convolution integral is recorded
    in cell_history_map.md.

    IMPORTANT
    ---------
    K_fourier is a quadratic construction in the canonical coefficient
    vector v. It represents the Fourier-side form of the quadratic
    K_v kernel.

    It must not be confused with coefficient-weighted constructions
    such as sum_v_F or sum_v_G.
    """

    r = mp.mpf(r)
    L = mp.mpf(L)
    N = len(v) - 1

    # --------------------------------------------------------
    # Diagonal terms
    # --------------------------------------------------------

    diag = mp.mpf(0)

    for m in range(0, N+1):
        diag += (
            v[m]
            * v[m]
            * C_mode(
                m,
                r,
                L,
            )
        )

    total = 2 * diag

    # --------------------------------------------------------
    # Off-diagonal terms
    # --------------------------------------------------------

    # S_m(r) is common to many terms, so evaluate it once.
    S = {
        m: S_mode(
            m,
            r,
            L,
        )
        for m in range(1, N+1)
    }

    off_diag = mp.mpf(0)
    for m in range(1, N+1):
        off_diag += v[m] * v[m] * S[m] / m
    total -= off_diag / mp.pi

    off_zero = mp.mpf(0)
    for m in range(1, N+1):
        off_zero += v[m] * S[m] / m
    total -= 2 * mp.sqrt(2) * v[0] * off_zero / mp.pi

    off = mp.mpf(0)
    for m in range(1, N):
        for n in range(m+1, N+1):
            off += v[m] * v[n] * ((m * S[m] - n * S[n]) / (n * n - m * m))
    total += 4 * off / mp.pi

    return total


# ============================================================
# ARCHIMEDEAN SOURCE
# ============================================================

def h_plus(r):
    """
    Archimedean source function:

        h_+(r)
          = Re psi(1/4 + i r/2) - log(pi).

    This is the source appearing in the explicit Archimedean
    functional

        A_arch(T)
          = 1/pi * integral_0^T
                h_+(r) K_fourier(v,r,L)
            dr.
    """

    r = mp.mpf(r)

    return (
        mp.re(
            mp.digamma(
                mp.mpf("0.25")
                + 1j * r / 2
            )
        )
        - mp.log(mp.pi)
    )


# ============================================================
# ARCHIMEDEAN SOURCE
# ============================================================
def archimedean_integral(T, v, L):
    return (
        mp.quad(
            lambda r:
                h_plus(r)
                * K_fourier(v, r, L),
            [0, T],
        )
        / mp.pi
    )


# ============================================================
# FOURIER WEIGHT
# ============================================================

def ghat(v, xi, L):
    """
    Fourier weight ghat_v(xi):

        ghat_v(xi) = pi K_v(1 - |xi| / Delta)

    for |xi| <= Delta.

    Since K_v is quadratic in v, ghat_v is also quadratic in v.
    """

    xi = mp.mpf(xi)

    Delta = compute_delta(L)

    if abs(xi) > Delta:
        return mp.mpf(0)

    omega = (
        1
        - abs(xi) / Delta
    )

    return (
        mp.pi
        * K_canonical(v, omega)
    )
