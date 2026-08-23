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
# The cache is deliberately generic.  It knows only:
#
#   namespace
#   parameters
#   results
#
# Mathematical validation belongs to the wrapper for the particular
# calculation being cached.
#
# Parameters and results must be JSON-compatible.
#
# The cache key is SHA-256 over a canonical JSON representation of
# the namespace and parameters.
# ============================================================

CACHE_SCHEMA_VERSION = 2

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

    payload = {
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "cache_key": digest,
        "namespace": str(namespace),
        "parameters": parameters,
        "results": results,
    }

    if timing is not None:
        payload["timing"] = timing

    # Atomic write.  The temporary file is in the same directory so
    # os/filesystem rename semantics remain atomic.
    temporary = path.with_suffix(".tmp")

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

    with path.open(
        "r",
        encoding="utf-8",
    ) as f:
        payload = json.load(f)

    # Validate the cache record itself rather than trusting the
    # filename.
    if payload.get("cache_schema_version") != CACHE_SCHEMA_VERSION:
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


def cache_get_or_compute(
    namespace,
    parameters,
    compute_fn,
    *,
    verbose=True,
):
    """
    Generic cached calculation.

    `compute_fn` must return a JSON-compatible result fragment.

    Returns:

        (results, metadata)

    The metadata contains timing information for both cache hits and
    misses.
    """
    digest, identity = cache_key(
        namespace,
        parameters,
    )

    lookup_start = time.perf_counter()

    try:
        results, metadata = cache_load(
            namespace,
            parameters,
        )

        lookup_elapsed = (
            time.perf_counter()
            - lookup_start
        )

        metadata = dict(metadata)

        metadata["lookup_seconds"] = (
            lookup_elapsed
        )
        metadata["total_seconds"] = (
            lookup_elapsed
        )

        if verbose:
            print()
            print(
                f"CACHE HIT [{namespace}]"
            )
            print(
                f"  key            = {digest}"
            )
            print(
                f"  lookup         = "
                f"{lookup_elapsed:.6f} s"
            )
            print(
                f"  total          = "
                f"{lookup_elapsed:.6f} s"
            )

        return results, metadata

    except FileNotFoundError:
        lookup_elapsed = (
            time.perf_counter()
            - lookup_start
        )

        if verbose:
            print()
            print(
                f"CACHE MISS [{namespace}]"
            )
            print(
                f"  key            = {digest}"
            )

    except ValueError as exc:
        lookup_elapsed = (
            time.perf_counter()
            - lookup_start
        )

        if verbose:
            print()
            print(
                f"CACHE INVALID [{namespace}]"
            )
            print(
                f"  key            = {digest}"
            )
            print(
                f"  reason         = {exc}"
            )
            print(
                "  recalculating"
            )

    compute_start = time.perf_counter()

    results = compute_fn()

    compute_elapsed = (
        time.perf_counter()
        - compute_start
    )

    save_start = time.perf_counter()

    save_timing = {
        "compute_seconds": compute_elapsed,
    }

    digest, path = cache_save(
        namespace,
        parameters,
        results,
        timing=save_timing,
    )

    save_elapsed = (
        time.perf_counter()
        - save_start
    )

    total_elapsed = (
        lookup_elapsed
        + compute_elapsed
        + save_elapsed
    )

    metadata = {
        "cache_hit": False,
        "cache_key": digest,
        "cache_path": path,
        "cache_record": identity,
        "lookup_seconds": lookup_elapsed,
        "compute_seconds": compute_elapsed,
        "save_seconds": save_elapsed,
        "total_seconds": total_elapsed,
        "stored_timing": save_timing,
    }

    if verbose:
        print()
        print(
            f"CACHE COMPUTED [{namespace}]"
        )
        print(
            f"  key            = {digest}"
        )
        print(
            f"  lookup         = "
            f"{lookup_elapsed:.6f} s"
        )
        print(
            f"  computation     = "
            f"{compute_elapsed:.6f} s"
        )
        print(
            f"  save            = "
            f"{save_elapsed:.6f} s"
        )
        print(
            f"  total           = "
            f"{total_elapsed:.6f} s"
        )

    return results, metadata


# ============================================================
# GROUND-STATE CACHE WRAPPER
# ============================================================

GROUND_STATE_OPERATOR_VERSION = (
    "cell.py-ground-state-v1"
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
    v_full,
    dps,
):
    """
    Convert arbitrary-precision numerical results into a
    JSON-compatible fragment.
    """
    digits = int(dps) + 10

    return {
        "lambda_min": mp.nstr(
            lambda_min,
            digits,
        ),
        "v_full": [
            mp.nstr(
                v_full[i, 0],
                digits,
            )
            for i in range(v_full.rows)
        ],
    }


def _ground_state_decode(
    results,
    N,
):
    """
    Reconstruct mpmath numerical objects from cached JSON.
    """
    lambda_min = mp.mpf(
        results["lambda_min"]
    )

    values = results["v_full"]

    expected = 2 * N + 1

    if len(values) != expected:
        raise ValueError(
            "cached ground-state vector has "
            f"length {len(values)}, expected {expected}"
        )

    v_full = mp.matrix(
        expected,
        1,
    )

    for i, value in enumerate(values):
        v_full[i, 0] = mp.mpf(value)

    return lambda_min, v_full


def _validate_ground_state_structure(
    lambda_min,
    v_full,
    N,
):
    """
    Cheap validation which requires no Galerkin matrix.

    This is suitable for the fast cache-hit path.
    """
    expected = 2 * N + 1

    if v_full.rows != expected:
        raise ValueError(
            "ground-state vector has incorrect row count"
        )

    if v_full.cols != 1:
        raise ValueError(
            "ground-state vector must be a column vector"
        )

    if not mp.isfinite(lambda_min):
        raise ValueError(
            "ground-state eigenvalue is not finite"
        )

    for i in range(expected):
        value = v_full[i, 0]

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
        mp.fdot(v_full, v_full)
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


def _validate_ground_state_full(
    Q,
    lambda_min,
    v_full,
    N,
):
    """
    Full mathematical validation against the requested operator.
    """
    structural = (
        _validate_ground_state_structure(
            lambda_min,
            v_full,
            N,
        )
    )

    residual = Q * v_full - (
        lambda_min * v_full
    )

    residual_norm = mp.sqrt(
        mp.fdot(
            residual,
            residual,
        )
    )

    tolerance = structural["tolerance"]

    if residual_norm > tolerance:
        raise ValueError(
            "ground-state eigenvector residual "
            "validation failed: "
            f"{mp.nstr(residual_norm, 10)} > "
            f"{mp.nstr(tolerance, 10)}"
        )

    v_canonical = full_to_canonical(
        v_full,
        N,
    )

    v_roundtrip = canonical_to_full(
        v_canonical,
        N,
    )

    roundtrip_error = mp.sqrt(
        mp.fdot(
            v_roundtrip - v_full,
            v_roundtrip - v_full,
        )
    )

    if roundtrip_error > tolerance:
        raise ValueError(
            "ground-state canonical/full "
            "round-trip validation failed: "
            f"{mp.nstr(roundtrip_error, 10)} > "
            f"{mp.nstr(tolerance, 10)}"
        )

    return {
        **structural,
        "residual_norm": residual_norm,
        "roundtrip_error": roundtrip_error,
    }


def get_ground_state(
    c,
    N,
    T,
    dps,
    *,
    cache=True,
    validation="fast",
    flint_bits=None,
    verbose=True,
):
    """
    Obtain the CvS ground state.

    Parameters
    ----------
    validation : {"fast", "full", "none"}
        fast:
            Validate cached results structurally without rebuilding Q.

        full:
            Rebuild Q and validate
            ||Qv - lambda v|| as well as the canonical/full
            coordinate transformation.

        none:
            Only validate the cache identity and JSON structure.

    Returns
    -------
    lambda_min, v_full, metadata
    """

    caller_dps = mp.mp.dps

    generation_dps = int(dps)

    if generation_dps <= 0:
        raise ValueError(
            "dps must be positive"
        )

    if validation not in (
        "fast",
        "full",
        "none",
    ):
        raise ValueError(
            "validation must be 'fast', 'full', or 'none'"
        )

    parameters = _ground_state_parameters(
        c=c,
        N=N,
        T=T,
        dps=dps,
        flint_bits=flint_bits,
    )

    namespace = "ground_state"

    # --------------------------------------------------------
    # FAST CACHE PATH
    #
    # Crucially, no Galerkin matrix is constructed here.
    # --------------------------------------------------------

    if cache:

        try:
            lookup_start = time.perf_counter()

            results, cache_meta = cache_load(
            	namespace,
            	parameters,
            )

            lookup_elapsed = (
            	time.perf_counter()
            	- lookup_start
            )

            decode_start = time.perf_counter()

            lambda_min, v_full = (
            	_ground_state_decode(
            		results,
            		N,
            	)
            )

            decode_elapsed = (
            	time.perf_counter()
            	- decode_start
            )

            validation_start = time.perf_counter()

            structural = (
            	None
            	if validation == "none"
            	else _validate_ground_state_structure(
            		lambda_min,
            		v_full,
            		N,
            	)
            )

            validation_elapsed = (
            	time.perf_counter()
            	- validation_start
            )

            cache_meta = dict(cache_meta)

            cache_meta["generation_dps"] = generation_dps
            cache_meta["caller_dps"] = caller_dps
            cache_meta["lookup_seconds"] = lookup_elapsed
            cache_meta["decode_seconds"] = decode_elapsed
            cache_meta["validation_seconds"] = validation_elapsed

            cache_meta["total_seconds"] = (
            	lookup_elapsed
            	+ decode_elapsed
            	+ validation_elapsed
            )

            cache_meta["validation_mode"] = validation
            cache_meta["structural_validation"] = structural

            if validation == "full":

            	Q_start = time.perf_counter()

            	Q = build_galerkin_matrix(
            		c=c,
            		N=N,
            		T=T,
            		dps=dps,
            		flint_bits=flint_bits,
            	)

            	Q_elapsed = (
            		time.perf_counter()
            		- Q_start
            	)

            	full_validation_start = time.perf_counter()

            	full_validation = (
            		_validate_ground_state_full(
            			Q,
            			lambda_min,
            			v_full,
            			N,
            		)
            	)

            	full_validation_elapsed = (
            		time.perf_counter()
            		- full_validation_start
            	)

            	cache_meta[
            		"operator_build_seconds"
            	] = Q_elapsed

            	cache_meta[
            		"full_validation_seconds"
            	] = full_validation_elapsed

            	cache_meta[
            		"total_seconds"
            	] = (
            		cache_meta["lookup_seconds"]
            		+ cache_meta["decode_seconds"]
            		+ cache_meta["validation_seconds"]
            		+ Q_elapsed
            		+ full_validation_elapsed
            	)

            	cache_meta[
            		"full_validation"
            	] = full_validation

            if verbose:
                print()
                print(
                    "GROUND STATE CACHE: HIT"
                )
                print(
                    f"  key        = "
                    f"{cache_meta['cache_key']}"
                )

                print(
                    f"  decode     = "
                    f"{cache_meta['decode_seconds']:.6f} s"
                )

                print(
                    f"  validation = "
                    f"{cache_meta['validation_seconds']:.6f} s"
                )

                if validation == "full":
                    print(
                        f"  Q build    = "
                        f"{cache_meta['operator_build_seconds']:.6f} s"
                    )

                print(
                    f"  total      = "
                    f"{cache_meta['total_seconds']:.6f} s"
                )

            return (
                lambda_min,
                v_full,
                cache_meta,
            )

        except FileNotFoundError:
            if verbose:
                print()
                print(
                    "GROUND STATE CACHE: MISS"
                )

        except ValueError as exc:
            if verbose:
                print()
                print(
                    "GROUND STATE CACHE: INVALID"
                )
                print(
                    f"  reason = {exc}"
                )
                print(
                    "  recalculating"
                )

    # --------------------------------------------------------
    # COLD CALCULATION
    # --------------------------------------------------------

    compute_start = time.perf_counter()

    mp.mp.dps = generation_dps

    Q = build_galerkin_matrix(
        c=c,
        N=N,
        T=T,
        dps=dps,
        flint_bits=flint_bits,
    )

    Q_build_elapsed = (
        time.perf_counter()
        - compute_start
    )

    eig_start = time.perf_counter()

    lambda_min, v_full = (
        compute_ground_state(Q)
    )

    eig_elapsed = (
        time.perf_counter()
        - eig_start
    )

    validation_start = time.perf_counter()

    structural = (
        _validate_ground_state_structure(
            lambda_min,
            v_full,
            N,
        )
    )

    full_validation = None

    if validation == "full":
        full_validation = (
            _validate_ground_state_full(
                Q,
                lambda_min,
                v_full,
                N,
            )
        )

    validation_elapsed = (
        time.perf_counter()
        - validation_start
    )

    results = _ground_state_encode(
        lambda_min,
        v_full,
        dps,
    )

    save_start = time.perf_counter()

    digest, path = cache_save(
        namespace,
        parameters,
        results,
        timing={
            "Q_build_seconds": Q_build_elapsed,
            "eigensolve_seconds": eig_elapsed,
            "validation_seconds": validation_elapsed,
        },
    )

    save_elapsed = (
        time.perf_counter()
        - save_start
    )

    total_elapsed = (
        Q_build_elapsed
        + eig_elapsed
        + validation_elapsed
        + save_elapsed
    )

    mp.mp.dps = caller_dps

    metadata = {
        "cache_hit": False,
        "cache_key": digest,
        "cache_path": path,
        "parameters": parameters,
        "generation_dps": generation_dps,
        "caller_dps": caller_dps,
        "validation_mode": validation,
        "Q_build_seconds": Q_build_elapsed,
        "eigensolve_seconds": eig_elapsed,
        "validation_seconds": validation_elapsed,
        "save_seconds": save_elapsed,
        "total_seconds": total_elapsed,
        "structural_validation": structural,
        "full_validation": full_validation,
    }

    if verbose:
        print()
        print(
            "GROUND STATE CACHE: COMPUTED"
        )
        print(
            f"  key        = {digest}"
        )
        print(
            f"  Q build    = "
            f"{Q_build_elapsed:.6f} s"
        )
        print(
            f"  eigensolve = "
            f"{eig_elapsed:.6f} s"
        )
        print(
            f"  validation = "
            f"{validation_elapsed:.6f} s"
        )
        print(
            f"  save       = "
            f"{save_elapsed:.6f} s"
        )
        print(
            f"  total      = "
            f"{total_elapsed:.6f} s"
        )

    return (
        lambda_min,
        v_full,
        metadata,
    )


# ============================================================
# CANONICAL FORENSIC GROUND-STATE CONFIGURATION
#
# All diagnostic cells investigating the current Cell-5
# discrepancy should use these parameters unless they are
# explicitly performing a separate convergence experiment.
# ============================================================

FORENSIC_GROUND_STATE = {
    "c": 13,
    "N": 8,
    "T": 400,
    "dps": 150,
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

def canonical_to_full(v, N):
    """
    Convert canonical real-even coordinates

        v = (v_0, v_1, ..., v_N)

    to the full symmetric coefficient vector

        u_{-N}, ..., u_0, ..., u_N

    with

        u_0 = v_0
        u_{+k} = u_{-k} = v_k / sqrt(2).
    """

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


def full_to_canonical(u, N):
    """
    Convert a symmetric full-space vector to canonical
    real-even coordinates.
    """

    v = mp.matrix(N + 1, 1)

    v[0] = u[N]

    for k in range(1, N + 1):

        v[k] = (
            mp.sqrt(2)
            * u[N + k]
        )

    return v


def canonical_norm(v):
    """
    Euclidean norm in canonical coordinates.
    """

    return mp.sqrt(
        mp.fdot(v, v)
    )


# ============================================================
# GROUND-STATE NORMALISATION
# ============================================================

def normalise_ground_state(eigvec, N):
    """
    Convert a repository ground-state eigenvector into the
    canonical real-even normalisation used by the cells.
    """

    coefficients = [
        mp.mpf(eigvec[i, 0])
        for i in range(eigvec.rows)
    ]

    norm = mp.sqrt(
        sum(x * x for x in coefficients)
    )

    coefficients = [
        x / norm
        for x in coefficients
    ]

    v = [coefficients[N]]

    for k in range(1, N + 1):

        v.append(
            mp.sqrt(2)
            * coefficients[N + k]
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


def F_vector(v, tau, L):
    """
    F_v(tau) for a canonical vector v.
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


def Fprime_vector(v, tau, L):
    """
    F'_v(tau) for a canonical vector v.
    """

    return sum(
        v[k]
        * Fprime_basis(k, tau, L)
        for k in range(len(v))
    )


# ============================================================
# COMPLEX G BASIS RESPONSE
# ============================================================

def G_complex_basis(k, z, L):
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


def G_complex(v, z, L):
    """
    Complex analytic response G_v(z).
    """

    return sum(
        v[k]
        * G_complex_basis(k, z, L)
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
    Canonical pole functional row.
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
    Trigonometric polynomial corresponding to canonical vector v.
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
# FOURIER WEIGHT
# ============================================================

def ghat(v, xi, L):
    """
    ghat(xi) = pi K(1 - |xi| / Delta)

    for |xi| <= Delta.
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
