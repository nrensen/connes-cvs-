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
    # filesystem rename semantics remain atomic.
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
    v_full,
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

    Decoding occurs at the caller's current mp.mp.dps. This is
    deliberate: the cached artefact may have been generated at a
    substantially higher precision than the current calculation
    requires.
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
    Cheap intrinsic validation performed on every retrieval.

    This validates the integrity and basic numerical structure of the
    cached result. It does NOT rebuild the Galerkin matrix.

    Mathematical correctness of the ground state is deliberately not
    part of ordinary cache retrieval. If that is ever in doubt, an
    explicit audit should be performed.
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

    Q = build_galerkin_matrix(
        c=c,
        N=N,
        T=T,
        dps=dps,
        flint_bits=flint_bits,
    )

    Q_build_elapsed = (
        time.perf_counter()
        - Q_start
    )

    eig_start = time.perf_counter()

    lambda_min, v_full = (
        compute_ground_state(Q)
    )

    eig_elapsed = (
        time.perf_counter()
        - eig_start
    )

    # Generation-time structural validation.
    validation_start = time.perf_counter()

    structural = (
        _validate_ground_state_structure(
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
    lambda_min, v_full, metadata
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

    # --------------------------------------------------------
    # STRUCTURAL VALIDATION
    #
    # Always performed at working precision.
    # --------------------------------------------------------

    validation_start = time.perf_counter()

    structural = (
        _validate_ground_state_structure(
            lambda_min,
            v_full,
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
        v_full,
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


def K_fourier(u, r, L):
    """
    Analytic Fourier-side representation of the quadratic K_v
    kernel.

    Parameters
    ----------
    u:
        Full symmetric Fourier coefficient vector

            u = (u_{-N}, u_{-N+1}, ..., u_N).

    r:
        Real spectral variable.

    L:
        Logarithmic interval length, L = log(c).

    Returns
    -------
    mp.mpf
        The value of the Fourier representation J_v(r) of the
        quadratic K_v kernel entering the explicit Archimedean
        calculation.

    Construction
    ------------
    The is evaluated as

        J_v(r)
          = sum_m 2 u_m^2 C_m(r)

            + 2 sum_{m<n}
                u_m u_n / pi
                * [S_n(r) - S_m(r)] / (m-n).

    The (m,n) and (n,m) terms are identical, so only m < n
    is evaluated. Each S_m(r) is also evaluated only once.

    IMPORTANT
    ---------
    K_fourier is a quadratic construction in the coefficient
    vector v. It represents the Fourier-side form of the
    quadratic K_v kernel.

    It must not be confused with coefficient-weighted
    constructions such as sum_v_F or sum_v_G.
    """

    r = mp.mpf(r)
    L = mp.mpf(L)

    N = int((len(u) - 1)/ 2)

    modes = list(
        range(-N, N + 1)
    )

    # S_m(r) is common to many terms, so evaluate it once.
    S = {
        m: S_mode(
            m,
            r,
            L,
        )
        for m in modes
    }

    total = mp.mpf(0)

    # --------------------------------------------------------
    # Diagonal terms
    # --------------------------------------------------------

    for m in modes:

        um = u[m + N]

        total += (
            2
            * um
            * um
            * C_mode(
                m,
                r,
                L,
            )
        )

    # --------------------------------------------------------
    # Off-diagonal terms
    #
    # The (m,n) and (n,m) contributions are identical, so only
    # m < n is required.
    # --------------------------------------------------------

    for i, m in enumerate(modes):

        um = u[i]

        for j in range(i + 1, len(modes)):

            n = modes[j]
            un = u[j]

            total += (
                2
                * um
                * un
                / mp.pi
                * (
                    S[n]
                    - S[m]
                )
                / mp.mpf(m - n)
            )

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
                h_+(r) K_fourier(u,r,L)
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
