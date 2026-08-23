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
from pathlib import Path

# ============================================================
# DEFAULT NUMERICAL PARAMETERS
# ============================================================

DEFAULT_DPS = 80

# ============================================================
# GROUND-STATE CACHE
# ============================================================
#
# Persistent cache for expensive ground-state calculations.
#
# The cache is deliberately self-describing and content-addressed.
# The hash incorporates all parameters which define the calculation,
# together with explicit schema/operator versions.
#
# Cached numerical values are stored as decimal strings so that JSON
# serialization does not introduce binary floating-point rounding.
#
# No existing mathematical calculation is changed by this machinery.
# ============================================================

GROUND_STATE_CACHE_VERSION = 1

# Increment this whenever the mathematical/operator construction
# changes in a way which can invalidate previously computed states.
GROUND_STATE_OPERATOR_VERSION = "cell.py-v1"

GROUND_STATE_CACHE_DIR = Path(__file__).resolve().parent / "ground_state_cache"


def _ground_state_parameter_record(
    c,
    N,
    T,
    dps,
    flint_bits=None,
):
    """
    Return the canonical parameter record defining a ground-state
    calculation.

    c is represented as a decimal string rather than a binary float.
    """
    if isinstance(c, float):
        raise TypeError(
            "Ground-state cache parameter c must not be a Python float; "
            "pass an integer, decimal string, or mp.mpf."
        )

    c_mp = mp.mpf(c)

    return {
        "cache_version": GROUND_STATE_CACHE_VERSION,
        "operator_version": GROUND_STATE_OPERATOR_VERSION,
        "c": mp.nstr(c_mp, max(50, int(dps) + 10)),
        "N": int(N),
        "T": int(T),
        "dps": int(dps),
        "flint_bits": (
            int(flint_bits)
            if flint_bits is not None
            else int(int(dps) * 3.5)
        ),
    }


def ground_state_cache_key(
    c,
    N,
    T,
    dps,
    flint_bits=None,
):
    """
    Return the SHA-256 cache key for a ground-state calculation.

    The returned tuple is:

        (hex_digest, parameter_record)
    """
    record = _ground_state_parameter_record(
        c=c,
        N=N,
        T=T,
        dps=dps,
        flint_bits=flint_bits,
    )

    payload = json.dumps(
        record,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")

    digest = hashlib.sha256(payload).hexdigest()

    return digest, record


def _ground_state_cache_path(digest):
    """
    Return the cache filename corresponding to a cache digest.
    """
    return GROUND_STATE_CACHE_DIR / f"{digest}.json"


def _matrix_residual(Q, lambda_min, v_full):
    """
    Euclidean residual ||Q v - lambda v||.
    """
    residual = Q * v_full - lambda_min * v_full

    return mp.sqrt(
        mp.fdot(residual, residual)
    )


def _validate_ground_state(
    Q,
    lambda_min,
    v_full,
    N,
    *,
    label="ground state",
):
    """
    Validate a ground state against its operator.

    Returns a diagnostic dictionary. Raises ValueError if the state
    is structurally invalid or numerically inconsistent.

    The tolerances are deliberately tied to the active arithmetic
    precision rather than being fixed decimal thresholds.
    """
    DIM = 2 * N + 1

    if not hasattr(v_full, "rows") or not hasattr(v_full, "cols"):
        raise ValueError(
            f"{label}: eigenvector is not an mpmath matrix"
        )

    if v_full.rows != DIM or v_full.cols != 1:
        raise ValueError(
            f"{label}: expected {DIM}x1 eigenvector, got "
            f"{v_full.rows}x{v_full.cols}"
        )

    if not mp.isfinite(lambda_min):
        raise ValueError(
            f"{label}: eigenvalue is not finite"
        )

    for i in range(DIM):
        if mp.im(v_full[i, 0]) != 0:
            raise ValueError(
                f"{label}: eigenvector has non-real entry at {i}"
            )
        if not mp.isfinite(v_full[i, 0]):
            raise ValueError(
                f"{label}: eigenvector has non-finite entry at {i}"
            )

    norm_full = mp.sqrt(
        mp.fdot(v_full, v_full)
    )

    if norm_full == 0:
        raise ValueError(
            f"{label}: eigenvector has zero norm"
        )

    norm_error = abs(norm_full - 1)

    # Reconstruct canonical coordinates and then return to full space.
    v_canonical = full_to_canonical(v_full, N)
    v_roundtrip = canonical_to_full(v_canonical, N)

    roundtrip_error = mp.sqrt(
        mp.fdot(
            v_roundtrip - v_full,
            v_roundtrip - v_full,
        )
    )

    residual = _matrix_residual(
        Q,
        lambda_min,
        v_full,
    )

    # These are deliberately generous relative to mp.eps. The purpose
    # is to reject corrupted/stale cache entries, not to establish a
    # new eigenvalue theorem.
    tol = mp.sqrt(mp.eps)

    if norm_error > tol:
        raise ValueError(
            f"{label}: norm validation failed: "
            f"{mp.nstr(norm_error, 8)} > {mp.nstr(tol, 8)}"
        )

    if roundtrip_error > tol:
        raise ValueError(
            f"{label}: canonical/full round-trip failed: "
            f"{mp.nstr(roundtrip_error, 8)} > {mp.nstr(tol, 8)}"
        )

    if residual > tol:
        raise ValueError(
            f"{label}: eigenvector residual failed: "
            f"{mp.nstr(residual, 8)} > {mp.nstr(tol, 8)}"
        )

    return {
        "norm_full": norm_full,
        "norm_error": norm_error,
        "roundtrip_error": roundtrip_error,
        "residual_norm": residual,
        "validation_tolerance": tol,
    }


def save_ground_state(
    c,
    N,
    T,
    dps,
    lambda_min,
    v_full,
    *,
    Q=None,
    flint_bits=None,
):
    """
    Validate and persist a ground state.

    If Q is supplied, validation includes the eigenvector residual
    against that exact operator. If Q is omitted, structural
    validation is still performed, but the operator residual cannot
    be checked.
    """
    digest, record = ground_state_cache_key(
        c=c,
        N=N,
        T=T,
        dps=dps,
        flint_bits=flint_bits,
    )

    if Q is not None:
        validation = _validate_ground_state(
            Q,
            lambda_min,
            v_full,
            N,
            label="new ground state",
        )
    else:
        DIM = 2 * N + 1

        if v_full.rows != DIM or v_full.cols != 1:
            raise ValueError(
                "new ground state: incorrect eigenvector dimensions"
            )

        norm_full = mp.sqrt(mp.fdot(v_full, v_full))

        validation = {
            "norm_full": norm_full,
            "norm_error": abs(norm_full - 1),
            "roundtrip_error": mp.mpf("nan"),
            "residual_norm": mp.mpf("nan"),
            "validation_tolerance": mp.sqrt(mp.eps),
        }

    payload = {
        "cache_key": digest,
        "parameters": record,
        "lambda_min": mp.nstr(lambda_min, int(dps) + 10),
        "v_full": [
            mp.nstr(v_full[i, 0], int(dps) + 10)
            for i in range(v_full.rows)
        ],
        "validation": {
            "norm_full": mp.nstr(
                validation["norm_full"],
                int(dps) + 10,
            ),
            "norm_error": mp.nstr(
                validation["norm_error"],
                int(dps) + 10,
            ),
            "roundtrip_error": mp.nstr(
                validation["roundtrip_error"],
                int(dps) + 10,
            ),
            "residual_norm": mp.nstr(
                validation["residual_norm"],
                int(dps) + 10,
            ),
            "validation_tolerance": mp.nstr(
                validation["validation_tolerance"],
                int(dps) + 10,
            ),
        },
    }

    GROUND_STATE_CACHE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    path = _ground_state_cache_path(digest)

    # Atomic write: write a temporary file in the same directory,
    # then replace the final path.
    temporary = path.with_suffix(".tmp")

    with temporary.open("w", encoding="utf-8") as f:
        json.dump(
            payload,
            f,
            indent=2,
            sort_keys=True,
        )
        f.write("\n")

    temporary.replace(path)

    return digest, path


def load_ground_state(
    c,
    N,
    T,
    dps,
    *,
    Q=None,
    flint_bits=None,
):
    """
    Load and validate a cached ground state.

    Returns:

        (lambda_min, v_full, validation, path)

    Raises FileNotFoundError if the cache entry does not exist.

    Raises ValueError if the cache entry exists but does not exactly
    correspond to the requested calculation or fails validation.
    """
    digest, record = ground_state_cache_key(
        c=c,
        N=N,
        T=T,
        dps=dps,
        flint_bits=flint_bits,
    )

    path = _ground_state_cache_path(digest)

    if not path.exists():
        raise FileNotFoundError(path)

    try:
        with path.open("r", encoding="utf-8") as f:
            payload = json.load(f)
    except Exception as exc:
        raise ValueError(
            f"Could not read ground-state cache {path}: {exc}"
        ) from exc

    if payload.get("cache_key") != digest:
        raise ValueError(
            "Ground-state cache key mismatch"
        )

    stored_parameters = payload.get("parameters")

    if stored_parameters != record:
        raise ValueError(
            "Ground-state cache parameters do not match request"
        )

    try:
        lambda_min = mp.mpf(
            payload["lambda_min"]
        )

        stored_vector = payload["v_full"]

        if len(stored_vector) != 2 * N + 1:
            raise ValueError(
                "cached eigenvector has incorrect dimension"
            )

        v_full = mp.matrix(
            [
                mp.mpf(value)
                for value in stored_vector
            ]
        )

        v_full = mp.matrix(
            2 * N + 1,
            1,
        )

        for i, value in enumerate(stored_vector):
            v_full[i, 0] = mp.mpf(value)

    except Exception as exc:
        raise ValueError(
            f"Invalid numerical data in ground-state cache {path}: {exc}"
        ) from exc

    validation = _validate_ground_state(
        Q,
        lambda_min,
        v_full,
        N,
        label="cached ground state",
    ) if Q is not None else None

    return (
        lambda_min,
        v_full,
        validation,
        path,
    )


def get_ground_state(
    c,
    N,
    T,
    dps,
    *,
    cache=True,
    validate=True,
    flint_bits=None,
    verbose=True,
):
    """
    Obtain the ground state, using the persistent cache when enabled.

    Returns:

        (lambda_min, v_full, metadata)

    The metadata dictionary records whether the result came from the
    cache and contains validation diagnostics.

    Existing mathematical routines are used unchanged:
        build_galerkin_matrix()
        compute_ground_state()
    """
    digest, record = ground_state_cache_key(
        c=c,
        N=N,
        T=T,
        dps=dps,
        flint_bits=flint_bits,
    )

    if cache:
        try:
            # Build Q first. This is intentionally conservative:
            # validation of a cache entry is against the exact operator
            # requested by this call.
            mp.mp.dps = int(dps)

            Q = build_galerkin_matrix(
                c=c,
                N=N,
                T=T,
                dps=dps,
                flint_bits=flint_bits,
            )

            lambda_min, v_full, validation, path = load_ground_state(
                c=c,
                N=N,
                T=T,
                dps=dps,
                Q=Q if validate else None,
                flint_bits=flint_bits,
            )

            if verbose:
                print("GROUND STATE CACHE: HIT")
                print(f"  key = {digest}")
                print(f"  path = {path}")

            return (
                lambda_min,
                v_full,
                {
                    "cache_hit": True,
                    "cache_key": digest,
                    "cache_path": path,
                    "parameters": record,
                    "validation": validation,
                },
            )

        except FileNotFoundError:
            if verbose:
                print("GROUND STATE CACHE: MISS")
                print(f"  key = {digest}")

        except ValueError as exc:
            if verbose:
                print("GROUND STATE CACHE: INVALID")
                print(f"  reason = {exc}")
                print("  recalculating ground state")

    # Cache miss, cache disabled, or invalid cache.
    mp.mp.dps = int(dps)

    Q = build_galerkin_matrix(
        c=c,
        N=N,
        T=T,
        dps=dps,
        flint_bits=flint_bits,
    )

    lambda_min, v_full = compute_ground_state(Q)

    validation = (
        _validate_ground_state(
            Q,
            lambda_min,
            v_full,
            N,
            label="fresh ground state",
        )
        if validate
        else None
    )

    cache_path = None

    if cache:
        _, cache_path = save_ground_state(
            c=c,
            N=N,
            T=T,
            dps=dps,
            lambda_min=lambda_min,
            v_full=v_full,
            Q=Q if validate else None,
            flint_bits=flint_bits,
        )

        if verbose:
            print("GROUND STATE CACHE: SAVED")
            print(f"  key = {digest}")
            print(f"  path = {cache_path}")

    return (
        lambda_min,
        v_full,
        {
            "cache_hit": False,
            "cache_key": digest,
            "cache_path": cache_path,
            "parameters": record,
            "validation": validation,
        },
    )

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
