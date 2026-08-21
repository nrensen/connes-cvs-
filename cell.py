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


# ============================================================
# DEFAULT NUMERICAL PARAMETERS
# ============================================================

DEFAULT_DPS = 80


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
