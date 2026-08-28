# ============================================================
# CELL 23 — OPTIMISED ANALYTIC ARCHIMEDEAN REDUCTION
#
# Purpose
# -------
# Cell 22 demonstrated that the y-integral in Cell 21 can be
# performed analytically, reducing the nested quadrature to a
# single r-integral.
#
# Cell 23 makes that same calculation more efficient and
# numerically better conditioned:
#
#   1. exploit (m,n) <-> (n,m) symmetry;
#   2. evaluate each Fourier-mode integral only once;
#   3. use sinc-style expressions for removable singularities.
#
# No mathematical change is intended relative to Cell 22.
#
# This is the candidate efficient implementation that will later
# be compared against:
#
#   Cell 22 @ 40 dps
#   Cell 21 @ 40 dps
#
# The latter remains the independent nested numerical control.
# ============================================================

import time

import mpmath as mp

from cell import (
    FORENSIC_GROUND_STATE,
    get_ground_state,
    canonical_to_full,
    full_to_canonical,
    compute_L,
)


# ============================================================
# PARAMETERS
# ============================================================

WORKING_DPS = 60
T = 60

mp.mp.dps = WORKING_DPS

c = FORENSIC_GROUND_STATE["c"]
N = FORENSIC_GROUND_STATE["N"]

L = compute_L(c)

DISPLAY_DIGITS = 120


def nstr(x):
    return mp.nstr(x, DISPLAY_DIGITS)


def elapsed(start):
    return time.perf_counter() - start


# ============================================================
# HEADER
# ============================================================

print("=" * 78)
print("CELL 23 — OPTIMISED ANALYTIC ARCHIMEDEAN REDUCTION")
print("=" * 78)

print()
print("Parameters:")
print(f"  c              = {c}")
print(f"  N              = {N}")
print(f"  T              = {T}")
print(f"  working_dps    = {WORKING_DPS}")

print()
print("Forensic ground state:")
print(f"  c              = {FORENSIC_GROUND_STATE['c']}")
print(f"  N              = {FORENSIC_GROUND_STATE['N']}")
print(f"  T              = {FORENSIC_GROUND_STATE['T']}")
print(f"  generation_dps = {FORENSIC_GROUND_STATE['dps']}")

print()
print("L =")
print(nstr(L))


# ============================================================
# 1. FORENSIC GROUND STATE
# ============================================================

print()
print("-" * 78)
print("1. FORENSIC GROUND STATE")
print("-" * 78)

ground_start = time.perf_counter()

lambda_forensic, u_star, ground_meta = get_ground_state(
    **FORENSIC_GROUND_STATE,
    verbose=True,
)

ground_elapsed = elapsed(ground_start)

u_star = mp.matrix(u_star)

print()
print("lambda_forensic =")
print(nstr(lambda_forensic))

print()
print("||u_star|| =")
print(
    nstr(
        mp.sqrt(
            mp.fdot(
                u_star,
                u_star,
            )
        )
    )
)

print()
print(
    f"ground-state retrieval elapsed = "
    f"{ground_elapsed:.6f} s"
)


# ============================================================
# 2. CANONICAL / FULL REPRESENTATIONS
# ============================================================

print()
print("-" * 78)
print("2. CANONICAL / FULL REPRESENTATIONS")
print("-" * 78)

v_star = full_to_canonical(
    u_star,
    N,
)

u = canonical_to_full(
    v_star,
    N,
)

print()
print("||v_star|| =")
print(
    nstr(
        mp.sqrt(
            mp.fdot(
                v_star,
                v_star,
            )
        )
    )
)

print()
print("||u - u_star|| =")
print(
    nstr(
        mp.sqrt(
            mp.fdot(
                u - u_star,
                u - u_star,
            )
        )
    )
)


# ============================================================
# 3. FOURIER MODE DATA
#
# Define
#
#     a_m = 2*pi*m/L.
#
# The full coefficients are real for the present ground state.
#
# We precompute a_m and the coefficient products needed by the
# symmetric double sum.
# ============================================================

modes = list(
    range(-N, N + 1)
)

a = {
    m: (
        2
        * mp.pi
        * m
        / L
    )
    for m in modes
}

coeff = {
    m: u[i]
    for i, m in enumerate(modes)
}


# ============================================================
# 4. STABLE ELEMENTARY FUNCTIONS
#
# The Cell-22 expressions contain:
#
#     1 - cos(x)
#
# and divisions by quantities which can approach zero.
#
# We use
#
#     1 - cos(x) = 2 sin^2(x/2)
#
# and
#
#     sin(z)/z
#
# explicitly to make the removable limits transparent.
# ============================================================


def sinc(z):
    """
    sin(z) / z, with the removable value 1 at z=0.
    """

    z = mp.mpf(z)

    if z == 0:
        return mp.mpf("1")

    return mp.sin(z) / z


def one_minus_cos(x):
    """
    1 - cos(x), evaluated through 2 sin^2(x/2).
    """

    x = mp.mpf(x)

    return (
        2
        * mp.sin(x / 2) ** 2
    )


# ============================================================
# 5. S_m(r)
#
# Cell 22 used
#
#   S(a,r)
#     = int_0^L sin(a*y) cos(r*y) dy
#
# with
#
#   a = 2*pi*m/L.
#
# Because a*L = 2*pi*m,
#
#   S(a,r)
#     =
#       1/2 [
#         (1-cos((a+r)L))/(a+r)
#         +
#         (1-cos((a-r)L))/(a-r)
#       ].
#
# We retain this exact algebraic form, but evaluate the
# numerators using sin^2.
#
# At a = 0 the result is exactly zero.
# ============================================================


def S_mode(m, r):

    am = a[m]

    if m == 0:
        return mp.mpf("0")

    k_plus = am + r
    k_minus = am - r

    plus = (
        one_minus_cos(
            k_plus * L
        )
        / k_plus
        if k_plus != 0
        else mp.mpf("0")
    )

    minus = (
        one_minus_cos(
            k_minus * L
        )
        / k_minus
        if k_minus != 0
        else mp.mpf("0")
    )

    return (
        plus + minus
    ) / 2


# ============================================================
# 6. W(k)
#
# Cell 22 used
#
#   W(k)
#     = (1-cos(kL))/(L*k^2)
#
# for
#
#   int_0^L (1-y/L) cos(k*y) dy.
#
# Rewrite:
#
#   W(k)
#     = L/2 * sinc(kL/2)^2.
#
# This makes the k -> 0 limit exactly
#
#   W(0) = L/2.
# ============================================================


def W(k):

    k = mp.mpf(k)

    return (
        L
        / 2
        * sinc(
            k * L / 2
        ) ** 2
    )


# ============================================================
# 7. C_m(r)
#
#     C(a,r)
#       = int_0^L
#           (1-y/L)
#           cos(a*y)
#           cos(r*y)
#         dy
#
#     = 1/2 [ W(a-r) + W(a+r) ].
#
# ============================================================


def C_mode(m, r):

    am = a[m]

    return (
        W(am - r)
        + W(am + r)
    ) / 2


# ============================================================
# 8. ANALYTIC J_v(r)
#
# Cell 22:
#
#   J_v(r)
#     = sum_{m,n} diagonal + off-diagonal.
#
# For m=n:
#
#   contribution = 2 u_m^2 C_m(r).
#
# For m != n:
#
#   contribution =
#
#     u_m u_n / pi
#       * [S_n(r) - S_m(r)]
#       / (m-n).
#
# The (m,n) and (n,m) terms are identical.
#
# Therefore:
#
#   off-diagonal
#     = 2 sum_{m<n}
#         u_m u_n / pi
#         * [S_n-S_m]/(m-n).
#
# This reduces the number of off-diagonal terms from
#
#     (2N+1)^2 - (2N+1)
#
# to
#
#     N_modes * (N_modes-1) / 2.
#
# For N=8:
#
#     272  -> 136.
#
# More importantly, each S_m(r) is calculated only once.
# ============================================================


def analytic_J(r):

    r = mp.mpf(r)

    # S_m(r) is common to many pair terms.
    S = {
        m: S_mode(m, r)
        for m in modes
    }

    total = mp.mpf("0")

    # --------------------------------------------------------
    # Diagonal
    # --------------------------------------------------------

    for m in modes:

        um = coeff[m]

        total += (
            2
            * um
            * um
            * C_mode(m, r)
        )

    # --------------------------------------------------------
    # Off-diagonal
    #
    # Only m < n is required because the two orientations
    # contribute identically.
    # --------------------------------------------------------

    for i, m in enumerate(modes):

        um = coeff[m]

        for n in modes[i + 1:]:

            un = coeff[n]

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
# 9. h_+(r)
# ============================================================


def h_plus(r):

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
# 10. EXPLICIT ARCHIMEDEAN FUNCTIONAL
#
# Only the outer r-integral remains numerical:
#
#   A_arch
#     = 1/pi int_0^T h_+(r) J_v(r) dr.
#
# ============================================================

print()
print("-" * 78)
print("3. OPTIMISED ANALYTIC ARCHIMEDEAN")
print("-" * 78)

print()
print(
    "Computing:"
)
print(
    "  J_v(r) from the symmetric finite Fourier sum"
)
print(
    "  A_arch = (1/pi) int_0^T h_+(r) J_v(r) dr"
)
print()
print(
    "No numerical y-integration is performed."
)

arch_start = time.perf_counter()

explicit_arch = (
    mp.quad(
        lambda r:
            h_plus(r)
            * analytic_J(r),
        [0, T],
    )
    / mp.pi
)

arch_elapsed = elapsed(
    arch_start
)

print()
print("Analytic Archimedean =")
print(nstr(explicit_arch))

print()
print(
    f"analytic Archimedean elapsed = "
    f"{arch_elapsed:.6f} s"
)


# ============================================================
# 11. FINAL RESULT
# ============================================================

print()
print("=" * 78)
print("CELL 23 COMPLETE")
print("=" * 78)

print()
print("Result:")
print(
    "  analytic Archimedean =",
    nstr(explicit_arch),
)

print()
print("Timing:")
print(
    f"  ground state         = "
    f"{ground_elapsed:.6f} s"
)
print(
    f"  analytic Archimedean = "
    f"{arch_elapsed:.6f} s"
)

print()
print("=" * 78)
print("CELL 23 COMPLETE")
print("=" * 78)
