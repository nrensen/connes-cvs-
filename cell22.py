# ============================================================
# CELL 22 — ANALYTIC ARCHIMEDEAN REDUCTION
#
# Purpose
# -------
# Replace the expensive nested numerical integration in Cell 21
# by analytically evaluating the y-integral.
#
# Cell 17 established:
#
#     D_v(omega) = pi K_v(omega)
#
# with
#
#     D_v(omega)
#       = sum_{m,n} u_m u_n
#           [sin(2*pi*m*omega) - sin(2*pi*n*omega)]
#           / (m-n)
#
# and the diagonal interpreted by continuity.
#
# With
#
#     omega = 1 - y/L
#
# and integer m,n, the y-integral reduces to elementary
# trigonometric integrals.
#
# The expensive nested quadrature of Cell 21 is therefore
# replaced by a single numerical integral over r.
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

WORKING_DPS = 20
T = 60

mp.mp.dps = WORKING_DPS

c = FORENSIC_GROUND_STATE["c"]
N = FORENSIC_GROUND_STATE["N"]

L = compute_L(c)

DISPLAY_DIGITS = 80


def nstr(x):
    return mp.nstr(x, DISPLAY_DIGITS)


def elapsed(start):
    return time.perf_counter() - start


# ============================================================
# HEADER
# ============================================================

print("=" * 78)
print("CELL 22 — ANALYTIC ARCHIMEDEAN REDUCTION")
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
# 2. CANONICAL VECTOR / FULL COEFFICIENTS
# ============================================================

v_star = full_to_canonical(
    u_star,
    N,
)

u = canonical_to_full(
    v_star,
    N,
)

print()
print("Canonical/full consistency:")
print(
    "  ||u - u_star|| =",
    nstr(
        mp.sqrt(
            mp.fdot(
                u - u_star,
                u - u_star,
            )
        )
    ),
)


# ============================================================
# 3. ANALYTIC ELEMENTARY INTEGRALS
# ============================================================
#
# We require:
#
#   S(a,r) =
#       int_0^L sin(a*y) cos(r*y) dy
#
# and
#
#   C(a,r) =
#       int_0^L (1-y/L) cos(a*y) cos(r*y) dy.
#
# The formulas are written with explicit limiting cases so that
# r = +/- a does not create a removable 0/0 singularity.
# ============================================================


def sinc_integral(k, L):
    """
    int_0^L sin(k*y) dy
    """
    k = mp.mpf(k)

    if abs(k) < mp.mpf("1e-30"):
        return mp.mpf("0")

    return (
        1 - mp.cos(k * L)
    ) / k


def cos_integral(k, L):
    """
    int_0^L cos(k*y) dy
    """
    k = mp.mpf(k)

    if abs(k) < mp.mpf("1e-30"):
        return L

    return (
        mp.sin(k * L)
        / k
    )


def sin_cos_integral(a, r):
    """
    S(a,r) =
        int_0^L sin(a*y) cos(r*y) dy.

    Using

        sin(a*y) cos(r*y)
          = 1/2 [
                sin((a+r)y)
                + sin((a-r)y)
            ].
    """

    return (
        sinc_integral(
            a + r,
            L,
        )
        +
        sinc_integral(
            a - r,
            L,
        )
    ) / 2


def linear_cos_cos_integral(a, r):
    """
    C(a,r) =
        int_0^L
            (1-y/L)
            cos(a*y)
            cos(r*y)
        dy.

    Using

        cos(a*y) cos(r*y)
          = 1/2 [
                cos((a-r)y)
                + cos((a+r)y)
            ]

    and

        int_0^L (1-y/L) cos(k*y) dy
          = (1-cos(kL))/(L*k^2)

    for k != 0, with limiting value L/2 at k=0.
    """

    def weighted_cos_integral(k):

        k = mp.mpf(k)

        if abs(k) < mp.mpf("1e-30"):
            return L / 2

        return (
            1 - mp.cos(k * L)
        ) / (
            L * k * k
        )

    return (
        weighted_cos_integral(a - r)
        +
        weighted_cos_integral(a + r)
    ) / 2


# ============================================================
# 4. ANALYTIC J_v(r)
# ============================================================
#
# Define
#
#     J_v(r)
#       = int_0^L K_v(1-y/L) cos(r*y) dy.
#
# For m != n:
#
#   K_mn(1-y/L)
#     =
#       u_m u_n / pi
#       *
#       [
#         -sin(a_m*y)
#         +sin(a_n*y)
#       ]
#       /(m-n)
#
# where
#
#     a_m = 2*pi*m/L.
#
# Hence:
#
#   J_mn(r)
#     =
#       u_m u_n / pi
#       *
#       [
#         -S(a_m,r)
#         +S(a_n,r)
#       ]
#       /(m-n).
#
# For m=n:
#
#   K_nn(1-y/L)
#     =
#       2(1-y/L) cos(a_n*y)
#
# and therefore:
#
#   J_nn(r)
#     =
#       2 C(a_n,r).
#
# ============================================================


def analytic_J(r):

    r = mp.mpf(r)

    total = mp.mpf("0")

    for i, m in enumerate(
        range(-N, N + 1)
    ):

        um = u[i]

        a_m = (
            2
            * mp.pi
            * m
            / L
        )

        for j, n in enumerate(
            range(-N, N + 1)
        ):

            un = u[j]

            if m == n:

                total += (
                    um
                    * un
                    * 2
                    * linear_cos_cos_integral(
                        a_m,
                        r,
                    )
                )

            else:

                a_n = (
                    2
                    * mp.pi
                    * n
                    / L
                )

                total += (
                    um
                    * un
                    / mp.pi
                    * (
                        -sin_cos_integral(
                            a_m,
                            r,
                        )
                        + sin_cos_integral(
                            a_n,
                            r,
                        )
                    )
                    / mp.mpf(m - n)
                )

    return total


# ============================================================
# 5. h_+(r)
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
# 6. EXPLICIT ARCHIMEDEAN FUNCTIONAL
#
# Now there is only ONE numerical quadrature:
#
#     A_arch =
#       1/pi int_0^T h_+(r) J_v(r) dr.
#
# There is no numerical integration over y.
# ============================================================

print()
print("-" * 78)
print("2. ANALYTIC ARCHIMEDEAN REDUCTION")
print("-" * 78)

print()
print(
    "Computing:"
)
print(
    "  J_v(r) analytically"
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
# 7. COMPARISON WITH CELL 21
# ============================================================
#
# Cell 21 at 20 dps gave:
#
#   -1.65903308749093566911612062081
#
# This value is deliberately recorded here as an external
# validation target rather than recalculating the expensive
# nested integral.
# ============================================================

CELL21_ARCH_20 = mp.mpf(
    "-1.65903308749093566911612062081"
)

difference = (
    explicit_arch
    - CELL21_ARCH_20
)

relative_difference = (
    abs(difference)
    / max(
        abs(explicit_arch),
        abs(CELL21_ARCH_20),
        mp.mpf("1"),
    )
)

print()
print("-" * 78)
print("3. COMPARISON WITH CELL 21")
print("-" * 78)

print()
print("Cell 21 direct K-fourier =")
print(nstr(CELL21_ARCH_20))

print()
print("Cell 22 analytic reduction =")
print(nstr(explicit_arch))

print()
print("difference =")
print(nstr(difference))

print()
print("relative difference =")
print(nstr(relative_difference))


# ============================================================
# 8. TIMING / FINAL SUMMARY
# ============================================================

print()
print("=" * 78)
print("CELL 22 COMPLETE")
print("=" * 78)

print()
print("Results:")
print(
    "  Cell 21 direct K-fourier =",
    nstr(CELL21_ARCH_20),
)
print(
    "  Cell 22 analytic         =",
    nstr(explicit_arch),
)
print(
    "  difference               =",
    nstr(difference),
)
print(
    "  relative difference      =",
    nstr(relative_difference),
)

print()
print("Timing:")
print(
    f"  ground state             = "
    f"{ground_elapsed:.6f} s"
)
print(
    f"  analytic Archimedean     = "
    f"{arch_elapsed:.6f} s"
)

print()
print("=" * 78)
print("CELL 22 COMPLETE")
print("=" * 78)
