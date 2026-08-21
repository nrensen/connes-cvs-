# ============================================================
# CELL_6B — ARCHIMEDEAN QUADRATURE / BREAKPOINT AUDIT
#
# Purpose:
#
# Determine whether the discrepancy seen in Cell 6 is caused
# by the way the oscillatory Archimedean integral is split.
#
# We compare three evaluations of
#
#     A_T(x) =
#         1/pi^2 integral_0^T
#             h_+(r) S(r,x,L) dr
#
# 1. CELL-6 STYLE:
#       split into fixed intervals of length 5.
#
# 2. BREAKPOINT STYLE:
#       split at alpha_x = 2*pi*x/L.
#
# 3. REPOSITORY STYLE:
#       integrate over [-T,T], with the repository's
#       positive/negative breakpoint structure.
#
# No Galerkin matrix is constructed.
# No eigenvector is calculated.
#
# ============================================================

import mpmath as mp

from cell import (
    DEFAULT_DPS,
    compute_L,
)

from connes_cvs.operator import (
    _re_S_and_dS_fused,
)


# ============================================================
# PARAMETERS
# ============================================================

mp.mp.dps = DEFAULT_DPS

c = 13
T = 40

L = compute_L(c)

print("\n" + "=" * 70)
print("CELL_6B — ARCHIMEDEAN QUADRATURE / BREAKPOINT AUDIT")
print("=" * 70)

print("\nParameters:")
print("c =", c)
print("T =", T)
print("dps =", mp.mp.dps)
print("L =", mp.nstr(L, 60))

print(
    "\n2*pi/L =",
    mp.nstr(
        2 * mp.pi / L,
        60
    )
)


# ============================================================
# 1. h_+(r)
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
# 2. S(r,x)
#
# Use the repository's analytic implementation here.
#
# Cell 6A established that this agrees with the independent
# direct y-integral to approximately 1e-80, so using it here
# isolates the OUTER quadrature question.
# ============================================================

def S_repo(r, x):

    r = mp.mpf(r)
    x = mp.mpf(x)

    S_value, dS_value = _re_S_and_dS_fused(
        r,
        x,
        L,
    )

    return mp.mpf(S_value)


# ============================================================
# 3. ARCHIMEDEAN INTEGRAND
# ============================================================

def arch_integrand(r, x):

    return (
        h_plus(r)
        * S_repo(r, x)
    )


# ============================================================
# 4. FIXED-STEP QUADRATURE
#
# This reproduces the strategy used in Cell 6:
#
#     0, 5, 10, 15, ...
#
# ============================================================

def integrate_fixed_step(x, T_local, step=5):

    x = mp.mpf(x)
    T_local = mp.mpf(T_local)
    step = mp.mpf(step)

    points = [mp.mpf(0)]

    r = step

    while r < T_local:

        points.append(r)

        r += step

    points.append(T_local)

    total = mp.mpf(0)

    for a, b in zip(
        points[:-1],
        points[1:]
    ):

        total += mp.quad(
            lambda r: arch_integrand(r, x),
            [a, b]
        )

    return total


# ============================================================
# 5. BREAKPOINT QUADRATURE
#
# The natural breakpoint is
#
#     alpha_x = 2*pi*x/L.
#
# On the positive half-line we therefore use
#
#     0 -> alpha_x -> T
#
# when alpha_x lies inside the interval.
#
# ============================================================

def integrate_breakpoint_positive(
    x,
    T_local,
):

    x = mp.mpf(x)
    T_local = mp.mpf(T_local)

    alpha_x = (
        2 * mp.pi * abs(x) / L
    )

    points = [mp.mpf(0)]

    if (
        alpha_x > 0
        and alpha_x < T_local
    ):
        points.append(alpha_x)

    points.append(T_local)

    total = mp.mpf(0)

    for a, b in zip(
        points[:-1],
        points[1:]
    ):

        total += mp.quad(
            lambda r: arch_integrand(r, x),
            [a, b]
        )

    return total


# ============================================================
# 6. REPOSITORY-STYLE FULL INTERVAL
#
# Because S and h_+ are even in r, this should agree with
#
#     2 * integral_0^T
#
# but we calculate it directly over [-T,T].
#
# We explicitly include the natural breakpoints
#
#     -alpha_x, 0, +alpha_x.
#
# ============================================================

def integrate_repository_style(
    x,
    T_local,
):

    x = mp.mpf(x)
    T_local = mp.mpf(T_local)

    alpha_x = (
        2 * mp.pi * abs(x) / L
    )

    points = [-T_local]

    if (
        alpha_x > 0
        and alpha_x < T_local
    ):
        points.append(-alpha_x)

    points.append(mp.mpf(0))

    if (
        alpha_x > 0
        and alpha_x < T_local
    ):
        points.append(alpha_x)

    points.append(T_local)

    # Ensure ordering.
    points = sorted(points)

    total = mp.mpf(0)

    for a, b in zip(
        points[:-1],
        points[1:]
    ):

        total += mp.quad(
            lambda r: arch_integrand(r, x),
            [a, b]
        )

    return total


# ============================================================
# 7. MAIN COMPARISON
#
# Test x values chosen specifically because their breakpoints
# occur at different locations relative to the 5-unit grid.
# ============================================================

x_values = [
    1,
    2,
    4,
    8,
]

print("\n" + "-" * 70)
print("QUADRATURE COMPARISON")
print("-" * 70)

results = {}


for x in x_values:

    alpha_x = (
        2 * mp.pi * x / L
    )

    print("\n")
    print("x =", x)

    print(
        "alpha_x =",
        mp.nstr(alpha_x, 50)
    )

    print(
        "alpha_x / 5 =",
        mp.nstr(
            alpha_x / 5,
            30
        )
    )

    # --------------------------------------------------------
    # A. Fixed-step positive integral
    # --------------------------------------------------------

    fixed_positive = integrate_fixed_step(
        x,
        T,
        step=5,
    )

    # --------------------------------------------------------
    # B. Breakpoint positive integral
    # --------------------------------------------------------

    breakpoint_positive = (
        integrate_breakpoint_positive(
            x,
            T,
        )
    )

    # --------------------------------------------------------
    # C. Full repository-style integral
    # --------------------------------------------------------

    repository_full = (
        integrate_repository_style(
            x,
            T,
        )
    )

    # --------------------------------------------------------
    # Convert to the actual psi_arch normalization:
    #
    #     1/pi^2 integral_0^T
    #
    # for the positive-half versions, and
    #
    #     1/(2*pi^2) integral_-T^T
    #
    # for the full version.
    # --------------------------------------------------------

    A_fixed = (
        fixed_positive
        / (mp.pi ** 2)
    )

    A_breakpoint = (
        breakpoint_positive
        / (mp.pi ** 2)
    )

    A_repository = (
        repository_full
        / (2 * mp.pi ** 2)
    )

    fixed_vs_breakpoint = (
        A_fixed
        - A_breakpoint
    )

    breakpoint_vs_repository = (
        A_breakpoint
        - A_repository
    )

    fixed_vs_repository = (
        A_fixed
        - A_repository
    )

    results[x] = (
        A_fixed,
        A_breakpoint,
        A_repository,
    )

    print("\n  Fixed-step Cell-6 style:")
    print(
        "    ",
        mp.nstr(
            A_fixed,
            60
        )
    )

    print("\n  Breakpoint positive:")
    print(
        "    ",
        mp.nstr(
            A_breakpoint,
            60
        )
    )

    print("\n  Repository-style full:")
    print(
        "    ",
        mp.nstr(
            A_repository,
            60
        )
    )

    print("\n  Fixed - breakpoint:")
    print(
        "    ",
        mp.nstr(
            fixed_vs_breakpoint,
            50
        )
    )

    print("\n  Breakpoint - repository:")
    print(
        "    ",
        mp.nstr(
            breakpoint_vs_repository,
            50
        )
    )

    print("\n  Fixed - repository:")
    print(
        "    ",
        mp.nstr(
            fixed_vs_repository,
            50
        )
    )


# ============================================================
# 8. PARITY CHECK
#
# This confirms the full integral really is twice the positive
# integral, independently of the breakpoint splitting.
# ============================================================

print("\n" + "-" * 70)
print("PARITY / FULL-INTEGRAL CHECK")
print("-" * 70)

max_full_positive_error = mp.mpf(0)

for x in x_values:

    positive = integrate_breakpoint_positive(
        x,
        T,
    )

    full = integrate_repository_style(
        x,
        T,
    )

    error = abs(
        full - 2 * positive
    )

    if error > max_full_positive_error:
        max_full_positive_error = error

    print("\nx =", x)

    print(
        "  full integral =",
        mp.nstr(
            full,
            50
        )
    )

    print(
        "  2 * positive  =",
        mp.nstr(
            2 * positive,
            50
        )
    )

    print(
        "  absolute error =",
        mp.nstr(
            error,
            30
        )
    )


# ============================================================
# 9. FINAL SUMMARY
# ============================================================

print("\n" + "=" * 70)
print("CELL_6B SUMMARY")
print("=" * 70)

print(
    "\nThe important quantities are:"
)

print(
    "\n  Fixed-step Cell-6 result"
    "\n  Breakpoint-split result"
    "\n  Repository-style full result"
)

print(
    "\nIf all three agree to high precision, the outer"
    "\nquadrature splitting is NOT responsible for the"
    "\nArchimedean discrepancy."
)

print(
    "\nIf Fixed differs substantially from Breakpoint,"
    "\nwhile Breakpoint agrees with Repository, then Cell 6's"
    "\n5-unit subdivision was inadequate."
)

print(
    "\nMaximum |full - 2*positive| ="
)

print(
    mp.nstr(
        max_full_positive_error,
        50
    )
)

print("\n" + "=" * 70)
print("END CELL_6B")
print("=" * 70)
