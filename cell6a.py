# ============================================================
# CELL_6A — POINTWISE ARCHIMEDEAN SOURCE CHECK
#
# Compare the independent direct y-integral used in Cell 6
# against the repository's analytic _re_S_and_dS_fused()
# implementation.
#
# This is deliberately a small, cheap diagnostic.
#
# It does NOT build the Galerkin matrix.
# It does NOT calculate any eigenvectors.
# It does NOT perform the expensive outer Archimedean integral.
#
# We compare:
#
#   Re S(r,x,L)
#
# and
#
#   Re dS/dx(r,x,L)
#
# at selected (r,x) points.
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
N = 8

L = compute_L(c)

alpha_factor = 2 * mp.pi / L


print("\n" + "=" * 70)
print("CELL_6A — POINTWISE ARCHIMEDEAN SOURCE CHECK")
print("=" * 70)

print("\nParameters:")
print("c =", c)
print("N =", N)
print("dps =", mp.mp.dps)
print("L =", mp.nstr(L, 60))
print(
    "2*pi/L =",
    mp.nstr(alpha_factor, 60)
)


# ============================================================
# 1. INDEPENDENT DIRECT DEFINITIONS
#
# These reproduce the definitions used in Cell 6.
# ============================================================

def S_direct(r, x):

    r = mp.mpf(r)
    x = mp.mpf(x)

    integrand = lambda y: (
        mp.sin(
            2
            * mp.pi
            * x
            * (1 - y / L)
        )
        * mp.cos(r * y)
    )

    return mp.quad(
        integrand,
        [0, L]
    )


def dS_direct(r, x):

    r = mp.mpf(r)
    x = mp.mpf(x)

    integrand = lambda y: (
        2
        * mp.pi
        * (1 - y / L)
        * mp.cos(
            2
            * mp.pi
            * x
            * (1 - y / L)
        )
        * mp.cos(r * y)
    )

    return mp.quad(
        integrand,
        [0, L]
    )


# ============================================================
# 2. REPOSITORY WRAPPER
#
# _re_S_and_dS_fused() returns the repository's analytic
# real S and derivative.
#
# We deliberately use the repository implementation here,
# rather than reproducing its algebra ourselves.
# ============================================================

def repository_S_and_dS(r, x):

    r = mp.mpf(r)
    x = mp.mpf(x)

    result = _re_S_and_dS_fused(
        r,
        x,
        L,
    )

    return (
        mp.mpf(result[0]),
        mp.mpf(result[1]),
    )


# ============================================================
# 3. TEST POINTS
#
# Include:
#
#   r = 0
#   ordinary interior points
#   the natural breakpoint alpha_x
#   points around the breakpoint
#   a relatively large r
#
# The breakpoint is
#
#   alpha_x = 2*pi*x/L.
#
# ============================================================

x_values = [
    0,
    1,
    2,
    3,
    8,
]

r_base_values = [
    mp.mpf("0"),
    mp.mpf("1"),
    mp.mpf("2"),
    mp.mpf("10"),
    mp.mpf("20"),
    mp.mpf("40"),
]


# ============================================================
# 4. NUMERICAL COMPARISON
# ============================================================

max_S_error = mp.mpf(0)
max_dS_error = mp.mpf(0)

max_S_location = None
max_dS_location = None


print("\n" + "-" * 70)
print("POINTWISE COMPARISON")
print("-" * 70)


for x in x_values:

    alpha_x = alpha_factor * x

    r_values = list(r_base_values)

    if x != 0:

        # Exact breakpoint
        r_values.append(alpha_x)

        # Points immediately around the breakpoint
        r_values.append(
            alpha_x - mp.mpf("0.1")
        )

        r_values.append(
            alpha_x + mp.mpf("0.1")
        )

        # Halfway to the breakpoint
        r_values.append(
            alpha_x / 2
        )

    # Remove duplicates while preserving numerical values.
    unique_r_values = []

    for r in r_values:

        if not any(
            r == existing
            for existing in unique_r_values
        ):
            unique_r_values.append(r)

    print("\n")
    print("x =", x)

    if x != 0:
        print(
            "alpha_x =",
            mp.nstr(alpha_x, 50)
        )

    for r in unique_r_values:

        direct_S = mp.re(
            S_direct(r, x)
        )

        direct_dS = mp.re(
            dS_direct(r, x)
        )

        repo_S, repo_dS = (
            repository_S_and_dS(r, x)
        )

        S_error = abs(
            direct_S - repo_S
        )

        dS_error = abs(
            direct_dS - repo_dS
        )

        if S_error > max_S_error:

            max_S_error = S_error
            max_S_location = (x, r)

        if dS_error > max_dS_error:

            max_dS_error = dS_error
            max_dS_location = (x, r)

        print("\n  r =", mp.nstr(r, 30))

        print(
            "    direct S  =",
            mp.nstr(direct_S, 50)
        )

        print(
            "    repo S    =",
            mp.nstr(repo_S, 50)
        )

        print(
            "    |error|   =",
            mp.nstr(S_error, 20)
        )

        print(
            "    direct dS =",
            mp.nstr(direct_dS, 50)
        )

        print(
            "    repo dS   =",
            mp.nstr(repo_dS, 50)
        )

        print(
            "    |error|   =",
            mp.nstr(dS_error, 20)
        )


# ============================================================
# 5. SUMMARY
# ============================================================

print("\n" + "-" * 70)
print("SUMMARY")
print("-" * 70)

print("\nMaximum |S_direct - S_repository| =")
print(
    mp.nstr(
        max_S_error,
        50
    )
)

print("\nLocation of maximum S error =")

if max_S_location is not None:

    print(
        "x =",
        max_S_location[0],
        " r =",
        mp.nstr(
            max_S_location[1],
            50
        )
    )


print("\nMaximum |dS_direct - dS_repository| =")
print(
    mp.nstr(
        max_dS_error,
        50
    )
)

print("\nLocation of maximum dS error =")

if max_dS_location is not None:

    print(
        "x =",
        max_dS_location[0],
        " r =",
        mp.nstr(
            max_dS_location[1],
            50
        )
    )


# ============================================================
# 6. PARITY CHECK
#
# Since the outer Archimedean integral is reduced using
# evenness in Cell 6, explicitly verify that the source has
# the expected parity.
#
# We compare S(r,x) with S(-r,x).
# ============================================================

print("\n" + "-" * 70)
print("PARITY CHECK")
print("-" * 70)

max_parity_error = mp.mpf(0)
max_parity_location = None

for x in x_values:

    for r in [
        mp.mpf("0.7"),
        mp.mpf("3.1"),
        mp.mpf("10.0"),
        mp.mpf("25.0"),
    ]:

        S_plus = mp.re(
            S_direct(r, x)
        )

        S_minus = mp.re(
            S_direct(-r, x)
        )

        error = abs(
            S_plus - S_minus
        )

        if error > max_parity_error:

            max_parity_error = error
            max_parity_location = (x, r)

print("\nMaximum |S(r,x) - S(-r,x)| =")
print(
    mp.nstr(
        max_parity_error,
        50
    )
)

if max_parity_location is not None:

    print("\nLocation:")
    print(
        "x =",
        max_parity_location[0],
        " r =",
        mp.nstr(
            max_parity_location[1],
            50
        )
    )


# ============================================================
# 7. END
# ============================================================

print("\n" + "=" * 70)
print("END CELL_6A")
print("=" * 70)
