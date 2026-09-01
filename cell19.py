# ============================================================
# cell19.py — LINEAR VS QUADRATIC HOMOGENEITY AUDIT
#
# Purpose:
#
#   Demonstrate directly that the historical Cell-5
#   Archimedean quantity is linear in v, whereas the
#   required Archimedean functional is quadratic in v.
#
# Therefore:
#
#     A_linear(a v) = a A_linear(v)
#
# whereas:
#
#     A_arch(a v) = a^2 A_arch(v)
#
# This is deliberately independent of the ground state.
#
# Historical h_plus and K_canonical definitions are reproduced
# locally from Cell 5 so that this audit has no dependency on
# their availability through cell.py.
# ============================================================

import mpmath as mp

from cell import sum_v_F


# ============================================================
# PARAMETERS
# ============================================================

mp.mp.dps = 50

c = mp.mpf("13")
N = 8
T = mp.mpf("10")
L = mp.log(c)


# ============================================================
# HISTORICAL CELL-5 h_plus
# ============================================================

def h_plus_cell5(r):
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
# HISTORICAL CELL-5 T_canonical
# ============================================================

def T_canonical(v, t):
    """
    Trigonometric polynomial T_v(t) in the full symmetric
    coefficient representation.
    """

    total = mp.mpc(0)

    # m = 0
    total += v[0]

    # positive / negative pairs
    for k in range(1, N + 1):

        uk = v[k] / mp.sqrt(2)

        total += (
            uk * mp.exp(2j * mp.pi * k * t)
            + uk * mp.exp(-2j * mp.pi * k * t)
        )

    return total


# ============================================================
# HISTORICAL CELL-5 K_canonical
# ============================================================

def K_canonical(v, omega):
    """
    Volterra sine-chord kernel

        K_v(omega)
          = 2 int_0^omega
              T_v(t) T_v(omega-t) dt.
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
        [0, omega],
    )


# ============================================================
# TEST VECTOR
# ============================================================

v = mp.matrix([
    mp.mpf("1.0"),
    mp.mpf("0.3"),
    mp.mpf("-0.2"),
    mp.mpf("0.15"),
    mp.mpf("0.1"),
    mp.mpf("-0.05"),
    mp.mpf("0.07"),
    mp.mpf("-0.04"),
    mp.mpf("0.02"),
])

v2 = 2 * v


# ============================================================
# HISTORICAL CELL-5 LINEAR FUNCTIONAL
# ============================================================

def A_linear(v):

    def integrand(r):

        return (
            h_plus_cell5(r)
            * mp.re(
                sum_v_F(
                    v,
                    r,
                    L
                )
            )
        )

    return mp.quad(
        integrand,
        [0, T],
    ) / mp.pi


# ============================================================
# REQUIRED QUADRATIC ARCHIMEDEAN FUNCTIONAL
# ============================================================

def K_inner(v, r):

    return mp.quad(
        lambda y:
            K_canonical(
                v,
                1 - y / L,
            )
            * mp.cos(r * y),
        [0, L],
    )


def A_quadratic(v):

    return mp.quad(
        lambda r:
            h_plus_cell5(r)
            * K_inner(v, r),
        [0, T],
    ) / mp.pi


# ============================================================
# COMPUTE
# ============================================================

print()
print("=" * 72)
print("CELL 19 — LINEAR VS QUADRATIC HOMOGENEITY AUDIT")
print("=" * 72)

print()
print("Parameters:")
print("  c   =", c)
print("  N   =", N)
print("  T   =", T)
print("  dps =", mp.mp.dps)
print("  L   =", mp.nstr(L, 40))


print()
print("-" * 72)
print("1. HISTORICAL CELL-5 LINEAR QUANTITY")
print("-" * 72)

A1 = A_linear(v)
A2 = A_linear(v2)

print()
print("A_linear(v) =")
print(mp.nstr(A1, 40))

print()
print("A_linear(2v) =")
print(mp.nstr(A2, 40))

print()
print("A_linear(2v) / A_linear(v) =")
print(mp.nstr(A2 / A1, 40))

print()
print("Expected = 2")


print()
print("-" * 72)
print("2. REQUIRED QUADRATIC ARCHIMEDEAN QUANTITY")
print("-" * 72)

Q1 = A_quadratic(v)
Q2 = A_quadratic(v2)

print()
print("A_arch(v) =")
print(mp.nstr(Q1, 40))

print()
print("A_arch(2v) =")
print(mp.nstr(Q2, 40))

print()
print("A_arch(2v) / A_arch(v) =")
print(mp.nstr(Q2 / Q1, 40))

print()
print("Expected = 4")


print()
print("-" * 72)
print("3. HOMOGENEITY ERRORS")
print("-" * 72)

linear_error = abs(
    A2 - 2 * A1
)

quadratic_error = abs(
    Q2 - 4 * Q1
)

print()
print("|A_linear(2v) - 2 A_linear(v)| =")
print(mp.nstr(linear_error, 30))

print()
print("|A_arch(2v) - 4 A_arch(v)| =")
print(mp.nstr(quadratic_error, 30))


print()
print("=" * 72)
print("CELL 19 COMPLETE")
print("=" * 72)
