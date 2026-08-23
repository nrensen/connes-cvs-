# ============================================================
# CELL 8 — EXPLICIT-FORMULA / CVs ARCHIMEDEAN DICTIONARY
#
# Purpose:
#
# Establish directly that the Archimedean source used by the
# repository is the Proposition-4.1 source
#
#   psi_arch(x)
#       = 1/(2*pi^2) integral h_+(r) S(r,x,L) dr
#
# where
#
#   S(r,x,L)
#       = integral_0^L
#           sin(2*pi*x*(1-y/L)) cos(r*y) dy.
#
# We independently evaluate S in two ways:
#
#   A. direct y quadrature
#   B. closed-form trigonometric expression
#
# Then compare both against the repository's kernel.
#
# Finally we integrate the resulting source and compare it
# against operator.py::psi_arch.
#
# This cell deliberately does NOT construct a Galerkin matrix.
# ============================================================

import mpmath as mp

from cell import compute_L
from connes_cvs.operator import (
    psi_arch,
    h_plus,
    _re_S_and_dS_fused,
)


# ============================================================
# PARAMETERS
# ============================================================

mp.mp.dps = 70

c = 13
T = 20

L = compute_L(c)

print("=" * 70)
print("CELL 8 — EXPLICIT-FORMULA / CvS ARCHIMEDEAN DICTIONARY")
print("=" * 70)
print()
print("Parameters:")
print(f"c = {c}")
print(f"T = {T}")
print(f"dps = {mp.mp.dps}")
print(f"L = {mp.nstr(L, 60)}")
print(f"2*pi/L = {mp.nstr(2*mp.pi/L, 60)}")
print()


# ============================================================
# 1. DIRECT DEFINITION OF S(r,x,L)
#
# S(r,x,L) =
#
#   integral_0^L
#       sin(2*pi*x*(1-y/L))
#       cos(r*y) dy
#
# This is the source kernel appearing in the explicit
# Archimedean formula.
# ============================================================

def S_direct(r, x):
    r = mp.mpf(r)
    x = mp.mpf(x)

    f = lambda y: (
        mp.sin(2 * mp.pi * x * (1 - y / L))
        * mp.cos(r * y)
    )

    return mp.quad(f, [0, L])


# ============================================================
# 2. CLOSED FORM FOR S
#
# Write
#
#   a = 2*pi*x/L
#
# Then
#
#   sin(2*pi*x - a*y) cos(r*y)
#
# is reduced using product-to-sum identities.
#
# We evaluate the resulting elementary expression directly.
# ============================================================

def _J(k):
    """
    Integral_0^L cos(k*y) dy.
    """
    k = mp.mpf(k)

    if k == 0:
        return L

    return mp.sin(k * L) / k


def S_closed(r, x):
    r = mp.mpf(r)
    x = mp.mpf(x)

    a = 2 * mp.pi * x / L
    theta = 2 * mp.pi * x

    # sin(theta - a*y) cos(r*y)
    #
    # = sin(theta) cos(a*y) cos(r*y)
    #   - cos(theta) sin(a*y) cos(r*y)
    #
    # First:
    #
    # cos(a y) cos(r y)
    #   = 1/2[cos((a-r)y) + cos((a+r)y)]
    #
    # Second:
    #
    # sin(a y) cos(r y)
    #   = 1/2[sin((a+r)y) + sin((a-r)y)]

    I1 = (
        _J(a - r)
        + _J(a + r)
    ) / 2

    def K(k):
        k = mp.mpf(k)

        if k == 0:
            return mp.mpf(0)

        return (1 - mp.cos(k * L)) / k

    I2 = (
        K(a + r)
        + K(a - r)
    ) / 2

    return mp.sin(theta) * I1 - mp.cos(theta) * I2


# ============================================================
# 3. REPOSITORY KERNEL
#
# operator.py computes Re(S_hat_x). For the real kernel used
# in the Archimedean integral this should coincide with S.
# ============================================================

def S_repository(r, x):
    re_S, _ = _re_S_and_dS_fused(
        mp.mpf(r),
        mp.mpf(x),
        L,
    )
    return re_S


# ============================================================
# 4. POINTWISE KERNEL COMPARISON
# ============================================================

print("-" * 70)
print("1. POINTWISE S(r,x,L) COMPARISON")
print("-" * 70)
print()

r_values = [
    mp.mpf("0"),
    mp.mpf("0.3"),
    mp.mpf("1"),
    mp.mpf("2"),
    mp.mpf("2.4496332798546520107426711685"),
    mp.mpf("3"),
    mp.mpf("5"),
]

x_values = [0, 1, 2, 4]

max_direct_closed = mp.mpf("0")
max_direct_repo = mp.mpf("0")
max_closed_repo = mp.mpf("0")

for x in x_values:
    for r in r_values:

        direct = S_direct(r, x)
        closed = S_closed(r, x)
        repo = S_repository(r, x)

        e1 = abs(direct - closed)
        e2 = abs(direct - repo)
        e3 = abs(closed - repo)

        max_direct_closed = max(max_direct_closed, e1)
        max_direct_repo = max(max_direct_repo, e2)
        max_closed_repo = max(max_closed_repo, e3)

        print(
            f"x={x:2d}, r={mp.nstr(r, 15)}"
        )
        print(
            f"  direct      = {mp.nstr(direct, 45)}"
        )
        print(
            f"  closed      = {mp.nstr(closed, 45)}"
        )
        print(
            f"  repository  = {mp.nstr(repo, 45)}"
        )
        print(
            f"  |direct-closed| = {mp.nstr(e1, 20)}"
        )
        print(
            f"  |direct-repo|   = {mp.nstr(e2, 20)}"
        )
        print(
            f"  |closed-repo|   = {mp.nstr(e3, 20)}"
        )
        print()


print("Maximum |direct - closed| =",
      mp.nstr(max_direct_closed, 30))
print()

print("Maximum |direct - repository| =",
      mp.nstr(max_direct_repo, 30))
print()

print("Maximum |closed - repository| =",
      mp.nstr(max_closed_repo, 30))
print()


# ============================================================
# 5. EVENNESS IN r
#
# The Archimedean multiplier h_+(r) is even, and the source
# kernel should also be even in r.
# ============================================================

print("-" * 70)
print("2. PARITY CHECK")
print("-" * 70)
print()

max_parity = mp.mpf("0")

for x in x_values:
    for r in [
        mp.mpf("0.3"),
        mp.mpf("1"),
        mp.mpf("2"),
        mp.mpf("5"),
    ]:

        e = abs(
            S_repository(r, x)
            - S_repository(-r, x)
        )

        max_parity = max(max_parity, e)

        print(
            f"x={x:2d}, r={mp.nstr(r, 12)}  "
            f"|S(r)-S(-r)| = {mp.nstr(e, 20)}"
        )

print()
print(
    "Maximum parity error =",
    mp.nstr(max_parity, 30),
)
print()


# ============================================================
# 6. ARCHIMEDEAN SOURCE FROM THE EXPLICIT FORMULA
#
# Independently construct
#
#   psi_explicit(x)
#       = 1/(2*pi^2)
#         integral_{-T}^T
#           h_+(r) S(r,x,L) dr.
#
# We use the CLOSED FORM S, not the repository implementation.
#
# This is the key dictionary test.
# ============================================================

def psi_explicit(x):
    x = mp.mpf(x)

    if x == 0:
        return mp.mpf("0")

    def integrand(r):
        return (
            h_plus(r, mp.mp.dps)
            * S_closed(r, x)
        )

    alpha = 2 * mp.pi * x / L

    points = [
        p for p in [-T, -alpha, 0, alpha, T]
        if -T <= p <= T
    ]

    points = sorted(set(points))

    total = mp.mpf("0")

    for a, b in zip(points[:-1], points[1:]):
        if a != b:
            total += mp.quad(integrand, [a, b])

    return total / (2 * mp.pi**2)


# ============================================================
# 7. COMPARE EXPLICIT SOURCE WITH REPOSITORY psi_arch
# ============================================================

print("-" * 70)
print("3. ARCHIMEDEAN SOURCE COMPARISON")
print("-" * 70)
print()

max_psi_error = mp.mpf("0")

for x in [0, 1, 2, 4, 8]:

    explicit = psi_explicit(x)

    repository = psi_arch(
        mp.mpf(x),
        L,
        T,
        mp.mp.dps,
    )

    err = abs(explicit - repository)

    max_psi_error = max(max_psi_error, err)

    print(f"x = {x}")
    print(
        f"  explicit  = {mp.nstr(explicit, 55)}"
    )
    print(
        f"  repository = {mp.nstr(repository, 55)}"
    )
    print(
        f"  |error|    = {mp.nstr(err, 30)}"
    )
    print()


print(
    "Maximum |psi_explicit - psi_repository| =",
    mp.nstr(max_psi_error, 30),
)
print()


# ============================================================
# 8. NORMALISATION / SIGN DIAGNOSTIC
#
# If the two source constructions agree directly, this section
# should show ratio 1.
#
# If not, it will tell us immediately whether the remaining
# discrepancy is:
#
#   * a constant factor
#   * a sign
#   * neither
# ============================================================

print("-" * 70)
print("4. NORMALISATION DIAGNOSTIC")
print("-" * 70)
print()

for x in [1, 2, 4, 8]:

    explicit = psi_explicit(x)
    repository = psi_arch(
        mp.mpf(x),
        L,
        T,
        mp.mp.dps,
    )

    print(f"x = {x}")

    if repository != 0:
        print(
            "  explicit / repository =",
            mp.nstr(explicit / repository, 40),
        )

    print(
        "  explicit - repository =",
        mp.nstr(explicit - repository, 30),
    )

    print()


# ============================================================
# 9. h_+ CHECK
#
# Reconfirm the factor-of-two identity established in Cell 7:
#
# h_+(r)
#   =
# 2 Re[d/ds log(pi^(-s/2) Gamma(s/2))]
# at s=1/2+ir.
# ============================================================

def h_completed(r):
    r = mp.mpf(r)

    s = mp.mpf("0.5") + 1j * r

    return mp.re(
        -mp.log(mp.pi) / 2
        + mp.digamma(s / 2) / 2
    )


print("-" * 70)
print("5. COMPLETED-ZETA ARCHIMEDEAN FACTOR")
print("-" * 70)
print()

max_h_error = mp.mpf("0")

for r in [
    mp.mpf("0"),
    mp.mpf("0.5"),
    mp.mpf("1"),
    mp.mpf("2"),
    mp.mpf("5"),
]:

    lhs = h_plus(r, mp.mp.dps)
    rhs = 2 * h_completed(r)

    err = abs(lhs - rhs)
    max_h_error = max(max_h_error, err)

    print(f"r = {mp.nstr(r, 15)}")
    print(
        "  h_+      =",
        mp.nstr(lhs, 50),
    )
    print(
        "  2*h_comp =",
        mp.nstr(rhs, 50),
    )
    print(
        "  error    =",
        mp.nstr(err, 20),
    )
    print()


# ============================================================
# SUMMARY
# ============================================================

print("=" * 70)
print("CELL 8 SUMMARY")
print("=" * 70)
print()

print(
    "Maximum |S_direct - S_closed| =",
    mp.nstr(max_direct_closed, 30),
)
print()

print(
    "Maximum |S_direct - S_repository| =",
    mp.nstr(max_direct_repo, 30),
)
print()

print(
    "Maximum |S_closed - S_repository| =",
    mp.nstr(max_closed_repo, 30),
)
print()

print(
    "Maximum parity error =",
    mp.nstr(max_parity, 30),
)
print()

print(
    "Maximum |psi_explicit - psi_repository| =",
    mp.nstr(max_psi_error, 30),
)
print()

print(
    "Maximum h_+ completed-factor error =",
    mp.nstr(max_h_error, 30),
)
print()

print("=" * 70)
print("END CELL 8")
print("=" * 70)
