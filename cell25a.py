# ============================================================
# CELL 25a
#
# Reproduce the HISTORICAL Cell-5 Archimedean test vector
# using the CURRENT canonical definitions in cell.py.
#
# This is deliberately NOT the Galerkin ground state.
#
# Historical cell5_corrected.py overwrites v_star with the
# maximiser of |F'(gamma_1)| subject to
#
#     P(v) = 0
#     F(gamma_1) = 0
#     F(gamma_2) = 0
#     ||v|| = 1
#
# This cell reconstructs that exact vector and then evaluates
# the new analytic Archimedean integral K_fourier against it.
#
# ============================================================

import time
import mpmath as mp

from cell import (
    compute_L,
    F_basis,
    Fprime_basis,
    pole_row,
    archimedean_integral,
)


# ============================================================
# NUMERICAL PARAMETERS
# ============================================================

c = 13
N = 8
WORKING_DPS = 80

m_test = 2

# Same broad T scan used for the Cell-25 investigation.
T_VALUES = [
    20,
    40,
    60,
    80,
    120,
    200,
    300,
    400,
    600,
    800,
    1200,
    1600,
    2000,
]

mp.mp.dps = WORKING_DPS

L = compute_L(c)


# ============================================================
# HISTORICAL REFERENCE VALUES
#
# These are the values printed by cell5_corrected.out.
# They are retained here only as forensic cross-checks.
# ============================================================

HISTORICAL_D_M2 = mp.mpf(
    "0.848726327096095877763669591987044657653672658277128165086134"
)

HISTORICAL_ARCH = {
    20: mp.mpf(
        "0.7125875193025273831377680272"
    ),
    40: mp.mpf(
        "0.7421874519938487471064025458"
    ),
}


# ============================================================
# BUILD THE HISTORICAL LINEAR-FUNCTIONAL ROWS
# ============================================================

print("============================================================")
print("CELL 25a — HISTORICAL CELL-5 MAXIMISER")
print("============================================================")
print()
print(f"c = {c}")
print(f"N = {N}")
print(f"m = {m_test}")
print(f"dps = {WORKING_DPS}")
print(f"L = {mp.nstr(L, 50)}")


# ------------------------------------------------------------
# Pole functional
#
# P(v) = pole_row . v
# ------------------------------------------------------------

pole = pole_row(N, L)


# ------------------------------------------------------------
# Zero constraints
#
# F(gamma_j) = zero_row(j) . v
# ------------------------------------------------------------

def zero_row(j):
    gamma_j = mp.im(mp.zetazero(j))

    return mp.matrix([
        F_basis(k, gamma_j, L)
        for k in range(N + 1)
    ])


gamma = [
    mp.im(mp.zetazero(j))
    for j in range(1, m_test + 1)
]


# ------------------------------------------------------------
# Derivative functional
#
# F'(gamma_1) = derivative_row . v
# ------------------------------------------------------------

derivative_row = mp.matrix([
    Fprime_basis(k, gamma[0], L)
    for k in range(N + 1)
])


# ============================================================
# CONSTRAINT MATRIX
#
# Exactly the construction used by historical constrained_D(m):
#
#     rows = [pole_row]
#     for j = 1,...,m:
#         rows.append(zero_row(j))
#
# ============================================================

rows = [pole]

for j in range(1, m_test + 1):
    rows.append(zero_row(j))

C = mp.matrix(
    len(rows),
    N + 1,
)

for i in range(len(rows)):
    for k in range(N + 1):
        C[i, k] = rows[i][k]


# ============================================================
# ORTHOGONAL PROJECTION
#
# Project derivative_row onto ker(C).
#
#     G = C C^T
#     Cd = C d
#     y = G^{-1} Cd
#     d_perp = d - C^T y
#
# This is exactly the historical Cell-5 construction.
# ============================================================

G = C * C.T

Cd = C * derivative_row

y = mp.lu_solve(
    G,
    Cd,
)

d_perp = (
    derivative_row
    - C.T * y
)

D_m2 = mp.sqrt(
    mp.fdot(
        d_perp,
        d_perp,
    )
)


# ============================================================
# NORMALISED HISTORICAL MAXIMISER
# ============================================================

v_maximiser_m2 = (
    d_perp / D_m2
)


# ============================================================
# FORENSIC CHECKS
# ============================================================

norm_v = mp.sqrt(
    mp.fdot(
        v_maximiser_m2,
        v_maximiser_m2,
    )
)

pole_residual = mp.fdot(
    pole,
    v_maximiser_m2,
)

zero_residuals = [
    mp.fdot(
        zero_row(j),
        v_maximiser_m2,
    )
    for j in range(1, m_test + 1)
]

derivative_value = mp.fdot(
    derivative_row,
    v_maximiser_m2,
)


print()
print("------------------------------------------------------------")
print("RECONSTRUCTED HISTORICAL VECTOR")
print("------------------------------------------------------------")

print()
print("D_m2 =")
print(mp.nstr(D_m2, 70))

print()
print("Historical Cell-5 D_m2 =")
print(mp.nstr(HISTORICAL_D_M2, 70))

print()
print("D_m2 difference =")
print(mp.nstr(
    D_m2 - HISTORICAL_D_M2,
    50,
))

print()
print("||v_maximiser_m2|| =")
print(mp.nstr(norm_v, 50))

print()
print("Constraint residuals:")
print("P(v) =")
print(mp.nstr(pole_residual, 50))

for j, residual in enumerate(zero_residuals, start=1):
    print(f"F(gamma{j}) =")
    print(mp.nstr(residual, 50))

print()
print("F'(gamma1) =")
print(mp.nstr(derivative_value, 60))


# ============================================================
# HARD FORENSIC GUARD
#
# If this fails, DO NOT proceed to the expensive Archimedean
# integration. It means we have failed to reproduce the
# historical Cell-5 vector.
# ============================================================

D_TOL = mp.mpf("1e-55")

if abs(D_m2 - HISTORICAL_D_M2) > D_TOL:
    raise RuntimeError(
        "Historical Cell-5 maximiser was not reproduced: "
        f"|D_m2 - D_reference| = "
        f"{mp.nstr(abs(D_m2 - HISTORICAL_D_M2), 20)}"
    )


# ============================================================
# ARCHIMEDEAN INTEGRAL SCAN
#
# This is now an apples-to-apples comparison with historical
# Cell 5:
#
#     same c
#     same N
#     same m = 2
#     same canonical maximising vector
#
# The only intended change is the Archimedean evaluation:
#
#     historical: explicit finite-T calculation
#     current:    analytic K_fourier representation
# ============================================================

print()
print("============================================================")
print("ARCHIMEDEAN INTEGRAL — HISTORICAL m=2 VECTOR")
print("============================================================")

print()
print(
    " T"
    "    A_arch(T)"
    "                                      "
    "historical difference"
)
print("------------------------------------------------------------")

arch_results = {}

total_start = time.perf_counter()

for T in T_VALUES:

    start = time.perf_counter()

    A = archimedean_integral(
        T,
        v_maximiser_m2,
        L,
    )

    elapsed = (
        time.perf_counter()
        - start
    )

    arch_results[T] = A

    if T in HISTORICAL_ARCH:
        difference = (
            A - HISTORICAL_ARCH[T]
        )

        historical_text = mp.nstr(
            difference,
            40,
        )
    else:
        historical_text = ""

    print(
        f"{T:4d}"
        f"  {mp.nstr(A, 45):<50}"
        f"  {historical_text}"
        f"    ({elapsed:.3f} s)"
    )


total_elapsed = (
    time.perf_counter()
    - total_start
)


# ============================================================
# HISTORICAL CROSS-CHECK
# ============================================================

print()
print("============================================================")
print("HISTORICAL CELL-5 CROSS-CHECK")
print("============================================================")

for T in (20, 40):

    A_new = arch_results[T]
    A_old = HISTORICAL_ARCH[T]

    difference = A_new - A_old

    print()
    print(f"T = {T}")
    print("current analytic A_arch(T) =")
    print(mp.nstr(A_new, 60))

    print("historical Cell-5 A_arch(T) =")
    print(mp.nstr(A_old, 60))

    print("difference =")
    print(mp.nstr(difference, 50))


print()
print("Total Archimedean scan time =")
print(f"{total_elapsed:.3f} s")
