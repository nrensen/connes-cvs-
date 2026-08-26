# ============================================================
# cell_refactor_pre.py
#
# Regression fingerprint for the cell.py mathematical-vocabulary
# refactor — PRE-REFACTOR version.
#
# PURPOSE
# -------
# Exercise the affected mathematical objects and their immediate
# neighbours, producing a deterministic numerical output.
#
# Protocol:
#
#     python cell_refactor_pre.py > cell_refactor_pre.out
#
# After the cell.py refactor, the corresponding post-refactor test
# will be:
#
#     python cell_refactor_post.py > cell_refactor_post.out
#
# The two output files can then be compared directly with a
# standard diff tool.
#
# IMPORTANT
# ---------
# This test does not introduce new mathematics.  It evaluates the
# existing definitions and independently checks the simple
# coefficient-weighted-sum identities for F, Fprime and G.
# ============================================================

import mpmath as mp

from cell import (
    compute_L,
    compute_delta,
    canonical_pairs,
    F_basis,
    sum_v_F,
    Fprime_basis,
    sum_v_Fprime,
    G_basis_complex,
    sum_v_G,
    pole_basis,
    pole_row,
    T_canonical,
    K_canonical,
    ghat,
)


# ============================================================
# TEST PARAMETERS
# ============================================================

mp.mp.dps = 80

C = mp.mpf("13")
L = compute_L(C)

# Small, deliberately nontrivial canonical coefficient vector.
V = [
    mp.mpf("0.73"),
    mp.mpf("-0.41"),
    mp.mpf("0.19"),
    mp.mpf("0.057"),
]

N = len(V) - 1

TAU = mp.mpf("2.375")
Z = mp.mpc("2.375", "0.61")
T = mp.mpf("0.217")
OMEGA = mp.mpf("0.63")

DELTA = compute_delta(L)
XI = DELTA * mp.mpf("0.37")


# ============================================================
# FORMATTING
# ============================================================

def fmt(x):
    """
    Deterministic high-precision representation of an mpmath
    scalar.
    """
    if isinstance(x, mp.mpc):
        return (
            f"({mp.nstr(mp.re(x), 70)}"
            f" + {mp.nstr(mp.im(x), 70)}j)"
        )

    return mp.nstr(x, 70)


def print_vector(values):
    for i, value in enumerate(values):
        print(f"  [{i}] = {fmt(value)}")


# ============================================================
# SEMANTIC CHECK
# ============================================================

def assert_close(lhs, rhs, label):
    """
    Verify a mathematical identity at substantially higher
    precision than the comparison tolerance.
    """
    error = abs(lhs - rhs)
    tolerance = mp.mpf("1e-60")

    if error > tolerance:
        raise AssertionError(
            f"{label} FAILED\n"
            f"  lhs   = {fmt(lhs)}\n"
            f"  rhs   = {fmt(rhs)}\n"
            f"  error = {fmt(error)}"
        )


# ============================================================
# HEADER
# ============================================================

print("=" * 78)
print("cell.py MATHEMATICAL-VOCABULARY REFACTOR REGRESSION")
print("PRE-REFACTOR")
print("=" * 78)

print()
print("TEST PARAMETERS")
print("----------------")
print(f"mp.dps = {mp.mp.dps}")
print(f"C      = {fmt(C)}")
print(f"N      = {N}")
print(f"L      = {fmt(L)}")
print(f"tau    = {fmt(TAU)}")
print(f"z      = {fmt(Z)}")
print(f"t      = {fmt(T)}")
print(f"omega  = {fmt(OMEGA)}")
print(f"delta  = {fmt(DELTA)}")
print(f"xi     = {fmt(XI)}")

print()
print("V")
print("-")
print_vector(V)


# ============================================================
# CANONICAL PAIRS
# ============================================================

print()
print("canonical_pairs(0)")
print("------------------")

for pair in canonical_pairs(0):
    print(
        f"  ({pair[0]}, {fmt(pair[1])})"
    )

print()
print("canonical_pairs(2)")
print("------------------")

for pair in canonical_pairs(2):
    print(
        f"  ({pair[0]}, {fmt(pair[1])})"
    )


# ============================================================
# F BASIS
# ============================================================

F_basis_values = [
    F_basis(k, TAU, L)
    for k in range(N + 1)
]

print()
print("F_basis(k, tau, L)")
print("------------------")
print_vector(F_basis_values)


# ============================================================
# F VECTOR / COEFFICIENT-WEIGHTED SUM
# ============================================================

F_vector_value = sum_v_F(V, TAU, L)

F_explicit = sum(
    V[k] * F_basis(k, TAU, L)
    for k in range(N + 1)
)

print()
print("F_vector(V, tau, L)")
print("-------------------")
print(f"  value = {fmt(F_vector_value)}")

print()
print("Explicit sum_k V[k] * F_basis(k, tau, L)")
print("------------------------------------------")
print(f"  value = {fmt(F_explicit)}")

assert_close(
    F_vector_value,
    F_explicit,
    "F coefficient-weighted-sum identity",
)

print("  CHECK = PASS")


# ============================================================
# FPRIME BASIS
# ============================================================

Fprime_basis_values = [
    Fprime_basis(k, TAU, L)
    for k in range(N + 1)
]

print()
print("Fprime_basis(k, tau, L)")
print("-----------------------")
print_vector(Fprime_basis_values)


# ============================================================
# FPRIME VECTOR / COEFFICIENT-WEIGHTED SUM
# ============================================================

Fprime_vector_value = sum_v_Fprime(V, TAU, L)

Fprime_explicit = sum(
    V[k] * Fprime_basis(k, TAU, L)
    for k in range(N + 1)
)

print()
print("Fprime_vector(V, tau, L)")
print("------------------------")
print(f"  value = {fmt(Fprime_vector_value)}")

print()
print("Explicit sum_k V[k] * Fprime_basis(k, tau, L)")
print("-----------------------------------------------")
print(f"  value = {fmt(Fprime_explicit)}")

assert_close(
    Fprime_vector_value,
    Fprime_explicit,
    "Fprime coefficient-weighted-sum identity",
)

print("  CHECK = PASS")


# ============================================================
# G BASIS
# ============================================================

G_complex_basis_values = [
    G_basis_complex(k, Z, L)
    for k in range(N + 1)
]

print()
print("G_complex_basis(k, z, L)")
print("------------------------")
print_vector(G_complex_basis_values)


# ============================================================
# G VECTOR / COEFFICIENT-WEIGHTED SUM
# ============================================================

G_complex_value = sum_v_G(V, Z, L)

G_explicit = sum(
    V[k] * G_basis_complex(k, Z, L)
    for k in range(N + 1)
)

print()
print("G_complex(V, z, L)")
print("------------------")
print(f"  value = {fmt(G_complex_value)}")

print()
print("Explicit sum_k V[k] * G_complex_basis(k, z, L)")
print("------------------------------------------------")
print(f"  value = {fmt(G_explicit)}")

assert_close(
    G_complex_value,
    G_explicit,
    "G coefficient-weighted-sum identity",
)

print("  CHECK = PASS")


# ============================================================
# POLE BASIS
# ============================================================

pole_basis_values = [
    pole_basis(k, L)
    for k in range(N + 1)
]

print()
print("pole_basis(k, L)")
print("----------------")
print_vector(pole_basis_values)


# ============================================================
# POLE ROW
# ============================================================

P_row = pole_row(N, L)

print()
print("pole_row(N, L)")
print("--------------")

for i in range(P_row.rows):
    print(
        f"  [{i}] = {fmt(P_row[i, 0])}"
    )


# ============================================================
# T
# ============================================================

T_value = T_canonical(V, T)

print()
print("T_canonical(V, t)")
print("-----------------")
print(f"  value = {fmt(T_value)}")


# ============================================================
# K
# ============================================================

K_value = K_canonical(V, OMEGA)

print()
print("K_canonical(V, omega)")
print("---------------------")
print(f"  value = {fmt(K_value)}")


# ============================================================
# GHAT
# ============================================================

ghat_value = ghat(V, XI, L)

print()
print("ghat(V, xi, L)")
print("---------------")
print(f"  value = {fmt(ghat_value)}")


# ============================================================
# COMPLETION
# ============================================================

print()
print("=" * 78)
print("ALL SEMANTIC CHECKS: PASS")
print("=" * 78)
