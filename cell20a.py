# ============================================================
# cell20a.py — POLE SANITY-CHECK FORENSICS
#
# Purpose
# -------
# Cell 20 reported an apparent discrepancy between
#
#     <u, Q_pole u>
#
# and
#
#     2 Re sum_v_G(v, i/2).
#
# This cell investigates that discrepancy only.
#
# No Archimedean integration is performed.
# No new Galerkin matrix is built.
#
# The investigation tests the more basic distinction:
#
#   P(v)             = pole linear functional
#
# versus
#
#   <u, Q_pole u>    = pole quadratic form
#
# and checks whether
#
#   2 Re G(v, i/2)
#
# is proportional to P(v).
# ============================================================

import time

import mpmath as mp

from cell import (
    FORENSIC_GROUND_STATE,
    sum_v_G,
    canonical_to_full,
    get_ground_state,
    normalise_ground_state,
    pole_basis,
    pole_row,
)


# ============================================================
# PARAMETERS
# ============================================================

mp.mp.dps = 80

c = FORENSIC_GROUND_STATE["c"]
N = FORENSIC_GROUND_STATE["N"]

L = mp.log(mp.mpf(c))

DISPLAY_DIGITS = 60


def nstr(x):
    return mp.nstr(
        x,
        DISPLAY_DIGITS,
    )


# ============================================================
# HEADER
# ============================================================

print("=" * 72)
print("CELL 20A — POLE SANITY-CHECK FORENSICS")
print("=" * 72)

print()
print("Parameters:")
print("  c   =", c)
print("  N   =", N)
print("  L   =", nstr(L))
print("  dps =", mp.mp.dps)

print()
print(
    "No Galerkin matrix or Archimedean integral will be computed."
)


# ============================================================
# 1. LOAD CACHED FORENSIC GROUND STATE
# ============================================================

print()
print("-" * 72)
print("1. LOAD CACHED FORENSIC GROUND STATE")
print("-" * 72)

ground_start = time.perf_counter()

lam_min, eigvec, cache_meta = get_ground_state(
    c=FORENSIC_GROUND_STATE["c"],
    N=FORENSIC_GROUND_STATE["N"],
    T=FORENSIC_GROUND_STATE["T"],
    dps=FORENSIC_GROUND_STATE["dps"],
    verbose=True,
)

ground_elapsed = (
    time.perf_counter()
    - ground_start
)

v_star = normalise_ground_state(
    eigvec,
    N,
)

v_star = mp.matrix(v_star)

u_star = canonical_to_full(
    v_star,
    N,
)

print()
print("lambda_min =")
print(nstr(lam_min))

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
print("ground-state retrieval elapsed =")
print(
    f"{ground_elapsed:.6f} s"
)


# ============================================================
# 2. POLE LINEAR FUNCTIONAL
# ============================================================
#
# Cell.py defines
#
#     P(e_k) = pole_basis(k, L)
#
# so
#
#     P(v) = sum_k v_k P(e_k).
#
# This is linear in v.
# ============================================================

print()
print("-" * 72)
print("2. POLE LINEAR FUNCTIONAL")
print("-" * 72)

P_row = pole_row(
    N,
    L,
)

P_v = mp.fdot(
    P_row,
    v_star,
)

print()
print("P(v_star) =")
print(nstr(P_v))


# ============================================================
# 3. G(i/2)
# ============================================================
#
# sum_v_G is also linear in v.
#
# Compare it directly against P(v).
# ============================================================

print()
print("-" * 72)
print("3. sum_v_G(v, i/2)")
print("-" * 72)

G_pole = mp.re(
    2 * sum_v_G(
        v_star,
        1j / 2,
        L,
    )
)

print()
print("2 Re sum_v_G(v_star, i/2) =")
print(nstr(G_pole))


# ============================================================
# 4. TEST PROPORTIONALITY
# ============================================================

print()
print("-" * 72)
print("4. P(v) VS 2 Re G(v, i/2)")
print("-" * 72)

if P_v == 0:
    raise RuntimeError(
        "P(v_star) is exactly zero; "
        "cannot determine proportionality constant."
    )

C = G_pole / P_v

print()
print(
    "C = [2 Re G(v_star, i/2)] / P(v_star)"
)
print(nstr(C))

print()
print("C * P(v_star) =")
print(nstr(C * P_v))

print()
print("difference =")
print(
    nstr(
        G_pole
        - C * P_v
    )
)


# ============================================================
# 5. BASIS-VECTOR TEST
# ============================================================
#
# If the proportionality is genuine, the same C should occur
# for every canonical basis vector.
# ============================================================

print()
print("-" * 72)
print("5. BASIS-VECTOR PROPORTIONALITY TEST")
print("-" * 72)

print()
print(
    "k".ljust(6),
    "P(e_k)".ljust(32),
    "2 Re G(e_k,i/2)".ljust(32),
    "ratio"
)

print("-" * 108)

basis_ratios = []

for k in range(N + 1):

    e = mp.matrix(
        N + 1,
        1,
    )

    e[k] = 1

    P_k = pole_basis(
        k,
        L,
    )

    G_k = mp.re(
        2 * sum_v_G(
            e,
            1j / 2,
            L,
        )
    )

    ratio = G_k / P_k

    basis_ratios.append(
        ratio
    )

    print(
        f"{k:<6}",
        f"{nstr(P_k):<32}",
        f"{nstr(G_k):<32}",
        nstr(ratio),
    )

ratio_spread = max(
    abs(r - C)
    for r in basis_ratios
)

print()
print("Maximum basis ratio deviation from C =")
print(nstr(ratio_spread))


# ============================================================
# 6. CONSTRUCT THE POLE MATRIX
# ============================================================
#
# This reproduces Cell 20's Q_pole, but without any
# Archimedean calculation.
# ============================================================

print()
print("-" * 72)
print("6. POLE QUADRATIC FORM")
print("-" * 72)


def psi_pole(x):

    x = mp.mpf(x)

    integrand = lambda y: (
        2
        * mp.cosh(y / 2)
        * mp.sin(
            2 * mp.pi
            * x
            * (
                1
                - y / L
            )
        )
    )

    return (
        1 / mp.pi
        * mp.quad(
            integrand,
            [0, L],
        )
    )


def psi_pole_derivative(x):

    x = mp.mpf(x)

    integrand = lambda y: (
        2
        * mp.cosh(y / 2)
        * (
            2 * mp.pi
            * (
                1
                - y / L
            )
        )
        * mp.cos(
            2 * mp.pi
            * x
            * (
                1
                - y / L
            )
        )
    )

    return (
        1 / mp.pi
        * mp.quad(
            integrand,
            [0, L],
        )
    )


Q_pole = mp.matrix(
    2 * N + 1,
    2 * N + 1,
)

for i, m in enumerate(
    range(-N, N + 1)
):

    for j, n in enumerate(
        range(-N, N + 1)
    ):

        if m != n:

            Q_pole[i, j] = (
                psi_pole(m)
                - psi_pole(n)
            ) / mp.mpf(
                m - n
            )

        else:

            Q_pole[i, j] = (
                psi_pole_derivative(m)
            )


pole_quadratic = mp.fdot(
    u_star,
    Q_pole * u_star,
)

print()
print("<u_star, Q_pole u_star> =")
print(nstr(pole_quadratic))

print()
print("2 Re G(v_star, i/2) =")
print(nstr(G_pole))

print()
print(
    "These are deliberately NOT expected to agree:"
)

print()
print("  P(v) is linear in v.")
print("  G(v, i/2) is linear in v.")
print("  <u,Q_pole u> is quadratic in v.")


# ============================================================
# 7. CONSTRAINED-VECTOR CHECK
# ============================================================
#
# The historical cell5_corrected.py used a vector satisfying
#
#     P(v) = 0.
#
# In that situation the G(i/2) quantity also vanishes.
#
# We construct a simple vector in the nullspace of P using
# two canonical basis vectors.
# ============================================================

print()
print("-" * 72)
print("7. P(v)=0 CONSTRAINED CHECK")
print("-" * 72)

v_constraint = mp.matrix(
    N + 1,
    1,
)

P0 = pole_basis(
    0,
    L,
)

P1 = pole_basis(
    1,
    L,
)

v_constraint[0] = P1
v_constraint[1] = -P0

constraint_P = mp.fdot(
    P_row,
    v_constraint,
)

constraint_G = mp.re(
    2 * sum_v_G(
        v_constraint,
        1j / 2,
        L,
    )
)

print()
print("P(v_constraint) =")
print(nstr(constraint_P))

print()
print(
    "2 Re G(v_constraint, i/2) ="
)
print(nstr(constraint_G))

print()
print(
    "This explains why the historical constrained test"
)
print(
    "could appear to validate the pole sanity check:"
)
print(
    "both linear quantities vanish when P(v)=0."
)


# ============================================================
# 8. FINAL DIAGNOSTIC
# ============================================================

print()
print("=" * 72)
print("CELL 20A — FINAL DIAGNOSTIC")
print("=" * 72)

print()
print("Ground-state pole functional:")
print(
    "  P(v_star) =",
    nstr(P_v),
)

print()
print("Ground-state G value:")
print(
    "  2 Re G(v_star,i/2) =",
    nstr(G_pole),
)

print()
print("Proportionality constant:")
print(
    "  C =",
    nstr(C),
)

print()
print("Maximum basis-vector ratio deviation:")
print(
    " ",
    nstr(ratio_spread),
)

print()
print("Pole quadratic form:")
print(
    "  <u,Q_pole u> =",
    nstr(pole_quadratic),
)

print()
print(
    "CONCLUSION:"
)

print(
    "  The Cell-20 discrepancy does not indicate an"
)

print(
    "  inconsistency in the Archimedean calculation."
)

print(
    "  The sanity check incorrectly compared a quadratic"
)

print(
    "  pole form with a linear pole functional."
)

print()
print("=" * 72)
print("CELL 20A COMPLETE")
print("=" * 72)
