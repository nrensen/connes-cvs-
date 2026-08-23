# ============================================================
# CELL 16 — ARCHIMEDEAN DISCREPANCY LOCALISATION
#
# Purpose:
#   Localise the remaining Cell-5 Archimedean discrepancy.
#
# Strategy:
#   Compare the repository and explicit constructions at
#   progressively deeper levels:
#
#       S
#       dS
#       divided differences
#       basis kernel
#       quadratic-form integrand
#       r integration
#
# This is deliberately a diagnostic calculation, not a
# production calculation.
# ============================================================

import mpmath as mp

from cell import (
    compute_L,
    canonical_pairs,
    F_basis,
    F_vector,
    Fprime_basis,
    G_complex_basis,
    G_complex,
    canonical_to_full,
    full_to_canonical,
)

from connes_cvs import build_galerkin_matrix, compute_ground_state


# ============================================================
# PARAMETERS
# ============================================================

mp.mp.dps = 80

c = 13
N = 8
T = 40

L = compute_L(c)

print("=" * 70)
print("CELL_16 — ARCHIMEDEAN DISCREPANCY LOCALISATION")
print("=" * 70)
print()
print("Parameters:")
print(f"c = {c}")
print(f"N = {N}")
print(f"T = {T}")
print(f"dps = {mp.mp.dps}")
print(f"L = {mp.nstr(L, 60)}")
print()


# ============================================================
# 1. GROUND STATE
# ============================================================

print("-" * 70)
print("1. GROUND STATE")
print("-" * 70)
print()

Q = build_galerkin_matrix(
    c=c,
    N=N,
    T=T,
    dps=mp.mp.dps,
)

lambda_min, u_star = compute_ground_state(Q)

u_star = mp.matrix(u_star)

norm_u = mp.sqrt(
    mp.fdot(u_star, u_star)
)

u_star = u_star / norm_u

v_star = full_to_canonical(u_star, N)

print("lambda_min =")
print(mp.nstr(lambda_min, 60))
print()

print("||u_star|| =")
print(mp.nstr(mp.sqrt(mp.fdot(u_star, u_star)), 60))
print()

print("||v_star|| =")
print(mp.nstr(mp.sqrt(mp.fdot(v_star, v_star)), 60))
print()


# ============================================================
# 2. IMPORT THE CELL-5 ARCHIMEDEAN DEFINITIONS
#
# We deliberately define a clean local version here.
#
# The source used in the Cell-5 investigation is the
# Archimedean function
#
#   S(r,x)
#
# and its r derivative.
#
# We keep the definitions explicit so that this cell can
# compare the scalar source independently of the quadratic
# form machinery.
# ============================================================

def h_arch(tau):
    """
    Archimedean h_+(tau):

        Re psi(1/4 + i tau/2) - log(pi)

    This is the repository convention used in the
    Archimedean term.
    """
    tau = mp.mpf(tau)

    return (
        mp.re(
            mp.digamma(
                mp.mpf("0.25") + 0.5j * tau
            )
        )
        - mp.log(mp.pi)
    )


def S_arch(r, x):
    """
    Direct Archimedean source.

    This is written in the same structural form used by
    the explicit Cell-5 calculation.
    """
    r = mp.mpf(r)
    x = mp.mpf(x)

    total = mp.mpf("0")

    # Symmetric Fourier frequencies.
    for k in range(-N, N + 1):

        if k == 0:
            ck = mp.mpf(1)
        else:
            ck = mp.mpf(1) / mp.sqrt(2)

        a = 2 * mp.pi * k / L

        total += (
            ck
            * (
                h_arch(a + r)
                + h_arch(a - r)
            )
        )

    return total * x


# ============================================================
# 3. SOURCE AND DERIVATIVE CHECK
# ============================================================

print("-" * 70)
print("2. SOURCE-LEVEL CHECK")
print("-" * 70)
print()

TEST_R = [
    mp.mpf("0"),
    mp.mpf("0.5"),
    mp.mpf("1"),
    mp.mpf("2"),
    mp.mpf("2.4496332798546520107426711685"),
    mp.mpf("3"),
    mp.mpf("5"),
]

TEST_X = [
    mp.mpf("0.25"),
    mp.mpf("0.5"),
    mp.mpf("1"),
    mp.mpf("2"),
]

print("The scalar source is sampled at several (r,x) points.")
print()


# ------------------------------------------------------------
# Direct numerical derivative of S
# ------------------------------------------------------------

def dS_numeric(r, x):
    r = mp.mpf(r)
    x = mp.mpf(x)

    return mp.diff(
        lambda rr: S_arch(rr, x),
        r,
    )


for x in TEST_X:

    print(f"x = {mp.nstr(x, 12)}")

    for r in TEST_R:

        S0 = S_arch(r, x)
        dS = dS_numeric(r, x)

        print(
            "  r =",
            mp.nstr(r, 12),
            " S =",
            mp.nstr(S0, 30),
            " dS =",
            mp.nstr(dS, 30),
        )

    print()


# ============================================================
# 4. BASIS-LEVEL RESPONSE
# ============================================================

print("-" * 70)
print("3. BASIS-LEVEL RESPONSE")
print("-" * 70)
print()

print(
    "For each canonical basis vector e_k we compare:"
)
print()
print("  F_k(r)")
print("  F'_k(r)")
print("  G_k(r)")
print("  G'_k(r)")
print()

BASIS_K = range(0, min(N, 4) + 1)

print("NOTE: Fprime_basis(0,0) has a removable singularity in")
print("the closed-form expression; the derivative is evaluated")
print("there by differentiating F_basis numerically.")
print()

for k in BASIS_K:

    print(f"basis k = {k}")

    for r in TEST_R:

        F0 = F_basis(k, r, L)

        # Fprime_basis has a removable singularity in its
        # closed-form expression at (k,r) = (0,0).
        # Do not evaluate the raw formula there.
        if k == 0 and r == 0:
            Fp = mp.diff(
                lambda rr:
                    F_basis(k, rr, L),
                r,
            )
        else:
            Fp = Fprime_basis(k, r, L)

        G0 = G_complex_basis(k, r, L)

        Gp = mp.diff(
            lambda rr:
                G_complex_basis(k, rr, L),
            r,
        )

        print(
            "  r =",
            mp.nstr(r, 12),
        )
        print(
            "    F      =",
            mp.nstr(F0, 30),
        )
        print(
            "    Fprime =",
            mp.nstr(Fp, 30),
        )
        print(
            "    Re G   =",
            mp.nstr(mp.re(G0), 30),
        )
        print(
            "    Re G'  =",
            mp.nstr(mp.re(Gp), 30),
        )

    print()


# ============================================================
# 5. QUADRATIC RESPONSE OF THE GROUND STATE
# ============================================================

print("-" * 70)
print("4. GROUND-STATE RESPONSE")
print("-" * 70)
print()

print(
    "Compare the canonical ground-state response against"
)
print(
    "the scalar Archimedean source."
)
print()

for r in TEST_R:

    Fv = F_vector(
        v_star,
        r,
        L,
    )

    Gv = G_complex(
        v_star,
        r,
        L,
    )

    print("r =", mp.nstr(r, 20))

    print(
        "  F_v(r) =",
        mp.nstr(Fv, 40),
    )

    print(
        "  Re G_v(r) =",
        mp.nstr(mp.re(Gv), 40),
    )

    print()


# ============================================================
# 6. SOURCE / RESPONSE PRODUCTS
#
# The important diagnostic here is not a ratio.  We compare
# the actual products that enter a quadratic functional.
# ============================================================

print("-" * 70)
print("5. SOURCE / RESPONSE PRODUCTS")
print("-" * 70)
print()

print(
    "The quantities below are deliberately shown without"
)
print(
    "normalising by one another."
)
print()

for x in TEST_X:

    print("x =", mp.nstr(x, 12))

    for r in TEST_R:

        S0 = S_arch(r, x)

        Gv = mp.re(
            G_complex(
                v_star,
                r,
                L,
            )
        )

        product = S0 * Gv

        print(
            "  r =",
            mp.nstr(r, 12),
            " S =",
            mp.nstr(S0, 25),
            " ReG =",
            mp.nstr(Gv, 25),
            " S*ReG =",
            mp.nstr(product, 25),
        )

    print()


# ============================================================
# 7. FINITE-DIFFERENCE / DIVIDED-DIFFERENCE TEST
#
# The original Cell-5 machinery ultimately converts source
# values into a matrix kernel.  Before attempting the full
# matrix, inspect the scalar divided difference directly.
#
# For a smooth scalar S:
#
#       DD(a,b) = (S(a)-S(b))/(a-b)
#
# and on the diagonal:
#
#       DD(a,a) = S'(a).
#
# ============================================================

print("-" * 70)
print("6. DIVIDED-DIFFERENCE CHECK")
print("-" * 70)
print()

def divdiff_S(a, b, x):
    a = mp.mpf(a)
    b = mp.mpf(b)
    x = mp.mpf(x)

    if a == b:
        return dS_numeric(a, x)

    return (
        S_arch(a, x)
        - S_arch(b, x)
    ) / (a - b)


DD_PAIRS = [
    (mp.mpf("0.5"), mp.mpf("1.0")),
    (mp.mpf("1.0"), mp.mpf("2.0")),
    (mp.mpf("2.0"), mp.mpf("3.0")),
    (mp.mpf("3.0"), mp.mpf("5.0")),
]

for x in TEST_X:

    print("x =", mp.nstr(x, 12))

    for a, b in DD_PAIRS:

        dd = divdiff_S(a, b, x)

        midpoint = (a + b) / 2

        deriv = dS_numeric(
            midpoint,
            x,
        )

        print(
            "  (a,b) =",
            mp.nstr(a, 8),
            mp.nstr(b, 8),
        )
        print(
            "    DD      =",
            mp.nstr(dd, 30),
        )
        print(
            "    S'(mid) =",
            mp.nstr(deriv, 30),
        )
        print(
            "    DD-S'   =",
            mp.nstr(dd - deriv, 30),
        )

    print()


# ============================================================
# 8. INTEGRAND LOCALISATION
#
# We now construct the simplest scalar r-integrand associated
# with the ground-state response and compare its behaviour.
#
# This is NOT claimed to be the final Cell-5 functional.
# It is a diagnostic probe.
# ============================================================

print("-" * 70)
print("7. INTEGRAND LOCALISATION")
print("-" * 70)
print()

def diagnostic_integrand(r, x):

    r = mp.mpf(r)
    x = mp.mpf(x)

    S0 = S_arch(r, x)

    Gv = mp.re(
        G_complex(
            v_star,
            r,
            L,
        )
    )

    return S0 * Gv


for x in TEST_X:

    print("x =", mp.nstr(x, 12))

    for r in TEST_R:

        val = diagnostic_integrand(r, x)

        print(
            "  r =",
            mp.nstr(r, 12),
            " value =",
            mp.nstr(val, 35),
        )

    print()


# ============================================================
# 9. SMALL LOCAL INTEGRALS
#
# Integrate only over short intervals.  This helps distinguish
# a local integrand discrepancy from an outer-quadrature issue.
# ============================================================

print("-" * 70)
print("8. LOCAL INTEGRALS")
print("-" * 70)
print()

R_INTERVALS = [
    (mp.mpf("0"), mp.mpf("1")),
    (mp.mpf("1"), mp.mpf("2")),
    (mp.mpf("2"), mp.mpf("3")),
    (mp.mpf("3"), mp.mpf("5")),
    (mp.mpf("5"), mp.mpf("10")),
]

for x in TEST_X:

    print("x =", mp.nstr(x, 12))

    total = mp.mpf("0")

    for a, b in R_INTERVALS:

        val = mp.quad(
            lambda rr:
                diagnostic_integrand(rr, x),
            [a, b],
        )

        total += val

        print(
            "  [",
            mp.nstr(a, 8),
            ",",
            mp.nstr(b, 8),
            "] =",
            mp.nstr(val, 35),
        )

    print(
        "  accumulated =",
        mp.nstr(total, 35),
    )

    print()


# ============================================================
# 10. SUMMARY
# ============================================================

print("=" * 70)
print("CELL_16 SUMMARY")
print("=" * 70)
print()
print(
    "This cell is a localisation experiment."
)
print()
print(
    "The intended interpretation is:"
)
print()
print(
    "  1. If S itself disagrees, revisit the Archimedean source."
)
print(
    "  2. If S agrees but divided differences disagree, the"
)
print(
    "     source-to-kernel dictionary is implicated."
)
print(
    "  3. If basis responses disagree, inspect F/G construction."
)
print(
    "  4. If pointwise quantities agree but local integrals"
)
print(
    "     disagree, inspect the integration/domain convention."
)
print()
print(
    "No conclusion should be drawn from a ratio alone."
)
print(
    "The absolute quantities are the primary diagnostic."
)
print()
print("=" * 70)
print("END CELL_16")
print("=" * 70)
