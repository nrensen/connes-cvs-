# ============================================================
# CELL 15 — CELL-5 DISCREPANCY CLOSURE AUDIT
#
# Purpose
# -------
# Close the historical Cell-5 discrepancy now that Cell 14 has
# established the full/canonical coordinate distinction.
#
# This cell deliberately does NOT reproduce the entire original
# Cell-5 calculation.
#
# Instead it performs four independent checks:
#
#   1. Construct the ground state in canonical coordinates.
#   2. Construct its full Fourier representation correctly.
#   3. Deliberately reproduce the OLD Cell-5 coordinate mistake.
#   4. Compare the resulting quadratic forms in the two
#      representations.
#
# The central question is:
#
#       Was the historical discrepancy caused by feeding a
#       full-space vector into a canonical-space functional?
#
# We also verify the result directly at the matrix level.
#
# ============================================================

from __future__ import annotations

import mpmath as mp

from cell import (
    build_galerkin_matrix,
    compute_ground_state,
    canonical_to_full,
    full_to_canonical,
    canonical_pairs,
)


# ============================================================
# PARAMETERS
# ============================================================

C = 13
N = 8
T = 40
DPS = 80

mp.mp.dps = DPS

L = mp.log(mp.mpf(C))


# ============================================================
# HELPERS
# ============================================================

def norm_sq(v):
    return mp.fsum(abs(v[k]) ** 2 for k in range(len(v)))


def max_vector_error(a, b):
    return max(
        abs(a[k] - b[k])
        for k in range(len(a))
    )


def canonical_embedding_matrix(N):
    """
    E maps canonical coordinates v in C^(N+1)
    into the symmetric full Fourier coordinates
    u in C^(2N+1).

    k=0:
        u_0 = v_0

    k>0:
        u_{+k} = u_{-k} = v_k/sqrt(2)
    """

    E = mp.matrix(2 * N + 1, N + 1)

    for k in range(N + 1):
        for m, coeff in canonical_pairs(k):
            E[m + N, k] = coeff

    return E


def canonical_gram_from_full(Q_full, N):
    """
    Pull a full-space quadratic form back into canonical
    coordinates:

        Q_can = E^T Q_full E

    where E is the canonical embedding.
    """

    E = canonical_embedding_matrix(N)

    return E.T * Q_full * E


def quadratic_form(Q, v):
    return mp.fdot(v, Q * v)


# ============================================================
# HEADER
# ============================================================

print()
print("=" * 70)
print("CELL 15 — CELL-5 DISCREPANCY CLOSURE AUDIT")
print("=" * 70)

print()
print("Parameters:")
print("c =", C)
print("N =", N)
print("T =", T)
print("dps =", DPS)
print("L =", mp.nstr(L, 60))


# ============================================================
# 1. BUILD GROUND STATE
# ============================================================

print()
print("-" * 70)
print("1. GROUND STATE")
print("-" * 70)

Q_full = build_galerkin_matrix(
    c=C,
    N=N,
    T=T,
    dps=DPS,
)

lambda_min, u_star = compute_ground_state(Q_full)

print()
print("lambda_min =")
print(mp.nstr(lambda_min, 70))

print()
print("full dimension =", len(u_star))

print()
print("||u_star||^2 =")
print(mp.nstr(norm_sq(u_star), 70))

print()
print("| ||u_star||^2 - 1 | =")
print(mp.nstr(abs(norm_sq(u_star) - 1), 30))


# ============================================================
# 2. FULL -> CANONICAL -> FULL
# ============================================================

print()
print("-" * 70)
print("2. FULL / CANONICAL REPRESENTATION")
print("-" * 70)

v_star = full_to_canonical(u_star, N)

u_roundtrip = canonical_to_full(v_star, N)

print()
print("canonical dimension =", len(v_star))

print()
print("||v_star||^2 =")
print(mp.nstr(norm_sq(v_star), 70))

print()
print("| ||v_star||^2 - 1 | =")
print(mp.nstr(abs(norm_sq(v_star) - 1), 30))

print()
print("max full -> canonical -> full error =")
print(
    mp.nstr(
        max_vector_error(u_star, u_roundtrip),
        40,
    )
)


# ============================================================
# 3. CANONICAL MATRIX REPRESENTATION
# ============================================================

print()
print("-" * 70)
print("3. CANONICAL MATRIX")
print("-" * 70)

Q_can = canonical_gram_from_full(Q_full, N)

lambda_can = quadratic_form(Q_can, v_star)

print()
print("v_star^T Q_can v_star =")
print(mp.nstr(lambda_can, 70))

print()
print("ground-state lambda_min =")
print(mp.nstr(lambda_min, 70))

print()
print("|canonical quadratic form - lambda_min| =")
print(
    mp.nstr(
        abs(lambda_can - lambda_min),
        40,
    )
)


# ============================================================
# 4. FULL-SPACE QUADRATIC FORM
# ============================================================

print()
print("-" * 70)
print("4. FULL-SPACE QUADRATIC FORM")
print("-" * 70)

lambda_full = quadratic_form(Q_full, u_star)

print()
print("u_star^T Q_full u_star =")
print(mp.nstr(lambda_full, 70))

print()
print("|full quadratic form - lambda_min| =")
print(
    mp.nstr(
        abs(lambda_full - lambda_min),
        40,
    )
)

print()
print("|full - canonical| =")
print(
    mp.nstr(
        abs(lambda_full - lambda_can),
        40,
    )
)


# ============================================================
# 5. DELIBERATELY REPRODUCE THE OLD MISTAKE
# ============================================================

print()
print("-" * 70)
print("5. DELIBERATELY REPRODUCE THE HISTORICAL MISTAKE")
print("-" * 70)

print()
print(
    "The historical mistake is to interpret the first N+1 "
    "entries of the full vector as canonical coefficients."
)

v_wrong = mp.matrix(N + 1, 1)

for k in range(N + 1):
    v_wrong[k] = u_star[k]

wrong_norm = norm_sq(v_wrong)

print()
print("wrong coefficient vector norm^2 =")
print(mp.nstr(wrong_norm, 70))

print()
print("correct canonical norm^2 =")
print(mp.nstr(norm_sq(v_star), 70))

print()
print("|wrong norm^2 - correct norm^2| =")
print(
    mp.nstr(
        abs(wrong_norm - norm_sq(v_star)),
        40,
    )
)


# ============================================================
# 6. WRONG CANONICAL QUADRATIC FORM
# ============================================================

print()
print("-" * 70)
print("6. WRONG CANONICAL QUADRATIC FORM")
print("-" * 70)

wrong_Q_value = quadratic_form(
    Q_can,
    v_wrong,
)

correct_Q_value = quadratic_form(
    Q_can,
    v_star,
)

print()
print("Q_can(v_wrong) =")
print(mp.nstr(wrong_Q_value, 70))

print()
print("Q_can(v_star) =")
print(mp.nstr(correct_Q_value, 70))

print()
print("|wrong - correct| =")
print(
    mp.nstr(
        abs(wrong_Q_value - correct_Q_value),
        40,
    )
)


# ============================================================
# 7. BASIS-BY-BASIS EMBEDDING CHECK
# ============================================================

print()
print("-" * 70)
print("7. BASIS-BY-BASIS EMBEDDING CHECK")
print("-" * 70)

E = canonical_embedding_matrix(N)

max_embedding_error = mp.mpf("0")

for k in range(N + 1):

    e = mp.matrix(N + 1, 1)
    e[k] = 1

    embedded = E * e

    expected = mp.matrix(2 * N + 1, 1)

    for m, coeff in canonical_pairs(k):
        expected[m + N] = coeff

    err = max_vector_error(
        embedded,
        expected,
    )

    max_embedding_error = max(
        max_embedding_error,
        err,
    )

    print(
        f"k = {k:2d}"
        f"   max embedding error = "
        f"{mp.nstr(err, 12)}"
    )

print()
print("maximum embedding error =")
print(mp.nstr(max_embedding_error, 30))


# ============================================================
# 8. RAYLEIGH QUOTIENT CHECK
# ============================================================

print()
print("-" * 70)
print("8. RAYLEIGH QUOTIENT CHECK")
print("-" * 70)

full_rayleigh = (
    quadratic_form(Q_full, u_star)
    / norm_sq(u_star)
)

canonical_rayleigh = (
    quadratic_form(Q_can, v_star)
    / norm_sq(v_star)
)

wrong_rayleigh = (
    quadratic_form(Q_can, v_wrong)
    / norm_sq(v_wrong)
)

print()
print("full Rayleigh quotient =")
print(mp.nstr(full_rayleigh, 70))

print()
print("canonical Rayleigh quotient =")
print(mp.nstr(canonical_rayleigh, 70))

print()
print("wrong-coordinate Rayleigh quotient =")
print(mp.nstr(wrong_rayleigh, 70))

print()
print("|full - canonical| =")
print(
    mp.nstr(
        abs(full_rayleigh - canonical_rayleigh),
        40,
    )
)


# ============================================================
# 9. CONVENTION DIAGNOSTIC
# ============================================================

print()
print("-" * 70)
print("9. CONVENTION DIAGNOSTIC")
print("-" * 70)

correct_conversion_error = max_vector_error(
    u_star,
    canonical_to_full(v_star, N),
)

wrong_conversion_error = max_vector_error(
    u_star,
    canonical_to_full(v_wrong, N),
)

print()
print("correct representation reconstruction error =")
print(mp.nstr(correct_conversion_error, 40))

print()
print("wrong representation reconstruction error =")
print(mp.nstr(wrong_conversion_error, 40))


# ============================================================
# 10. SUMMARY
# ============================================================

print()
print("=" * 70)
print("CELL 15 SUMMARY")
print("=" * 70)

print()
print("Ground-state normalization error:")
print(
    mp.nstr(
        abs(norm_sq(u_star) - 1),
        20,
    )
)

print()
print("Canonical normalization error:")
print(
    mp.nstr(
        abs(norm_sq(v_star) - 1),
        20,
    )
)

print()
print("Full/canonical round-trip error:")
print(
    mp.nstr(
        correct_conversion_error,
        20,
    )
)

print()
print("Full-vs-canonical quadratic-form difference:")
print(
    mp.nstr(
        abs(lambda_full - lambda_can),
        20,
    )
)

print()
print("Wrong-coordinate norm^2:")
print(
    mp.nstr(
        wrong_norm,
        30,
    )
)

print()
print("Correct-coordinate norm^2:")
print(
    mp.nstr(
        norm_sq(v_star),
        30,
    )
)

print()
print("Wrong-coordinate Rayleigh quotient:")
print(
    mp.nstr(
        wrong_rayleigh,
        30,
    )
)

print()
print("Correct-coordinate Rayleigh quotient:")
print(
    mp.nstr(
        canonical_rayleigh,
        30,
    )
)

print()
print("=" * 70)
print("END CELL 15")
print("=" * 70)
