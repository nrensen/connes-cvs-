# ============================================================
# CELL 14 — CORRECTED FULL/CANONICAL PARSEVAL AUDIT
#
# Purpose:
#   Confirm that the apparent Parseval discrepancy in Cell 13
#   was caused by treating the full (2N+1)-vector as though it
#   were the canonical (N+1)-vector.
#
# This cell:
#   1. Builds Q using the current cell.py implementation.
#   2. Computes the full ground-state eigenvector.
#   3. Converts full -> canonical correctly.
#   4. Converts canonical -> full and checks round-trip accuracy.
#   5. Checks normalization in both representations.
#   6. Reconstructs f(t) from the canonical coefficients.
#   7. Computes its L2 norm directly.
#   8. Computes Fourier coefficients directly from f(t).
#   9. Checks Parseval.
#
# No source/Weil dictionary is tested here.
# ============================================================

from __future__ import annotations

import mpmath as mp

from cell import (
    build_galerkin_matrix,
    compute_ground_state,
    full_to_canonical,
    canonical_to_full,
)


# ------------------------------------------------------------
# PARAMETERS
# ------------------------------------------------------------

C = 13
N = 8
T = 40
DPS = 80

mp.mp.dps = DPS

L = mp.log(mp.mpf(C))
OMEGA = 2 * mp.pi / L


# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------

def max_abs_vector_error(a, b):
    assert len(a) == len(b)
    return max(abs(a[i] - b[i]) for i in range(len(a)))


def l2_norm_sq_vector(v):
    return mp.fsum(abs(x) ** 2 for x in v)


def f_from_canonical(v, t):
    """
    Real trigonometric reconstruction from canonical coefficients

        v = (v_0, ..., v_N)

    with orthonormal basis

        phi_0(t) = 1/sqrt(L)
        phi_k(t) = sqrt(2/L) cos(2*pi*k*t/L), k >= 1.
    """
    total = v[0] / mp.sqrt(L)

    for k in range(1, N + 1):
        total += (
            v[k]
            * mp.sqrt(2 / L)
            * mp.cos(OMEGA * k * t)
        )

    return total


def fourier_coefficient_from_f(v, k):
    """
    Compute the complex Fourier coefficient

        u_k = (1/sqrt(L)) integral_0^L f(t)
              exp(-2*pi*i*k*t/L) dt

    directly from the reconstructed f.

    This uses the full Fourier-normalized coefficient convention.
    """
    integrand = lambda t: (
        f_from_canonical(v, t)
        * mp.exp(-2j * mp.pi * k * t / L)
        / mp.sqrt(L)
    )

    return mp.quad(integrand, [0, L])


def direct_function_norm_sq(v):
    """
    Direct continuous L2 norm:

        integral_0^L |f(t)|^2 dt.
    """
    integrand = lambda t: abs(f_from_canonical(v, t)) ** 2
    return mp.quad(integrand, [0, L])


# ------------------------------------------------------------
# HEADER
# ------------------------------------------------------------

print("=" * 70)
print("CELL 14 — CORRECTED FULL/CANONICAL PARSEVAL AUDIT")
print("=" * 70)

print()
print("Parameters:")
print(f"c = {C}")
print(f"N = {N}")
print(f"T = {T}")
print(f"dps = {DPS}")
print(f"L = {mp.nstr(L, 60)}")
print(f"2*pi/L = {mp.nstr(OMEGA, 60)}")


# ------------------------------------------------------------
# 1. BUILD OPERATOR
# ------------------------------------------------------------

print()
print("-" * 70)
print("1. GROUND STATE")
print("-" * 70)

Q = build_galerkin_matrix(
    c=C,
    N=N,
    T=T,
    dps=DPS,
)

lambda_min, full_star = compute_ground_state(Q)

print()
print("lambda_min =")
print(mp.nstr(lambda_min, 70))

print()
print("full dimension =", len(full_star))

full_norm_sq = l2_norm_sq_vector(full_star)

print()
print("||full_star||^2 =")
print(mp.nstr(full_norm_sq, 70))

print()
print("| ||full_star||^2 - 1 | =")
print(mp.nstr(abs(full_norm_sq - 1), 30))


# ------------------------------------------------------------
# 2. CORRECT FULL -> CANONICAL CONVERSION
# ------------------------------------------------------------

print()
print("-" * 70)
print("2. CORRECT FULL -> CANONICAL CONVERSION")
print("-" * 70)

v_star = full_to_canonical(full_star, N)

print()
print("canonical dimension =", len(v_star))

canonical_norm_sq = l2_norm_sq_vector(v_star)

print()
print("||v_star||^2 =")
print(mp.nstr(canonical_norm_sq, 70))

print()
print("| ||v_star||^2 - 1 | =")
print(mp.nstr(abs(canonical_norm_sq - 1), 30))

print()
print("| ||full_star||^2 - ||v_star||^2 | =")
print(mp.nstr(abs(full_norm_sq - canonical_norm_sq), 30))


# ------------------------------------------------------------
# 3. CANONICAL -> FULL ROUND TRIP
# ------------------------------------------------------------

print()
print("-" * 70)
print("3. CANONICAL -> FULL ROUND TRIP")
print("-" * 70)

full_roundtrip = canonical_to_full(v_star, N)

roundtrip_error = max_abs_vector_error(
    full_star,
    full_roundtrip,
)

print()
print("max |full_star - canonical_to_full(full_to_canonical(full_star))| =")
print(mp.nstr(roundtrip_error, 50))


# ------------------------------------------------------------
# 4. SHOW THE COORDINATES
# ------------------------------------------------------------

print()
print("-" * 70)
print("4. COORDINATE CHECK")
print("-" * 70)

print()
print(" k       canonical v_k                  full u_+k                 ")
print()

for k in range(N + 1):
    if k == 0:
        u_plus = full_roundtrip[N]
    else:
        u_plus = full_roundtrip[N + k]

    print(
        f"{k:2d}   "
        f"{mp.nstr(v_star[k], 30):>32}   "
        f"{mp.nstr(u_plus, 30):>32}"
    )


# ------------------------------------------------------------
# 5. RECONSTRUCT f(t)
# ------------------------------------------------------------

print()
print("-" * 70)
print("5. FUNCTION RECONSTRUCTION")
print("-" * 70)

sample_points = [
    mp.mpf("0"),
    L / 8,
    L / 4,
    3 * L / 8,
    L / 2,
    5 * L / 8,
    3 * L / 4,
    7 * L / 8,
]

print()
for t in sample_points:
    print(
        "t =",
        mp.nstr(t, 25),
        "   f(t) =",
        mp.nstr(f_from_canonical(v_star, t), 40),
    )


# ------------------------------------------------------------
# 6. DIRECT CONTINUOUS L2 NORM
# ------------------------------------------------------------

print()
print("-" * 70)
print("6. DIRECT CONTINUOUS L2 NORM")
print("-" * 70)

function_norm_sq = direct_function_norm_sq(v_star)

print()
print("integral_0^L |f(t)|^2 dt =")
print(mp.nstr(function_norm_sq, 70))

print()
print("difference from canonical norm =")
print(mp.nstr(abs(function_norm_sq - canonical_norm_sq), 40))

print()
print("difference from 1 =")
print(mp.nstr(abs(function_norm_sq - 1), 40))


# ------------------------------------------------------------
# 7. DIRECT FOURIER COEFFICIENTS
# ------------------------------------------------------------

print()
print("-" * 70)
print("7. DIRECT FOURIER COEFFICIENTS")
print("-" * 70)

print()
print("Comparing direct Fourier coefficients against")
print("canonical_to_full(v_star, N).")
print()

fourier_errors = []

for k in range(-N, N + 1):

    u_expected = full_roundtrip[k + N]

    u_direct = fourier_coefficient_from_f(v_star, k)

    err = abs(u_direct - u_expected)
    fourier_errors.append(err)

    print(
        f"k = {k:2d}"
        f"   expected = {mp.nstr(u_expected, 30)}"
        f"   direct = {mp.nstr(u_direct, 30)}"
        f"   |error| = {mp.nstr(err, 12)}"
    )

max_fourier_error = max(fourier_errors)

print()
print("Maximum Fourier coefficient error =")
print(mp.nstr(max_fourier_error, 40))


# ------------------------------------------------------------
# 8. PARSEVAL IN FULL FOURIER COORDINATES
# ------------------------------------------------------------

print()
print("-" * 70)
print("8. PARSEVAL — FULL FOURIER COORDINATES")
print("-" * 70)

u_direct = []

for k in range(-N, N + 1):
    u_direct.append(
        fourier_coefficient_from_f(v_star, k)
    )

fourier_norm_sq = mp.fsum(abs(u) ** 2 for u in u_direct)

print()
print("sum_{k=-N}^N |u_k|^2 =")
print(mp.nstr(fourier_norm_sq, 70))

print()
print("continuous norm =")
print(mp.nstr(function_norm_sq, 70))

print()
print("|Fourier norm - continuous norm| =")
print(mp.nstr(abs(fourier_norm_sq - function_norm_sq), 40))

print()
print("|Fourier norm - 1| =")
print(mp.nstr(abs(fourier_norm_sq - 1), 40))


# ------------------------------------------------------------
# 9. PARSEVAL IN CANONICAL REAL COORDINATES
# ------------------------------------------------------------

print()
print("-" * 70)
print("9. PARSEVAL — CANONICAL REAL COORDINATES")
print("-" * 70)

canonical_parseval = (
    abs(v_star[0]) ** 2
    + mp.fsum(abs(v_star[k]) ** 2 for k in range(1, N + 1))
)

print()
print("sum_{k=0}^N |v_k|^2 =")
print(mp.nstr(canonical_parseval, 70))

print()
print("continuous norm =")
print(mp.nstr(function_norm_sq, 70))

print()
print("|canonical norm - continuous norm| =")
print(mp.nstr(abs(canonical_parseval - function_norm_sq), 40))


# ------------------------------------------------------------
# 10. THE OLD CELL-13 ERROR REPRODUCED
#
# This is included deliberately. It demonstrates that if one
# incorrectly treats the full vector as canonical, one recovers
# the anomalous norm from Cell 13.
# ------------------------------------------------------------

print()
print("-" * 70)
print("10. REPRODUCE THE CELL-13 COORDINATE MISTAKE")
print("-" * 70)

v_wrong = [
    full_star[k]
    for k in range(N + 1)
]

wrong_norm_sq = l2_norm_sq_vector(v_wrong)

print()
print("Incorrectly treating full_star[0:N+1] as canonical:")
print()
print("wrong coefficient norm =")
print(mp.nstr(wrong_norm_sq, 70))

print()
print("correct canonical norm =")
print(mp.nstr(canonical_norm_sq, 70))

print()
print("|wrong - correct| =")
print(mp.nstr(abs(wrong_norm_sq - canonical_norm_sq), 40))


# ------------------------------------------------------------
# 11. SUMMARY
# ------------------------------------------------------------

print()
print("=" * 70)
print("CELL 14 SUMMARY")
print("=" * 70)

print()
print("Ground-state full-vector normalization error:")
print(
    mp.nstr(abs(full_norm_sq - 1), 20)
)

print()
print("Full -> canonical norm error:")
print(
    mp.nstr(abs(canonical_norm_sq - 1), 20)
)

print()
print("Full -> canonical -> full round-trip error:")
print(
    mp.nstr(roundtrip_error, 20)
)

print()
print("Direct Fourier coefficient maximum error:")
print(
    mp.nstr(max_fourier_error, 20)
)

print()
print("Continuous-vs-Fourier Parseval error:")
print(
    mp.nstr(abs(fourier_norm_sq - function_norm_sq), 20)
)

print()
print("Canonical-vs-continuous Parseval error:")
print(
    mp.nstr(abs(canonical_parseval - function_norm_sq), 20)
)

print()
print("Incorrect Cell-13 coefficient norm:")
print(
    mp.nstr(wrong_norm_sq, 30)
)

print()
print("Correct coefficient norm:")
print(
    mp.nstr(canonical_norm_sq, 30)
)

print()
print("=" * 70)
print("END CELL 14")
print("=" * 70)
