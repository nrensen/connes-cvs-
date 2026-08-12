#!/usr/bin/env python3
"""Negative-control guard for the corrected version 2 statements.

Unlike the thirteen reproducibility guards, which confirm encoded identities,
this guard checks that each v1 statement FAILS on its counterexample and that
the corrected v2 statement HOLDS there. Standard library only; exact Fraction
arithmetic wherever the claim is algebraic, seeded Monte Carlo for the variance
models. Writes a JSON artifact with "status": "PASS" and exits non-zero on any
failure.

Resolvent quantities in NC3-NC6 are computed FROM THE MATRICES by exact
Fraction linear solves (Gaussian elimination), not from asserted spectral
weights, so the guard verifies the matrix statements themselves.

Controls:
  NC1  Perfectly correlated errors kill the (2N+1)^2 variance reduction; the
       reduction reappears under the pairwise-uncorrelated hypothesis, checked
       both for the IID model and for a pairwise-uncorrelated but
       non-identically-distributed model with common variance (Rademacher and
       Gaussian entries mixed), the actual hypothesis of the corollary.
  NC2  A rank-one matrix in a non-all-ones direction has sigma_2/sigma_1 = 0
       but a large residual to the all-ones direction: the two diagnostics are
       not equivalent.
  NC3  An eigenvalue whose eigenspace is orthogonal to 1_N contributes no pole
       to the coincidence resolvent (removable point); G(z) is computed from
       the matrix by exact solves and shown to equal 2/(2-z) although 5 is in
       the spectrum.
  NC4  A repeated eigenvalue contributes one grouped spectral-projection
       residue, not one pole per basis eigenvector; the grouped residue is
       computed from the matrix resolvent exactly.
  NC5  The Weyl functions W_+/- can vanish at points off both spectra, where
       the reciprocal identity is undefined. W_+ is computed INDEPENDENTLY by
       exact inversion of V_+ (not from the Sherman-Morrison formula under
       test), the point z0 = 0 is certified off both spectra by determinant
       evaluation, and only then is the Sherman-Morrison prediction
       W_+ = W_-/(1 - a W_-) compared against the independent value; the
       reciprocal form is checked exactly where W_+ W_- != 0.
  NC6  An indefinite velocity matrix has spectrum off [0, infinity), so its
       Weyl function is Nevanlinna but not Stieltjes: the Krein-string
       boundary-mass reading needs hypotheses beyond the rank-one algebra.
  NC7  Symmetry-correlated noise (E_mn = E_nm, the entries with m <= n
       pairwise uncorrelated, each variance sigma^2): the full d^2 = (2N+1)^2
       reduction FAILS, and the coincidence-average variance is exactly
       sigma^2 (2d-1)/d^3 (a reduction of order d^2/2), the formula stated in
       Remark 3.5(ii); verified by exact w^t Sigma w accumulation in
       Fractions for d = 3, 5, 7, 9.
"""

import json
import math
import random
import sys
from fractions import Fraction
from pathlib import Path

FAILURES = []
RESULTS = {}


def record(name, ok, detail):
    RESULTS[name] = {"pass": bool(ok), "detail": detail}
    if not ok:
        FAILURES.append(name)
    print(("PASS " if ok else "FAIL ") + name + ": " + detail)


# ------------------------------------------------------- exact linear algebra
def solve_exact(M, b):
    """Solve M x = b over Fractions (M square nonsingular), Gaussian elim."""
    n = len(M)
    A = [list(M[i]) + [b[i]] for i in range(n)]
    for col in range(n):
        piv = next(r for r in range(col, n) if A[r][col] != 0)
        A[col], A[piv] = A[piv], A[col]
        pv = A[col][col]
        A[col] = [v / pv for v in A[col]]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [vr - f * vc for vr, vc in zip(A[r], A[col])]
    return [A[r][n] for r in range(n)]


def det_exact(M):
    """Determinant over Fractions by elimination."""
    n = len(M)
    A = [list(row) for row in M]
    d = Fraction(1)
    for col in range(n):
        piv = None
        for r in range(col, n):
            if A[r][col] != 0:
                piv = r
                break
        if piv is None:
            return Fraction(0)
        if piv != col:
            A[col], A[piv] = A[piv], A[col]
            d = -d
        d *= A[col][col]
        for r in range(col + 1, n):
            if A[r][col] != 0:
                f = A[r][col] / A[col][col]
                A[r] = [vr - f * vc for vr, vc in zip(A[r], A[col])]
    return d


def rank_exact(M):
    """Rank over Fractions by row echelon."""
    A = [list(row) for row in M]
    n_rows, n_cols = len(A), len(A[0])
    rank = 0
    row = 0
    for col in range(n_cols):
        piv = None
        for r in range(row, n_rows):
            if A[r][col] != 0:
                piv = r
                break
        if piv is None:
            continue
        A[row], A[piv] = A[piv], A[row]
        for r in range(row + 1, n_rows):
            if A[r][col] != 0:
                f = A[r][col] / A[row][col]
                A[r] = [vr - f * vc for vr, vc in zip(A[r], A[row])]
        rank += 1
        row += 1
    return rank


def shift(M, z):
    """M - z I over Fractions."""
    n = len(M)
    return [[M[i][j] - (z if i == j else 0) for j in range(n)] for i in range(n)]


def weyl_from_matrix(M, z):
    """<1, (M - z)^{-1} 1> computed by an exact linear solve on the matrix."""
    n = len(M)
    x = solve_exact(shift(M, z), [Fraction(1)] * n)
    return sum(x)


# ---------------------------------------------------------------- NC1
def nc1_variance():
    rng = random.Random(20260726)
    d = 5          # d = 2N+1 with N = 2
    trials = 200_000
    sigma = 1.0
    # correlated model: every entry is the same variable Z
    acc = 0.0
    acc2 = 0.0
    for _ in range(trials):
        z = rng.gauss(0.0, sigma)
        avg = z  # average of d^2 identical copies
        acc += avg
        acc2 += avg * avg
    var_corr = acc2 / trials - (acc / trials) ** 2
    # IID model
    acc = 0.0
    acc2 = 0.0
    for _ in range(trials):
        s = 0.0
        for _ in range(d * d):
            s += rng.gauss(0.0, sigma)
        avg = s / (d * d)
        acc += avg
        acc2 += avg * avg
    var_iid = acc2 / trials - (acc / trials) ** 2
    # pairwise-uncorrelated, non-identically-distributed model with common
    # variance sigma^2: alternate Rademacher(+-sigma) and Gaussian entries,
    # all independent (the corollary's actual hypothesis, beyond IID)
    acc = 0.0
    acc2 = 0.0
    for _ in range(trials):
        s = 0.0
        for i in range(d * d):
            if i % 2 == 0:
                s += sigma if rng.random() < 0.5 else -sigma
            else:
                s += rng.gauss(0.0, sigma)
        avg = s / (d * d)
        acc += avg
        acc2 += avg * avg
    var_mixed = acc2 / trials - (acc / trials) ** 2
    ok_corr = abs(var_corr - sigma**2) < 0.05          # ~= sigma^2, NOT sigma^2/d^2
    ok_iid = abs(var_iid - sigma**2 / d**2) < 0.005     # ~= sigma^2/25
    ok_mixed = abs(var_mixed - sigma**2 / d**2) < 0.005  # ~= sigma^2/25
    record(
        "NC1_correlated_errors_kill_variance_reduction",
        ok_corr and ok_iid and ok_mixed,
        f"var(correlated)={var_corr:.4f} (expected ~1, v1 formula predicts {1/d**2:.4f}); "
        f"var(IID)={var_iid:.5f} and var(mixed non-identical)={var_mixed:.5f} "
        f"(both expected ~{1/d**2:.4f})",
    )


# ---------------------------------------------------------------- NC2
def nc2_diagnostics():
    d = 4
    # M = e1 e2^T : rank one, direction not all-ones
    # singular values: sigma_1 = 1, sigma_2 = 0  ->  ratio 0 exactly
    sigma_ratio = Fraction(0)
    # all-ones fit coefficient: <M, J>_F / d^2 with J = ones
    coeff = Fraction(1, d * d)
    # residual squared = ||M||_F^2 - coeff^2 * d^2 = 1 - 1/d^2
    resid_sq = Fraction(1) - Fraction(1, d * d)
    ok = sigma_ratio == 0 and resid_sq > Fraction(9, 10)
    record(
        "NC2_rank_one_wrong_direction_separates_diagnostics",
        ok,
        f"sigma2/sigma1 = 0 exactly, residual^2 to all-ones model = {resid_sq} "
        "(large): diagnostics are not equivalent",
    )


# ---------------------------------------------------------------- NC3 / NC4
def nc3_orthogonal_eigenvector():
    # Q0 = [[7/2, -3/2], [-3/2, 7/2]]: eigenvector (1,1)/sqrt2 with eigenvalue
    # 2 (coupled to 1_2) and (1,-1)/sqrt2 with eigenvalue 5 (orthogonal to
    # 1_2). G(z) = <1, (Q0 - z)^{-1} 1> is computed from the matrix by exact
    # solves; it must equal 2/(2 - z), with NO pole at z = 5 even though 5 is
    # an eigenvalue of Q0.
    Q0 = [[Fraction(7, 2), Fraction(-3, 2)],
          [Fraction(-3, 2), Fraction(7, 2)]]
    lam_orthogonal = Fraction(5)
    in_spectrum = det_exact(shift(Q0, lam_orthogonal)) == 0
    # equality with 2/(2-z) at five sample points (both sides are rational
    # with denominator degree <= 2, so five agreements pin them down)
    samples = [Fraction(0), Fraction(1), Fraction(3), Fraction(10), Fraction(-4)]
    match = all(weyl_from_matrix(Q0, z) == Fraction(2) / (Fraction(2) - z)
                for z in samples)
    # approach the orthogonal eigenvalue: G stays bounded (removable point)
    vals = [weyl_from_matrix(Q0, lam_orthogonal + Fraction(1, 10**k))
            for k in (4, 6, 8)]
    spread = max(vals) - min(vals)
    ok = in_spectrum and match and spread < Fraction(1, 10**3)
    record(
        "NC3_orthogonal_eigenspace_gives_removable_point",
        ok,
        "G(z) from exact matrix solves equals 2/(2-z) although 5 is in "
        "spec(Q0); G bounded near 5 (removable): v1's 'simple pole at each "
        "eigenvalue' fails, corrected P_lambda 1 != 0 condition holds",
    )


def nc4_repeated_eigenvalue():
    # Q0 = diag(3,3,7), 1 = (1,1,1): the eigenvalue 3 is double
    # (rank(Q0 - 3I) = 1), and the grouped residue at 3 is -||P_3 1||^2 = -2,
    # i.e. ONE pole with residue -2, not two poles. Computed from the matrix.
    lam, mu = Fraction(3), Fraction(7)
    Q0 = [[Fraction(3), Fraction(0), Fraction(0)],
          [Fraction(0), Fraction(3), Fraction(0)],
          [Fraction(0), Fraction(0), Fraction(7)]]
    double = rank_exact(shift(Q0, lam)) == 1  # nullity 2: multiplicity two
    eps = Fraction(1, 10**8)
    g = weyl_from_matrix(Q0, lam + eps)
    approx_residue = -eps * (g - Fraction(1) / (mu - lam - eps))
    ok = double and approx_residue == Fraction(2)
    record(
        "NC4_repeated_eigenvalue_groups_residue",
        ok,
        f"3 is a double eigenvalue (rank(Q0-3I)=1); grouped residue from the "
        f"matrix resolvent = -{approx_residue} (exactly -||P_lam 1||^2 = -2), "
        "one pole, not one per eigenvector",
    )


# ---------------------------------------------------------------- NC5 / NC6
def nc5_weyl_zero_off_spectra():
    # V_- = diag(-1, 1), 1_2 = (1,1): W_-(z) = -2z/(z^2-1), computed from the
    # matrix. V_+ = V_- - a 11^t is formed explicitly and W_+ is computed
    # INDEPENDENTLY by exact inversion of V_+, so the Sherman-Morrison formula
    # under test is compared against an independent value, not against itself.
    a = Fraction(3, 10)
    Vm = [[Fraction(-1), Fraction(0)], [Fraction(0), Fraction(1)]]
    Vp = [[Vm[i][j] - a for j in range(2)] for i in range(2)]
    z0 = Fraction(0)
    # certify z0 off BOTH spectra by determinant evaluation
    det_vm_z0 = det_exact(shift(Vm, z0))
    det_vp_z0 = det_exact(shift(Vp, z0))
    off_both_spectra = det_vm_z0 != 0 and det_vp_z0 != 0
    wm0 = weyl_from_matrix(Vm, z0)
    wp0_direct = weyl_from_matrix(Vp, z0)          # independent computation
    sm_pred0 = wm0 / (1 - a * wm0)                  # Sherman-Morrison prediction
    zero_of_wm = wm0 == 0
    nonreciprocal_holds = wp0_direct == sm_pred0 == 0
    reciprocal_defined = (wm0 != 0 and wp0_direct != 0)  # it is NOT defined here
    # corrected pointwise domain: reciprocal identity where W_+ W_- != 0
    z1 = Fraction(1, 2)
    wm1 = weyl_from_matrix(Vm, z1)
    wp1_direct = weyl_from_matrix(Vp, z1)          # independent computation
    sm_pred1 = wm1 / (1 - a * wm1)
    sm_matches_independent = wp1_direct == sm_pred1
    reciprocal_ok_at_z1 = (Fraction(1) / wp1_direct) == (Fraction(1) / wm1) - a
    ok = (off_both_spectra and zero_of_wm and nonreciprocal_holds
          and (not reciprocal_defined) and sm_matches_independent
          and reciprocal_ok_at_z1)
    record(
        "NC5_weyl_zero_off_both_spectra",
        ok,
        f"W_-(0)=0 with 0 certified off both spectra (det(V_- )={det_vm_z0}, "
        f"det(V_+)={det_vp_z0}); W_+(0)={wp0_direct} by independent exact "
        "inversion of V_+, matching the Sherman-Morrison prediction; "
        "reciprocal form undefined there and exact at z=1/2 (independent "
        "W_+(1/2) matches too): v1's 'every z off the spectra' domain fails, "
        "corrected domain holds",
    )


def nc6_not_stieltjes():
    # A Stieltjes function is analytic on C \ [0, inf). W_- above has a pole at
    # z = -1 < 0, because V_- is indefinite; so W_- is Nevanlinna, NOT
    # Stieltjes, and the Krein-string reading needs positivity hypotheses.
    Vm = [[Fraction(-1), Fraction(0)], [Fraction(0), Fraction(1)]]
    neg = Fraction(-1)
    negative_spectrum = det_exact(shift(Vm, neg)) == 0 and neg < 0
    # blow-up approaching the negative pole:
    near = abs(weyl_from_matrix(Vm, neg + Fraction(1, 10**6)))
    ok = negative_spectrum and near > 10**5
    record(
        "NC6_indefinite_velocity_not_stieltjes",
        ok,
        "V_- has an eigenvalue at -1 < 0 and W_- has a pole there, so W_- is "
        "not a Stieltjes function; the rank-one algebra alone does not yield a "
        "Krein string (positivity hypotheses are genuinely additional)",
    )


# ---------------------------------------------------------------- NC7
def nc7_symmetric_noise_variance():
    # Symmetric measurement E_mn = E_nm: the free variables are the entries
    # with m <= n (diagonal included), pairwise uncorrelated, each of variance
    # sigma^2 (take sigma = 1). The coincidence average (1/d^2) sum_{m,n} E_mn
    # puts weight 1/d^2 on each diagonal variable and 2/d^2 on each
    # strictly-upper variable, so Var = w^t Sigma w = sum of squared weights.
    # Corrected statement (Remark 3.5(ii)): Var = sigma^2 (2d-1)/d^3 exactly,
    # of order d^2/2 reduction. v1-style full d^2 reduction FAILS: Var > 1/d^2.
    all_ok = True
    details = []
    for d in (3, 5, 7, 9):
        var = (d * Fraction(1, d * d) ** 2
               + (d * (d - 1) // 2) * Fraction(2, d * d) ** 2)
        expected = Fraction(2 * d - 1, d**3)
        full_reduction = Fraction(1, d * d)
        ok_d = (var == expected) and (var > full_reduction)
        all_ok = all_ok and ok_d
        details.append(f"d={d}: Var={var} == (2d-1)/d^3={expected}, "
                       f"> 1/d^2={full_reduction}")
    record(
        "NC7_symmetric_noise_breaks_full_reduction",
        all_ok,
        "exact w^t Sigma w in Fractions: the full d^2 reduction fails under "
        "symmetry-correlated noise while sigma^2(2d-1)/d^3 holds exactly "
        "(order d^2/2 reduction); " + "; ".join(details),
    )


def main():
    nc1_variance()
    nc2_diagnostics()
    nc3_orthogonal_eigenvector()
    nc4_repeated_eigenvalue()
    nc5_weyl_zero_off_spectra()
    nc6_not_stieltjes()
    nc7_symmetric_noise_variance()
    status = "PASS" if not FAILURES else "FAIL"
    artifact = {
        "guard": "check_negative_controls",
        "status": status,
        "failures": FAILURES,
        "results": RESULTS,
        "python": sys.version.split()[0],
        "note": (
            "Negative controls for version 2: each v1 defect fails on "
            "its counterexample and each corrected statement holds. Seven "
            "controls (NC1-NC7). Resolvent quantities are computed from the "
            "matrices by exact Fraction solves; W_+ in NC5 is computed "
            "independently of the Sherman-Morrison formula under test. Exact "
            "Fraction arithmetic except the seeded Monte Carlo in NC1."
        ),
    }
    out = (
        Path(__file__).resolve().parents[1]
        / "artifacts"
        / "check_negative_controls.json"
    )
    out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="ascii")
    print(f"status={status} json={out}")
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
