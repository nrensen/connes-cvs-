#!/usr/bin/env python3
"""
================================================================================
CELL 76 — SCALE-INVARIANT STIELTJES-RESIDUE REPAIR AUDIT, NORMALIZED ZERO
         RESIDUALS, RECONSTRUCTION CONVERGENCE FLOOR, AND BOUNDARY WEIGHT
         RATIO SPECTRAL FACTORIZATION
================================================================================

Strategic Roadmap Milestone M22 (Paper NR2 Section 8.24 & Proposition 8.25):
----------------------------------------------------------------------------
Following the execution of cell 75, this cell executes the definitive scale-
invariant audit and repair of the Stieltjes-residue product formula for boundary
weights d_k^2 in the finite-rank Connes-van Suijlekom Galerkin truncation.

BACKGROUND & ROOT CAUSE ANALYSIS:
---------------------------------
In cell 75, the Stieltjes zero finder employed an absolute function-value stopping
criterion:
    if abs(fc) <= 1e-70: return c
For higher modes this was benign, but in the lowest interlacing interval
(E_0, E_1) at N=24, the natural physical scale of the regularized bracket function
    f_0(z) = (z - E_0)(E_1 - z) G_d(z)
is of order:
    |f_0(z)| ~ d_0^2 (E_1 - E_0) ~ 1.3e-40 * 1e-34 ~ 1e-74.
Because 1e-74 << 1e-70, the root finder terminated prematurely upon entering the
bracket, returning a spurious value for z_0^*. Because z_0^* enters every residue
product as a factor |E_k - z_0^*|, this error contaminated all reconstructed
boundary weights d_k^2.

THE FOUR AUDIT TESTS OF CELL 76:
--------------------------------
1. Test A — Scale-Invariant Pure Bracketed Bisection Root Solver:
   Replaces false-position/Pegasus with pure bracketed bisection testing strictly
   the algebraic sign of f_j(c). Uses exclusively a relative bracket-width stopping
   criterion:
       (b - a) / (E_{j+1} - E_j) <= 1e-55.
   Zero absolute function-value stopping thresholds are applied, ensuring absolute
   scale invariance across all intervals from 10^{-74} to 10^2.

2. Test B — Independent Scale-Normalized Residual Verification:
   Evaluates G_d(z_j^*) directly at the computed root and measures the scale-
   normalized residual:
       Res_norm(z_j^*) = |G_d(z_j^*)| / sum_{k=0}^N [ d_k^2 / |E_k - z_j^*| ].
   This provides an audit-proof, scale-free metric of root cancellation.

3. Test C — Exact Residue Reconstruction of Boundary Weights d_k^2:
   Re-runs the exact Stieltjes residue product formula:
       d_{k, prod}^2 = (2N + 1) * [ prod_{j=0}^{N-1} |E_k - z_j^*| ]
                                 / [ prod_{l != k} |E_k - E_l| ]
   and audits the relative error |d_{k, prod}^2 / d_k^2 - 1| down to the numerical
   precision floor across all modes k in {0, ..., N}.

4. Test D — Boundary Weight Growth Ratio alpha_j Spectral Factorization:
   Audits the consecutive ratio alpha_j = d_{j+1}^2 / d_j^2 both ways:
       (i)  Directly from eigenvectors: alpha_j^{eigen} = d_{j+1}^2 / d_j^2,
       (ii) From the Stieltjes product formula (Proposition 8.25):
            alpha_j^{prod} = zeta_j * (Pi_j^{zeros} * Pi_j^{evals}),
            where zeta_j = (E_{j+1} - z_j^*) / (z_j^* - E_j) is the isolated local zero
            ratio, Pi_j^{zeros} = prod_{l != j} |E_{j+1} - z_l^*| / |E_j - z_l^*|,
            and Pi_j^{evals} = prod_{l != j, j+1} |E_j - E_l| / |E_{j+1} - E_l|.
   Compares zeta_j with the odd-eigenvalue interlacing ratio R_j / L_j and confirms
   the bridge to the exact three-factor pole asymmetry identity:
       H_{j+1}(mu_j) / H_j(mu_j) = alpha_j * (L_j / R_j)^2.

OUTPUT CONSTRAINTS:
-------------------
- Dispassionate, dry, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 76 EXECUTION COMPLETE.
================================================================================
"""

from __future__ import annotations

import time
import mpmath as mp

from connes_cvs import build_galerkin_matrix
from cell import get_galerkin_matrix

# -----------------------------------------------------------------------------
# Precision and Parameter Configuration
# -----------------------------------------------------------------------------

# Set working precision to 70 decimal digits for safe 50-digit verification
mp.mp.dps = 70

C_PARAM = 13
L_PARAM = mp.log(C_PARAM)
T_PARAM = 400
GROUND_DPS = 50

N_LIST = [8, 12, 16, 20, 24]
BARRIER_V_STAR = mp.mpf('1.0')

# Relative bracket tolerance for pure bisection root solver
BISECTION_REL_TOL = mp.mpf('1e-55')
BISECTION_MAX_ITERS = 250


# -----------------------------------------------------------------------------
# Parity Basis & Operator Construction
# -----------------------------------------------------------------------------

def full_parity_basis(N: int) -> tuple[mp.matrix, mp.matrix]:
    """Construct orthonormal basis matrices E (even) and O (odd) for R^{2N+1}."""
    dim = 2 * N + 1
    E = mp.matrix(dim, N + 1)
    O = mp.matrix(dim, N)

    E[N, 0] = mp.mpf(1)

    inv_sqrt2 = 1 / mp.sqrt(2)
    for m in range(1, N + 1):
        E[N + m, m] = inv_sqrt2
        E[N - m, m] = inv_sqrt2
        O[N + m, m - 1] = inv_sqrt2
        O[N - m, m - 1] = -inv_sqrt2

    return E, O


def solve_parity_eigensystems(
    Q_full: mp.matrix, N: int
) -> tuple[
    mp.mpf, mp.matrix, mp.matrix,
    list[mp.mpf], mp.matrix,
    list[mp.mpf], mp.matrix
]:
    """Compute even and odd spectra and eigenvectors."""
    E, O = full_parity_basis(N)

    Q_even = E.T * Q_full * E
    Q_even = mp.mpf('0.5') * (Q_even + Q_even.T)

    Q_odd = O.T * Q_full * O
    Q_odd = mp.mpf('0.5') * (Q_odd + Q_odd.T)

    evals_e, V_e = mp.eigsy(Q_even)
    evals_o, V_o = mp.eigsy(Q_odd)

    # Sort even
    idx_e = sorted(range(N + 1), key=lambda i: evals_e[i])
    sorted_evals_e = [evals_e[i] for i in idx_e]
    sorted_V_e = mp.matrix(N + 1, N + 1)
    for col_idx, orig_col in enumerate(idx_e):
        for row_idx in range(N + 1):
            sorted_V_e[row_idx, col_idx] = V_e[row_idx, orig_col]

    # Sort odd
    idx_o = sorted(range(N), key=lambda i: evals_o[i])
    sorted_evals_o = [evals_o[i] for i in idx_o]
    sorted_V_o = mp.matrix(N, N)
    for col_idx, orig_col in enumerate(idx_o):
        for row_idx in range(N):
            sorted_V_o[row_idx, col_idx] = V_o[row_idx, orig_col]

    lam_0 = sorted_evals_e[0]
    return lam_0, E, O, sorted_evals_e, sorted_V_e, sorted_evals_o, sorted_V_o


# -----------------------------------------------------------------------------
# Regularized Stieltjes Function & Scale-Invariant Root Finder
# -----------------------------------------------------------------------------

def regularized_stieltjes_bracket_func(
    z: mp.mpf, j: int, evals_e: list[mp.mpf], d_k_sq: list[mp.mpf]
) -> mp.mpf:
    """
    Evaluates f_j(z) = (z - E_j)(E_{j+1} - z) * G_d(z) on [E_j, E_{j+1}].
    This function is regular, continuous, strictly monotonic, with:
        f_j(E_j)     = - d_j^2 * (E_{j+1} - E_j) < 0
        f_j(E_{j+1}) = + d_{j+1}^2 * (E_{j+1} - E_j) > 0.
    """
    E_j = evals_e[j]
    E_jp1 = evals_e[j + 1]

    # Leading two terms:
    # (z - E_j)(E_{j+1} - z) * [ d_j^2 / (E_j - z) + d_{j+1}^2 / (E_{j+1} - z) ]
    # = - d_j^2 * (E_{j+1} - z) + d_{j+1}^2 * (z - E_j)
    term_local = - d_k_sq[j] * (E_jp1 - z) + d_k_sq[j + 1] * (z - E_j)

    # Remaining terms:
    factor_outer = (z - E_j) * (E_jp1 - z)
    term_outer = mp.mpf(0)
    for k in range(len(evals_e)):
        if k != j and k != j + 1:
            term_outer += d_k_sq[k] / (evals_e[k] - z)

    return term_local + factor_outer * term_outer


def find_stieltjes_zero_bisection(
    j: int,
    evals_e: list[mp.mpf],
    d_k_sq: list[mp.mpf],
    rel_tol: mp.mpf = BISECTION_REL_TOL,
    max_iter: int = BISECTION_MAX_ITERS,
) -> tuple[mp.mpf, int]:
    """
    Test A: Finds the unique zero z_j^* in (E_j, E_{j+1}) of G_d(z) using
    pure bracketed bisection on the regularized bracket function f_j(z).

    The stopping criterion is strictly based on relative bracket width:
        (b - a) / (E_{j+1} - E_j) <= rel_tol.
    No absolute function-value threshold is applied, ensuring complete scale
    invariance regardless of interval magnitude.
    """
    E_j = evals_e[j]
    E_jp1 = evals_e[j + 1]
    interval_len = E_jp1 - E_j

    a = E_j
    b = E_jp1

    iters = 0
    while (b - a) / interval_len > rel_tol and iters < max_iter:
        iters += 1
        c = mp.mpf('0.5') * (a + b)
        fc = regularized_stieltjes_bracket_func(c, j, evals_e, d_k_sq)
        if fc < 0:
            a = c
        elif fc > 0:
            b = c
        else:
            return c, iters

    c_final = mp.mpf('0.5') * (a + b)
    return c_final, iters


def evaluate_scale_normalized_residual(
    z_val: mp.mpf, evals_e: list[mp.mpf], d_k_sq: list[mp.mpf]
) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    """
    Test B: Computes direct G_d(z) and scale-normalized residual:
        Res_norm(z) = |G_d(z)| / sum_{k=0}^N [ d_k^2 / |E_k - z| ].
    """
    gd_num = mp.mpf(0)
    denom_scale = mp.mpf(0)
    for k in range(len(evals_e)):
        diff = evals_e[k] - z_val
        abs_diff = abs(diff)
        gd_num += d_k_sq[k] / diff
        denom_scale += d_k_sq[k] / abs_diff

    res_norm = abs(gd_num) / denom_scale if denom_scale > 0 else mp.mpf(0)
    return gd_num, denom_scale, res_norm


# -----------------------------------------------------------------------------
# Residue Reconstruction & Spectral Factorization
# -----------------------------------------------------------------------------

def reconstruct_boundary_weights_product(
    evals_e: list[mp.mpf], zeros_z: list[mp.mpf], N: int
) -> list[mp.mpf]:
    """
    Test C: Evaluates the exact residue product formula for each d_k^2:
        d_{k, prod}^2 = (2N + 1) * [ prod_{j=0}^{N-1} |E_k - z_j^*| ]
                                  / [ prod_{l != k} |E_k - E_l| ].
    """
    d_k_sq_prod: list[mp.mpf] = []
    norm_factor = mp.mpf(2 * N + 1)

    for k in range(N + 1):
        E_k = evals_e[k]

        # Numerator: prod_{j=0}^{N-1} |E_k - z_j^*|
        prod_num = mp.mpf(1)
        for j in range(N):
            prod_num *= abs(E_k - zeros_z[j])

        # Denominator: prod_{l != k} |E_k - E_l|
        prod_den = mp.mpf(1)
        for l in range(N + 1):
            if l != k:
                prod_den *= abs(E_k - evals_e[l])

        val_prod = norm_factor * (prod_num / prod_den)
        d_k_sq_prod.append(val_prod)

    return d_k_sq_prod


def factorize_boundary_ratio(
    j: int,
    evals_e: list[mp.mpf],
    zeros_z: list[mp.mpf],
    d_k_sq: list[mp.mpf],
    evals_o: list[mp.mpf],
) -> dict:
    """
    Test D: Evaluates consecutive boundary weight ratio alpha_j = d_{j+1}^2 / d_j^2
    both directly and via the Stieltjes zero product factorization:
        alpha_j^{prod} = zeta_j * (Pi_j^{zeros} * Pi_j^{evals}),
    where:
        zeta_j        = (E_{j+1} - z_j^*) / (z_j^* - E_j),
        Pi_j^{zeros}  = prod_{l != j} |E_{j+1} - z_l^*| / |E_j - z_l^*|,
        Pi_j^{evals}  = prod_{l != j, j+1} |E_j - E_l| / |E_{j+1} - E_l|.
    """
    N = len(zeros_z)
    E_j = evals_e[j]
    E_jp1 = evals_e[j + 1]
    mu_j = evals_o[j]

    # Direct eigenvector ratio
    alpha_eigen = d_k_sq[j + 1] / d_k_sq[j]

    # Local zero factor
    z_j = zeros_z[j]
    zeta_j = (E_jp1 - z_j) / (z_j - E_j)

    # Outer zero product: prod_{l != j} |E_{j+1} - z_l^*| / |E_j - z_l^*|
    pi_zeros = mp.mpf(1)
    for l in range(N):
        if l != j:
            pi_zeros *= abs(E_jp1 - zeros_z[l]) / abs(E_j - zeros_z[l])

    # Outer eigenvalue product: prod_{l != j, j+1} |E_j - E_l| / |E_{j+1} - E_l|
    pi_evals = mp.mpf(1)
    for l in range(N + 1):
        if l != j and l != j + 1:
            pi_evals *= abs(E_j - evals_e[l]) / abs(E_jp1 - evals_e[l])

    pi_outer = pi_zeros * pi_evals
    alpha_prod = zeta_j * pi_outer
    rel_err_alpha = abs(alpha_prod / alpha_eigen - 1)

    # Interlacing gap ratio for odd eigenvalue mu_j
    L_j = mu_j - E_j
    R_j = E_jp1 - mu_j
    ratio_RL = R_j / L_j
    inv_ratio_RL_sq = (L_j / R_j) ** 2

    # Comparison zeta_j vs R_j / L_j
    zeta_vs_gap_ratio = zeta_j / ratio_RL

    return {
        'j': j,
        'alpha_eigen': alpha_eigen,
        'alpha_prod': alpha_prod,
        'rel_err_alpha': rel_err_alpha,
        'zeta_j': zeta_j,
        'ratio_RL': ratio_RL,
        'inv_ratio_RL_sq': inv_ratio_RL_sq,
        'pi_zeros': pi_zeros,
        'pi_evals': pi_evals,
        'pi_outer': pi_outer,
        'zeta_vs_gap_ratio': zeta_vs_gap_ratio,
        'L_j': L_j,
        'R_j': R_j,
    }


def compute_H_pole_decomposition(
    mu_val: mp.mpf, lam_0: mp.mpf, evals_e: list[mp.mpf], d_k_sq: list[mp.mpf]
) -> tuple[mp.mpf, list[mp.mpf]]:
    """Compute H(mu) and its individual pole contributions H_k(mu)."""
    N = len(evals_e) - 1
    d0_sq = d_k_sq[0]

    H_terms = [d0_sq]  # H_0 = D_0^2
    H_total = d0_sq

    shift = mu_val - lam_0
    shift_sq = shift ** 2

    for k in range(1, N + 1):
        gap = evals_e[k] - mu_val
        term_k = d_k_sq[k] * shift_sq / (gap ** 2)
        H_terms.append(term_k)
        H_total += term_k

    return H_total, H_terms


# -----------------------------------------------------------------------------
# Main Execution Function
# -----------------------------------------------------------------------------

def run_cell76() -> None:
    print("=" * 105)
    print("CELL 76 — SCALE-INVARIANT STIELTJES-RESIDUE REPAIR AUDIT, NORMALIZED ZERO RESIDUALS,")
    print("         RECONSTRUCTION CONVERGENCE FLOOR, AND BOUNDARY WEIGHT RATIO SPECTRAL FACTORIZATION")
    print("=" * 105)
    print(f"Parameters: c = {C_PARAM}, L = {float(L_PARAM):.14f}, T = {T_PARAM}, dps = {mp.mp.dps}, V_* = {float(BARRIER_V_STAR)}")
    print(f"Root Finder: Pure bracketed bisection, stopping tol = {mp.nstr(BISECTION_REL_TOL, 4)} (bracket width / gap)")
    print()

    synthesis_table1_records: list[dict] = []
    synthesis_table2_records: list[dict] = []
    synthesis_table3_records: list[dict] = []
    synthesis_table4_records: list[dict] = []

    for N in N_LIST:
        t0 = time.time()

        # Retrieve cached Galerkin matrix
        Q_full, _ = get_galerkin_matrix(
            c=C_PARAM,
            N=N,
            T=T_PARAM,
            dps=GROUND_DPS,
            verbose=False,
        )

        # Solve parity eigensystems
        lam_0, E, O, evals_e, V_e, evals_o, V_o = solve_parity_eigensystems(Q_full, N)

        # Boundary vector d in R^{N+1}: d_even = (1, sqrt(2), ..., sqrt(2))^T
        d_even = mp.matrix(N + 1, 1)
        d_even[0, 0] = mp.mpf(1)
        for m in range(1, N + 1):
            d_even[m, 0] = mp.sqrt(2)

        # Boundary overlaps in even sector: d_k = <u_k^{even}, d>
        d_k_list: list[mp.mpf] = []
        d_k_sq_list: list[mp.mpf] = []
        for k in range(N + 1):
            d_val = sum(d_even[m, 0] * V_e[m, k] for m in range(N + 1))
            d_k_list.append(d_val)
            d_k_sq_list.append(d_val ** 2)

        norm_d_sq = sum(d_k_sq_list)
        norm_expected = mp.mpf(2 * N + 1)
        norm_err = abs(norm_d_sq - norm_expected)

        print(f"--- DIMENSION N = {N:2d} (Elapsed: {time.time() - t0:.2f}s) ---")
        print(f"  Total Norm ||d||^2:  Exact = {int(norm_expected)}, Computed = {mp.nstr(norm_d_sq, 12)}, Residual = {mp.nstr(norm_err, 4)}")
        print(f"  Even Ground State E_0: {mp.nstr(evals_e[0], 12)}, Odd Ground State mu_0: {mp.nstr(evals_o[0], 12)}")

        # ---------------------------------------------------------------------
        # Tests A & B: Pure Bisection Root Finding and Scale-Normalized Residuals
        # ---------------------------------------------------------------------
        zeros_z: list[mp.mpf] = []
        bisection_iters: list[int] = []
        res_norm_list: list[mp.mpf] = []
        disp_norm_list: list[mp.mpf] = []

        print(f"  [Tests A & B] Pure Bisection Stieltjes Zeros z_j^* on (E_j, E_{{j+1}}):")
        print(f"    {'j':>2s} | {'E_j':>14s} | {'z_j^*':>14s} | {'E_{j+1}':>14s} | {'Iters':>5s} | {'Res_norm(z_j^*)':>15s} | {'|z_j^* - mu_j|/Gap':>18s}")
        print("    " + "-" * 95)

        for j in range(N):
            E_j = evals_e[j]
            E_jp1 = evals_e[j + 1]
            gap_j = E_jp1 - E_j
            mu_j = evals_o[j]

            # Test A: Pure bisection root finding
            z_j, iters = find_stieltjes_zero_bisection(j, evals_e, d_k_sq_list)
            zeros_z.append(z_j)
            bisection_iters.append(iters)

            # Test B: Scale-normalized residual
            _, _, res_norm = evaluate_scale_normalized_residual(z_j, evals_e, d_k_sq_list)
            res_norm_list.append(res_norm)

            # Alignment with odd eigenvalue mu_j
            disp_norm = abs(z_j - mu_j) / gap_j
            disp_norm_list.append(disp_norm)

            if j < 4 or j == N - 1:
                print(f"    {j:2d} | {mp.nstr(E_j, 8):>14s} | {mp.nstr(z_j, 8):>14s} | {mp.nstr(E_jp1, 8):>14s} | {iters:5d} | {mp.nstr(res_norm, 6):>15s} | {mp.nstr(disp_norm, 6):>18s}")
            elif j == 4:
                print(f"    {'...':>2s} | {'...':>14s} | {'...':>14s} | {'...':>14s} | {'...':>5s} | {'...':>15s} | {'...':>18s}")

        max_res_norm = max(res_norm_list)
        max_iters = max(bisection_iters)
        print(f"    Max Bisection Iterations: {max_iters}, Max Scale-Normalized Residual: {mp.nstr(max_res_norm, 4)}")

        # ---------------------------------------------------------------------
        # Test C: Exact Residue Reconstruction of Boundary Weights d_k^2
        # ---------------------------------------------------------------------
        d_k_sq_prod = reconstruct_boundary_weights_product(evals_e, zeros_z, N)
        rel_errors_d_sq = [abs(d_k_sq_prod[k] / d_k_sq_list[k] - 1) for k in range(N + 1)]
        max_rel_err_d_sq = max(rel_errors_d_sq)

        print(f"  [Test C] Stieltjes Residue Product Reconstruction of Boundary Weights d_k^2:")
        print(f"    {'k':>2s} | {'d_k^2 (Eigen)':>15s} | {'d_{k, prod}^2 (Residue)':>24s} | {'Rel Error':>15s}")
        print("    " + "-" * 62)

        for k in range(min(5, N + 1)):
            print(f"    {k:2d} | {mp.nstr(d_k_sq_list[k], 8):>15s} | {mp.nstr(d_k_sq_prod[k], 8):>24s} | {mp.nstr(rel_errors_d_sq[k], 6):>15s}")

        if N >= 5:
            print(f"    {'...':>2s} | {'...':>15s} | {'...':>24s} | {'...':>15s}")
            last_k = N
            print(f"    {last_k:2d} | {mp.nstr(d_k_sq_list[last_k], 8):>15s} | {mp.nstr(d_k_sq_prod[last_k], 8):>24s} | {mp.nstr(rel_errors_d_sq[last_k], 6):>15s}")

        print(f"    Max Relative Reconstruction Error across all k in {{0, ..., {N}}}: {mp.nstr(max_rel_err_d_sq, 4)}")

        # ---------------------------------------------------------------------
        # Test D: Boundary Weight Growth Ratio alpha_j Spectral Factorization
        # ---------------------------------------------------------------------
        print(f"  [Test D] Consecutive Weight Ratio alpha_j = d_{{j+1}}^2 / d_j^2 Factorization:")
        print(f"    {'j':>2s} | {'alpha_j (Eigen)':>15s} | {'alpha_j (Prod)':>15s} | {'Rel Error':>12s} | {'zeta_j':>14s} | {'R_j / L_j':>14s} | {'zeta/(R/L)':>12s}")
        print("    " + "-" * 95)

        fact_records_this_N: list[dict] = []
        for j in range(min(4, N)):
            fact = factorize_boundary_ratio(j, evals_e, zeros_z, d_k_sq_list, evals_o)
            fact_records_this_N.append(fact)
            print(f"    {j:2d} | {mp.nstr(fact['alpha_eigen'], 7):>15s} | {mp.nstr(fact['alpha_prod'], 7):>15s} | {mp.nstr(fact['rel_err_alpha'], 4):>12s} | {mp.nstr(fact['zeta_j'], 7):>14s} | {mp.nstr(fact['ratio_RL'], 7):>14s} | {mp.nstr(fact['zeta_vs_gap_ratio'], 6):>12s}")

        # Evaluate pole asymmetry identity for mode j = 2
        fact2 = fact_records_this_N[2]
        H_total_2, H_terms_2 = compute_H_pole_decomposition(evals_o[2], lam_0, evals_e, d_k_sq_list)
        H_left_2 = H_terms_2[2]   # pole E_2
        H_right_2 = H_terms_2[3]  # pole E_3
        ratio_H_right_left = H_right_2 / H_left_2
        pred_asym = fact2['alpha_eigen'] * fact2['inv_ratio_RL_sq']
        asym_res = abs(ratio_H_right_left - pred_asym)
        fidelity_2 = (H_left_2 + H_right_2) / H_total_2 * 100

        print(f"  [Pole Asymmetry Bridge for j = 2]:")
        print(f"    alpha_2 = {mp.nstr(fact2['alpha_eigen'], 8)},  (L_2 / R_2)^2 = {mp.nstr(fact2['inv_ratio_RL_sq'], 8)}")
        print(f"    Predicted Asymmetry alpha_2 * (L_2 / R_2)^2 = {mp.nstr(pred_asym, 8)}")
        print(f"    Actual H_3(mu_2) / H_2(mu_2)               = {mp.nstr(ratio_H_right_left, 8)}")
        print(f"    Asymmetry Identity Residual: {mp.nstr(asym_res, 4)}")
        print(f"    Two-Pole Cluster Fidelity:   {float(fidelity_2):.4f}%")
        print()

        # Store synthesis records
        synthesis_table1_records.append({
            'N': N,
            'max_iters': max_iters,
            'max_res_norm': max_res_norm,
            'res_norm_0': res_norm_list[0],
            'res_norm_2': res_norm_list[2],
            'disp_norm_2': disp_norm_list[2],
            'disp_norm_3': disp_norm_list[3] if N > 3 else mp.mpf(0),
        })

        synthesis_table2_records.append({
            'N': N,
            'rel_err_d0': rel_errors_d_sq[0],
            'rel_err_d1': rel_errors_d_sq[1],
            'rel_err_d2': rel_errors_d_sq[2],
            'rel_err_d3': rel_errors_d_sq[3],
            'max_rel_err': max_rel_err_d_sq,
            'norm_err': norm_err,
        })

        synthesis_table3_records.append({
            'N': N,
            'alpha_eigen': fact2['alpha_eigen'],
            'alpha_prod': fact2['alpha_prod'],
            'rel_err_alpha': fact2['rel_err_alpha'],
            'zeta_2': fact2['zeta_j'],
            'ratio_RL': fact2['ratio_RL'],
            'pi_outer': fact2['pi_outer'],
            'zeta_vs_gap_ratio': fact2['zeta_vs_gap_ratio'],
        })

        synthesis_table4_records.append({
            'N': N,
            'alpha_2': fact2['alpha_eigen'],
            'inv_ratio_RL_sq': fact2['inv_ratio_RL_sq'],
            'pred_asym': pred_asym,
            'actual_asym': ratio_H_right_left,
            'asym_res': asym_res,
            'fidelity': fidelity_2,
        })

    # -------------------------------------------------------------------------
    # Multi-Dimension Synthesis Tables
    # -------------------------------------------------------------------------
    print("=" * 105)
    print("SYNTHESIS TABLE 1: SCALE-INVARIANT ZERO FINDER & NORMALIZED RESIDUALS (TESTS A & B)")
    print("=" * 105)
    print(f"{'N':>3s} | {'MaxIters':>8s} | {'Max Res_norm':>14s} | {'Res_norm(z_0^*)':>16s} | {'Res_norm(z_2^*)':>16s} | {'|z_2^*-mu_2|/Gap':>18s} | {'|z_3^*-mu_3|/Gap':>18s}")
    print("-" * 105)
    for rec in synthesis_table1_records:
        print(f"{rec['N']:3d} | {rec['max_iters']:8d} | {mp.nstr(rec['max_res_norm'], 4):>14s} | {mp.nstr(rec['res_norm_0'], 4):>16s} | {mp.nstr(rec['res_norm_2'], 4):>16s} | {mp.nstr(rec['disp_norm_2'], 6):>18s} | {mp.nstr(rec['disp_norm_3'], 6):>18s}")
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 2: EXACT RESIDUE RECONSTRUCTION OF BOUNDARY WEIGHTS d_k^2 (TEST C)")
    print("=" * 105)
    print(f"{'N':>3s} | {'RelErr(d_0^2)':>15s} | {'RelErr(d_1^2)':>15s} | {'RelErr(d_2^2)':>15s} | {'RelErr(d_3^2)':>15s} | {'MaxRelErr (All)':>16s} | {'Norm Error':>15s}")
    print("-" * 105)
    for rec in synthesis_table2_records:
        print(f"{rec['N']:3d} | {mp.nstr(rec['rel_err_d0'], 4):>15s} | {mp.nstr(rec['rel_err_d1'], 4):>15s} | {mp.nstr(rec['rel_err_d2'], 4):>15s} | {mp.nstr(rec['rel_err_d3'], 4):>15s} | {mp.nstr(rec['max_rel_err'], 4):>16s} | {mp.nstr(rec['norm_err'], 4):>15s}")
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 3: BOUNDARY WEIGHT RATIO alpha_2 FACTORIZATION & SPECTRAL GAP STRUCTURE (TEST D)")
    print("=" * 105)
    print(f"{'N':>3s} | {'alpha_2 (Eigen)':>15s} | {'alpha_2 (Prod)':>15s} | {'Rel Error':>12s} | {'zeta_2':>14s} | {'R_2 / L_2':>14s} | {'Pi_2^{outer}':>14s} | {'zeta/(R/L)':>12s}")
    print("-" * 105)
    for rec in synthesis_table3_records:
        print(f"{rec['N']:3d} | {mp.nstr(rec['alpha_eigen'], 6):>15s} | {mp.nstr(rec['alpha_prod'], 6):>15s} | {mp.nstr(rec['rel_err_alpha'], 4):>12s} | {mp.nstr(rec['zeta_2'], 6):>14s} | {mp.nstr(rec['ratio_RL'], 6):>14s} | {mp.nstr(rec['pi_outer'], 6):>14s} | {mp.nstr(rec['zeta_vs_gap_ratio'], 6):>12s}")
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 4: BRIDGE TO POLE ASYMMETRY BALANCE FOR MODE j = 2 (TEST D)")
    print("=" * 105)
    print(f"{'N':>3s} | {'alpha_2':>14s} | {'(L_2 / R_2)^2':>15s} | {'alpha_2 * (L/R)^2':>18s} | {'H_3 / H_2':>12s} | {'Asymmetry Res':>15s} | {'Fidelity F(%)':>14s}")
    print("-" * 105)
    for rec in synthesis_table4_records:
        print(f"{rec['N']:3d} | {mp.nstr(rec['alpha_2'], 6):>14s} | {mp.nstr(rec['inv_ratio_RL_sq'], 6):>15s} | {mp.nstr(rec['pred_asym'], 6):>18s} | {mp.nstr(rec['actual_asym'], 6):>12s} | {mp.nstr(rec['asym_res'], 4):>15s} | {float(rec['fidelity']):.4f}%")
    print()

    print("=" * 80)
    print("CELL 76 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell76()
