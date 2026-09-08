#!/usr/bin/env python3
"""
================================================================================
CELL 77 — STIELTJES OUTER PRODUCT FACTOR Pi_j, PAIRWISE GAP LOCALIZATION,
         CHARACTERISTIC-POLYNOMIAL RATIO DECOMPOSITION, AND ASYMMETRY BALANCE
================================================================================

Strategic Roadmap Milestone M23 (Paper 4B Section 8.25):
---------------------------------------------------------
Following the certified completion of Milestone M22 in cell 76, this cell
investigates the analytical structure, mode-by-mode localization, and
characteristic-polynomial representation of the Stieltjes outer product factor:
    Pi_j = prod_{l != j} |E_{j+1} - z_l^*| / |E_j - z_l^*|
         * prod_{l != j, j+1} |E_j - E_l| / |E_{j+1} - E_l|.

BACKGROUND & ANALYTICAL MOTIVATION:
-----------------------------------
In cell 76, the boundary weight growth ratio alpha_j = d_{j+1}^2 / d_j^2 was
certified to factorize into:
    alpha_j = zeta_j * Pi_j,
where zeta_j = (E_{j+1} - z_j^*) / (z_j^* - E_j) is the local zero coordinate ratio.
Remarkably, at N=24, j=2:
    alpha_2 = 92102.75,   zeta_2 = 386869.2,   Pi_2 = 0.238072 = O(1).
While alpha_2 and zeta_2 grow into massive numbers, the outer factor Pi_2
remains a benign O(1) constant.

THE FOUR INVESTIGATIVE TESTS OF CELL 77:
----------------------------------------
1. Test A — Multi-Modal Audit of Pi_j Across Bound Modes j in {0, 1, 2, 3}:
   Evaluates Pi_j across tested dimensions N in {8, 12, 16, 20, 24} for all low
   bound modes, auditing whether Pi_j is uniformly bounded above and below:
       0 < c_min <= Pi_j(N) <= C_max < infinity.

2. Test B — Pairwise Mode Cancellation & Spatial Gap Localization:
   For each mode l != j, j+1, defines the pairwise cancellation factor:
       omega_{j, l} = [ |E_{j+1} - z_l^*| / |E_j - z_l^*| ]
                    * [ |E_j - E_l| / |E_{j+1} - E_l| ].
   Because z_l^* in (E_l, E_{l+1}) tracks E_l, as |E_j - E_l| >> Delta E_j,
   omega_{j, l} -> 1 rapidly. Audits |omega_{j, l} - 1| across l to verify that
   Pi_j is strictly localized to nearest neighbor modes.

3. Test C — Characteristic-Polynomial & Zero-Polynomial Ratio Identity:
   Audits the exact polynomial representations:
       P_{even}'(E_j) = prod_{l != j} (E_j - E_l),
       P_{zero}(E_j)  = prod_{m=0}^{N-1} (E_j - z_m^*).
   Proves and verifies that:
       alpha_j = |P_{zero}(E_{j+1}) / P_{zero}(E_j)| * |P_{even}'(E_j) / P_{even}'(E_{j+1})|,
       Pi_j    = (1 / zeta_j) * |P_{zero}(E_{j+1}) / P_{zero}(E_j)| * |P_{even}'(E_j) / P_{even}'(E_{j+1})|,
   completely eliminating the local coordinate from the total ratio.

4. Test D — Decomposition of the Master Asymmetry Balance:
   Decomposes the three-factor pole asymmetry ratio:
       H_{j+1}(mu_j) / H_j(mu_j) = zeta_j * Pi_j * (L_j / R_j)^2
                                  = [ zeta_j / (R_j / L_j) ] * Pi_j * (L_j / R_j).
   Evaluates each factor to reveal how the O(1) balance (1.528 at N=24) is achieved.

OUTPUT CONSTRAINTS:
-------------------
- Dispassionate, dry, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 77 EXECUTION COMPLETE.
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

# Working precision of 70 decimal digits
mp.mp.dps = 70

C_PARAM = 13
L_PARAM = mp.log(C_PARAM)
T_PARAM = 400
GROUND_DPS = 50

N_LIST = [8, 12, 16, 20, 24]
BARRIER_V_STAR = mp.mpf('1.0')

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

    idx_e = sorted(range(N + 1), key=lambda i: evals_e[i])
    sorted_evals_e = [evals_e[i] for i in idx_e]
    sorted_V_e = mp.matrix(N + 1, N + 1)
    for col_idx, orig_col in enumerate(idx_e):
        for row_idx in range(N + 1):
            sorted_V_e[row_idx, col_idx] = V_e[row_idx, orig_col]

    idx_o = sorted(range(N), key=lambda i: evals_o[i])
    sorted_evals_o = [evals_o[i] for i in idx_o]
    sorted_V_o = mp.matrix(N, N)
    for col_idx, orig_col in enumerate(idx_o):
        for row_idx in range(N):
            sorted_V_o[row_idx, col_idx] = V_o[row_idx, orig_col]

    lam_0 = sorted_evals_e[0]
    return lam_0, E, O, sorted_evals_e, sorted_V_e, sorted_evals_o, sorted_V_o


# -----------------------------------------------------------------------------
# Scale-Invariant Root Finder
# -----------------------------------------------------------------------------

def regularized_stieltjes_bracket_func(
    z: mp.mpf, j: int, evals_e: list[mp.mpf], d_k_sq: list[mp.mpf]
) -> mp.mpf:
    """Evaluates f_j(z) = (z - E_j)(E_{j+1} - z) * G_d(z) on [E_j, E_{j+1}]."""
    E_j = evals_e[j]
    E_jp1 = evals_e[j + 1]

    term_local = - d_k_sq[j] * (E_jp1 - z) + d_k_sq[j + 1] * (z - E_j)
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
    """Finds zero z_j^* using pure bracketed bisection."""
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

    return mp.mpf('0.5') * (a + b), iters


# -----------------------------------------------------------------------------
# Outer Factor Pi_j, Pairwise Modes, and Polynomial Ratios
# -----------------------------------------------------------------------------

def evaluate_outer_product_detailed(
    j: int,
    evals_e: list[mp.mpf],
    zeros_z: list[mp.mpf],
    d_k_sq: list[mp.mpf],
    evals_o: list[mp.mpf],
) -> dict:
    """
    Computes Pi_j, its pairwise cancellation decomposition omega_{j, l},
    and its characteristic-polynomial ratio representation.
    """
    N = len(zeros_z)
    E_j = evals_e[j]
    E_jp1 = evals_e[j + 1]
    mu_j = evals_o[j]

    # Local zero factor
    z_j = zeros_z[j]
    zeta_j = (E_jp1 - z_j) / (z_j - E_j)

    # Direct eigenvector ratio
    alpha_eigen = d_k_sq[j + 1] / d_k_sq[j]

    # 1. Zero product factor: prod_{l != j} |E_{j+1} - z_l^*| / |E_j - z_l^*|
    pi_zeros = mp.mpf(1)
    for l in range(N):
        if l != j:
            pi_zeros *= abs(E_jp1 - zeros_z[l]) / abs(E_j - zeros_z[l])

    # 2. Eigenvalue product factor: prod_{l != j, j+1} |E_j - E_l| / |E_{j+1} - E_l|
    pi_evals = mp.mpf(1)
    for l in range(N + 1):
        if l != j and l != j + 1:
            pi_evals *= abs(E_j - evals_e[l]) / abs(E_jp1 - evals_e[l])

    pi_outer = pi_zeros * pi_evals
    alpha_prod = zeta_j * pi_outer

    # 3. Pairwise mode factors omega_{j, l} for l != j, j+1
    # For l in {0, ..., N-1} (except j, j+1):
    #   omega_{j, l} = (|E_{j+1} - z_l^*| / |E_j - z_l^*|) * (|E_j - E_l| / |E_{j+1} - E_l|)
    # Note: for l = N, there is only an eigenvalue factor: |E_j - E_N| / |E_{j+1} - E_N|
    # And for l = j+1, there is only a zero factor: |E_{j+1} - z_{j+1}^*| / |E_j - z_{j+1}^*|
    pairwise_omega: list[dict] = []
    for l in range(N):
        if l != j and l != j + 1:
            ratio_zero = abs(E_jp1 - zeros_z[l]) / abs(E_j - zeros_z[l])
            ratio_eval = abs(E_j - evals_e[l]) / abs(E_jp1 - evals_e[l])
            omega_l = ratio_zero * ratio_eval
            deviation = abs(omega_l - mp.mpf(1))
            pairwise_omega.append({
                'l': l,
                'omega_l': omega_l,
                'deviation': deviation,
                'dist_index': l - j,
            })

    # 4. Characteristic-Polynomial Derivative and Value Ratios
    # P_{even}'(E_k) = prod_{l != k} (E_k - E_l)
    p_even_prime_Ej = mp.mpf(1)
    for l in range(N + 1):
        if l != j:
            p_even_prime_Ej *= abs(E_j - evals_e[l])

    p_even_prime_Ejp1 = mp.mpf(1)
    for l in range(N + 1):
        if l != j + 1:
            p_even_prime_Ejp1 *= abs(E_jp1 - evals_e[l])

    ratio_p_even_prime = p_even_prime_Ej / p_even_prime_Ejp1

    # P_{zero}(E_k) = prod_{m=0}^{N-1} (E_k - z_m^*)
    p_zero_Ej = mp.mpf(1)
    for m in range(N):
        p_zero_Ej *= abs(E_j - zeros_z[m])

    p_zero_Ejp1 = mp.mpf(1)
    for m in range(N):
        p_zero_Ejp1 *= abs(E_jp1 - zeros_z[m])

    ratio_p_zero = p_zero_Ejp1 / p_zero_Ej

    alpha_poly = ratio_p_zero * ratio_p_even_prime
    pi_poly = (mp.mpf(1) / zeta_j) * alpha_poly

    # 5. Master Asymmetry Balance
    L_j = mu_j - E_j
    R_j = E_jp1 - mu_j
    ratio_RL = R_j / L_j
    inv_ratio_RL = L_j / R_j
    inv_ratio_RL_sq = inv_ratio_RL ** 2

    # Ratio coordinate vs gap ratio
    zeta_over_gap = zeta_j / ratio_RL

    # Balanced factor: [zeta_j / (R_j / L_j)] * (L_j / R_j) * Pi_j = zeta_j * Pi_j * (L_j/R_j)^2
    pred_asym = zeta_j * pi_outer * inv_ratio_RL_sq

    return {
        'j': j,
        'N': N,
        'alpha_eigen': alpha_eigen,
        'alpha_prod': alpha_prod,
        'alpha_poly': alpha_poly,
        'zeta_j': zeta_j,
        'pi_zeros': pi_zeros,
        'pi_evals': pi_evals,
        'pi_outer': pi_outer,
        'pi_poly': pi_poly,
        'ratio_p_even_prime': ratio_p_even_prime,
        'ratio_p_zero': ratio_p_zero,
        'pairwise_omega': pairwise_omega,
        'L_j': L_j,
        'R_j': R_j,
        'ratio_RL': ratio_RL,
        'inv_ratio_RL_sq': inv_ratio_RL_sq,
        'zeta_over_gap': zeta_over_gap,
        'pred_asym': pred_asym,
    }


def compute_H_pole_decomposition(
    mu_val: mp.mpf, lam_0: mp.mpf, evals_e: list[mp.mpf], d_k_sq: list[mp.mpf]
) -> tuple[mp.mpf, list[mp.mpf]]:
    """Compute H(mu) and its individual pole contributions H_k(mu)."""
    N = len(evals_e) - 1
    d0_sq = d_k_sq[0]

    H_terms = [d0_sq]
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

def run_cell77() -> None:
    print("=" * 105)
    print("CELL 77 — STIELTJES OUTER PRODUCT FACTOR Pi_j, PAIRWISE GAP LOCALIZATION,")
    print("         CHARACTERISTIC-POLYNOMIAL RATIO DECOMPOSITION, AND ASYMMETRY BALANCE")
    print("=" * 105)
    print(f"Parameters: c = {C_PARAM}, L = {float(L_PARAM):.14f}, T = {T_PARAM}, dps = {mp.mp.dps}, V_* = {float(BARRIER_V_STAR)}")
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

        # Find Stieltjes zeros z_j^*
        zeros_z: list[mp.mpf] = []
        for j in range(N):
            z_j, _ = find_stieltjes_zero_bisection(j, evals_e, d_k_sq_list)
            zeros_z.append(z_j)

        print(f"--- DIMENSION N = {N:2d} (Elapsed: {time.time() - t0:.2f}s) ---")

        # ---------------------------------------------------------------------
        # Test A: Multi-Modal Audit of Pi_j
        # ---------------------------------------------------------------------
        print("  [Test A] Outer Stieltjes Factor Pi_j Across Bound Modes:")
        print(f"    {'j':>2s} | {'alpha_j':>14s} | {'zeta_j':>14s} | {'Pi_j':>14s} | {'Pi_zeros':>14s} | {'Pi_evals':>14s}")
        print("    " + "-" * 85)

        records_modes: list[dict] = []
        for j in range(min(4, N)):
            rec_j = evaluate_outer_product_detailed(j, evals_e, zeros_z, d_k_sq_list, evals_o)
            records_modes.append(rec_j)
            print(f"    {j:2d} | {mp.nstr(rec_j['alpha_eigen'], 6):>14s} | {mp.nstr(rec_j['zeta_j'], 6):>14s} | {mp.nstr(rec_j['pi_outer'], 6):>14s} | {mp.nstr(rec_j['pi_zeros'], 6):>14s} | {mp.nstr(rec_j['pi_evals'], 6):>14s}")

        synthesis_table1_records.append({
            'N': N,
            'pi_0': records_modes[0]['pi_outer'],
            'pi_1': records_modes[1]['pi_outer'],
            'pi_2': records_modes[2]['pi_outer'],
            'pi_3': records_modes[3]['pi_outer'] if N > 3 else mp.mpf(0),
        })

        # ---------------------------------------------------------------------
        # Test B: Pairwise Mode Cancellation & Localization for Mode j = 2
        # ---------------------------------------------------------------------
        rec_2 = records_modes[2]
        print(f"  [Test B] Pairwise Cancellation Factors omega_{{2, l}} for Mode j = 2:")
        print(f"    {'l':>2s} | {'l - j':>5s} | {'omega_{2, l}':>16s} | {'|omega - 1|':>14s}")
        print("    " + "-" * 45)

        max_dev_far = mp.mpf(0)
        for pw in rec_2['pairwise_omega']:
            l_idx = pw['l']
            dist = pw['dist_index']
            if abs(dist) >= 2 and pw['deviation'] > max_dev_far:
                max_dev_far = pw['deviation']

            if abs(dist) <= 3 or l_idx == N - 1:
                print(f"    {l_idx:2d} | {dist:5d} | {mp.nstr(pw['omega_l'], 8):>16s} | {mp.nstr(pw['deviation'], 4):>14s}")
            elif abs(dist) == 4:
                print(f"    {'...':>2s} | {'...':>5s} | {'...':>16s} | {'...':>14s}")

        synthesis_table2_records.append({
            'N': N,
            'omega_l0': rec_2['pairwise_omega'][0]['omega_l'],  # l = 0 (dist = -2)
            'omega_l1': (rec_2['pairwise_omega'][1]['omega_l'] if len(rec_2['pairwise_omega']) > 1 else mp.mpf(1)),
            'max_dev_far': max_dev_far,
        })

        # ---------------------------------------------------------------------
        # Test C: Characteristic Polynomial & Zero Polynomial Representation
        # ---------------------------------------------------------------------
        err_poly_alpha = abs(rec_2['alpha_poly'] / rec_2['alpha_eigen'] - 1)
        err_poly_pi = abs(rec_2['pi_poly'] / rec_2['pi_outer'] - 1)

        print(f"  [Test C] Polynomial Ratio Identity for Mode j = 2:")
        print(f"    |P_zero(E_3) / P_zero(E_2)|          = {mp.nstr(rec_2['ratio_p_zero'], 8)}")
        print(f"    |P_even'(E_2) / P_even'(E_3)|        = {mp.nstr(rec_2['ratio_p_even_prime'], 8)}")
        print(f"    alpha_2 (via Polynomials)           = {mp.nstr(rec_2['alpha_poly'], 8)}")
        print(f"    alpha_2 (Direct Eigenvalues)         = {mp.nstr(rec_2['alpha_eigen'], 8)}")
        print(f"    Polynomial Identity Relative Error   = {mp.nstr(err_poly_alpha, 4)}")
        print(f"    Pi_2 (via Polynomials)              = {mp.nstr(rec_2['pi_poly'], 8)}")
        print(f"    Pi_2 (Direct Outer Product)          = {mp.nstr(rec_2['pi_outer'], 8)}")
        print(f"    Pi_2 Identity Relative Error         = {mp.nstr(err_poly_pi, 4)}")

        synthesis_table3_records.append({
            'N': N,
            'ratio_p_zero': rec_2['ratio_p_zero'],
            'ratio_p_even_prime': rec_2['ratio_p_even_prime'],
            'err_poly_alpha': err_poly_alpha,
            'err_poly_pi': err_poly_pi,
        })

        # ---------------------------------------------------------------------
        # Test D: Master Asymmetry Balance Decomposition for Mode j = 2
        # ---------------------------------------------------------------------
        H_total_2, H_terms_2 = compute_H_pole_decomposition(evals_o[2], lam_0, evals_e, d_k_sq_list)
        actual_asym = H_terms_2[3] / H_terms_2[2]
        pred_asym = rec_2['pred_asym']
        asym_res = abs(actual_asym - pred_asym)

        print(f"  [Test D] Master Asymmetry Balance for Mode j = 2:")
        print(f"    zeta_2 = {mp.nstr(rec_2['zeta_j'], 8)},  R_2 / L_2 = {mp.nstr(rec_2['ratio_RL'], 8)},  zeta / (R/L) = {mp.nstr(rec_2['zeta_over_gap'], 8)}")
        print(f"    Pi_2   = {mp.nstr(rec_2['pi_outer'], 8)},  (L_2 / R_2)^2 = {mp.nstr(rec_2['inv_ratio_RL_sq'], 8)}")
        print(f"    Asymmetry H_3 / H_2: Predicted = {mp.nstr(pred_asym, 8)}, Actual = {mp.nstr(actual_asym, 8)}, Residual = {mp.nstr(asym_res, 4)}")
        print()

        synthesis_table4_records.append({
            'N': N,
            'zeta_2': rec_2['zeta_j'],
            'ratio_RL': rec_2['ratio_RL'],
            'zeta_over_gap': rec_2['zeta_over_gap'],
            'pi_2': rec_2['pi_outer'],
            'inv_ratio_RL_sq': rec_2['inv_ratio_RL_sq'],
            'pred_asym': pred_asym,
            'actual_asym': actual_asym,
            'asym_res': asym_res,
        })

    # -------------------------------------------------------------------------
    # Multi-Dimension Synthesis Tables
    # -------------------------------------------------------------------------
    print("=" * 105)
    print("SYNTHESIS TABLE 1: OUTER STIELTJES FACTOR Pi_j ACROSS MODES (TEST A)")
    print("=" * 105)
    print(f"{'N':>3s} | {'Pi_0':>16s} | {'Pi_1':>16s} | {'Pi_2':>16s} | {'Pi_3':>16s}")
    print("-" * 105)
    for rec in synthesis_table1_records:
        print(f"{rec['N']:3d} | {mp.nstr(rec['pi_0'], 6):>16s} | {mp.nstr(rec['pi_1'], 6):>16s} | {mp.nstr(rec['pi_2'], 6):>16s} | {mp.nstr(rec['pi_3'], 6):>16s}")
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 2: PAIRWISE LOCALIZATION FOR MODE j = 2 (TEST B)")
    print("=" * 105)
    print(f"{'N':>3s} | {'omega_{2, 0} (l=0)':>22s} | {'omega_{2, 1} (l=1)':>22s} | {'Max |omega - 1| (|dist|>=2)':>30s}")
    print("-" * 105)
    for rec in synthesis_table2_records:
        print(f"{rec['N']:3d} | {mp.nstr(rec['omega_l0'], 8):>22s} | {mp.nstr(rec['omega_l1'], 8):>22s} | {mp.nstr(rec['max_dev_far'], 4):>30s}")
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 3: CHARACTERISTIC-POLYNOMIAL RATIO IDENTITY FOR MODE j = 2 (TEST C)")
    print("=" * 105)
    col_even_prime = "|P_even'(E2)/P_even'(E3)|"
    print(f"{'N':>3s} | {'|P_zero(E3)/P_zero(E2)|':>25s} | {col_even_prime:>25s} | {'RelErr(alpha)':>16s} | {'RelErr(Pi)':>16s}")
    print("-" * 105)
    for rec in synthesis_table3_records:
        print(f"{rec['N']:3d} | {mp.nstr(rec['ratio_p_zero'], 8):>25s} | {mp.nstr(rec['ratio_p_even_prime'], 8):>25s} | {mp.nstr(rec['err_poly_alpha'], 4):>16s} | {mp.nstr(rec['err_poly_pi'], 4):>16s}")
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 4: DECOMPOSITION OF THE MASTER ASYMMETRY BALANCE FOR MODE j = 2 (TEST D)")
    print("=" * 105)
    print(f"{'N':>3s} | {'zeta_2':>12s} | {'R_2 / L_2':>10s} | {'zeta/(R/L)':>12s} | {'Pi_2':>10s} | {'(L_2/R_2)^2':>14s} | {'H_3 / H_2':>12s} | {'Residual':>12s}")
    print("-" * 105)
    for rec in synthesis_table4_records:
        print(f"{rec['N']:3d} | {mp.nstr(rec['zeta_2'], 6):>12s} | {mp.nstr(rec['ratio_RL'], 5):>10s} | {mp.nstr(rec['zeta_over_gap'], 6):>12s} | {mp.nstr(rec['pi_2'], 6):>10s} | {mp.nstr(rec['inv_ratio_RL_sq'], 6):>14s} | {mp.nstr(rec['actual_asym'], 6):>12s} | {mp.nstr(rec['asym_res'], 4):>12s}")
    print()

    print("=" * 80)
    print("CELL 77 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell77()
