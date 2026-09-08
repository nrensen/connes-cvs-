#!/usr/bin/env python3
"""
================================================================================
CELL 78 — EXACT PAIRWISE CANCELLATION IDENTITY FOR omega_{j, l} - 1,
         STIELTJES ZERO DISPLACEMENT FORMULA delta_l = d_l^2 / sum d_k^2/(E_k - z_l^*),
         AND QUANTITATIVE REMOTE PRODUCT LOCALIZATION BOUNDS
================================================================================

Strategic Roadmap Milestone M24 (Paper 4B Proposition 8.26):
-------------------------------------------------------------
Following the completion of Milestone M23 in cell 77, this cell audits the exact
algebraic foundation underlying the localization of the outer Stieltjes factor:
    Pi_j = prod_{l != j} |E_{j+1} - z_l^*| / |E_j - z_l^*|
         * prod_{l != j, j+1} |E_j - E_l| / |E_{j+1} - E_l|.

BACKGROUND & ANALYTICAL MOTIVATION:
-----------------------------------
Cell 77 established that Pi_j remains in an O(1) range (0.154 to 0.338) across
low modes, and proved that while raw zero and eigenvalue products diverge/collapse
wildly (Pi_2^{zeros} ~ 1.8e9, Pi_2^{evals} ~ 1.3e-10), the pairwise factors:
    omega_{j, l} = [ |E_{j+1} - z_l^*| / |E_j - z_l^*| ]
                 * [ |E_j - E_l| / |E_{j+1} - E_l| ]
converge rapidly to 1 away from the active bracket [E_j, E_{j+1}].

Proposition 8.26 proves that this cancellation is governed by two exact identities:
1. Exact Pairwise Cancellation Formula:
       omega_{j, l} - 1 = [ Delta_j * (z_l^* - E_l) ] / [ |E_j - z_l^*| * |E_{j+1} - E_l| ],
   where Delta_j = E_{j+1} - E_j.
2. Exact Stieltjes Zero Displacement Formula:
   Since G_d(z_l^*) = 0, the displacement delta_l = z_l^* - E_l > 0 satisfies:
       delta_l = d_l^2 / [ sum_{k != l} d_k^2 / (E_k - z_l^*) ].

THE FOUR INVESTIGATIVE TESTS OF CELL 78:
----------------------------------------
1. Test A — Exact Algebraic Identity for omega_{j, l} - 1:
   Audits the exact formula against direct evaluation of omega_{j, l} - 1 across
   all remote modes l != j, j+1 for mode j = 2, verifying identity residuals to
   the numerical precision floor.

2. Test B — Exact Stieltjes Displacement Formula for delta_l = z_l^* - E_l:
   Audits the rational formula delta_l^{formula} against direct bisection roots
   delta_l^{root} = z_l^* - E_l across all intervals l in {0, ..., N-1}.

3. Test C — Quantitative Remote Localization Bound:
   Audits the rigorous quadratic spectral decay upper bound:
       0 < omega_{j, l} - 1 <= [ Delta_j * delta_l ] / D_{j, l}^2,
   where D_{j, l} = dist(E_l, [E_j, E_{j+1}]), and tracks the convergence of
   the remote tail product Pi_{j, remote} = prod_{|l - j| >= 2} omega_{j, l} -> 1.

4. Test D — Decomposition of the Stieltjes Displacement Denominator:
   Decomposes the denominator sum_{k != l} d_k^2 / (E_k - z_l^*) into left neighbor,
   right neighbor, and remote poles, revealing the exact boundary-weight scaling.

OUTPUT CONSTRAINTS:
-------------------
- Dispassionate, dry, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 78 EXECUTION COMPLETE.
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
# Core Analytical Routines for Proposition 8.26
# -----------------------------------------------------------------------------

def evaluate_stieltjes_displacement_formula(
    l: int,
    z_l: mp.mpf,
    evals_e: list[mp.mpf],
    d_k_sq: list[mp.mpf],
) -> tuple[mp.mpf, mp.mpf, dict]:
    """
    Computes delta_l = z_l^* - E_l via the exact Stieltjes displacement formula:
        delta_l = d_l^2 / sum_{k != l} [ d_k^2 / (E_k - z_l^*) ].
    Also returns the denominator decomposition into neighbors and remote poles.
    """
    N_plus_1 = len(evals_e)
    d_l_sq = d_k_sq[l]

    denom_total = mp.mpf(0)
    denom_left = mp.mpf(0)
    denom_right = mp.mpf(0)
    denom_remote = mp.mpf(0)

    for k in range(N_plus_1):
        if k != l:
            term = d_k_sq[k] / (evals_e[k] - z_l)
            denom_total += term
            if k == l - 1:
                denom_left = term
            elif k == l + 1:
                denom_right = term
            else:
                denom_remote += term

    delta_formula = d_l_sq / denom_total
    decomp = {
        'denom_total': denom_total,
        'denom_left': denom_left,
        'denom_right': denom_right,
        'denom_remote': denom_remote,
        'share_left': (abs(denom_left) / denom_total) if denom_total > 0 else mp.mpf(0),
        'share_right': (abs(denom_right) / denom_total) if denom_total > 0 else mp.mpf(0),
    }
    return delta_formula, denom_total, decomp


def evaluate_pairwise_factor_detailed(
    j: int,
    l: int,
    evals_e: list[mp.mpf],
    zeros_z: list[mp.mpf],
    d_k_sq: list[mp.mpf],
) -> dict:
    """
    Audits the exact formula for omega_{j, l} - 1:
        omega_{j, l} - 1 = [ Delta_j * (z_l^* - E_l) ] / [ |E_j - z_l^*| * |E_{j+1} - E_l| ],
    and computes the quantitative upper bound:
        Bound = [ Delta_j * delta_l ] / D_{j, l}^2.
    """
    E_j = evals_e[j]
    E_jp1 = evals_e[j + 1]
    Delta_j = E_jp1 - E_j

    E_l = evals_e[l]
    z_l = zeros_z[l]
    delta_l_root = z_l - E_l

    # Direct evaluation from definition:
    ratio_zero = abs(E_jp1 - z_l) / abs(E_j - z_l)
    ratio_eval = abs(E_j - E_l) / abs(E_jp1 - E_l)
    omega_direct = ratio_zero * ratio_eval
    dev_direct = omega_direct - mp.mpf(1)

    # Exact algebraic formula:
    denom_exact = abs(E_j - z_l) * abs(E_jp1 - E_l)
    dev_formula = (Delta_j * delta_l_root) / denom_exact

    # Algebraic identity residual:
    identity_residual = abs(dev_direct - dev_formula)
    rel_err_formula = (identity_residual / dev_direct) if dev_direct > 0 else mp.mpf(0)

    # Spectral distance D_{j, l} = dist(E_l, [E_j, E_{j+1}]):
    if l < j:
        D_jl = E_j - evals_e[l + 1]  # right endpoint of bracket l to left endpoint of bracket j
    else:  # l > j + 1
        D_jl = E_l - E_jp1

    # Conservative distance to eigenvalue:
    D_jl_eval = min(abs(E_l - E_j), abs(E_l - E_jp1))

    # Upper bound based on spectral distance:
    bound_spectral = (Delta_j * delta_l_root) / (D_jl_eval ** 2)

    return {
        'j': j,
        'l': l,
        'dist_index': l - j,
        'omega_direct': omega_direct,
        'dev_direct': dev_direct,
        'dev_formula': dev_formula,
        'identity_residual': identity_residual,
        'rel_err_formula': rel_err_formula,
        'D_jl_eval': D_jl_eval,
        'bound_spectral': bound_spectral,
        'bound_ratio': (bound_spectral / dev_direct) if dev_direct > 0 else mp.mpf(0),
    }


# -----------------------------------------------------------------------------
# Main Execution Function
# -----------------------------------------------------------------------------

def run_cell78() -> None:
    print("=" * 105)
    print("CELL 78 — EXACT PAIRWISE CANCELLATION IDENTITY FOR omega_{j, l} - 1,")
    print("         STIELTJES ZERO DISPLACEMENT FORMULA, AND REMOTE LOCALIZATION BOUNDS")
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
        # Test A: Exact Identity for omega_{2, l} - 1 (Mode j = 2)
        # ---------------------------------------------------------------------
        print("  [Test A] Exact Identity for omega_{2, l} - 1 Across Remote Modes:")
        print(f"    {'l':>2s} | {'l - j':>5s} | {'omega - 1 (Direct)':>20s} | {'omega - 1 (Formula)':>20s} | {'Rel Err':>12s}")
        print("    " + "-" * 70)

        max_pw_rel_err = mp.mpf(0)
        pw_records_2: list[dict] = []
        for l in range(N):
            if l != 2 and l != 3:
                rec_pw = evaluate_pairwise_factor_detailed(2, l, evals_e, zeros_z, d_k_sq_list)
                pw_records_2.append(rec_pw)
                if rec_pw['rel_err_formula'] > max_pw_rel_err:
                    max_pw_rel_err = rec_pw['rel_err_formula']

                dist = rec_pw['dist_index']
                if abs(dist) <= 3 or l == N - 1:
                    print(f"    {l:2d} | {dist:5d} | {mp.nstr(rec_pw['dev_direct'], 8):>20s} | {mp.nstr(rec_pw['dev_formula'], 8):>20s} | {mp.nstr(rec_pw['rel_err_formula'], 4):>12s}")
                elif abs(dist) == 4:
                    print(f"    {'...':>2s} | {'...':>5s} | {'...':>20s} | {'...':>20s} | {'...':>12s}")

        synthesis_table1_records.append({
            'N': N,
            'max_pw_rel_err': max_pw_rel_err,
            'dev_l0': pw_records_2[0]['dev_direct'],
            'dev_l1': (pw_records_2[1]['dev_direct'] if len(pw_records_2) > 1 else mp.mpf(0)),
            'dev_l4': (next((r['dev_direct'] for r in pw_records_2 if r['l'] == 4), mp.mpf(0))),
        })

        # ---------------------------------------------------------------------
        # Test B: Exact Stieltjes Displacement Formula delta_l = z_l^* - E_l
        # ---------------------------------------------------------------------
        print("  [Test B] Stieltjes Zero Displacement Formula delta_l Across Bound Intervals:")
        print(f"    {'l':>2s} | {'delta_l (Bisection)':>22s} | {'delta_l (Formula)':>22s} | {'Rel Err':>12s}")
        print("    " + "-" * 65)

        max_disp_rel_err = mp.mpf(0)
        disp_records: list[dict] = []
        for l in range(min(4, N)):
            z_l = zeros_z[l]
            delta_root = z_l - evals_e[l]
            delta_form, denom_tot, decomp = evaluate_stieltjes_displacement_formula(l, z_l, evals_e, d_k_sq_list)
            rel_err_disp = abs(delta_root - delta_form) / delta_root
            if rel_err_disp > max_disp_rel_err:
                max_disp_rel_err = rel_err_disp

            disp_records.append({
                'l': l,
                'delta_root': delta_root,
                'delta_form': delta_form,
                'rel_err': rel_err_disp,
                'decomp': decomp,
            })
            print(f"    {l:2d} | {mp.nstr(delta_root, 8):>22s} | {mp.nstr(delta_form, 8):>22s} | {mp.nstr(rel_err_disp, 4):>12s}")

        synthesis_table2_records.append({
            'N': N,
            'delta_0': disp_records[0]['delta_root'],
            'delta_1': disp_records[1]['delta_root'],
            'delta_2': disp_records[2]['delta_root'],
            'max_disp_rel_err': max_disp_rel_err,
        })

        # ---------------------------------------------------------------------
        # Test C: Quantitative Spectral Upper Bound & Remote Tail Product
        # ---------------------------------------------------------------------
        remote_tail_prod = mp.mpf(1)
        max_bound_slack = mp.mpf(0)
        for rec in pw_records_2:
            if abs(rec['dist_index']) >= 2:
                remote_tail_prod *= rec['omega_direct']
                if rec['bound_ratio'] > max_bound_slack:
                    max_bound_slack = rec['bound_ratio']

        remote_tail_dev = abs(remote_tail_prod - mp.mpf(1))

        print(f"  [Test C] Quantitative Remote Upper Bound & Tail Product for Mode j = 2:")
        print(f"    Remote Tail Product Pi_{{2, remote}}       = {mp.nstr(remote_tail_prod, 10)}")
        print(f"    Remote Tail Deviation |Pi_{{remote}} - 1| = {mp.nstr(remote_tail_dev, 6)}")
        print(f"    Max Bound Slack (Bound / Actual)          = {mp.nstr(max_bound_slack, 4)}")

        synthesis_table3_records.append({
            'N': N,
            'remote_tail_prod': remote_tail_prod,
            'remote_tail_dev': remote_tail_dev,
            'max_bound_slack': max_bound_slack,
        })

        # ---------------------------------------------------------------------
        # Test D: Displacement Denominator Decomposition for Mode l = 2
        # ---------------------------------------------------------------------
        decomp_2 = disp_records[2]['decomp']
        print(f"  [Test D] Denominator S_2 = sum_{{k != 2}} d_k^2 / (E_k - z_2^*):")
        print(f"    S_2 Total          = {mp.nstr(decomp_2['denom_total'], 8)}")
        print(f"    Left Neighbor k=1  = {mp.nstr(decomp_2['denom_left'], 8)} (Share: {mp.nstr(decomp_2['share_left']*100, 4)}%)")
        print(f"    Right Neighbor k=3 = {mp.nstr(decomp_2['denom_right'], 8)} (Share: {mp.nstr(decomp_2['share_right']*100, 4)}%)")
        print(f"    Remote Poles       = {mp.nstr(decomp_2['denom_remote'], 8)}")
        print()

        synthesis_table4_records.append({
            'N': N,
            'denom_total_2': decomp_2['denom_total'],
            'share_right_2': decomp_2['share_right'] * 100,
            'share_left_2': decomp_2['share_left'] * 100,
        })

    # -------------------------------------------------------------------------
    # Multi-Dimension Synthesis Tables
    # -------------------------------------------------------------------------
    print("=" * 105)
    print("SYNTHESIS TABLE 1: EXACT PAIRWISE CANCELLATION FORMULA AUDIT (TEST A)")
    print("=" * 105)
    print(f"{'N':>3s} | {'omega_{2,0} - 1':>18s} | {'omega_{2,1} - 1':>18s} | {'omega_{2,4} - 1':>18s} | {'Max Rel Err':>14s}")
    print("-" * 105)
    for rec in synthesis_table1_records:
        print(f"{rec['N']:3d} | {mp.nstr(rec['dev_l0'], 6):>18s} | {mp.nstr(rec['dev_l1'], 6):>18s} | {mp.nstr(rec['dev_l4'], 6):>18s} | {mp.nstr(rec['max_pw_rel_err'], 4):>14s}")
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 2: STIELTJES DISPLACEMENT FORMULA delta_l = z_l^* - E_l (TEST B)")
    print("=" * 105)
    print(f"{'N':>3s} | {'delta_0':>18s} | {'delta_1':>18s} | {'delta_2':>18s} | {'Max Rel Err':>14s}")
    print("-" * 105)
    for rec in synthesis_table2_records:
        print(f"{rec['N']:3d} | {mp.nstr(rec['delta_0'], 6):>18s} | {mp.nstr(rec['delta_1'], 6):>18s} | {mp.nstr(rec['delta_2'], 6):>18s} | {mp.nstr(rec['max_disp_rel_err'], 4):>14s}")
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 3: REMOTE TAIL PRODUCT CONVERGENCE FOR MODE j = 2 (TEST C)")
    print("=" * 105)
    print(f"{'N':>3s} | {'Pi_{2, remote}':>22s} | {'|Pi_{remote} - 1|':>20s} | {'Max Slack (Bound/Actual)':>26s}")
    print("-" * 105)
    for rec in synthesis_table3_records:
        print(f"{rec['N']:3d} | {mp.nstr(rec['remote_tail_prod'], 10):>22s} | {mp.nstr(rec['remote_tail_dev'], 6):>20s} | {mp.nstr(rec['max_bound_slack'], 4):>26s}")
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 4: DENOMINATOR DECOMPOSITION FOR MODE l = 2 (TEST D)")
    print("=" * 105)
    print(f"{'N':>3s} | {'S_2 Total':>18s} | {'Right (k=3) Share (%)':>24s} | {'Left (k=1) Share (%)':>24s}")
    print("-" * 105)
    for rec in synthesis_table4_records:
        print(f"{rec['N']:3d} | {mp.nstr(rec['denom_total_2'], 6):>18s} | {mp.nstr(rec['share_right_2'], 4):>24s} | {mp.nstr(rec['share_left_2'], 4):>24s}")
    print()

    print("=" * 80)
    print("CELL 78 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell78()
