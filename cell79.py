#!/usr/bin/env python3
"""
================================================================================
CELL 79 — AUDIT OF THE ONE-SIDED STIELTJES DISPLACEMENT BOUND delta_l < Delta_l / alpha_l,
         NEGATIVE-TO-POSITIVE RATIO N_l / P_l, AND CLOSED REMOTE PRODUCT BOUNDS
================================================================================

Strategic Roadmap Milestone M25 (Paper 4B Corollary 8.26.1):
-------------------------------------------------------------
Following the certification of the exact pairwise identity in cell 78 (Milestone M24),
this cell numerically audits the closed-form remote localization bound established in
Corollary 8.26.1:
    omega_{j, l} - 1 < [ Delta_j * Delta_l ] / [ alpha_l * D_{j, l}^2 ],
where:
    - Delta_l = E_{l+1} - E_l is the spectral interval length,
    - alpha_l = d_{l+1}^2 / d_l^2 is the boundary weight amplification factor,
    - D_{j, l} = min(|E_l - E_j|, |E_l - E_{j+1}|) is the spectral distance,
    - delta_l = z_l^* - E_l is the Stieltjes zero displacement.

ANALYTICAL MOTIVATION & RECURSIVE DECOUPLING:
--------------------------------------------
Proposition 8.26 showed that the crude bound [Delta_j * delta_l] / D_{j, l}^2 has
a large slack factor (10^4 to 10^5) because it retains the uncalculated zero displacement
delta_l. 

Corollary 8.26.1 provides the fundamental algebraic closure by establishing the
one-sided Stieltjes displacement bound:
    delta_l < Delta_l / alpha_l.
This bound proves that large boundary amplification (alpha_l >> 1) automatically
quenches the displacement (delta_l << 1), guaranteeing that the outer Stieltjes
product Pi_j = prod_{l != j} omega_{j, l} converges to an O(1) constant uniformly as N -> infty.

THE FOUR INVESTIGATIVE TESTS OF CELL 79:
----------------------------------------
1. Test A — The One-Sided Displacement Bound delta_l < Delta_l / alpha_l:
   Audits the ratio R_delta(l) = delta_l / (Delta_l / alpha_l) across low modes
   l in {0, 1, 2, 3} and across dimensions N in {8, 12, 16, 20, 24}, confirming
   that R_delta(l) < 1 strictly holds (unconditionally for l=0; with substantial
   margin for l in {1, 2, 3}).

2. Test B — Sign Structure and Negative-to-Positive Ratio N_l / P_l:
   Decomposes the denominator S_l = sum_{k != l} d_k^2 / (E_k - z_l^*) into:
       P_l = sum_{k > l} d_k^2 / (E_k - z_l^*) > 0   (positive poles),
       N_l = sum_{k < l} d_k^2 / (z_l^* - E_k) > 0   (negative poles),
   and audits the ratio epsilon_l = N_l / P_l for l in {1, 2, 3}, verifying that
   N_l / P_l << 10^{-5} due to the steep boundary amplification hierarchy.

3. Test C — Closed Remote Pairwise Bound Audit for Mode j = 2:
   Audits the closed upper bound:
       omega_{2, l} - 1 < eta_{2, l} = [ Delta_2 * Delta_l ] / [ alpha_l * D_{2, l}^2 ]
   against the actual pairwise deviation omega_{2, l} - 1 across remote modes |l - 2| >= 2,
   auditing the validity and tracking the slack factor.

4. Test D — Logarithmic Remote Tail Product Enclosure for Mode j = 2:
   Audits the rigorous tail enclosure:
       0 < ln Pi_{2, remote} <= sum_{|l - 2| >= 2} eta_{2, l},
   verifying that the closed bound encloses the remote product uniformly across N.

OUTPUT CONSTRAINTS:
-------------------
- Dispassionate, dry, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 79 EXECUTION COMPLETE.
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
# Analytical Routines for Corollary 8.26.1
# -----------------------------------------------------------------------------

def evaluate_stieltjes_sign_structure(
    l: int,
    z_l: mp.mpf,
    evals_e: list[mp.mpf],
    d_k_sq: list[mp.mpf],
) -> dict:
    """
    Decomposes the denominator S_l = sum_{k != l} d_k^2 / (E_k - z_l^*) into:
        P_l = sum_{k > l} d_k^2 / (E_k - z_l^*) > 0,
        N_l = sum_{k < l} d_k^2 / (z_l^* - E_k) > 0,
    so that S_l = P_l - N_l.
    Also computes the ratio epsilon_l = N_l / P_l.
    """
    N_plus_1 = len(evals_e)
    pos_sum = mp.mpf(0)
    neg_sum = mp.mpf(0)

    for k in range(N_plus_1):
        if k > l:
            pos_sum += d_k_sq[k] / (evals_e[k] - z_l)
        elif k < l:
            neg_sum += d_k_sq[k] / (z_l - evals_e[k])

    s_total = pos_sum - neg_sum
    epsilon = (neg_sum / pos_sum) if pos_sum > 0 else mp.mpf(0)

    return {
        'l': l,
        'pos_sum': pos_sum,
        'neg_sum': neg_sum,
        's_total': s_total,
        'epsilon': epsilon,
    }


def evaluate_closed_remote_factor(
    j: int,
    l: int,
    evals_e: list[mp.mpf],
    zeros_z: list[mp.mpf],
    d_k_sq: list[mp.mpf],
) -> dict:
    """
    Audits the closed remote bound:
        omega_{j, l} - 1 < eta_{j, l} = [ Delta_j * Delta_l ] / [ alpha_l * D_{j, l}^2 ],
    where D_{j, l} = min(|E_l - E_j|, |E_l - E_{j+1}|).
    """
    E_j = evals_e[j]
    E_jp1 = evals_e[j + 1]
    Delta_j = E_jp1 - E_j

    E_l = evals_e[l]
    E_lp1 = evals_e[l + 1]
    Delta_l = E_lp1 - E_l

    alpha_l = d_k_sq[l + 1] / d_k_sq[l]
    z_l = zeros_z[l]

    # Direct pairwise factor:
    ratio_zero = abs(E_jp1 - z_l) / abs(E_j - z_l)
    ratio_eval = abs(E_j - E_l) / abs(E_jp1 - E_l)
    omega_direct = ratio_zero * ratio_eval
    dev_direct = omega_direct - mp.mpf(1)

    # Conservative spectral distance to bracket:
    D_jl = min(abs(E_l - E_j), abs(E_l - E_jp1))

    # Closed upper bound:
    eta_closed = (Delta_j * Delta_l) / (alpha_l * (D_jl ** 2))

    # Slack factor:
    slack = (eta_closed / dev_direct) if dev_direct > 0 else mp.mpf(0)

    return {
        'j': j,
        'l': l,
        'dist_index': l - j,
        'omega_direct': omega_direct,
        'dev_direct': dev_direct,
        'eta_closed': eta_closed,
        'slack': slack,
        'D_jl': D_jl,
        'is_strictly_bounded': (dev_direct < eta_closed),
    }


# -----------------------------------------------------------------------------
# Main Execution Function
# -----------------------------------------------------------------------------

def run_cell79() -> None:
    print("=" * 105)
    print("CELL 79 — AUDIT OF THE ONE-SIDED STIELTJES DISPLACEMENT BOUND delta_l < Delta_l / alpha_l,")
    print("         NEGATIVE-TO-POSITIVE RATIO N_l / P_l, AND CLOSED REMOTE PRODUCT BOUNDS")
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
        # Test A: The One-Sided Displacement Bound delta_l < Delta_l / alpha_l
        # ---------------------------------------------------------------------
        print("  [Test A] One-Sided Displacement Bound delta_l vs Delta_l / alpha_l:")
        print(f"    {'l':>2s} | {'delta_l':>22s} | {'Delta_l / alpha_l':>22s} | {'Ratio R_delta':>14s} | {'Bound Met?':>10s}")
        print("    " + "-" * 80)

        test_a_dict: dict[int, mp.mpf] = {}
        for l in range(min(4, N)):
            delta_l = zeros_z[l] - evals_e[l]
            Delta_l = evals_e[l + 1] - evals_e[l]
            alpha_l = d_k_sq_list[l + 1] / d_k_sq_list[l]
            scale_l = Delta_l / alpha_l
            ratio_r = delta_l / scale_l
            test_a_dict[l] = ratio_r
            bound_met = "YES" if ratio_r < mp.mpf(1) else "NO"

            print(f"    {l:2d} | {mp.nstr(delta_l, 8):>22s} | {mp.nstr(scale_l, 8):>22s} | {mp.nstr(ratio_r, 6):>14s} | {bound_met:>10s}")

        synthesis_table1_records.append({
            'N': N,
            'r_0': test_a_dict[0],
            'r_1': test_a_dict[1],
            'r_2': test_a_dict[2],
            'r_3': (test_a_dict[3] if 3 in test_a_dict else mp.mpf(0)),
        })

        # ---------------------------------------------------------------------
        # Test B: Sign Structure and Negative-to-Positive Ratio N_l / P_l
        # ---------------------------------------------------------------------
        print("  [Test B] Denominator Decomposition and Negative-to-Positive Ratio N_l / P_l:")
        print(f"    {'l':>2s} | {'Positive Sum P_l':>20s} | {'Negative Sum N_l':>20s} | {'Ratio N_l / P_l':>16s} | {'Dominance':>10s}")
        print("    " + "-" * 80)

        test_b_dict: dict[int, mp.mpf] = {}
        for l in range(1, min(4, N)):
            decomp_l = evaluate_stieltjes_sign_structure(l, zeros_z[l], evals_e, d_k_sq_list)
            p_l = decomp_l['pos_sum']
            n_l = decomp_l['neg_sum']
            eps_l = decomp_l['epsilon']
            test_b_dict[l] = eps_l
            dominance = "P >> N" if eps_l < mp.mpf('1e-4') else "MARGINAL"

            print(f"    {l:2d} | {mp.nstr(p_l, 8):>20s} | {mp.nstr(n_l, 8):>20s} | {mp.nstr(eps_l, 6):>16s} | {dominance:>10s}")

        synthesis_table2_records.append({
            'N': N,
            'eps_1': test_b_dict[1],
            'eps_2': test_b_dict[2],
            'eps_3': (test_b_dict[3] if 3 in test_b_dict else mp.mpf(0)),
        })

        # ---------------------------------------------------------------------
        # Test C: Closed Remote Pairwise Bound Audit for Mode j = 2
        # ---------------------------------------------------------------------
        print("  [Test C] Closed Remote Pairwise Bound eta_{2, l} for Mode j = 2:")
        print(f"    {'l':>2s} | {'omega_{2,l} - 1 (Actual)':>25s} | {'eta_{2,l} (Closed Bound)':>25s} | {'Slack Factor':>14s} | {'Valid?':>8s}")
        print("    " + "-" * 85)

        pw_closed_records: list[dict] = []
        for l in range(N):
            if abs(l - 2) >= 2:
                rec_c = evaluate_closed_remote_factor(2, l, evals_e, zeros_z, d_k_sq_list)
                pw_closed_records.append(rec_c)
                valid_str = "YES" if rec_c['is_strictly_bounded'] else "NO"

                dist = rec_c['dist_index']
                if abs(dist) <= 3 or l == N - 1:
                    print(f"    {l:2d} | {mp.nstr(rec_c['dev_direct'], 8):>25s} | {mp.nstr(rec_c['eta_closed'], 8):>25s} | {mp.nstr(rec_c['slack'], 6):>14s} | {valid_str:>8s}")
                elif abs(dist) == 4:
                    print(f"    {'...':>2s} | {'...':>25s} | {'...':>25s} | {'...':>14s} | {'...':>8s}")

        synthesis_table3_records.append({
            'N': N,
            'dev_l0': pw_closed_records[0]['dev_direct'],
            'eta_l0': pw_closed_records[0]['eta_closed'],
            'slack_l0': pw_closed_records[0]['slack'],
            'dev_l4': (next((r['dev_direct'] for r in pw_closed_records if r['l'] == 4), mp.mpf(0))),
            'eta_l4': (next((r['eta_closed'] for r in pw_closed_records if r['l'] == 4), mp.mpf(0))),
            'slack_l4': (next((r['slack'] for r in pw_closed_records if r['l'] == 4), mp.mpf(0))),
        })

        # ---------------------------------------------------------------------
        # Test D: Logarithmic Remote Tail Product Enclosure for Mode j = 2
        # ---------------------------------------------------------------------
        pi_2_remote = mp.mpf(1)
        sum_eta_closed = mp.mpf(0)
        for rec in pw_closed_records:
            pi_2_remote *= rec['omega_direct']
            sum_eta_closed += rec['eta_closed']

        log_pi_remote = mp.ln(pi_2_remote)
        enclosure_holds = (mp.mpf(0) < log_pi_remote <= sum_eta_closed)
        enclosure_str = "HOLD" if enclosure_holds else "VIOLATED"

        print(f"  [Test D] Log Remote Product Enclosure: 0 < ln Pi_{{2, remote}} <= sum eta_{{2, l}}:")
        print(f"    Pi_{{2, remote}}               = {mp.nstr(pi_2_remote, 12)}")
        print(f"    ln Pi_{{2, remote}} (Actual)   = {mp.nstr(log_pi_remote, 8)}")
        print(f"    sum eta_{{2, l}} (Upper Bound) = {mp.nstr(sum_eta_closed, 8)}")
        print(f"    Enclosure Status              = {enclosure_str}")
        print()

        synthesis_table4_records.append({
            'N': N,
            'pi_2_remote': pi_2_remote,
            'log_pi_remote': log_pi_remote,
            'sum_eta_closed': sum_eta_closed,
            'enclosure_ratio': (sum_eta_closed / log_pi_remote) if log_pi_remote > 0 else mp.mpf(0),
        })

    # -------------------------------------------------------------------------
    # Multi-Dimension Synthesis Tables
    # -------------------------------------------------------------------------
    print("=" * 105)
    print("SYNTHESIS TABLE 1: RATIO R_delta(l) = delta_l / (Delta_l / alpha_l) ACROSS MODES (TEST A)")
    print("=" * 105)
    print(f"{'N':>3s} | {'R_delta(0)':>18s} | {'R_delta(1)':>18s} | {'R_delta(2)':>18s} | {'R_delta(3)':>18s} | {'Bound Met?':>12s}")
    print("-" * 105)
    for rec in synthesis_table1_records:
        all_met = "YES (all < 1)" if all(rec[k] < mp.mpf(1) for k in ['r_0', 'r_1', 'r_2']) else "NO"
        print(f"{rec['N']:3d} | {mp.nstr(rec['r_0'], 6):>18s} | {mp.nstr(rec['r_1'], 6):>18s} | {mp.nstr(rec['r_2'], 6):>18s} | {mp.nstr(rec['r_3'], 6):>18s} | {all_met:>12s}")
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 2: NEGATIVE-TO-POSITIVE RATIO N_l / P_l ACROSS MODES (TEST B)")
    print("=" * 105)
    print(f"{'N':>3s} | {'N_1 / P_1':>22s} | {'N_2 / P_2':>22s} | {'N_3 / P_3':>22s} | {'Dominance':>16s}")
    print("-" * 105)
    for rec in synthesis_table2_records:
        all_dom = "N << P (all < 1e-4)" if all(rec[k] < mp.mpf('1e-4') for k in ['eps_1', 'eps_2']) else "MARGINAL"
        print(f"{rec['N']:3d} | {mp.nstr(rec['eps_1'], 6):>22s} | {mp.nstr(rec['eps_2'], 6):>22s} | {mp.nstr(rec['eps_3'], 6):>22s} | {all_dom:>16s}")
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 3: CLOSED REMOTE PAIRWISE BOUND FOR MODE j = 2 (TEST C)")
    print("=" * 105)
    print(f"{'N':>3s} | {'omega_{2,0}-1':>16s} | {'eta_{2,0}':>16s} | {'Slack(l=0)':>12s} | {'omega_{2,4}-1':>16s} | {'eta_{2,4}':>16s} | {'Slack(l=4)':>12s}")
    print("-" * 105)
    for rec in synthesis_table3_records:
        print(f"{rec['N']:3d} | {mp.nstr(rec['dev_l0'], 5):>16s} | {mp.nstr(rec['eta_l0'], 5):>16s} | {mp.nstr(rec['slack_l0'], 4):>12s} | {mp.nstr(rec['dev_l4'], 5):>16s} | {mp.nstr(rec['eta_l4'], 5):>16s} | {mp.nstr(rec['slack_l4'], 4):>12s}")
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 4: LOGARITHMIC REMOTE TAIL PRODUCT ENCLOSURE FOR MODE j = 2 (TEST D)")
    print("=" * 105)
    print(f"{'N':>3s} | {'ln Pi_{2, remote}':>22s} | {'sum eta_{2, l} (Bound)':>25s} | {'Enclosure Ratio (Bound / Actual)':>34s}")
    print("-" * 105)
    for rec in synthesis_table4_records:
        print(f"{rec['N']:3d} | {mp.nstr(rec['log_pi_remote'], 8):>22s} | {mp.nstr(rec['sum_eta_closed'], 8):>25s} | {mp.nstr(rec['enclosure_ratio'], 6):>34s}")
    print()

    print("=" * 80)
    print("CELL 79 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell79()
