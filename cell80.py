#!/usr/bin/env python3
"""
================================================================================
CELL 80 — SPECTRAL-WIDE AUDIT OF THE STIELTJES SIGN RATIO epsilon_l = N_l / P_l,
         UPPER-EDGE MODE FORENSICS, CORRECTED CLOSED REMOTE PAIRWISE BOUNDS,
         AND BOUNDARY-WEIGHT LADDER GEOMETRIC DECAY
================================================================================

Strategic Roadmap Milestone M26 (Paper 4B Section 8.25 & Corollary 8.26.1):
----------------------------------------------------------------------------
Following the discovery in cell 79 (Milestone M25) that the bare closed remote bound
eta_{2, l} = [ Delta_2 * Delta_l ] / [ alpha_l * D_{2, l}^2 ] failed at the upper
spectral edge (l = N-1), this cell audits the complete spectral profile of the
Stieltjes sign ratio:
    epsilon_l = N_l / P_l,
where:
    P_l = sum_{k > l} d_k^2 / (E_k - z_l^*) > 0,
    N_l = sum_{k < l} d_k^2 / (z_l^* - E_k) > 0,
and tests the mathematically rigorous closed remote bound:
    omega_{j, l} - 1 < eta_{j, l}^{corr} = [ Delta_j * Delta_l ] / [ alpha_l * (1 - epsilon_l) * D_{j, l}^2 ].

ANALYTICAL MOTIVATION:
----------------------
1. Corollary 8.26.1 established that for l >= 1, the displacement bound is:
       delta_l < Delta_l / [ alpha_l * (1 - epsilon_l) ].
   While delta_l < Delta_l / alpha_l held empirically for low modes l in {0, 1, 2, 3},
   omitting the factor (1 - epsilon_l)^{-1} in the remote pairwise bound caused
   failure at the upper boundary l = N-1 (l = 15 at N=16, l = 23 at N=24).

2. This cell investigates whether:
   (a) epsilon_l < 1 holds uniformly across all modes l in {1, ..., N-1},
   (b) The corrected bound eta_{2, l}^{corr} restores strict validity across all modes,
   (c) The negative sum N_l is suppressed by geometric decay in the boundary-weight ladder:
       d_k^2 / d_l^2 <= C * q^{l - k}  (for k < l).

THE FOUR INVESTIGATIVE TESTS OF CELL 80:
----------------------------------------
1. Test A — Full Spectral Sweep of the Sign Ratio epsilon_l = N_l / P_l:
   Evaluates epsilon_l across all l in {1, ..., N-1} for N in {8, 12, 16, 20, 24},
   tracking max_l epsilon_l and verifying whether epsilon_l <= epsilon_* < 1 uniformly.

2. Test B — Upper-Edge Mode Forensics for l = N - 1:
   Decomposes the top interval l = N-1, auditing d_{N-1}^2, d_N^2, alpha_{N-1},
   Delta_{N-1}, delta_{N-1}, P_{N-1}, N_{N-1}, epsilon_{N-1}, and the correction
   factor (1 - epsilon_{N-1})^{-1}.

3. Test C — Corrected Closed Remote Pairwise Bound Audit for Mode j = 2:
   Audits the corrected bound eta_{2, l}^{corr} against actual omega_{2, l} - 1
   across all remote modes |l - 2| >= 2, verifying whether validity is restored at
   l = 15 (N=16) and l = 23 (N=24).

4. Test D — Boundary-Weight Ladder Geometric Decay Audit:
   Evaluates the consecutive weight ratios r_l = d_{l-1}^2 / d_l^2 and effective
   geometric decay factors q_l = max_{k < l} (d_k^2 / d_l^2)^{1 / (l - k)} across
   modes, testing whether boundary weights supply geometric suppression of N_l.

OUTPUT CONSTRAINTS:
-------------------
- Dispassionate, dry, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 80 EXECUTION COMPLETE.
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
# Analytical Routines for Milestone M26
# -----------------------------------------------------------------------------

def evaluate_sign_ratio_full(
    l: int,
    z_l: mp.mpf,
    evals_e: list[mp.mpf],
    d_k_sq: list[mp.mpf],
) -> dict:
    """
    Computes positive sum P_l, negative sum N_l, and ratio epsilon_l = N_l / P_l.
    """
    N_plus_1 = len(evals_e)
    pos_sum = mp.mpf(0)
    neg_sum = mp.mpf(0)

    for k in range(N_plus_1):
        if k > l:
            pos_sum += d_k_sq[k] / (evals_e[k] - z_l)
        elif k < l:
            neg_sum += d_k_sq[k] / (z_l - evals_e[k])

    epsilon = (neg_sum / pos_sum) if pos_sum > 0 else mp.mpf(0)
    correction_factor = (1 / (mp.mpf(1) - epsilon)) if epsilon < mp.mpf(1) else mp.mpf('inf')

    return {
        'l': l,
        'pos_sum': pos_sum,
        'neg_sum': neg_sum,
        'epsilon': epsilon,
        'correction_factor': correction_factor,
    }


def evaluate_corrected_closed_remote_factor(
    j: int,
    l: int,
    evals_e: list[mp.mpf],
    zeros_z: list[mp.mpf],
    d_k_sq: list[mp.mpf],
    epsilon_l: mp.mpf,
) -> dict:
    """
    Audits the corrected closed remote bound:
        omega_{j, l} - 1 < eta_{j, l}^{corr} = [ Delta_j * Delta_l ] / [ alpha_l * (1 - epsilon_l) * D_{j, l}^2 ],
    and compares with the bare candidate bound:
        eta_{j, l}^{bare} = [ Delta_j * Delta_l ] / [ alpha_l * D_{j, l}^2 ].
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

    # Conservative spectral distance:
    D_jl = min(abs(E_l - E_j), abs(E_l - E_jp1))

    # Bare candidate bound:
    eta_bare = (Delta_j * Delta_l) / (alpha_l * (D_jl ** 2))

    # Corrected bound with (1 - epsilon_l)^{-1}:
    if epsilon_l < mp.mpf(1):
        corr_factor = 1 / (mp.mpf(1) - epsilon_l)
        eta_corr = eta_bare * corr_factor
    else:
        corr_factor = mp.mpf('inf')
        eta_corr = mp.mpf('inf')

    is_bare_valid = (dev_direct < eta_bare)
    is_corr_valid = (dev_direct < eta_corr)

    return {
        'j': j,
        'l': l,
        'dist_index': l - j,
        'omega_direct': omega_direct,
        'dev_direct': dev_direct,
        'eta_bare': eta_bare,
        'eta_corr': eta_corr,
        'is_bare_valid': is_bare_valid,
        'is_corr_valid': is_corr_valid,
        'corr_factor': corr_factor,
        'slack_corr': (eta_corr / dev_direct) if dev_direct > 0 else mp.mpf(0),
    }


def evaluate_ladder_geometric_decay(
    l: int,
    d_k_sq: list[mp.mpf],
) -> dict:
    """
    Audits the boundary-weight ratio ladder for mode l:
        r_l = d_{l-1}^2 / d_l^2,
        q_l = max_{k < l} (d_k^2 / d_l^2)^{1 / (l - k)}.
    """
    if l == 0:
        return {'l': 0, 'r_l': mp.mpf(0), 'q_l': mp.mpf(0)}

    d_l_sq = d_k_sq[l]
    r_l = d_k_sq[l - 1] / d_l_sq

    max_root_rate = mp.mpf(0)
    for k in range(l):
        step = l - k
        ratio = d_k_sq[k] / d_l_sq
        root_rate = mp.power(ratio, mp.mpf(1) / step)
        if root_rate > max_root_rate:
            max_root_rate = root_rate

    return {
        'l': l,
        'r_l': r_l,
        'q_l': max_root_rate,
    }


# -----------------------------------------------------------------------------
# Main Execution Function
# -----------------------------------------------------------------------------

def run_cell80() -> None:
    print("=" * 105)
    print("CELL 80 — SPECTRAL-WIDE AUDIT OF SIGN RATIO epsilon_l = N_l / P_l, UPPER-EDGE FORENSICS,")
    print("         CORRECTED CLOSED REMOTE PAIRWISE BOUNDS, AND WEIGHT LADDER GEOMETRIC DECAY")
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
        # Test A: Full Spectral Sweep of the Sign Ratio epsilon_l = N_l / P_l
        # ---------------------------------------------------------------------
        print("  [Test A] Full Spectral Sweep of Sign Ratio epsilon_l across l in {1, ..., N-1}:")
        print(f"    {'l':>2s} | {'Positive Sum P_l':>20s} | {'Negative Sum N_l':>20s} | {'epsilon_l = N_l/P_l':>20s} | {'(1 - eps)^{-1}':>16s}")
        print("    " + "-" * 88)

        eps_list: list[mp.mpf] = []
        eps_records: list[dict] = []
        max_eps = mp.mpf(0)
        max_eps_l = -1

        for l in range(1, N):
            rec_eps = evaluate_sign_ratio_full(l, zeros_z[l], evals_e, d_k_sq_list)
            eps_records.append(rec_eps)
            eps_val = rec_eps['epsilon']
            eps_list.append(eps_val)
            if eps_val > max_eps:
                max_eps = eps_val
                max_eps_l = l

            # Print representative modes (first 2, middle, and last 2)
            if l <= 2 or l == N // 2 or l >= N - 2:
                print(f"    {l:2d} | {mp.nstr(rec_eps['pos_sum'], 6):>20s} | {mp.nstr(rec_eps['neg_sum'], 6):>20s} | {mp.nstr(eps_val, 6):>20s} | {mp.nstr(rec_eps['correction_factor'], 6):>16s}")
            elif l == 3 and N > 6:
                print(f"    {'...':>2s} | {'...':>20s} | {'...':>20s} | {'...':>20s} | {'...':>16s}")

        synthesis_table1_records.append({
            'N': N,
            'max_eps': max_eps,
            'max_eps_l': max_eps_l,
            'eps_1': eps_records[0]['epsilon'],
            'eps_mid': eps_records[N // 2 - 1]['epsilon'],
            'eps_last': eps_records[-1]['epsilon'],
            'is_all_less_than_one': all(e < mp.mpf(1) for e in eps_list),
        })

        # ---------------------------------------------------------------------
        # Test B: Upper-Edge Mode Forensics for l = N - 1
        # ---------------------------------------------------------------------
        l_top = N - 1
        rec_top_eps = eps_records[-1]
        delta_top = zeros_z[l_top] - evals_e[l_top]
        Delta_top = evals_e[l_top + 1] - evals_e[l_top]
        alpha_top = d_k_sq_list[l_top + 1] / d_k_sq_list[l_top]
        corr_scale_top = Delta_top / (alpha_top * (mp.mpf(1) - rec_top_eps['epsilon']))
        bare_scale_top = Delta_top / alpha_top
        ratio_r_corr_top = delta_top / corr_scale_top
        ratio_r_bare_top = delta_top / bare_scale_top

        print(f"  [Test B] Upper-Edge Mode Forensics for l = N - 1 = {l_top}:")
        print(f"    Delta_{{N-1}} = E_N - E_{{N-1}}      = {mp.nstr(Delta_top, 8)}")
        print(f"    delta_{{N-1}} = z_{{N-1}}^* - E_{{N-1}}  = {mp.nstr(delta_top, 8)}")
        print(f"    d_{{N-1}}^2, d_N^2                = {mp.nstr(d_k_sq_list[l_top], 6)}, {mp.nstr(d_k_sq_list[N], 6)}")
        print(f"    alpha_{{N-1}} = d_N^2 / d_{{N-1}}^2     = {mp.nstr(alpha_top, 8)}")
        print(f"    epsilon_{{N-1}} = N_{{N-1}} / P_{{N-1}} = {mp.nstr(rec_top_eps['epsilon'], 8)}")
        print(f"    (1 - epsilon_{{N-1}})^{{-1}}          = {mp.nstr(rec_top_eps['correction_factor'], 8)}")
        print(f"    delta_{{N-1}} / (Delta / alpha)      = {mp.nstr(ratio_r_bare_top, 6)}  (Bare ratio)")
        print(f"    delta_{{N-1}} / [Delta / (alpha(1-eps))] = {mp.nstr(ratio_r_corr_top, 6)}  (Corrected ratio)")

        synthesis_table2_records.append({
            'N': N,
            'l_top': l_top,
            'alpha_top': alpha_top,
            'eps_top': rec_top_eps['epsilon'],
            'corr_factor_top': rec_top_eps['correction_factor'],
            'ratio_bare_top': ratio_r_bare_top,
            'ratio_corr_top': ratio_r_corr_top,
        })

        # ---------------------------------------------------------------------
        # Test C: Corrected Closed Remote Pairwise Bound Audit for Mode j = 2
        # ---------------------------------------------------------------------
        print(f"  [Test C] Corrected Closed Remote Pairwise Bound Audit for Mode j = 2:")
        print(f"    {'l':>2s} | {'omega - 1 (Actual)':>20s} | {'eta (Bare)':>18s} | {'eta (Corrected)':>18s} | {'Bare?':>6s} | {'Corr?':>6s} | {'Corr Slack':>12s}")
        print("    " + "-" * 95)

        pw_corr_records: list[dict] = []
        all_bare_valid = True
        all_corr_valid = True

        for l in range(N):
            if abs(l - 2) >= 2:
                eps_val = eps_records[l - 1]['epsilon'] if l > 0 else mp.mpf(0)
                rec_c = evaluate_corrected_closed_remote_factor(2, l, evals_e, zeros_z, d_k_sq_list, eps_val)
                pw_corr_records.append(rec_c)

                if not rec_c['is_bare_valid']:
                    all_bare_valid = False
                if not rec_c['is_corr_valid']:
                    all_corr_valid = False

                bare_str = "YES" if rec_c['is_bare_valid'] else "NO"
                corr_str = "YES" if rec_c['is_corr_valid'] else "NO"

                dist = rec_c['dist_index']
                if l == 0 or abs(dist) <= 3 or l == N - 1:
                    print(f"    {l:2d} | {mp.nstr(rec_c['dev_direct'], 6):>20s} | {mp.nstr(rec_c['eta_bare'], 6):>18s} | {mp.nstr(rec_c['eta_corr'], 6):>18s} | {bare_str:>6s} | {corr_str:>6s} | {mp.nstr(rec_c['slack_corr'], 4):>12s}")
                elif abs(dist) == 4:
                    print(f"    {'...':>2s} | {'...':>20s} | {'...':>18s} | {'...':>18s} | {'...':>6s} | {'...':>6s} | {'...':>12s}")

        rec_edge = pw_corr_records[-1]
        synthesis_table3_records.append({
            'N': N,
            'l_edge': rec_edge['l'],
            'dev_edge': rec_edge['dev_direct'],
            'eta_bare_edge': rec_edge['eta_bare'],
            'eta_corr_edge': rec_edge['eta_corr'],
            'is_bare_valid_edge': rec_edge['is_bare_valid'],
            'is_corr_valid_edge': rec_edge['is_corr_valid'],
            'all_corr_valid': all_corr_valid,
        })

        # ---------------------------------------------------------------------
        # Test D: Boundary-Weight Ladder Geometric Decay Audit
        # ---------------------------------------------------------------------
        print("  [Test D] Boundary-Weight Ladder Ratio r_l and Geometric Decay Rate q_l:")
        print(f"    {'l':>2s} | {'d_l^2':>18s} | {'r_l = d_{l-1}^2 / d_l^2':>24s} | {'Max Rate q_l':>16s}")
        print("    " + "-" * 68)

        ladder_records: list[dict] = []
        max_ladder_q = mp.mpf(0)
        for l in range(1, N + 1):
            rec_l = evaluate_ladder_geometric_decay(l, d_k_sq_list)
            ladder_records.append(rec_l)
            if rec_l['q_l'] > max_ladder_q:
                max_ladder_q = rec_l['q_l']

            if l <= 3 or l == (N + 1) // 2 or l >= N:
                print(f"    {l:2d} | {mp.nstr(d_k_sq_list[l], 6):>18s} | {mp.nstr(rec_l['r_l'], 6):>24s} | {mp.nstr(rec_l['q_l'], 6):>16s}")
            elif l == 4 and N > 6:
                print(f"    {'...':>2s} | {'...':>18s} | {'...':>24s} | {'...':>16s}")

        print()

        synthesis_table4_records.append({
            'N': N,
            'r_1': ladder_records[0]['r_l'],
            'r_2': ladder_records[1]['r_l'],
            'r_mid': ladder_records[N // 2 - 1]['r_l'],
            'r_top': ladder_records[-1]['r_l'],
            'max_ladder_q': max_ladder_q,
        })

    # -------------------------------------------------------------------------
    # Multi-Dimension Synthesis Tables
    # -------------------------------------------------------------------------
    print("=" * 105)
    print("SYNTHESIS TABLE 1: SPECTRAL-WIDE SIGN RATIO epsilon_l = N_l / P_l AUDIT (TEST A)")
    print("=" * 105)
    print(f"{'N':>3s} | {'epsilon_1':>16s} | {'epsilon_{mid}':>16s} | {'epsilon_{N-1}':>16s} | {'max_l epsilon_l':>18s} | {'at l':>5s} | {'All eps < 1?':>14s}")
    print("-" * 105)
    for rec in synthesis_table1_records:
        all_less_one = "YES" if rec['is_all_less_than_one'] else "NO"
        print(f"{rec['N']:3d} | {mp.nstr(rec['eps_1'], 4):>16s} | {mp.nstr(rec['eps_mid'], 4):>16s} | {mp.nstr(rec['eps_last'], 4):>16s} | {mp.nstr(rec['max_eps'], 4):>18s} | {rec['max_eps_l']:5d} | {all_less_one:>14s}")
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 2: UPPER-EDGE MODE FORENSICS FOR l = N - 1 (TEST B)")
    print("=" * 105)
    print(f"{'N':>3s} | {'alpha_{N-1}':>16s} | {'epsilon_{N-1}':>16s} | {'(1 - eps)^{-1}':>16s} | {'Bare Ratio':>16s} | {'Corr Ratio':>16s}")
    print("-" * 105)
    for rec in synthesis_table2_records:
        print(f"{rec['N']:3d} | {mp.nstr(rec['alpha_top'], 4):>16s} | {mp.nstr(rec['eps_top'], 4):>16s} | {mp.nstr(rec['corr_factor_top'], 4):>16s} | {mp.nstr(rec['ratio_bare_top'], 4):>16s} | {mp.nstr(rec['ratio_corr_top'], 4):>16s}")
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 3: UPPER-EDGE CLOSED REMOTE BOUND VALIDITY FOR MODE j = 2 (TEST C)")
    print("=" * 105)
    print(f"{'N':>3s} | {'l_edge':>6s} | {'omega - 1 (Actual)':>20s} | {'eta (Bare)':>18s} | {'eta (Corr)':>18s} | {'Bare Valid?':>12s} | {'Corr Valid?':>12s}")
    print("-" * 105)
    for rec in synthesis_table3_records:
        bare_valid_str = "YES" if rec['is_bare_valid_edge'] else "NO (FAILED)"
        corr_valid_str = "YES (RESTORED)" if rec['is_corr_valid_edge'] else "NO"
        print(f"{rec['N']:3d} | {rec['l_edge']:6d} | {mp.nstr(rec['dev_edge'], 4):>20s} | {mp.nstr(rec['eta_bare_edge'], 4):>18s} | {mp.nstr(rec['eta_corr_edge'], 4):>18s} | {bare_valid_str:>12s} | {corr_valid_str:>12s}")
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 4: BOUNDARY-WEIGHT LADDER GEOMETRIC DECAY RATIOS (TEST D)")
    print("=" * 105)
    print(f"{'N':>3s} | {'r_1 = d_0^2/d_1^2':>18s} | {'r_2 = d_1^2/d_2^2':>18s} | {'r_{mid}':>16s} | {'r_N':>16s} | {'Max Rate q_l':>16s}")
    print("-" * 105)
    for rec in synthesis_table4_records:
        print(f"{rec['N']:3d} | {mp.nstr(rec['r_1'], 4):>18s} | {mp.nstr(rec['r_2'], 4):>18s} | {mp.nstr(rec['r_mid'], 4):>16s} | {mp.nstr(rec['r_top'], 4):>16s} | {mp.nstr(rec['max_ladder_q'], 4):>16s}")
    print()

    print("=" * 80)
    print("CELL 80 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell80()
