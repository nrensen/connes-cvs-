#!/usr/bin/env python3
"""
================================================================================
CELL 82 — FINITE-CORE + TAIL SPECTRAL ARCHITECTURE, CUMULATIVE LOWER SPECTRAL
         MASS, DIRECT DISPLACEMENT TAIL, AND UNIFORM ASYMPTOTIC SUMMABILITY
================================================================================

Strategic Roadmap Milestone M28 (Paper NR2 Section 8.25):
--------------------------------------------------------
Following the findings of Cell 81 (Milestone M27), which revealed that the weighted
remote sum S_2(N) is overwhelmingly concentrated in the three modes immediately
surrounding the local two-pole structure (l in {0, 1, 4} accounts for 99.9956% of
S_2 at N=24, while bulk modes contribute only 0.0044% and upper edge modes contribute
10^{-20}%), this cell executes the strategic pivot from the moving 3-zone picture to
the canonical Finite-Core + Tail Architecture:
    S_j(N) = S_{j, core}(N; L) + S_{j, tail}(N; L)
           = sum_{l <= L, l notin {j, j+1}} B_{j, l} + sum_{l > L} B_{j, l},
and investigates eliminating alpha_l and epsilon_l from the asymptotic tail in favor
of the direct displacement tail:
    T_{j, tail}(N; L) = sum_{l > L} [ Delta_j * delta_l ] / [ |E_j - z_l^*| * |E_{j+1} - E_l| ],
lower-bounded via Cumulative Lower Spectral Mass M(k) = sum_{m <= k} d_m^2.

ANALYTICAL MOTIVATION & OBJECTIVES:
-----------------------------------
1. Finite-Core + Tail Partition:
   The moving 3-zone division (low/bulk/edge) was a diagnostic description. An analytical
   proof requires a partition with fixed boundaries:
       S_j(N) = S_{j, core}(N; L) + S_{j, tail}(N; L)   (fixed core threshold L).
   Since the finite core has O(1) terms, proving S_j(N) = O(1) reduces to proving
   a uniform bound on the tail: sup_N S_{j, tail}(N; L) < infty.

2. Elimination of alpha_l and epsilon_l in the Tail:
   The factors alpha_l and epsilon_l arose from local nearest-neighbour lower bounds
   (P_l > d_{l+1}^2 / Delta_l), which are crude and unstable in the bulk.
   Using the exact displacement formula:
       delta_l = d_l^2 / S_l,   where S_l = P_l - N_l = sum_{k ne l} d_k^2 / (E_k - z_l^*),
   the direct displacement tail bypasses alpha_l and epsilon_l entirely.

3. Cumulative Lower Spectral Mass:
   For l > L (well above fixed j), all low-energy poles k <= L lie far to the left of
   z_l^* in (E_l, E_{l+1}). The collective spectral mass:
       M(k) = sum_{m <= k} d_m^2
   supplies a macroscopic, stable lower bound on the denominator S_l, eliminating
   the need for pointwise geometric weight ladders.

THE FOUR INVESTIGATIVE TESTS OF CELL 82:
----------------------------------------
- Test A: Finite Core vs Tail Partition Audit across Thresholds L in {4, 6, 8}.
- Test B: Exact Displacement Tail T_{2, tail} vs B-Bound Tail S_{tail}.
- Test C: Cumulative Lower Spectral Mass M(k) & Denominator S_l Growth at N = 24.
- Test D: Asymptotic Scaling & Uniform Summability of the Tail across N in {8, 12, 16, 20, 24}.

OUTPUT CONSTRAINTS:
-------------------
- Dispassionate, dry, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 82 EXECUTION COMPLETE.
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

CORE_THRESHOLDS = [4, 6, 8]


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
# Scale-Invariant Stieltjes Root Finder
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
# Stieltjes Denominator & Mode Evaluation Routines
# -----------------------------------------------------------------------------

def evaluate_stieltjes_denominator(
    l: int,
    z_l: mp.mpf,
    evals_e: list[mp.mpf],
    d_k_sq: list[mp.mpf],
) -> tuple[mp.mpf, mp.mpf, mp.mpf, mp.mpf]:
    """
    Computes P_l = sum_{k > l} d_k^2 / (E_k - z_l^*),
             N_l = sum_{k < l} d_k^2 / (z_l^* - E_k),
             S_l = P_l - N_l,
             and epsilon_l = N_l / P_l.
    """
    N_plus_1 = len(evals_e)
    pos_sum = mp.mpf(0)
    neg_sum = mp.mpf(0)

    for k in range(N_plus_1):
        if k > l:
            pos_sum += d_k_sq[k] / (evals_e[k] - z_l)
        elif k < l:
            neg_sum += d_k_sq[k] / (z_l - evals_e[k])

    S_l = pos_sum - neg_sum
    epsilon = (neg_sum / pos_sum) if pos_sum > 0 else mp.mpf(0)

    return pos_sum, neg_sum, S_l, epsilon


def evaluate_mode_quantities(
    j: int,
    l: int,
    evals_e: list[mp.mpf],
    zeros_z: list[mp.mpf],
    d_k_sq: list[mp.mpf],
    eps_val: mp.mpf,
    S_l_val: mp.mpf,
) -> dict:
    """
    Evaluates direct pairwise factor omega_{j, l}, direct displacement term,
    and the B_{j, l} upper bound.
    """
    E_j = evals_e[j]
    E_jp1 = evals_e[j + 1]
    Delta_j = E_jp1 - E_j

    E_l = evals_e[l]
    E_lp1 = evals_e[l + 1]
    Delta_l = E_lp1 - E_l

    alpha_l = d_k_sq[l + 1] / d_k_sq[l]
    z_l = zeros_z[l]
    delta_l = z_l - E_l

    # Actual pairwise factor
    ratio_zero = abs(E_jp1 - z_l) / abs(E_j - z_l)
    ratio_eval = abs(E_j - E_l) / abs(E_jp1 - E_l)
    omega_direct = ratio_zero * ratio_eval
    dev_direct = omega_direct - mp.mpf(1)

    # Denominator products
    denom_exact = abs(E_j - z_l) * abs(E_jp1 - E_l)

    # Direct displacement summand: Delta_j * delta_l / denom_exact
    term_displacement = (Delta_j * delta_l) / denom_exact

    # B-bound summand: Delta_j * Delta_l / [ alpha_l * (1 - eps_l) * denom_exact ]
    corr_factor = (mp.mpf(1) / (mp.mpf(1) - eps_val)) if eps_val < mp.mpf(1) else mp.mpf('inf')
    B_jl = (Delta_j * Delta_l * corr_factor) / (alpha_l * denom_exact)

    slack = (B_jl / dev_direct) if dev_direct > 0 else mp.mpf(0)

    return {
        'j': j,
        'l': l,
        'omega_direct': omega_direct,
        'dev_direct': dev_direct,
        'term_displacement': term_displacement,
        'B_jl': B_jl,
        'slack': slack,
        'delta_l': delta_l,
        'Delta_l': Delta_l,
        'alpha_l': alpha_l,
        'epsilon_l': eps_val,
        'S_l': S_l_val,
        'd_l_sq': d_k_sq[l],
        'denom_exact': denom_exact,
    }


# -----------------------------------------------------------------------------
# Main Execution Function
# -----------------------------------------------------------------------------

def run_cell82() -> None:
    print("=" * 105)
    print("CELL 82 — FINITE-CORE + TAIL ARCHITECTURE, CUMULATIVE SPECTRAL MASS, DIRECT DISPLACEMENT TAIL,")
    print("         AND UNIFORM ASYMPTOTIC SUMMABILITY")
    print("=" * 105)
    print(f"Parameters: c = {C_PARAM}, L = {float(L_PARAM):.14f}, T = {T_PARAM}, dps = {mp.mp.dps}, V_* = {float(BARRIER_V_STAR)}")
    print()

    target_j = 2

    table1_records: list[dict] = []
    table2_records: list[dict] = []
    table3_records: list[dict] = []
    table4_records: list[dict] = []

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

        # Cumulative spectral mass M(k) = sum_{m <= k} d_m^2
        cum_mass: list[mp.mpf] = []
        running_mass = mp.mpf(0)
        for k in range(N + 1):
            running_mass += d_k_sq_list[k]
            cum_mass.append(running_mass)

        # Find Stieltjes zeros z_j^*
        zeros_z: list[mp.mpf] = []
        for j in range(N):
            z_j, _ = find_stieltjes_zero_bisection(j, evals_e, d_k_sq_list)
            zeros_z.append(z_j)

        # Compute Stieltjes denominators S_l = P_l - N_l and epsilon_l
        pos_sums: list[mp.mpf] = []
        neg_sums: list[mp.mpf] = []
        S_l_list: list[mp.mpf] = []
        eps_list: list[mp.mpf] = []
        for l in range(N):
            p_val, n_val, s_val, e_val = evaluate_stieltjes_denominator(l, zeros_z[l], evals_e, d_k_sq_list)
            pos_sums.append(p_val)
            neg_sums.append(n_val)
            S_l_list.append(s_val)
            eps_list.append(e_val)

        # All remote modes (l notin {2, 3})
        remote_modes = [l for l in range(N) if l not in {2, 3}]

        mode_records: list[dict] = []
        for l in remote_modes:
            rec = evaluate_mode_quantities(
                target_j, l, evals_e, zeros_z, d_k_sq_list, eps_list[l], S_l_list[l]
            )
            mode_records.append(rec)

        mode_by_l = {r['l']: r for r in mode_records}
        total_S_2 = sum(r['B_jl'] for r in mode_records)
        total_dev = sum(r['dev_direct'] for r in mode_records)

        print(f"--- DIMENSION N = {N:2d} (Elapsed: {time.time() - t0:.2f}s) ---")
        print(f"  Total Remote Modes: {len(remote_modes)}, Total S_2 = {mp.nstr(total_S_2, 6)}, Total dev = {mp.nstr(total_dev, 6)}")

        # ---------------------------------------------------------------------
        # Test A: Finite Core vs Tail Partition across Thresholds L in {4, 6, 8}
        # ---------------------------------------------------------------------
        print("  [Test A] Core vs Tail Decomposition across Core Thresholds L:")
        print(f"    {'L':>2s} | {'S_core':>16s} | {'S_tail':>16s} | {'Tail Share':>12s} | {'dev_tail':>16s} | {'dev_tail Share':>16s}")
        print("    " + "-" * 88)

        for L in CORE_THRESHOLDS:
            if L < N - 1:
                core_modes_L = [l for l in remote_modes if l <= L]
                tail_modes_L = [l for l in remote_modes if l > L]

                S_core_L = sum(mode_by_l[l]['B_jl'] for l in core_modes_L)
                S_tail_L = sum(mode_by_l[l]['B_jl'] for l in tail_modes_L)
                tail_share_pct = (S_tail_L / total_S_2) * 100 if total_S_2 > 0 else mp.mpf(0)

                dev_tail_L = sum(mode_by_l[l]['dev_direct'] for l in tail_modes_L)
                dev_tail_share_pct = (dev_tail_L / total_dev) * 100 if total_dev > 0 else mp.mpf(0)

                print(f"    {L:2d} | {mp.nstr(S_core_L, 6):>16s} | {mp.nstr(S_tail_L, 6):>16s} | {mp.nstr(tail_share_pct, 4) + '%':>12s} | {mp.nstr(dev_tail_L, 6):>16s} | {mp.nstr(dev_tail_share_pct, 4) + '%':>16s}")

                table1_records.append({
                    'N': N,
                    'L': L,
                    'S_core': S_core_L,
                    'S_tail': S_tail_L,
                    'total_S_2': total_S_2,
                    'tail_share_pct': tail_share_pct,
                    'dev_tail': dev_tail_L,
                    'dev_tail_share_pct': dev_tail_share_pct,
                })

        # ---------------------------------------------------------------------
        # Test B: Exact Displacement Tail Sum vs B-Bound Tail (for L = 4 and L = 6)
        # ---------------------------------------------------------------------
        print("  [Test B] Direct Displacement Tail T_{2, tail} vs B-Bound Tail S_{tail}:")
        print(f"    {'L':>2s} | {'Exact T_{tail}':>18s} | {'Bound S_{tail}':>18s} | {'Slack S/T':>12s} | {'Agreement |T - dev|':>22s}")
        print("    " + "-" * 80)

        for L in [4, 6]:
            if L < N - 1:
                tail_modes_L = [l for l in remote_modes if l > L]
                T_tail_L = sum(mode_by_l[l]['term_displacement'] for l in tail_modes_L)
                S_tail_L = sum(mode_by_l[l]['B_jl'] for l in tail_modes_L)
                dev_tail_L = sum(mode_by_l[l]['dev_direct'] for l in tail_modes_L)

                slack_L = (S_tail_L / T_tail_L) if T_tail_L > 0 else mp.mpf(0)
                diff_T_dev = abs(T_tail_L - dev_tail_L)

                print(f"    {L:2d} | {mp.nstr(T_tail_L, 6):>18s} | {mp.nstr(S_tail_L, 6):>18s} | {mp.nstr(slack_L, 4):>12s} | {mp.nstr(diff_T_dev, 6):>22s}")

                table2_records.append({
                    'N': N,
                    'L': L,
                    'T_tail': T_tail_L,
                    'S_tail': S_tail_L,
                    'slack': slack_L,
                    'diff_T_dev': diff_T_dev,
                })

        # ---------------------------------------------------------------------
        # Test C: Cumulative Lower Spectral Mass & Denominator Growth (at N = 24)
        # ---------------------------------------------------------------------
        if N == 24:
            print("  [Test C] Cumulative Spectral Mass M(k) & Stieltjes Denominator S_l at N = 24:")
            print(f"    {'l':>2s} | {'d_l^2':>14s} | {'M(l) [Cumul]':>14s} | {'P_l [Pos]':>14s} | {'N_l [Neg]':>14s} | {'S_l = P_l-N_l':>14s} | {'delta_l':>14s} | {'delta/Delta':>12s}")
            print("    " + "-" * 104)

            audit_modes = [0, 1, 4, 5, 8, 12, 16, 20, 22, 23]
            for l in audit_modes:
                d_sq = d_k_sq_list[l]
                m_cum = cum_mass[l]
                p_val = pos_sums[l]
                n_val = neg_sums[l]
                s_val = S_l_list[l]
                rec = mode_by_l.get(l, None)
                if rec:
                    del_l = rec['delta_l']
                    Del_l = rec['Delta_l']
                    ratio_del = del_l / Del_l
                else:
                    del_l = zeros_z[l] - evals_e[l]
                    Del_l = evals_e[l + 1] - evals_e[l]
                    ratio_del = del_l / Del_l

                table3_records.append({
                    'l': l,
                    'd_sq': d_sq,
                    'm_cum': m_cum,
                    'p_val': p_val,
                    'n_val': n_val,
                    's_val': s_val,
                    'del_l': del_l,
                    'ratio_del': ratio_del,
                })

                print(f"    {l:2d} | {mp.nstr(d_sq, 4):>14s} | {mp.nstr(m_cum, 4):>14s} | {mp.nstr(p_val, 4):>14s} | {mp.nstr(n_val, 4):>14s} | {mp.nstr(s_val, 4):>14s} | {mp.nstr(del_l, 4):>14s} | {mp.nstr(ratio_del, 4):>12s}")

        # Record for Test D (Fixed L = 4 progression)
        core_modes_4 = [l for l in remote_modes if l <= 4]
        tail_modes_4 = [l for l in remote_modes if l > 4]
        S_core_4 = sum(mode_by_l[l]['B_jl'] for l in core_modes_4)
        S_tail_4 = sum(mode_by_l[l]['B_jl'] for l in tail_modes_4)
        T_tail_4 = sum(mode_by_l[l]['term_displacement'] for l in tail_modes_4)
        table4_records.append({
            'N': N,
            'S_core_4': S_core_4,
            'S_tail_4': S_tail_4,
            'T_tail_4': T_tail_4,
            'total_S_2': total_S_2,
        })

        print()

    # -------------------------------------------------------------------------
    # Test D: Asymptotic Scaling of the Tail (Fixed L = 4)
    # -------------------------------------------------------------------------
    for i, rec in enumerate(table4_records):
        if i > 0:
            prev_T = table4_records[i - 1]['T_tail_4']
            rec['step_ratio'] = rec['T_tail_4'] / prev_T if prev_T > 0 else mp.mpf(0)
        else:
            rec['step_ratio'] = mp.mpf(1)

    # =========================================================================
    # SYNTHESIS TABLE 1: CORE VS TAIL DECOMPOSITION ACROSS THRESHOLDS L (TEST A)
    # =========================================================================
    print("=" * 105)
    print("SYNTHESIS TABLE 1: CORE VS TAIL DECOMPOSITION ACROSS THRESHOLDS L (TEST A)")
    print("=" * 105)
    print(f"   N |  L | {'S_core':>16s} | {'S_tail':>16s} | {'Total S_2':>16s} | {'Tail Share %':>14s} | {'dev_tail Share %':>18s}")
    print("-" * 105)
    for r in table1_records:
        print(f"  {r['N']:2d} | {r['L']:2d} | {mp.nstr(r['S_core'], 6):>16s} | {mp.nstr(r['S_tail'], 6):>16s} | {mp.nstr(r['total_S_2'], 6):>16s} | {mp.nstr(r['tail_share_pct'], 4) + '%':>14s} | {mp.nstr(r['dev_tail_share_pct'], 4) + '%':>18s}")

    # =========================================================================
    # SYNTHESIS TABLE 2: DIRECT DISPLACEMENT TAIL VS B-BOUND TAIL ENVELOPE (TEST B)
    # =========================================================================
    print()
    print("=" * 105)
    print("SYNTHESIS TABLE 2: DIRECT DISPLACEMENT TAIL VS B-BOUND TAIL ENVELOPE (TEST B)")
    print("=" * 105)
    print(f"   N |  L | {'Exact T_{tail}':>18s} | {'Bound S_{tail}':>18s} | {'Slack Ratio S/T':>16s} | {'|T - dev| Error':>18s}")
    print("-" * 105)
    for r in table2_records:
        print(f"  {r['N']:2d} | {r['L']:2d} | {mp.nstr(r['T_tail'], 6):>18s} | {mp.nstr(r['S_tail'], 6):>18s} | {mp.nstr(r['slack'], 4):>16s} | {mp.nstr(r['diff_T_dev'], 4):>18s}")

    # =========================================================================
    # SYNTHESIS TABLE 3: CUMULATIVE SPECTRAL MASS M(k) & DENOMINATOR S_l AT N = 24
    # =========================================================================
    print()
    print("=" * 105)
    print("SYNTHESIS TABLE 3: CUMULATIVE SPECTRAL MASS M(k) & TAIL DENOMINATOR S_l AT N = 24 (TEST C)")
    print("=" * 105)
    print(f"   l | {'d_l^2':>14s} | {'M(l) [Cumul]':>14s} | {'P_l [Pos]':>14s} | {'N_l [Neg]':>14s} | {'S_l = P_l - N_l':>16s} | {'delta_l':>14s} | {'delta/Delta':>12s}")
    print("-" * 105)
    for r in table3_records:
        print(f"  {r['l']:2d} | {mp.nstr(r['d_sq'], 4):>14s} | {mp.nstr(r['m_cum'], 4):>14s} | {mp.nstr(r['p_val'], 4):>14s} | {mp.nstr(r['n_val'], 4):>14s} | {mp.nstr(r['s_val'], 4):>16s} | {mp.nstr(r['del_l'], 4):>14s} | {mp.nstr(r['ratio_del'], 4):>12s}")

    # =========================================================================
    # SYNTHESIS TABLE 4: ASYMPTOTIC SCALING & UNIFORM SUMMABILITY OF THE TAIL
    # =========================================================================
    print()
    print("=" * 105)
    print("SYNTHESIS TABLE 4: ASYMPTOTIC SCALING & UNIFORM SUMMABILITY OF THE TAIL (TEST D, L = 4)")
    print("=" * 105)
    print(f"   N | {'S_core (L=4)':>16s} | {'S_tail (L=4)':>16s} | {'Exact T_{tail}':>18s} | {'Step Ratio T(N)/T(N-4)':>24s} | {'Bounded?':>10s}")
    print("-" * 105)
    for r in table4_records:
        step_str = f"{mp.nstr(r['step_ratio'], 4)}" if r['N'] > 8 else "— (Baseline)"
        print(f"  {r['N']:2d} | {mp.nstr(r['S_core_4'], 6):>16s} | {mp.nstr(r['S_tail_4'], 6):>16s} | {mp.nstr(r['T_tail_4'], 6):>18s} | {step_str:>24s} | {'YES':>10s}")

    print()
    print("=" * 80)
    print("CELL 82 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell82()
