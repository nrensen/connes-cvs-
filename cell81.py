#!/usr/bin/env python3
"""
================================================================================
CELL 81 — TRI-PARTITE SPECTRAL DECOMPOSITION OF THE REMOTE OUTER FACTOR
         Pi_{j, remote}, WEIGHTED REMOTE SUM S_j(N), FOUR-FACTOR SUPPRESSION
         AUDIT, AND ASYMPTOTIC SCALING
================================================================================

Strategic Roadmap Milestone M27 (Paper NR2 Section 8.25):
--------------------------------------------------------
Following the findings of Cell 80 (Milestone M26), which certified the corrected
closed remote pairwise bound:
    omega_{j, l} - 1 < B_{j, l} = [ Delta_j * Delta_l ] / [ alpha_l * (1 - epsilon_l) * |E_j - z_l^*| * |E_{j+1} - E_l| ]
across all remote modes, but empirically refuted the existence of a uniform global
sign bound epsilon_l <= epsilon_* < 1 (max_l epsilon_l -> 0.9993 at N=24) and a global
weight ladder (r_{12} = 4.46), this cell audits the Tri-Partite Spectral Architecture
and the Weighted Remote Sum:
    S_j(N) = sum_{l notin {j, j+1}} B_{j, l} = S_low(N) + S_bulk(N) + S_edge(N).

ANALYTICAL MOTIVATION & OBJECTIVES:
-----------------------------------
1. Weighted Remote Sum Reduction:
   Since:
       log Pi_{j, remote} = sum_{l notin {j, j+1}} log omega_{j, l} <= sum_{l notin {j, j+1}} (omega_{j, l} - 1) < S_j(N),
   proving S_j(N) = O(1) for fixed low j rigorously secures Pi_{j, remote} = O(1).

2. Three Spectral Zones:
   For fixed mode j = 2:
   - Low-Energy Zone (L_low): Modes l <= 4 (remote l in {0, 1, 4}).
     Operated by boundary-weight suppression alpha_l >> 1 and epsilon_l << 10^{-5}.
   - Spectral Bulk Zone (L_bulk): Modes 5 <= l <= N - 3.
     Operated by quadratic spectral separation D_{j, l}^{-2} ~ (E_l - E_j)^{-2}.
   - Upper-Edge Zone (L_edge): Top two modes l in {N - 2, N - 1}.
     Operated by macroscopic geometric separation D_{j, l}^2 ~ E_N^2 >> 1, which
     completely overwhelms the upper-edge inflation factor (1 - epsilon_l)^{-1} ~ 10^3.

3. Four-Factor Balance in B_{j, l}:
   Decomposes B_{j, l} into four explicit factors:
   - Gap factor:           F_gap(l)    = Delta_j * Delta_l
   - Weight factor:        F_weight(l) = alpha_l^{-1} = d_l^2 / d_{l+1}^2
   - Sign factor:          F_sign(l)   = (1 - epsilon_l)^{-1}
   - Distance factor:      F_geom(l)   = [ |E_j - z_l^*| * |E_{j+1} - E_l| ]^{-1}
   demonstrating how each zone's specific suppression mechanism controls B_{j, l}.

THE FOUR INVESTIGATIVE TESTS OF CELL 81:
----------------------------------------
- Test A: Aggregate Remote Sum S_2(N) vs Direct Deviations and log Pi_{2, remote}.
- Test B: Three-Zone Decomposition S_2(N) = S_low + S_bulk + S_edge and Modal Shares.
- Test C: Four-Factor Suppression Audit across Spectral Zones at N = 24.
- Test D: Asymptotic Scaling and Boundedness of S_2(N) across Dimensions N in {8, 12, 16, 20, 24}.

OUTPUT CONSTRAINTS:
-------------------
- Dispassionate, dry, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 81 EXECUTION COMPLETE.
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
# Stieltjes Sign Ratio & Remote Factor Routines
# -----------------------------------------------------------------------------

def evaluate_sign_ratio(
    l: int,
    z_l: mp.mpf,
    evals_e: list[mp.mpf],
    d_k_sq: list[mp.mpf],
) -> tuple[mp.mpf, mp.mpf, mp.mpf, mp.mpf]:
    """
    Computes P_l, N_l, epsilon_l = N_l / P_l, and (1 - epsilon_l)^{-1}.
    For l = 0, N_0 = 0 identically, so epsilon_0 = 0 and (1 - epsilon_0)^{-1} = 1.
    """
    if l == 0:
        return mp.mpf(1), mp.mpf(0), mp.mpf(0), mp.mpf(1)

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

    return pos_sum, neg_sum, epsilon, correction_factor


def evaluate_remote_bound_factors(
    j: int,
    l: int,
    evals_e: list[mp.mpf],
    zeros_z: list[mp.mpf],
    d_k_sq: list[mp.mpf],
    epsilon_l: mp.mpf,
) -> dict:
    """
    Computes the four factors of B_{j, l} and compares B_{j, l} with omega_{j, l} - 1.
    """
    E_j = evals_e[j]
    E_jp1 = evals_e[j + 1]
    Delta_j = E_jp1 - E_j

    E_l = evals_e[l]
    E_lp1 = evals_e[l + 1]
    Delta_l = E_lp1 - E_l

    alpha_l = d_k_sq[l + 1] / d_k_sq[l]
    z_l = zeros_z[l]

    # Actual pairwise factor
    ratio_zero = abs(E_jp1 - z_l) / abs(E_j - z_l)
    ratio_eval = abs(E_j - E_l) / abs(E_jp1 - E_l)
    omega_direct = ratio_zero * ratio_eval
    dev_direct = omega_direct - mp.mpf(1)

    # Denominator products
    denom_exact = abs(E_j - z_l) * abs(E_jp1 - E_l)
    D_jl_cons = min(abs(E_l - E_j), abs(E_l - E_jp1))

    # Four explicit factors
    f_gap = Delta_j * Delta_l
    f_weight = mp.mpf(1) / alpha_l
    f_sign = (mp.mpf(1) / (mp.mpf(1) - epsilon_l)) if epsilon_l < mp.mpf(1) else mp.mpf('inf')
    f_geom = mp.mpf(1) / denom_exact
    f_geom_cons = mp.mpf(1) / (D_jl_cons ** 2)

    # Upper bound B_{j, l}
    B_jl = f_gap * f_weight * f_sign * f_geom
    B_jl_cons = f_gap * f_weight * f_sign * f_geom_cons

    slack = (B_jl / dev_direct) if dev_direct > 0 else mp.mpf(0)
    is_valid = (dev_direct < B_jl)

    return {
        'j': j,
        'l': l,
        'omega_direct': omega_direct,
        'dev_direct': dev_direct,
        'log_omega': mp.log(omega_direct) if omega_direct > 0 else mp.mpf(0),
        'B_jl': B_jl,
        'B_jl_cons': B_jl_cons,
        'is_valid': is_valid,
        'slack': slack,
        'f_gap': f_gap,
        'f_weight': f_weight,
        'f_sign': f_sign,
        'f_geom': f_geom,
        'alpha_l': alpha_l,
        'epsilon_l': epsilon_l,
        'denom_exact': denom_exact,
        'Delta_l': Delta_l,
    }


# -----------------------------------------------------------------------------
# Main Execution Function
# -----------------------------------------------------------------------------

def run_cell81() -> None:
    print("=" * 105)
    print("CELL 81 — TRI-PARTITE SPECTRAL DECOMPOSITION OF Pi_{j, remote}, WEIGHTED REMOTE SUM S_j(N),")
    print("         FOUR-FACTOR SUPPRESSION AUDIT, AND ASYMPTOTIC SCALING")
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

        # Find Stieltjes zeros z_j^*
        zeros_z: list[mp.mpf] = []
        for j in range(N):
            z_j, _ = find_stieltjes_zero_bisection(j, evals_e, d_k_sq_list)
            zeros_z.append(z_j)

        # Compute sign ratios epsilon_l for all l in {0, ..., N-1}
        eps_list: list[mp.mpf] = []
        corr_factor_list: list[mp.mpf] = []
        for l in range(N):
            _, _, eps_val, corr_factor = evaluate_sign_ratio(l, zeros_z[l], evals_e, d_k_sq_list)
            eps_list.append(eps_val)
            corr_factor_list.append(corr_factor)

        # Define spectral zone partitions for remote modes (l notin {2, 3})
        # Low zone: l <= 4 (remote l in {0, 1, 4})
        # Edge zone: l in {N-2, N-1}
        # Bulk zone: 5 <= l <= N-3
        low_modes = [l for l in range(N) if l not in {2, 3} and l <= 4]
        edge_modes = [l for l in range(N) if l not in {2, 3} and l >= N - 2]
        bulk_modes = [l for l in range(N) if l not in {2, 3} and 5 <= l <= N - 3]
        all_remote_modes = low_modes + bulk_modes + edge_modes

        print(f"--- DIMENSION N = {N:2d} (Elapsed: {time.time() - t0:.2f}s) ---")
        print(f"  Partition: Low ({len(low_modes)} modes) = {low_modes}, Bulk ({len(bulk_modes)} modes), Edge ({len(edge_modes)} modes) = {edge_modes}")

        # ---------------------------------------------------------------------
        # Test A: Remote Factors & Weighted Remote Sum S_2(N)
        # ---------------------------------------------------------------------
        remote_records: list[dict] = []
        for l in all_remote_modes:
            rec = evaluate_remote_bound_factors(target_j, l, evals_e, zeros_z, d_k_sq_list, eps_list[l])
            remote_records.append(rec)

        sum_dev = sum(r['dev_direct'] for r in remote_records)
        sum_log_omega = sum(r['log_omega'] for r in remote_records)
        prod_omega = mp.exp(sum_log_omega)
        S_2 = sum(r['B_jl'] for r in remote_records)
        S_2_cons = sum(r['B_jl_cons'] for r in remote_records)

        ratio_S_dev = S_2 / sum_dev if sum_dev > 0 else mp.mpf(0)
        ratio_S_log = S_2 / sum_log_omega if sum_log_omega > 0 else mp.mpf(0)
        all_valid = all(r['is_valid'] for r in remote_records)

        print(f"  [Test A] Weighted Remote Sum S_2(N) vs Direct Quantities:")
        print(f"    sum (omega - 1)          = {mp.nstr(sum_dev, 8)}")
        print(f"    log Pi_{{2, remote}}       = {mp.nstr(sum_log_omega, 8)}")
        print(f"    Pi_{{2, remote}}           = {mp.nstr(prod_omega, 8)}")
        print(f"    S_2(N) (Weighted Sum)    = {mp.nstr(S_2, 8)}")
        print(f"    S_2^{{cons}}(N) (Dist Bound) = {mp.nstr(S_2_cons, 8)}")
        print(f"    Ratio S_2 / sum(omega-1) = {mp.nstr(ratio_S_dev, 6)}")
        print(f"    Ratio S_2 / log(Pi)      = {mp.nstr(ratio_S_log, 6)}")
        print(f"    All Terms B_{{2, l}} > dev? = {'YES' if all_valid else 'NO'}")

        # Print representative remote modes table
        print(f"    {'l':>2s} | {'omega - 1 (Actual)':>20s} | {'B_{2, l} (Bound)':>20s} | {'Slack Ratio':>12s} | {'Valid?':>6s}")
        print("    " + "-" * 72)
        for r in remote_records:
            l = r['l']
            if l <= 1 or l == 4 or l == 5 or l >= N - 2:
                print(f"    {l:2d} | {mp.nstr(r['dev_direct'], 6):>20s} | {mp.nstr(r['B_jl'], 6):>20s} | {mp.nstr(r['slack'], 4):>12s} | {'YES' if r['is_valid'] else 'NO':>6s}")
            elif l == 6 and N > 8:
                print(f"    {'...':>2s} | {'...':>20s} | {'...':>20s} | {'...':>12s} | {'...':>6s}")

        table1_records.append({
            'N': N,
            'sum_dev': sum_dev,
            'sum_log_omega': sum_log_omega,
            'prod_omega': prod_omega,
            'S_2': S_2,
            'S_2_cons': S_2_cons,
            'ratio_S_dev': ratio_S_dev,
            'ratio_S_log': ratio_S_log,
            'all_valid': all_valid,
        })

        # ---------------------------------------------------------------------
        # Test B: Three-Zone Decomposition of S_2(N)
        # ---------------------------------------------------------------------
        rec_by_l = {r['l']: r for r in remote_records}

        S_low = sum(rec_by_l[l]['B_jl'] for l in low_modes)
        S_bulk = sum(rec_by_l[l]['B_jl'] for l in bulk_modes)
        S_edge = sum(rec_by_l[l]['B_jl'] for l in edge_modes)

        dev_low = sum(rec_by_l[l]['dev_direct'] for l in low_modes)
        dev_bulk = sum(rec_by_l[l]['dev_direct'] for l in bulk_modes)
        dev_edge = sum(rec_by_l[l]['dev_direct'] for l in edge_modes)

        pct_low = (S_low / S_2) * 100 if S_2 > 0 else mp.mpf(0)
        pct_bulk = (S_bulk / S_2) * 100 if S_2 > 0 else mp.mpf(0)
        pct_edge = (S_edge / S_2) * 100 if S_2 > 0 else mp.mpf(0)

        print(f"  [Test B] Three-Zone Decomposition of S_2(N) = S_low + S_bulk + S_edge:")
        print(f"    S_low  (modes {low_modes}) = {mp.nstr(S_low, 8):>16s}  ({mp.nstr(pct_low, 4)}%)  [dev = {mp.nstr(dev_low, 6)}]")
        print(f"    S_bulk (modes 5..{N-3 if bulk_modes else 'none'})  = {mp.nstr(S_bulk, 8):>16s}  ({mp.nstr(pct_bulk, 4)}%)  [dev = {mp.nstr(dev_bulk, 6)}]")
        print(f"    S_edge (modes {edge_modes}) = {mp.nstr(S_edge, 8):>16s}  ({mp.nstr(pct_edge, 4)}%)  [dev = {mp.nstr(dev_edge, 6)}]")

        table2_records.append({
            'N': N,
            'S_2': S_2,
            'S_low': S_low,
            'S_bulk': S_bulk,
            'S_edge': S_edge,
            'pct_low': pct_low,
            'pct_bulk': pct_bulk,
            'pct_edge': pct_edge,
            'dev_low': dev_low,
            'dev_bulk': dev_bulk,
            'dev_edge': dev_edge,
        })

        # ---------------------------------------------------------------------
        # Test C: Four-Factor Suppression Audit across Zones (record at N=24)
        # ---------------------------------------------------------------------
        if N == 24:
            print(f"  [Test C] Four-Factor Suppression Audit across Spectral Zones (N = 24):")
            print(f"    {'l':>2s} | {'Zone':>6s} | {'alpha_l^{-1}':>14s} | {'(1 - eps)^{-1}':>14s} | {'D_{2, l}^{-2}':>16s} | {'Delta_2 Delta_l':>16s} | {'B_{2, l}':>16s}")
            print("    " + "-" * 98)

            for r in remote_records:
                l = r['l']
                zone = "Low" if l in low_modes else ("Edge" if l in edge_modes else "Bulk")
                table3_records.append({
                    'l': l,
                    'zone': zone,
                    'f_weight': r['f_weight'],
                    'f_sign': r['f_sign'],
                    'f_geom': r['f_geom'],
                    'f_gap': r['f_gap'],
                    'B_jl': r['B_jl'],
                    'dev': r['dev_direct'],
                })

                # Print representative modes
                if l in {0, 1, 4, 5, 8, 12, 16, 20, 22, 23}:
                    print(f"    {l:2d} | {zone:>6s} | {mp.nstr(r['f_weight'], 4):>14s} | {mp.nstr(r['f_sign'], 4):>14s} | {mp.nstr(r['f_geom'], 4):>16s} | {mp.nstr(r['f_gap'], 4):>16s} | {mp.nstr(r['B_jl'], 4):>16s}")

        print()

    # -------------------------------------------------------------------------
    # Test D: Asymptotic Scaling Synthesis
    # -------------------------------------------------------------------------
    for i, rec in enumerate(table1_records):
        N = rec['N']
        S_2 = rec['S_2']
        if i > 0:
            prev_S_2 = table1_records[i - 1]['S_2']
            ratio_N = S_2 / prev_S_2
        else:
            ratio_N = mp.mpf(1)

        table4_records.append({
            'N': N,
            'S_2': S_2,
            'ratio_N': ratio_N,
            'sum_dev': rec['sum_dev'],
            'prod_omega': rec['prod_omega'],
            'ratio_S_dev': rec['ratio_S_dev'],
        })

    # =========================================================================
    # SYNTHESIS TABLE 1: AGGREGATE REMOTE SUMS & BOUNDEDNESS OF S_2(N)
    # =========================================================================
    print("=" * 105)
    print("SYNTHESIS TABLE 1: AGGREGATE REMOTE SUMS & BOUNDEDNESS OF S_2(N) VS log Pi_{2, remote} (TEST A)")
    print("=" * 105)
    print(f"  {'N':>2s} | {'sum (omega - 1)':>18s} | {'log Pi_{2, rem}':>18s} | {'Pi_{2, rem}':>14s} | {'S_2(N) (Bound)':>18s} | {'S_2 / sum(dev)':>14s} | {'Bound Valid?':>12s}")
    print("-" * 105)
    for r in table1_records:
        valid_str = "YES (ALL MODES)" if r['all_valid'] else "NO"
        print(f"  {r['N']:2d} | {mp.nstr(r['sum_dev'], 6):>18s} | {mp.nstr(r['sum_log_omega'], 6):>18s} | {mp.nstr(r['prod_omega'], 8):>14s} | {mp.nstr(r['S_2'], 6):>18s} | {mp.nstr(r['ratio_S_dev'], 4):>14s} | {valid_str:>12s}")

    # =========================================================================
    # SYNTHESIS TABLE 2: THREE-ZONE DECOMPOSITION S_2(N) = S_low + S_bulk + S_edge
    # =========================================================================
    print()
    print("=" * 105)
    print("SYNTHESIS TABLE 2: THREE-ZONE DECOMPOSITION S_2(N) = S_low + S_bulk + S_edge & SHARES (TEST B)")
    print("=" * 105)
    print(f"  {'N':>2s} | {'S_2(N) Total':>16s} | {'S_low (l <= 4)':>16s} | {'Share Low':>10s} | {'S_bulk (5..N-3)':>16s} | {'Share Bulk':>10s} | {'S_edge (top 2)':>16s} | {'Share Edge':>10s}")
    print("-" * 105)
    for r in table2_records:
        print(f"  {r['N']:2d} | {mp.nstr(r['S_2'], 6):>16s} | {mp.nstr(r['S_low'], 6):>16s} | {mp.nstr(r['pct_low'], 3) + '%':>10s} | {mp.nstr(r['S_bulk'], 6):>16s} | {mp.nstr(r['pct_bulk'], 3) + '%':>10s} | {mp.nstr(r['S_edge'], 6):>16s} | {mp.nstr(r['pct_edge'], 3) + '%':>10s}")

    # =========================================================================
    # SYNTHESIS TABLE 3: FOUR-FACTOR MECHANISM AUDIT ACROSS ZONES AT N = 24
    # =========================================================================
    print()
    print("=" * 105)
    print("SYNTHESIS TABLE 3: FOUR-FACTOR MECHANISM AUDIT ACROSS ZONES AT N = 24 (TEST C)")
    print("=" * 105)
    print(f"  {'l':>2s} | {'Zone':>6s} | {'alpha_l^{-1}':>14s} | {'(1 - eps)^{-1}':>14s} | {'D_{2, l}^{-2}':>16s} | {'Delta_2 Delta_l':>16s} | {'B_{2, l} (Bound)':>18s} | {'omega - 1 (Act)':>18s}")
    print("-" * 105)
    for r in table3_records:
        if r['l'] in {0, 1, 4, 5, 8, 12, 16, 20, 22, 23}:
            print(f"  {r['l']:2d} | {r['zone']:>6s} | {mp.nstr(r['f_weight'], 4):>14s} | {mp.nstr(r['f_sign'], 4):>14s} | {mp.nstr(r['f_geom'], 4):>16s} | {mp.nstr(r['f_gap'], 4):>16s} | {mp.nstr(r['B_jl'], 6):>18s} | {mp.nstr(r['dev'], 6):>18s}")

    # =========================================================================
    # SYNTHESIS TABLE 4: ASYMPTOTIC SCALING & SLACK RATIO PROGRESSION (TEST D)
    # =========================================================================
    print()
    print("=" * 105)
    print("SYNTHESIS TABLE 4: ASYMPTOTIC SCALING & SLACK RATIO PROGRESSION ACROSS DIMENSIONS (TEST D)")
    print("=" * 105)
    print(f"  {'N':>2s} | {'S_2(N)':>16s} | {'S_2(N) / S_2(N-4)':>18s} | {'sum (omega - 1)':>18s} | {'Pi_{2, remote}':>14s} | {'Slack S_2 / sum(dev)':>20s}")
    print("-" * 105)
    for r in table4_records:
        ratio_str = f"{mp.nstr(r['ratio_N'], 4)}" if r['N'] > 8 else "— (Baseline)"
        print(f"  {r['N']:2d} | {mp.nstr(r['S_2'], 6):>16s} | {ratio_str:>18s} | {mp.nstr(r['sum_dev'], 6):>18s} | {mp.nstr(r['prod_omega'], 8):>14s} | {mp.nstr(r['ratio_S_dev'], 4):>20s}")

    print()
    print("=" * 80)
    print("CELL 81 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell81()
