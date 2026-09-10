#!/usr/bin/env python3
"""
================================================================================
CELL 74 — AUDIT OF LOCAL TWO-POLE CLUSTERING ARCHITECTURE, INTERLACING
         NEIGHBOUR GAPS, AND CONSECUTIVE MODAL GROWTH LADDER C_j = H(mu_{j+1})/H(mu_j)
================================================================================

PURPOSE:
--------
Perform a dedicated numerical and structural audit of the local two-pole clustering
architecture of the positive regularized Stieltjes function H(mu) (Paper NR2 §8.24,
Milestone M21).

From Cell 73, the nearest even pole E_2 accounts for only 39.56% of H(mu_2) at N = 24,
while the adjacent higher pole E_3 accounts for 60.44%. Together, the bracketing
pair {E_2, E_3} accounts for 99.9973% of H(mu_2), with all other poles contributing
less than 6.3e-33. Similarly, {E_3, E_4} accounts for 99.996% of H(mu_3).

This script audits the precise mechanics of this two-pole cluster:

1. Exact Interlacing Bracketing & Spectral Distances:
   Audit the strict spectral interlacing:
       E_j < mu_j < E_{j+1}   (j = 1, 2, 3)
   and evaluate:
       delta_L(j) = mu_j - E_j        (left spectral distance, nearest pole),
       delta_R(j) = E_{j+1} - mu_j    (right spectral distance, dominant pole),
       rho_delta(j) = delta_R / delta_L (distance asymmetry ratio).

2. Boundary Weight Ladder:
   Tabulate d_k^2 for k in {0, 1, 2, 3, 4, 5} and compute the boundary weight
   amplification ratio:
       alpha_k = d_{k+1}^2 / d_k^2.
   Demonstrate how boundary weight growth alpha_2 compensates for the larger
   right distance (delta_R)^2 to produce pole asymmetry:
       A(j) = H_{j+1}(mu_j) / H_j(mu_j) = alpha_j * (delta_L / delta_R)^2.

3. Two-Pole Cluster Decomposition & Fidelity:
   Evaluate the two-pole cluster sum:
       H_{j, j+1}(mu_j) = H_j(mu_j) + H_{j+1}(mu_j)
   and audit:
       - Cluster fidelity:  F_{2-pole}(j) = [H_{j, j+1}(mu_j) / H(mu_j)] * 100%,
       - Outer residual:    R_{outer}(j) = H(mu_j) - H_{j, j+1}(mu_j).

4. Consecutive Modewise Growth Ladder:
   Evaluate the ladder ratios between consecutive odd modes:
       C_1(N) = H(mu_2) / H(mu_1),
       C_2(N) = H(mu_3) / H(mu_2),
   and fit consecutive logarithmic rates:
       sigma_{C1} = Delta log(C_1) / Delta N,
       sigma_{C2} = Delta log(C_2) / Delta N.

OUTPUT CONSTRAINTS:
-------------------
- Dispassionate, dry, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 74 EXECUTION COMPLETE.
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

mp.mp.dps = 50

C_PARAM = 13
L_PARAM = mp.log(C_PARAM)
T_PARAM = 400
GROUND_DPS = 50

N_LIST = [8, 12, 16, 20, 24]
BARRIER_V_STAR = mp.mpf('1.0')


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
# Main Execution Function
# -----------------------------------------------------------------------------

def run_cell74() -> None:
    print("=" * 105)
    print("CELL 74 — AUDIT OF LOCAL TWO-POLE CLUSTERING ARCHITECTURE, INTERLACING")
    print("         NEIGHBOUR GAPS, AND CONSECUTIVE MODAL GROWTH LADDER C_j = H(mu_{j+1})/H(mu_j)")
    print("=" * 105)
    print(f"Parameters: c = {C_PARAM}, L = {float(L_PARAM):.14f}, T = {T_PARAM}, dps = {mp.mp.dps}, V_* = {float(BARRIER_V_STAR)}")
    print()

    synthesis_records: list[dict] = []

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

        # Ground state c in R^{N+1}
        c_even = V_e[:, 0]

        # Boundary vector d in R^{N+1}: d_even = (1, sqrt(2), ..., sqrt(2))^T
        d_even = mp.matrix(N + 1, 1)
        d_even[0, 0] = mp.mpf(1)
        for m in range(1, N + 1):
            d_even[m, 0] = mp.sqrt(2)

        D_0 = sum(d_even[m, 0] * c_even[m, 0] for m in range(N + 1))
        D_0_sq = D_0 ** 2

        # Boundary overlaps in even sector: d_k = <u_k^{even}, d>
        d_k_list: list[mp.mpf] = []
        d_k_sq_list: list[mp.mpf] = []
        for k in range(N + 1):
            d_val = sum(d_even[m, 0] * V_e[m, k] for m in range(N + 1))
            d_k_list.append(d_val)
            d_k_sq_list.append(d_val ** 2)

        # Boundary weight amplification ratios: alpha_k = d_{k+1}^2 / d_k^2
        alpha_list: list[mp.mpf] = []
        for k in range(min(5, N)):
            if d_k_sq_list[k] > 0:
                alpha_list.append(d_k_sq_list[k + 1] / d_k_sq_list[k])
            else:
                alpha_list.append(mp.mpf(0))

        # Full term-by-term pole calculations for modes j = 1, 2, 3
        mode_data: dict[int, dict] = {}
        for j in [1, 2, 3]:
            if j >= N:
                continue

            mu_j = evals_o[j]
            gap_j = mu_j - lam_0
            gap_j_sq = gap_j ** 2

            # All pole terms k in 0..N
            H_k_terms = [D_0_sq]
            for k in range(1, N + 1):
                E_k = evals_e[k]
                denom = (E_k - mu_j) ** 2
                term = d_k_sq_list[k] * (gap_j_sq / denom)
                H_k_terms.append(term)

            H_total = sum(H_k_terms)

            # Left and right bracketing neighbours: E_j < mu_j < E_{j+1}
            E_left = evals_e[j]
            E_right = evals_e[j + 1] if j + 1 <= N else mp.mpf(0)

            delta_L = mu_j - E_left
            delta_R = E_right - mu_j if j + 1 <= N else mp.mpf(0)
            rho_delta = delta_R / delta_L if delta_L > 0 else mp.mpf(0)

            # Interlacing check: E_j < mu_j < E_{j+1}
            interlaced = (E_left < mu_j < E_right) if j + 1 <= N else (E_left < mu_j)

            # Two-pole cluster sum
            H_left = H_k_terms[j]
            H_right = H_k_terms[j + 1] if j + 1 <= N else mp.mpf(0)
            H_2pole = H_left + H_right

            # Cluster fidelity and outer residual
            fidelity = (H_2pole / H_total) * 100 if H_total > 0 else mp.mpf(0)
            outer_res = H_total - H_2pole

            # Pole asymmetry ratio
            asymmetry = H_right / H_left if H_left > 0 else mp.mpf(0)

            mode_data[j] = {
                "mu_j": mu_j,
                "gap_j": gap_j,
                "E_left": E_left,
                "E_right": E_right,
                "delta_L": delta_L,
                "delta_R": delta_R,
                "rho_delta": rho_delta,
                "interlaced": interlaced,
                "H_left": H_left,
                "H_right": H_right,
                "H_2pole": H_2pole,
                "H_total": H_total,
                "fidelity": fidelity,
                "outer_res": outer_res,
                "asymmetry": asymmetry,
                "share_left": (H_left / H_total) * 100 if H_total > 0 else mp.mpf(0),
                "share_right": (H_right / H_total) * 100 if H_total > 0 else mp.mpf(0),
            }

        # Consecutive ladder ratios: C_1 = H(mu_2)/H(mu_1), C_2 = H(mu_3)/H(mu_2)
        H_1 = mode_data[1]["H_total"]
        H_2 = mode_data[2]["H_total"]
        H_3 = mode_data[3]["H_total"]

        C_1 = H_2 / H_1
        C_2 = H_3 / H_2
        R_31 = H_3 / H_1

        elapsed = time.time() - t0

        record = {
            "N": N,
            "elapsed": elapsed,
            "D_0_sq": D_0_sq,
            "d_k_sq_list": d_k_sq_list,
            "alpha_list": alpha_list,
            "mode_data": mode_data,
            "C_1": C_1,
            "C_2": C_2,
            "R_31": R_31,
        }
        synthesis_records.append(record)

        print(f"--- DIMENSION N = {N:2d} (Elapsed: {elapsed:.2f}s) ---")
        m2 = mode_data[2]
        print(f"  Mode 2 Interlacing:  E_2 < mu_2 < E_3 = {m2['interlaced']}")
        print(f"    Left dist delta_L  = mu_2 - E_2     = {mp.nstr(m2['delta_L'], 9)}")
        print(f"    Right dist delta_R = E_3 - mu_2     = {mp.nstr(m2['delta_R'], 9)}")
        print(f"    Dist ratio delta_R/delta_L          = {float(m2['rho_delta']):.2f}")
        print(f"    Left pole H_2 (share)               = {mp.nstr(m2['H_left'], 7)} ({float(m2['share_left']):.2f}%)")
        print(f"    Right pole H_3 (share)              = {mp.nstr(m2['H_right'], 7)} ({float(m2['share_right']):.2f}%)")
        print(f"    Two-pole fidelity H_{{2,3}} / H_tot  = {float(m2['fidelity']):.4f}%")
        print(f"    Outer residual H_tot - H_{{2,3}}     = {mp.nstr(m2['outer_res'], 4)}")
        print(f"    Pole asymmetry H_3 / H_2            = {float(m2['asymmetry']):.4f}")
        print(f"  Ladder ratio C_1 = H(mu_2)/H(mu_1)    = {mp.nstr(C_1, 9)}")
        print(f"  Ladder ratio C_2 = H(mu_3)/H(mu_2)    = {mp.nstr(C_2, 9)}")
        print()

    # -------------------------------------------------------------------------
    # Synthesis Tables
    # -------------------------------------------------------------------------

    print("=" * 105)
    print("SYNTHESIS TABLE 1: EXACT INTERLACING BRACKETING AND SPECTRAL DISTANCES")
    print("=" * 105)
    header1 = (
        f"{'N':>3} | {'Mode':>4} | {'Interlaced?':>11} | {'delta_L = mu_j - E_j':>22} | "
        f"{'delta_R = E_{j+1} - mu_j':>24} | {'delta_R / delta_L':>17}"
    )
    print(header1)
    print("-" * len(header1))
    for rec in synthesis_records:
        N = rec["N"]
        for j in [2, 3]:
            mj = rec["mode_data"][j]
            int_str = "YES" if mj["interlaced"] else "NO"
            print(
                f"{N:3d} | {j:4d} | {int_str:>11} | {mp.nstr(mj['delta_L'], 9):>22} | "
                f"{mp.nstr(mj['delta_R'], 9):>24} | {float(mj['rho_delta']):17.4f}"
            )
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 2: BOUNDARY WEIGHT LADDER d_k^2 AND AMPLIFICATION RATIOS alpha_k")
    print("=" * 105)
    header2 = (
        f"{'N':>3} | {'d_0^2':>12} | {'d_1^2':>12} | {'d_2^2':>12} | {'d_3^2':>12} | "
        f"{'alpha_1 = d_2^2/d_1^2':>20} | {'alpha_2 = d_3^2/d_2^2':>20}"
    )
    print(header2)
    print("-" * len(header2))
    for rec in synthesis_records:
        N = rec["N"]
        dsq = rec["d_k_sq_list"]
        a_1 = dsq[2] / dsq[1] if dsq[1] > 0 else mp.mpf(0)
        a_2 = dsq[3] / dsq[2] if dsq[2] > 0 else mp.mpf(0)
        print(
            f"{N:3d} | {mp.nstr(dsq[0], 6):>12} | {mp.nstr(dsq[1], 6):>12} | {mp.nstr(dsq[2], 6):>12} | "
            f"{mp.nstr(dsq[3], 6):>12} | {float(a_1):20.2f} | {float(a_2):20.2f}"
        )
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 3: TWO-POLE CLUSTER DECOMPOSITION AND FIDELITY FOR MODE j = 2")
    print("=" * 105)
    header3 = (
        f"{'N':>3} | {'H_2 (left)':>15} | {'H_3 (right)':>15} | {'H_{2,3} (sum)':>15} | "
        f"{'H_total':>15} | {'Fidelity(%)':>12} | {'H_3/H_2':>10}"
    )
    print(header3)
    print("-" * len(header3))
    for rec in synthesis_records:
        N = rec["N"]
        m2 = rec["mode_data"][2]
        print(
            f"{N:3d} | {mp.nstr(m2['H_left'], 9):>15} | {mp.nstr(m2['H_right'], 9):>15} | "
            f"{mp.nstr(m2['H_2pole'], 9):>15} | {mp.nstr(m2['H_total'], 9):>15} | "
            f"{float(m2['fidelity']):12.6f} | {float(m2['asymmetry']):10.4f}"
        )
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 4: TWO-POLE CLUSTER DECOMPOSITION AND FIDELITY FOR MODE j = 3")
    print("=" * 105)
    header4 = (
        f"{'N':>3} | {'H_3 (left)':>15} | {'H_4 (right)':>15} | {'H_{3,4} (sum)':>15} | "
        f"{'H_total':>15} | {'Fidelity(%)':>12} | {'H_4/H_3':>10}"
    )
    print(header4)
    print("-" * len(header4))
    for rec in synthesis_records:
        N = rec["N"]
        m3 = rec["mode_data"][3]
        print(
            f"{N:3d} | {mp.nstr(m3['H_left'], 9):>15} | {mp.nstr(m3['H_right'], 9):>15} | "
            f"{mp.nstr(m3['H_2pole'], 9):>15} | {mp.nstr(m3['H_total'], 9):>15} | "
            f"{float(m3['fidelity']):12.6f} | {float(m3['asymmetry']):10.4f}"
        )
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 5: CONSECUTIVE MODAL GROWTH LADDER C_j = H(mu_{j+1})/H(mu_j)")
    print("=" * 105)
    header5 = (
        f"{'N':>3} | {'H(mu_1)':>15} | {'H(mu_2)':>15} | {'C_1 = H2/H1':>15} | "
        f"{'H(mu_3)':>15} | {'C_2 = H3/H2':>15} | {'R_31 = H3/H1':>15}"
    )
    print(header5)
    print("-" * len(header5))
    for rec in synthesis_records:
        N = rec["N"]
        H_1 = rec["mode_data"][1]["H_total"]
        H_2 = rec["mode_data"][2]["H_total"]
        H_3 = rec["mode_data"][3]["H_total"]
        C_1 = rec["C_1"]
        C_2 = rec["C_2"]
        R_31 = rec["R_31"]
        print(
            f"{N:3d} | {mp.nstr(H_1, 9):>15} | {mp.nstr(H_2, 9):>15} | {mp.nstr(C_1, 9):>15} | "
            f"{mp.nstr(H_3, 9):>15} | {mp.nstr(C_2, 9):>15} | {mp.nstr(R_31, 9):>15}"
        )
    print()

    print("=" * 105)
    print("SYNTHESIS TABLE 6: CONSECUTIVE LOGARITHMIC RATES FOR LADDER GROWTH RATIOS C_1 AND C_2")
    print("=" * 105)
    header6 = (
        f"{'Interval':>14} | {'Delta N':>7} | {'C_1(start)':>15} | {'C_1(end)':>15} | "
        f"{'sigma_{C1}':>12} | {'C_2(start)':>15} | {'C_2(end)':>15} | {'sigma_{C2}':>12}"
    )
    print(header6)
    print("-" * len(header6))
    for idx in range(len(synthesis_records) - 1):
        r_a = synthesis_records[idx]
        r_b = synthesis_records[idx + 1]
        n_a, n_b = r_a["N"], r_b["N"]
        dn = n_b - n_a

        c1_a = r_a["C_1"]
        c1_b = r_b["C_1"]
        sigma_c1 = (mp.log(c1_b) - mp.log(c1_a)) / dn

        c2_a = r_a["C_2"]
        c2_b = r_b["C_2"]
        sigma_c2 = (mp.log(c2_b) - mp.log(c2_a)) / dn

        interval_str = f"N = {n_a:2d} -> {n_b:2d}"
        print(
            f"{interval_str:>14} | {dn:7d} | {mp.nstr(c1_a, 9):>15} | {mp.nstr(c1_b, 9):>15} | "
            f"{float(sigma_c1):12.5f} | {mp.nstr(c2_a, 9):>15} | {mp.nstr(c2_b, 9):>15} | "
            f"{float(sigma_c2):12.5f}"
        )
    print()

    print("=" * 105)
    print("CELL 74 EXECUTION COMPLETE")
    print("=" * 105)


if __name__ == "__main__":
    run_cell74()
