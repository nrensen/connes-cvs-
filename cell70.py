#!/usr/bin/env python3
"""
================================================================================
CELL 70 — AUDIT OF FIRST-MODE SPECTRAL CONCENTRATION, RELATIVE TAIL RATIO eps_N,
AND MODE-BY-MODE SPECTRAL FILTERING IN THE COORDINATE WAVEPACKET
================================================================================

PURPOSE:
--------
Investigate the exact finite-N relative tail ratio and first-mode spectral
concentration of the excited coordinate wavepacket (Paper NR2 §8.22, Milestone M17):
    v_exc = P_{perp u_0} K c
across discrete Galerkin dimensions N in {8, 12, 16, 20, 24}:

1. Exact Relative Tail Ratio & D_0^2 Cancellation:
   Verify the exact finite-N identity (Proposition 8.22, Eq. 8.22.2):
       eps_N = (sum_{j >= 2} b_{0j}^2) / b_{01}^2
             = sum_{j=2}^{N-1} (a_j^2 / a_1^2) * ((mu_1 - lam_0) / (mu_j - lam_0))^2,
   confirming that the global boundary-layer amplitude D_0^2 cancels identically,
   leaving eps_N governed purely by overlap ratios O_j = a_j^2 / a_1^2 and automatic
   spectral suppression factors S_j = ((mu_1 - lam_0) / (mu_j - lam_0))^2.

2. Residual Equivalence:
   Verify the exact finite-N relation (Proposition 8.22, Eq. 8.22.11):
       ||v_exc||^2 = b_{01}^2 * (1 + eps_N),
   checking identity residual | ||v_exc||^2 - b_{01}^2 * (1 + eps_N) | < 10^{-45}.

3. First-Mode Spectral Concentration:
   Measure the normalized first-mode concentration ratio:
       C_1(N) = b_{01}^2 / ||v_exc||^2 = 1 / (1 + eps_N),
   testing whether C_1(N) -> 1 monotonically (exceeding 99.9995% at N = 24),
   corresponding to weak measure concentration dnu_v / ||v_exc||^2 -> delta_{mu_1}.

4. Two-Sector Tail Decomposition:
   Partition eps_N = eps_N^{bound} + eps_N^{high} at the above-barrier threshold
   V_* = 1.0 (mu_j < V_* for bound modes, mu_j >= V_* for high-energy modes):
       eps_N^{bound} = sum_{2 <= j, mu_j < V_*} T_j,
       eps_N^{high}  = sum_{j, mu_j >= V_*} T_j,
   testing whether eps_N^{high} is exponentially decoupled and negligible compared
   to the bound-state tail.

5. Modewise Spectral Filtering Audit:
   Audit individual tail terms T_j = (b_{0j} / b_{01})^2, overlap ratios O_j = a_j^2 / a_1^2,
   and spectral suppression factors S_j = ((mu_1 - lam_0) / (mu_j - lam_0))^2 for
   modes j in {2, 3, 4, 5}.

OUTPUT CONSTRAINTS:
-------------------
- Strictly dispassionate, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 70 EXECUTION COMPLETE.
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
# Main Diagnostic Execution
# -----------------------------------------------------------------------------

def run_cell70() -> None:
    print("=" * 96)
    print("CELL 70 — AUDIT OF FIRST-MODE SPECTRAL CONCENTRATION, RELATIVE TAIL RATIO eps_N,")
    print("         AND MODE-BY-MODE SPECTRAL FILTERING IN THE COORDINATE WAVEPACKET")
    print("=" * 96)
    print(
        f"Parameters: c = {C_PARAM}, L = {mp.nstr(L_PARAM, 15)}, "
        f"T = {T_PARAM}, dps = {mp.mp.dps}, V_* = {mp.nstr(BARRIER_V_STAR, 4)}"
    )
    print()

    synthesis_records: list[dict] = []

    for N in N_LIST:
        t0 = time.time()

        # Build Galerkin operator (cache-accelerated via cell.py)
        Q_full, _ = get_galerkin_matrix(
            c=C_PARAM,
            N=N,
            T=T_PARAM,
            dps=GROUND_DPS,
            verbose=False,
        )

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

        # Source vector psi = K Q e_0 on R^{2N+1}, mapped to H_odd via O
        # psi_full[N + m] = m * Q_full[N + m, N], psi_full[N - m] = -m * Q_full[N - m, N]
        # In odd basis: psi_odd[m-1] = sqrt(2) * m * Q_full[N + m, N]
        psi_odd = mp.matrix(N, 1)
        for m in range(1, N + 1):
            psi_odd[m - 1, 0] = mp.sqrt(2) * m * Q_full[N + m, N]

        # Mode overlaps in odd sector: a_j = <psi, u_j>
        a_j_list: list[mp.mpf] = []
        for j in range(N):
            a_val = sum(psi_odd[m - 1, 0] * V_o[m - 1, j] for m in range(1, N + 1))
            a_j_list.append(a_val)

        # Coordinate transition dipoles b_{0j} = <c, K u_j>
        b_0j_list: list[mp.mpf] = []
        for j in range(N):
            b_j = sum(c_even[m, 0] * m * V_o[m - 1, j] for m in range(1, N + 1))
            b_0j_list.append(b_j)

        # Mode 1 quantities
        mu_1 = evals_o[1]
        gap_1 = mu_1 - lam_0
        a_1 = a_j_list[1]
        b_01 = b_0j_list[1]
        b_01_sq = b_01 ** 2

        # Excited wavepacket residual norm: ||v_exc||^2 = sum_{j=1}^{N-1} b_{0j}^2
        v_exc_sq = sum(b_0j_list[j] ** 2 for j in range(1, N))

        # Modewise tail terms T_j = b_{0j}^2 / b_{01}^2
        # Overlap ratios O_j = a_j^2 / a_1^2
        # Spectral suppression factors S_j = ((mu_1 - lam_0) / (mu_j - lam_0))^2
        modewise_data: list[dict] = []
        eps_N_direct = mp.mpf(0)
        eps_N_formula = mp.mpf(0)
        eps_N_bound = mp.mpf(0)
        eps_N_high = mp.mpf(0)
        max_id_err = mp.mpf(0)
        n_bound = 0

        for j in range(N):
            if evals_o[j] < BARRIER_V_STAR:
                n_bound += 1

        for j in range(2, N):
            mu_j = evals_o[j]
            gap_j = mu_j - lam_0
            a_j = a_j_list[j]
            b_j = b_0j_list[j]

            t_direct = (b_j ** 2) / b_01_sq
            o_j = (a_j ** 2) / (a_1 ** 2)
            s_j = (gap_1 / gap_j) ** 2
            t_formula = o_j * s_j

            id_err = abs(t_direct - t_formula)
            if id_err > max_id_err:
                max_id_err = id_err

            eps_N_direct += t_direct
            eps_N_formula += t_formula

            if mu_j < BARRIER_V_STAR:
                eps_N_bound += t_direct
            else:
                eps_N_high += t_direct

            if j <= 5:
                modewise_data.append({
                    "j": j,
                    "mu_j": mu_j,
                    "gap_j": gap_j,
                    "T_j": t_direct,
                    "O_j": o_j,
                    "S_j": s_j,
                })

        # Residual equivalence: ||v_exc||^2 = b_{01}^2 * (1 + eps_N)
        v_exc_reconstructed = b_01_sq * (1 + eps_N_direct)
        equiv_residual = abs(v_exc_sq - v_exc_reconstructed)

        # First-mode spectral concentration: C_1 = b_{01}^2 / ||v_exc||^2
        c_1_ratio = b_01_sq / v_exc_sq
        c_1_pct = c_1_ratio * 100

        elapsed = time.time() - t0

        record = {
            "N": N,
            "elapsed": elapsed,
            "D_0_sq": D_0_sq,
            "gap_1": gap_1,
            "b_01_sq": b_01_sq,
            "v_exc_sq": v_exc_sq,
            "eps_N": eps_N_direct,
            "eps_N_bound": eps_N_bound,
            "eps_N_high": eps_N_high,
            "n_bound": n_bound,
            "c_1_ratio": c_1_ratio,
            "c_1_pct": c_1_pct,
            "equiv_residual": equiv_residual,
            "max_id_err": max_id_err,
            "modewise": modewise_data,
        }
        synthesis_records.append(record)

        print(f"--- DIMENSION N = {N:2d} (Elapsed: {elapsed:.2f}s) ---")
        print(f"  Excited gap:         mu_1 - lam_0           = {mp.nstr(gap_1, 9)}")
        print(f"  Dipole squared:      b_{{01}}^2               = {mp.nstr(b_01_sq, 9)}")
        print(f"  Excited norm:        ||v_exc||^2            = {mp.nstr(v_exc_sq, 9)}")
        print(f"  Tail ratio:          eps_N = sum_{{j>=2}} T_j = {mp.nstr(eps_N_direct, 9)}")
        print(f"  Bound tail:          eps_N^{{bound}} (< V_*)  = {mp.nstr(eps_N_bound, 9)}")
        print(f"  High-energy tail:    eps_N^{{high}} (>= V_*)  = {mp.nstr(eps_N_high, 9)}")
        print(f"  First-mode share:    C_1 = b_{{01}}^2/||v||^2  = {mp.nstr(c_1_pct, 10)}%")
        print(f"  Equiv residual:      |||v||^2 - b^2(1+eps)| = {mp.nstr(equiv_residual, 4)}")
        print(f"  Max D_0 cancel err:  |T_dir - O_j * S_j|    = {mp.nstr(max_id_err, 4)}")

        for m_rec in modewise_data:
            j = m_rec["j"]
            print(
                f"    Mode j = {j}: T_{j} = {mp.nstr(m_rec['T_j'], 8)}  "
                f"[O_{j} = {mp.nstr(m_rec['O_j'], 8)}, S_{j} = {mp.nstr(m_rec['S_j'], 8)}]"
            )
        print()

    # -------------------------------------------------------------------------
    # Multi-Dimension Synthesis Tables
    # -------------------------------------------------------------------------
    print("=" * 110)
    print("SYNTHESIS TABLE 1: RELATIVE TAIL RATIO eps_N, SECTOR DECOMPOSITION, AND FIRST-MODE CONCENTRATION")
    print("=" * 110)
    header1 = (
        f"{'N':>3} | {'||v_exc||^2':>14} | {'b_{01}^2':>14} | {'eps_N':>13} | "
        f"{'eps_N^{bound}':>13} | {'eps_N^{high}':>13} | {'C_1 (%)':>11} | {'Equiv Res':>11}"
    )
    print(header1)
    print("-" * len(header1))

    for rec in synthesis_records:
        row = (
            f"{rec['N']:3d} | "
            f"{mp.nstr(rec['v_exc_sq'], 8):>14} | "
            f"{mp.nstr(rec['b_01_sq'], 8):>14} | "
            f"{mp.nstr(rec['eps_N'], 7):>13} | "
            f"{mp.nstr(rec['eps_N_bound'], 7):>13} | "
            f"{mp.nstr(rec['eps_N_high'], 7):>13} | "
            f"{mp.nstr(rec['c_1_pct'], 8):>11} | "
            f"{mp.nstr(rec['equiv_residual'], 4):>11}"
        )
        print(row)

    print()
    print("=" * 110)
    print("SYNTHESIS TABLE 2: MODE-BY-MODE TAIL RATIOS T_j AND SPECTRAL SUPPRESSION FACTORS S_j (j = 2, 3, 4, 5)")
    print("=" * 110)
    header2 = (
        f"{'N':>3} | {'T_2':>12} | {'S_2':>12} | {'T_3':>12} | {'S_3':>12} | "
        f"{'T_4':>12} | {'S_4':>12} | {'T_5':>12} | {'S_5':>12}"
    )
    print(header2)
    print("-" * len(header2))

    for rec in synthesis_records:
        m_map = {m["j"]: m for m in rec["modewise"]}
        t2_str = mp.nstr(m_map[2]["T_j"], 6) if 2 in m_map else "N/A"
        s2_str = mp.nstr(m_map[2]["S_j"], 6) if 2 in m_map else "N/A"
        t3_str = mp.nstr(m_map[3]["T_j"], 6) if 3 in m_map else "N/A"
        s3_str = mp.nstr(m_map[3]["S_j"], 6) if 3 in m_map else "N/A"
        t4_str = mp.nstr(m_map[4]["T_j"], 6) if 4 in m_map else "N/A"
        s4_str = mp.nstr(m_map[4]["S_j"], 6) if 4 in m_map else "N/A"
        t5_str = mp.nstr(m_map[5]["T_j"], 6) if 5 in m_map else "N/A"
        s5_str = mp.nstr(m_map[5]["S_j"], 6) if 5 in m_map else "N/A"

        row = (
            f"{rec['N']:3d} | "
            f"{t2_str:>12} | {s2_str:>12} | "
            f"{t3_str:>12} | {s3_str:>12} | "
            f"{t4_str:>12} | {s4_str:>12} | "
            f"{t5_str:>12} | {s5_str:>12}"
        )
        print(row)

    print()
    print("=" * 80)
    print("CELL 70 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell70()
