#!/usr/bin/env python3
"""
================================================================================
CELL 71 — AUDIT OF THE EXACT STIELTJES DERIVATIVE OVERLAP FACTORIZATION
a_j^2 / a_1^2 = (||K u_j||^2 / ||K u_1||^2) * (G_d'(mu_1) / G_d'(mu_j))
AND THREE-FACTOR BOUND-STATE TAIL FILTERING T_j = K_j * G_j * S_j
================================================================================

PURPOSE:
--------
Investigate the exact finite-N Stieltjes derivative factorization of odd source
overlaps and modewise excited tail ratios (Paper NR2 §8.23, Milestone M18):
    a_j^2 / a_1^2 = K_j * G_j
across discrete Galerkin dimensions N in {8, 12, 16, 20, 24}:

1. Exact Stieltjes Derivative Representation:
   Verify the exact finite-N identity for every excited odd mode j >= 1 (Eq. 8.23.1):
       ||K u_j||^2 = a_j^2 * G_d'(mu_j) = a_j^2 * sum_{k=0}^N d_k^2 / (mu_j - E_k)^2,
   checking identity residual | a_j^2 - ||K u_j||^2 / G_d'(mu_j) | < 10^{-40}.

2. Overlap Ratio Factorization:
   Verify the exact finite-N factorization (Eq. 8.23.2):
       a_j^2 / a_1^2 = (||K u_j||^2 / ||K u_1||^2) * (G_d'(mu_1) / G_d'(mu_j))
                     = K_j * G_j,
   confirming complete decoupling from the boundary-layer amplitude D_0.

3. Three-Factor Bound-State Tail Decomposition:
   Audit each individual modewise tail term T_j = b_{0j}^2 / b_{01}^2 via:
       T_j = K_j * G_j * S_j,
   where:
       K_j = ||K u_j||^2 / ||K u_1||^2    (coordinate kinetic energy ratio),
       G_j = G_d'(mu_1) / G_d'(mu_j)      (inverse Stieltjes derivative ratio),
       S_j = ((mu_1 - lam_0)/(mu_j-lam_0))^2 (squared spectral gap suppression).

4. Bound-State Ladder Audit:
   Audit modes j in {2, 3, 4, 5} beneath the effective barrier V_* = 1.0, testing
   whether K_j = O(1), whether G_j remains bounded, and demonstrating that the
   exponential suppression factor S_j dominates to enforce T_j -> 0.

OUTPUT CONSTRAINTS:
-------------------
- Strictly dispassionate, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 71 EXECUTION COMPLETE.
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

def run_cell71() -> None:
    print("=" * 110)
    print("CELL 71 — AUDIT OF THE EXACT STIELTJES DERIVATIVE OVERLAP FACTORIZATION")
    print("         a_j^2 / a_1^2 = (||Ku_j||^2 / ||Ku_1||^2) * (G_d'(mu_1) / G_d'(mu_j))")
    print("         AND THREE-FACTOR BOUND-STATE TAIL FILTERING T_j = K_j * G_j * S_j")
    print("=" * 110)
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

        # Boundary overlaps in even sector: d_k = <u_k^{even}, d>
        d_k_list: list[mp.mpf] = []
        for k in range(N + 1):
            d_val = sum(d_even[m, 0] * V_e[m, k] for m in range(N + 1))
            d_k_list.append(d_val)

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

        # Coordinate kinetic energies: ||K u_j||^2 = sum_{m=1}^N m^2 (u_{j, m-1})^2
        ke_list: list[mp.mpf] = []
        for j in range(N):
            ke_val = sum((m ** 2) * (V_o[m - 1, j] ** 2) for m in range(1, N + 1))
            ke_list.append(ke_val)

        # Stieltjes derivatives: G_d'(mu_j) = sum_{k=0}^N d_k^2 / (mu_j - E_k)^2
        g_prime_list: list[mp.mpf] = []
        stieltjes_residuals: list[mp.mpf] = []
        for j in range(N):
            mu_j = evals_o[j]
            gp = sum((d_k_list[k] ** 2) / ((mu_j - evals_e[k]) ** 2) for k in range(N + 1))
            g_prime_list.append(gp)

            # Check identity: ||K u_j||^2 = a_j^2 * G_d'(mu_j)
            stieltjes_res = abs(ke_list[j] - (a_j_list[j] ** 2) * gp)
            stieltjes_residuals.append(stieltjes_res)

        # Mode 1 reference quantities
        mu_1 = evals_o[1]
        gap_1 = mu_1 - lam_0
        a_1_sq = a_j_list[1] ** 2
        ke_1 = ke_list[1]
        gp_1 = g_prime_list[1]
        b_01_sq = b_0j_list[1] ** 2

        # Detailed modewise breakdown for bound modes j in {2, 3, 4, 5}
        modewise_audit: list[dict] = []
        max_fact_err = mp.mpf(0)
        max_tail_err = mp.mpf(0)

        for j in range(2, N):
            mu_j = evals_o[j]
            gap_j = mu_j - lam_0
            a_j_sq = a_j_list[j] ** 2
            ke_j = ke_list[j]
            gp_j = g_prime_list[j]
            b_j_sq = b_0j_list[j] ** 2

            # Direct ratios
            overlap_ratio_direct = a_j_sq / a_1_sq
            tail_term_direct = b_j_sq / b_01_sq

            # Factorization components
            K_j = ke_j / ke_1
            G_j = gp_1 / gp_j
            S_j = (gap_1 / gap_j) ** 2

            overlap_ratio_fact = K_j * G_j
            tail_term_fact = K_j * G_j * S_j

            fact_err = abs(overlap_ratio_direct - overlap_ratio_fact)
            tail_err = abs(tail_term_direct - tail_term_fact)

            if fact_err > max_fact_err:
                max_fact_err = fact_err
            if tail_err > max_tail_err:
                max_tail_err = tail_err

            if j <= 5:
                modewise_audit.append({
                    "j": j,
                    "mu_j": mu_j,
                    "gap_j": gap_j,
                    "ke_j": ke_j,
                    "gp_j": gp_j,
                    "O_j_direct": overlap_ratio_direct,
                    "O_j_fact": overlap_ratio_fact,
                    "K_j": K_j,
                    "G_j": G_j,
                    "S_j": S_j,
                    "T_j_direct": tail_term_direct,
                    "T_j_fact": tail_term_fact,
                    "fact_err": fact_err,
                    "tail_err": tail_err,
                })

        elapsed = time.time() - t0

        record = {
            "N": N,
            "elapsed": elapsed,
            "D_0_sq": D_0_sq,
            "gap_1": gap_1,
            "ke_1": ke_1,
            "gp_1": gp_1,
            "b_01_sq": b_01_sq,
            "max_stieltjes_res": max(stieltjes_residuals[1:]),
            "max_fact_err": max_fact_err,
            "max_tail_err": max_tail_err,
            "modewise": modewise_audit,
        }
        synthesis_records.append(record)

        print(f"--- DIMENSION N = {N:2d} (Elapsed: {elapsed:.2f}s) ---")
        print(f"  Excited gap:         mu_1 - lam_0         = {mp.nstr(gap_1, 9)}")
        print(f"  Mode 1 KE:           ||K u_1||^2          = {mp.nstr(ke_1, 9)}")
        print(f"  Mode 1 Stieltjes:    G_d'(mu_1)           = {mp.nstr(gp_1, 9)}")
        print(f"  Max Stieltjes res:   |||Ku_j||^2 - a^2*G'| = {mp.nstr(record['max_stieltjes_res'], 4)}")
        print(f"  Max factor err:      |a_j^2/a_1^2 - K*G|  = {mp.nstr(max_fact_err, 4)}")
        print(f"  Max tail factor err: |T_j - K*G*S|        = {mp.nstr(max_tail_err, 4)}")

        for m_rec in modewise_audit:
            j = m_rec["j"]
            print(
                f"    Mode j = {j}: O_{j} = {mp.nstr(m_rec['O_j_direct'], 7)} = "
                f"K_{j} ({mp.nstr(m_rec['K_j'], 5)}) * G_{j} ({mp.nstr(m_rec['G_j'], 5)})  |  "
                f"S_{j} = {mp.nstr(m_rec['S_j'], 6)}  ==>  T_{j} = {mp.nstr(m_rec['T_j_direct'], 6)}"
            )
        print()

    # -------------------------------------------------------------------------
    # Multi-Dimension Synthesis Tables
    # -------------------------------------------------------------------------
    print("=" * 115)
    print("SYNTHESIS TABLE 1: OVERLAP RATIO FACTORIZATION a_j^2 / a_1^2 = K_j * G_j FOR BOUND MODES (j = 2, 3)")
    print("=" * 115)
    header1 = (
        f"{'N':>3} | {'O_2 = a_2^2/a_1^2':>18} | {'K_2 = ||Ku_2||^2/||Ku_1||^2':>27} | "
        f"{'G_2 = G\'(1)/G\'(2)':>20} | {'O_3 = a_3^2/a_1^2':>18} | {'K_3':>10} | {'G_3':>10}"
    )
    print(header1)
    print("-" * len(header1))

    for rec in synthesis_records:
        m_map = {m["j"]: m for m in rec["modewise"]}
        o2_str = mp.nstr(m_map[2]["O_j_direct"], 6) if 2 in m_map else "N/A"
        k2_str = mp.nstr(m_map[2]["K_j"], 6) if 2 in m_map else "N/A"
        g2_str = mp.nstr(m_map[2]["G_j"], 6) if 2 in m_map else "N/A"
        o3_str = mp.nstr(m_map[3]["O_j_direct"], 6) if 3 in m_map else "N/A"
        k3_str = mp.nstr(m_map[3]["K_j"], 6) if 3 in m_map else "N/A"
        g3_str = mp.nstr(m_map[3]["G_j"], 6) if 3 in m_map else "N/A"

        row = (
            f"{rec['N']:3d} | "
            f"{o2_str:>18} | "
            f"{k2_str:>27} | "
            f"{g2_str:>20} | "
            f"{o3_str:>18} | "
            f"{k3_str:>10} | "
            f"{g3_str:>10}"
        )
        print(row)

    print()
    print("=" * 115)
    print("SYNTHESIS TABLE 2: THREE-FACTOR TAIL TERM DECOMPOSITION T_j = K_j * G_j * S_j (j = 2, 3)")
    print("=" * 115)
    header2 = (
        f"{'N':>3} | {'T_2 (direct)':>14} | {'K_2 * G_2 * S_2':>16} | {'S_2 (gap ratio)':>16} | "
        f"{'T_3 (direct)':>14} | {'K_3 * G_3 * S_3':>16} | {'S_3 (gap ratio)':>16}"
    )
    print(header2)
    print("-" * len(header2))

    for rec in synthesis_records:
        m_map = {m["j"]: m for m in rec["modewise"]}
        t2_dir = mp.nstr(m_map[2]["T_j_direct"], 6) if 2 in m_map else "N/A"
        t2_fact = mp.nstr(m_map[2]["T_j_fact"], 6) if 2 in m_map else "N/A"
        s2_str = mp.nstr(m_map[2]["S_j"], 6) if 2 in m_map else "N/A"
        t3_dir = mp.nstr(m_map[3]["T_j_direct"], 6) if 3 in m_map else "N/A"
        t3_fact = mp.nstr(m_map[3]["T_j_fact"], 6) if 3 in m_map else "N/A"
        s3_str = mp.nstr(m_map[3]["S_j"], 6) if 3 in m_map else "N/A"

        row = (
            f"{rec['N']:3d} | "
            f"{t2_dir:>14} | "
            f"{t2_fact:>16} | "
            f"{s2_str:>16} | "
            f"{t3_dir:>14} | "
            f"{t3_fact:>16} | "
            f"{s3_str:>16}"
        )
        print(row)

    print()
    print("=" * 80)
    print("CELL 71 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell71()
