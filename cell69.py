#!/usr/bin/env python3
"""
================================================================================
CELL 69 — AUDIT OF THE EXACT TRANSITION DIPOLE FACTORIZATION R_gap = b_{01}^2 / R_1
AND EVEN-RESOLVENT CANCELLATION IN THE FIRST RELATIVE TUNNELING GAP
================================================================================

PURPOSE:
--------
Investigate the exact finite-N algebraic reduction of the first relative tunneling
gap (Paper 4B §8.19, Milestone M14):
    R_gap(N) = D_0^2 / (mu_1 - lambda_0)
across discrete Galerkin dimensions N in {8, 12, 16, 20, 24}:

1. Exact Dipole Factorization Identity:
   Verify the exact identity (Proposition 8.19, Eq. 8.19.2):
       R_gap(N) = b_{01}^2 / R_1,
   where:
       b_{01} = <u_0^{even}, K u_1> = <c, K u_1> = - (D_0 * a_1) / (mu_1 - lambda_0),
       R_1    = a_1^2 / (mu_1 - lambda_0),
   checking identity residual |R_gap - b_{01}^2 / R_1| < 10^{-45}.

2. Transmission Ratio Stability:
   Measure the first excited odd transmission ratio R_1 = a_1^2 / (mu_1 - lambda_0)
   across N, confirming that R_1 remains strictly O(1) (~ 5.13), so that R_gap(N) -> 0
   is driven identically by the decay of the transition dipole b_{01}^2 -> 0.

3. Coordinate Wavepacket Mode-1 Concentration:
   Compute the excited coordinate wavepacket residual norm:
       ||v_exc||^2 = ||P_{perp u_0} K c||^2 = sum_{j=1}^{N-1} b_{0j}^2 = D_0^2 M_{2,exc}.
   Measure the concentration ratio b_{01}^2 / ||v_exc||^2, testing whether mode 1
   carries > 99.999% of the excited wavepacket residual.

4. Even-Sector Resolvent Cancellation:
   Verify Parseval's identity on K u_1 in the even eigenbasis (Eq. 8.19.5):
       ||K u_1||^2 = b_{01}^2 + a_1^2 * sum_{k=1}^N d_k^2 / (mu_1 - E_k)^2,
   revealing that b_{01}^2 ~ 10^{-6} emerges from a near-perfect cancellation between
   the O(1) coordinate derivative norm ||K u_1||^2 ~ 4.5 and the O(1) excited even
   resolvent energy E_{even, 1} = a_1^2 * ||(Q_even - mu_1 I)^{-1} P_{perp c} d||^2.

OUTPUT CONSTRAINTS:
-------------------
- Strictly dispassionate, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 69 EXECUTION COMPLETE.
================================================================================
"""

from __future__ import annotations

import time
import mpmath as mp

from connes_cvs import build_galerkin_matrix

# -----------------------------------------------------------------------------
# Precision and Parameter Configuration
# -----------------------------------------------------------------------------

mp.mp.dps = 50

C_PARAM = 13
L_PARAM = mp.log(C_PARAM)
T_PARAM = 400
GROUND_DPS = 50

N_LIST = [8, 12, 16, 20, 24]


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

def run_cell69() -> None:
    print("=" * 96)
    print("CELL 69 — AUDIT OF THE EXACT TRANSITION DIPOLE FACTORIZATION R_gap = b_{01}^2 / R_1")
    print("         AND EVEN-RESOLVENT CANCELLATION IN THE FIRST RELATIVE TUNNELING GAP")
    print("=" * 96)
    print(
        f"Parameters: c = {C_PARAM}, L = {mp.nstr(L_PARAM, 15)}, "
        f"T = {T_PARAM}, dps = {mp.mp.dps}"
    )
    print()

    synthesis_records: list[dict] = []

    for N in N_LIST:
        t0 = time.time()

        # Build Galerkin operator
        Q_full = build_galerkin_matrix(N, C_PARAM, L_PARAM, T_PARAM, dps=GROUND_DPS)
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

        # First excited odd mode j = 1 (mode 0 is odd ground state)
        mu_1 = evals_o[1]
        gap_1 = mu_1 - lam_0
        a_1 = a_j_list[1]
        R_1 = (a_1 ** 2) / gap_1

        # Direct relative gap
        R_gap = D_0_sq / gap_1

        # Coordinate transition dipole b_{01} = <c, K u_1>
        # (K_{eo} u_1)[m] = m * u_{1, m} for m >= 1; 0 at m = 0.
        b_01 = sum(c_even[m, 0] * m * V_o[m - 1, 1] for m in range(1, N + 1))
        b_01_sq = b_01 ** 2

        # Commutator prediction: b_{01}^{comm} = - (D_0 * a_1) / gap_1
        b_01_comm = - (D_0 * a_1) / gap_1
        b_01_comm_err = abs(b_01 - b_01_comm)

        # Factorized relative gap: R_gap^{fact} = b_{01}^2 / R_1
        R_gap_fact = b_01_sq / R_1
        fact_error = abs(R_gap - R_gap_fact)

        # Excited coordinate wavepacket residual norm ||v_exc||^2 = sum_{j=1}^{N-1} b_{0j}^2
        b_0j_sq_list: list[mp.mpf] = []
        for j in range(1, N):
            b_j = sum(c_even[m, 0] * m * V_o[m - 1, j] for m in range(1, N + 1))
            b_0j_sq_list.append(b_j ** 2)

        v_exc_sq = sum(b_0j_sq_list)
        mode1_share = (b_01_sq / v_exc_sq) * 100 if v_exc_sq > 0 else mp.mpf(0)

        # Discrete coordinate kinetic energy of u_1: ||K u_1||^2 = sum_{m=1}^N m^2 (V_o[m-1, 1])^2
        ke_1 = sum((m ** 2) * (V_o[m - 1, 1] ** 2) for m in range(1, N + 1))

        # Excited even-resolvent energy: E_{even, 1} = a_1^2 * sum_{k=1}^N d_k^2 / (mu_1 - E_k)^2
        even_resolvent_sum = sum((d_k_list[k] ** 2) / ((mu_1 - evals_e[k]) ** 2) for k in range(1, N + 1))
        E_even_1 = (a_1 ** 2) * even_resolvent_sum

        # Cancellation residual: |ke_1 - E_even_1 - b_01_sq|
        cancellation_residual = abs(ke_1 - E_even_1 - b_01_sq)

        elapsed = time.time() - t0

        record = {
            "N": N,
            "elapsed": elapsed,
            "D_0_sq": D_0_sq,
            "gap_1": gap_1,
            "a_1": a_1,
            "R_1": R_1,
            "R_gap": R_gap,
            "b_01": b_01,
            "b_01_sq": b_01_sq,
            "b_01_comm_err": b_01_comm_err,
            "fact_error": fact_error,
            "v_exc_sq": v_exc_sq,
            "mode1_share": mode1_share,
            "ke_1": ke_1,
            "E_even_1": E_even_1,
            "cancellation_residual": cancellation_residual,
        }
        synthesis_records.append(record)

        print(f"--- DIMENSION N = {N:2d} (Elapsed: {elapsed:.2f}s) ---")
        print(f"  Tunneling gap:       mu_1 - lam_0         = {mp.nstr(gap_1, 9)}")
        print(f"  Odd overlap:         a_1                  = {mp.nstr(a_1, 9)}")
        print(f"  Transmission ratio:  R_1 = a_1^2 / gap_1  = {mp.nstr(R_1, 9)}")
        print(f"  Relative gap:        R_gap                = {mp.nstr(R_gap, 9)}")
        print(f"  Transition dipole:   b_{{01}}               = {mp.nstr(b_01, 9)}")
        print(f"  Dipole commutator:   |b_{{01}} - b_{{comm}}|   = {mp.nstr(b_01_comm_err, 4)}")
        print(f"  Dipole squared:      b_{{01}}^2             = {mp.nstr(b_01_sq, 9)}")
        print(f"  Factorized R_gap:    b_{{01}}^2 / R_1       = {mp.nstr(R_gap_fact, 9)}")
        print(f"  Factorization error: |R_gap - b^2/R_1|    = {mp.nstr(fact_error, 4)}")
        print(f"  Excited norm:        ||v_exc||^2          = {mp.nstr(v_exc_sq, 9)}")
        print(f"  Mode 1 share:        b_{{01}}^2 / ||v_exc||^2 = {mp.nstr(mode1_share, 8)}%")
        print(f"  Coordinate KE:       ||K u_1||^2          = {mp.nstr(ke_1, 9)}")
        print(f"  Even resolvent term: E_{{even, 1}}          = {mp.nstr(E_even_1, 9)}")
        print(f"  Cancellation error:  |||Ku_1||^2 - E - b^2| = {mp.nstr(cancellation_residual, 4)}")
        print()

    # -------------------------------------------------------------------------
    # Synthesis Tables
    # -------------------------------------------------------------------------
    print("=" * 116)
    print("SYNTHESIS TABLE 1: EXACT RELATIVE GAP FACTORIZATION R_gap = b_{01}^2 / R_1")
    print("=" * 116)
    print(f"{'N':>4} | {'R_gap':>13} | {'b_{01}^2':>13} | {'R_1':>10} | {'b_{01}^2 / R_1':>14} | {'Fact Error':>12} | {'Comm Error':>12}")
    print("-" * 116)
    for r in synthesis_records:
        print(
            f"{r['N']:4d} | "
            f"{mp.nstr(r['R_gap'], 6):>13} | "
            f"{mp.nstr(r['b_01_sq'], 6):>13} | "
            f"{mp.nstr(r['R_1'], 5):>10} | "
            f"{mp.nstr(r['b_01_sq'] / r['R_1'], 6):>14} | "
            f"{mp.nstr(r['fact_error'], 4):>12} | "
            f"{mp.nstr(r['b_01_comm_err'], 4):>12}"
        )
    print("=" * 116)
    print()

    print("=" * 116)
    print("SYNTHESIS TABLE 2: EVEN-RESOLVENT CANCELLATION & MODE-1 WAVEPACKET DOMINANCE")
    print("=" * 116)
    print(f"{'N':>4} | {'||K u_1||^2':>13} | {'E_{even, 1}':>13} | {'b_{01}^2':>13} | {'Cancel Residual':>16} | {'Mode-1 Share':>14}")
    print("-" * 116)
    for r in synthesis_records:
        print(
            f"{r['N']:4d} | "
            f"{mp.nstr(r['ke_1'], 6):>13} | "
            f"{mp.nstr(r['E_even_1'], 6):>13} | "
            f"{mp.nstr(r['b_01_sq'], 6):>13} | "
            f"{mp.nstr(r['cancellation_residual'], 4):>16} | "
            f"{mp.nstr(r['mode1_share'], 6):>13}%"
        )
    print("=" * 116)
    print()

    print("=" * 80)
    print("CELL 69 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell69()
