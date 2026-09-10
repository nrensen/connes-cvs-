#!/usr/bin/env python3
"""
================================================================================
CELL 92 — BOUNDARY-GAP STRESS TEST AT EXTENDED DIMENSIONS (N <= 56),
          RITZ RESIDUALS r_j, AND 11|12 DECOUPLING PERSISTENCE
================================================================================

Strategic Roadmap Milestone M38 (Paper NR2 Section 8.25 & Remark 8.34):
----------------------------------------------------------------------
Following the execution of Cell 91 and the reviewer's diagnostic:
  1. The Boundary-Gap Invariance:
     While individual eigenvalues continue to drift downward (E_{12}: 1.96 -> 0.5509),
     the spectral gap between modes 11 and 12 remains remarkably stable:
         g_{11}(N) = E_{12}^{(N)} - E_{11}^{(N)}:
         0.5388 (N=28), 0.5502 (N=32), 0.5693 (N=36), 0.5612 (N=40), 0.5416 (N=44).
     The central asymptotic object is the boundary gap H_{gap}(J), not an individual
     eigenvalue floor.
  2. Decoupling of Modes 11 and 12:
     Modes 11 and 12 initially underwent heavy finite-N character exchange
     (|<u_{11}^{20}, u_{12}^{24}>| = 0.865), but decoupled sharply at higher sizes
     (cross-overlap dropped to 0.0016, diagonal overlaps reached 0.9993 and 0.9744).
  3. Crossing of the 0.01 Threshold by Mode 11:
     At N = 44, E_{11} fell to 0.009226 < 0.01, shifting N(E < 0.01) from 11 to 12.
     Meanwhile, the descent ratio E_{11}(N)/E_{11}(N-4) rose from 0.33 to 0.82,
     confirming dramatic deceleration while leaving both E_{11} -> 0 and
     E_{11} -> E_{11}^{(infty)} > 0 open.
  4. Strategic Shift to Milestone M38:
     Push dimensions further to N in {36, 40, 44, 48, 52, 56} to stress-test the
     persistence of g_{11} ~ 0.54 - 0.57, evaluate embedded Ritz residuals r_j to
     prepare for Kato-Temple / Davis-Kahan perturbation proofs, and verify decoupling
     persistence.

THE FOUR INVESTIGATIVE TESTS OF CELL 92:
----------------------------------------
- Test A: Boundary-Gap Stress Test across Extended Dimensions (N in {36, ..., 56})
- Test B: Boundary Ritz Residuals r_j and Davis-Kahan Perturbative Bounds
- Test C: Decoupling Persistence: Individual & Cross-Overlaps for {u_{11}, u_{12}}
- Test D: Trajectory and Deceleration Audit of Mode 11

OUTPUT CONSTRAINTS:
-------------------
- Dispassionate, dry, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 92 EXECUTION COMPLETE.
================================================================================
"""

from __future__ import annotations

import time
import mpmath as mp

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

# Extended dimensions sweep up to N = 56
N_LIST = [36, 40, 44, 48, 52, 56]
BOUNDARY_MODES = [10, 11, 12, 13, 14]

FINE_THRESHOLDS = [
    mp.mpf('0.001'),
    mp.mpf('0.002'),
    mp.mpf('0.005'),
    mp.mpf('0.008'),
    mp.mpf('0.010'),
    mp.mpf('0.020'),
    mp.mpf('0.050'),
    mp.mpf('0.100'),
    mp.mpf('0.500'),
    mp.mpf('1.000'),
]


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


def solve_parity_eigensystem(
    Q_full: mp.matrix, N: int
) -> tuple[mp.matrix, list[mp.mpf], mp.matrix]:
    """Compute even spectrum, eigenvectors, and projected matrix."""
    E, _ = full_parity_basis(N)

    Q_even = E.T * Q_full * E
    Q_even = mp.mpf('0.5') * (Q_even + Q_even.T)

    evals_e, V_e = mp.eigsy(Q_even)

    idx_e = sorted(range(N + 1), key=lambda i: evals_e[i])
    sorted_evals_e = [evals_e[i] for i in idx_e]
    sorted_V_e = mp.matrix(N + 1, N + 1)
    for col_idx, orig_col in enumerate(idx_e):
        for row_idx in range(N + 1):
            sorted_V_e[row_idx, col_idx] = V_e[row_idx, orig_col]

    return Q_even, sorted_evals_e, sorted_V_e


# -----------------------------------------------------------------------------
# Ritz Residuals & Overlap Calculations
# -----------------------------------------------------------------------------

def compute_embedded_ritz_residual(
    Q_even_high: mp.matrix, N_low: int, N_high: int,
    u_col: mp.matrix, col_idx: int, eval_low: mp.mpf
) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    """
    Compute Ritz residual r_j = ||Q_{N_high} u_tilde - eval_low * u_tilde||_2
    where u_tilde in R^{N_high + 1} is zero-padded from u_col[:, col_idx].
    Also returns the spectral isolation distance delta_j to other evals of Q_{high}.
    """
    dim_high = N_high + 1
    # Zero-padded eigenvector
    u_tilde = mp.matrix(dim_high, 1)
    for m in range(N_low + 1):
        u_tilde[m, 0] = u_col[m, col_idx]

    # Matrix-vector product
    Qu = Q_even_high * u_tilde

    # Residual vector: Qu - eval_low * u_tilde
    res_sq = mp.mpf('0')
    for m in range(dim_high):
        diff = Qu[m, 0] - eval_low * u_tilde[m, 0]
        res_sq += diff * diff

    res_norm = mp.sqrt(max(mp.mpf('0'), res_sq))
    return res_norm


def compute_individual_overlap(
    V_low: mp.matrix, N_low: int, col_low: int,
    V_high: mp.matrix, col_high: int
) -> mp.mpf:
    """Compute inner product between zero-padded column of V_low and column of V_high."""
    val = mp.mpf('0')
    for m in range(N_low + 1):
        val += V_low[m, col_low] * V_high[m, col_high]
    return val


# -----------------------------------------------------------------------------
# Main Execution Function
# -----------------------------------------------------------------------------

def run_cell92() -> None:
    print("=" * 135)
    print("CELL 92 — BOUNDARY-GAP STRESS TEST AT EXTENDED DIMENSIONS (N <= 56),")
    print("          RITZ RESIDUALS r_j, AND 11|12 DECOUPLING PERSISTENCE")
    print("=" * 135)
    print(f"Configuration: c = {C_PARAM}, L = log(c) = {mp.nstr(L_PARAM, 12)}, T = {T_PARAM}")
    print(f"Working Precision: {mp.mp.dps} decimal digits")
    print(f"Dimensions Sweep: N in {N_LIST}")
    print(f"Boundary Modes Focus: j in {BOUNDARY_MODES}")
    print("=" * 135)

    all_data: dict[int, dict] = {}
    start_total_time = time.time()

    for N in N_LIST:
        step_start = time.time()
        print(f"\n>>> Processing Dimension N = {N} (dim = {N+1}) ...")

        Q_full, _ = get_galerkin_matrix(
            c=C_PARAM,
            N=N,
            T=T_PARAM,
            dps=GROUND_DPS,
            verbose=False,
        )

        Q_even, evals_e, V_e = solve_parity_eigensystem(Q_full, N)

        elapsed = time.time() - step_start
        e11_str = mp.nstr(evals_e[11], 6)
        e12_str = mp.nstr(evals_e[12], 6)
        g11_str = mp.nstr(evals_e[12] - evals_e[11], 6)
        print(f"   Completed N = {N} in {elapsed:.2f}s | E_11 = {e11_str} | E_12 = {e12_str} | g_11 = {g11_str}")

        all_data[N] = {
            'N': N,
            'Q_even': Q_even,
            'evals_e': evals_e,
            'V_e': V_e,
        }

    total_elapsed = time.time() - start_total_time
    print(f"\nAll dimensions processed in {total_elapsed:.2f}s.")

    step_pairs = [(N_LIST[i], N_LIST[i + 1]) for i in range(len(N_LIST) - 1)]

    # =========================================================================
    # TEST A: Boundary-Gap Stress Test across Extended Dimensions
    # =========================================================================
    print("\n" + "=" * 135)
    print("TEST A: Boundary-Gap Stress Test across Extended Dimensions N in {36, 40, 44, 48, 52, 56}")
    print("Direct Audit of the Persistence of the Boundary Gap g_{11} = E_{12} - E_{11} ~ 0.54 - 0.57")
    print("=" * 135)

    mode_headers = [f"E_{j}" for j in BOUNDARY_MODES]
    header_str = " | ".join(f"{h:>12}" for h in mode_headers)
    print(f"{'N':>3} | {header_str}")
    print("-" * 135)

    for N in N_LIST:
        evals = all_data[N]['evals_e']
        row_vals = [mp.nstr(evals[j], 7) for j in BOUNDARY_MODES]
        row_str = " | ".join(f"{v:>12}" for v in row_vals)
        print(f"{N:>3} | {row_str}")

    print("-" * 135)
    print("Boundary Spectral Gaps g_j = E_{j+1} - E_j across Dimensions:")
    print("-" * 135)

    gap_indices = [10, 11, 12, 13]
    gap_headers = [f"g_{j} (E_{j+1}-E_{j})" for j in gap_indices]
    gap_hdr_str = " | ".join(f"{gh:>18}" for gh in gap_headers)
    print(f"{'N':>3} | {gap_hdr_str}")
    print("-" * 135)

    for N in N_LIST:
        evals = all_data[N]['evals_e']
        gap_vals = [mp.nstr(evals[j + 1] - evals[j], 7) for j in gap_indices]
        gap_str = " | ".join(f"{gv:>18}" for gv in gap_vals)
        print(f"{N:>3} | {gap_str}")

    print("-" * 135)
    print("Step Decrements Delta E_j = E_j^{(N)} - E_j^{(N-4)} across Dimension Steps:")
    print("-" * 135)

    for n_low, n_high in step_pairs:
        evals_low = all_data[n_low]['evals_e']
        evals_high = all_data[n_high]['evals_e']
        diffs = [mp.nstr(evals_high[j] - evals_low[j], 7) for j in BOUNDARY_MODES]
        diffs_str = " | ".join(f"{d:>12}" for d in diffs)
        print(f"{f'{n_low}->{n_high}':>7} | {diffs_str}")

    # =========================================================================
    # TEST B: Ritz Residuals of Embedded Boundary Modes
    # =========================================================================
    print("\n" + "=" * 135)
    print("TEST B: Ritz Residuals r_j^{(N)} = ||Q_{N+4} u_tilde_j^{(N)} - E_j^{(N)} u_tilde_j^{(N)}||_2 for Boundary Modes")
    print("Auditing Operator Invariance of Embedded Eigenvectors & Perturbative Isolation Gaps")
    print("=" * 135)

    res_headers = [f"r_{j} (res)" for j in [10, 11, 12, 13]]
    res_hdr_str = " | ".join(f"{rh:>16}" for rh in res_headers)
    print(f"{'Step (N -> N+4)':>16} | {res_hdr_str}")
    print("-" * 135)

    for n_low, n_high in step_pairs:
        q_high = all_data[n_high]['Q_even']
        v_low = all_data[n_low]['V_e']
        evals_low = all_data[n_low]['evals_e']

        residuals = []
        for j in [10, 11, 12, 13]:
            r_val = compute_embedded_ritz_residual(q_high, n_low, n_high, v_low, j, evals_low[j])
            residuals.append(mp.nstr(r_val, 7))

        res_str = " | ".join(f"{r:>16}" for r in residuals)
        print(f"{f'{n_low} -> {n_high}':>16} | {res_str}")

    print("-" * 135)
    print("Perturbative Spectral Isolation Distance delta_j = min_{k != j} |E_j^{(N)} - E_k^{(N+4)}| and Tilt Ratio r_j / delta_j:")
    print("-" * 135)
    print(f"{'Step (N -> N+4)':>16} | {'delta_{11}':>14} | {'r_{11}/delta_{11}':>18} | {'delta_{12}':>14} | {'r_{12}/delta_{12}':>18}")
    print("-" * 135)

    for n_low, n_high in step_pairs:
        q_high = all_data[n_high]['Q_even']
        v_low = all_data[n_low]['V_e']
        evals_low = all_data[n_low]['evals_e']
        evals_high = all_data[n_high]['evals_e']

        r11 = compute_embedded_ritz_residual(q_high, n_low, n_high, v_low, 11, evals_low[11])
        r12 = compute_embedded_ritz_residual(q_high, n_low, n_high, v_low, 12, evals_low[12])

        delta11 = min(abs(evals_low[11] - evals_high[k]) for k in range(len(evals_high)) if k != 11)
        delta12 = min(abs(evals_low[12] - evals_high[k]) for k in range(len(evals_high)) if k != 12)

        ratio11 = r11 / delta11 if delta11 > 0 else mp.mpf('inf')
        ratio12 = r12 / delta12 if delta12 > 0 else mp.mpf('inf')

        print(f"{f'{n_low} -> {n_high}':>16} | {mp.nstr(delta11, 6):>14} | {mp.nstr(ratio11, 6):>18} | {mp.nstr(delta12, 6):>14} | {mp.nstr(ratio12, 6):>18}")

    # =========================================================================
    # TEST C: Decoupling Persistence: Individual & Cross-Overlaps
    # =========================================================================
    print("\n" + "=" * 135)
    print("TEST C: Decoupling Persistence: Diagonal & Cross-Overlaps between Modes 11 and 12")
    print("Testing Whether the Mode Decoupling Observed at N = 44 Persists through N = 56")
    print("=" * 135)

    print(f"{'Step (N -> N+4)':>16} | {'|<u_{11}, u_{11}>|':>18} | {'|<u_{11}, u_{12}>|':>18} | {'|<u_{12}, u_{11}>|':>18} | {'|<u_{12}, u_{12}>|':>18}")
    print("-" * 135)

    for n_low, n_high in step_pairs:
        v_low = all_data[n_low]['V_e']
        v_high = all_data[n_high]['V_e']

        s_11_11 = abs(compute_individual_overlap(v_low, n_low, 11, v_high, 11))
        s_11_12 = abs(compute_individual_overlap(v_low, n_low, 11, v_high, 12))
        s_12_11 = abs(compute_individual_overlap(v_low, n_low, 12, v_high, 11))
        s_12_12 = abs(compute_individual_overlap(v_low, n_low, 12, v_high, 12))

        print(f"{f'{n_low} -> {n_high}':>16} | {mp.nstr(s_11_11, 8):>18} | {mp.nstr(s_11_12, 8):>18} | {mp.nstr(s_12_11, 8):>18} | {mp.nstr(s_12_12, 8):>18}")

    # =========================================================================
    # TEST D: Trajectory and Deceleration Audit of Mode 11
    # =========================================================================
    print("\n" + "=" * 135)
    print("TEST D: Fine-Grained Mode Counting N(E; N) and Trajectory / Deceleration Audit of Mode 11")
    print("Tracking E_{11}(N), Step Decrements, and Ratios E_{11}(N) / E_{11}(N-4)")
    print("=" * 135)

    thresh_headers = [f"E < {mp.nstr(t, 3)}" for t in FINE_THRESHOLDS]
    thresh_hdr_str = " | ".join(f"{th:>11}" for th in thresh_headers)
    print(f"{'N':>3} | {thresh_hdr_str}")
    print("-" * 135)

    for N in N_LIST:
        evals = all_data[N]['evals_e']
        counts = [sum(1 for ev in evals if ev < t) for t in FINE_THRESHOLDS]
        counts_str = " | ".join(f"{c:>11}" for c in counts)
        print(f"{N:>3} | {counts_str}")

    print("-" * 135)
    print("Trajectory and Step Ratios of Mode 11 (E_{11}(N) / E_{11}(N-4)):")
    print("-" * 135)

    for i in range(1, len(N_LIST)):
        n_prev = N_LIST[i - 1]
        n_curr = N_LIST[i]
        e11_prev = all_data[n_prev]['evals_e'][11]
        e11_curr = all_data[n_curr]['evals_e'][11]
        ratio = e11_curr / e11_prev
        diff = e11_curr - e11_prev
        print(f"Step {n_prev} -> {n_curr}: E_11 = {mp.nstr(e11_curr, 8)} | Ratio = {mp.nstr(ratio, 6)} | Decrement = {mp.nstr(diff, 6)}")

    print("\n" + "=" * 80)
    print("CELL 92 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell92()
