#!/usr/bin/env python3
"""
================================================================================
CELL 91 — BOUNDARY CLUSTER DYNAMICS, INDIVIDUAL EIGENVECTOR OVERLAPS,
          BOUNDARY GAPS, AND PROJECTOR ENCLOSURE
================================================================================

Strategic Roadmap Milestone M37 (Paper NR2 Section 8.25 & Remark 8.33):
----------------------------------------------------------------------
Following the execution of Cell 90 and the reviewer's diagnostic:
  1. Convergence of K=10 Subspace:
     ||Delta P_{10}||_{op} fell monotonically from 0.8868 -> 0.0224 (N = 36 -> 40),
     confirming that the 11-dimensional deep tunneling subspace U_{10} stabilizes.
  2. Nonlinear Transition of K=12 into Convergence:
     ||Delta P_{12}||_{op} dropped from 0.9963 -> 0.6174 -> 0.1662 -> 0.1217.
  3. Persistent Instability of K=13:
     ||Delta P_{13}||_{op} remained stuck in [0.9892, 1.0000], confirming that mode 13
     remains unsettled and securing K=12 as the natural spectral boundary.
  4. Emergence of E_{12} as the Natural Positive Continuum Base:
     E_{12} stabilized at 0.572442 (recent decrements: -0.076, -0.003, -0.013),
     separated from mode 11 by a massive spectral gap:
         g_{11} = E_{12} - E_{11} = 0.57244 - 0.01122 = 0.56122.
     Meanwhile E_{13} continued drifting downwards (0.7816 -> 0.7261, Delta = -0.0555).
  5. Strategic Shift to Milestone M37:
     Zoom in on the boundary cluster: modes 10, 11, 12, 13, 14 across N in {20, ..., 44}.
     - Test A: Record boundary energies E_{10} ... E_{14}, step decrements Delta E_j,
               and boundary gaps g_{10} ... g_{13}.
     - Test B: Audit embedded projector differences ||Delta P_K||_{op} and ||Delta P_K||_F
               for K in {10, 11, 12, 13} to confirm whether P_{11} and P_{12} converge
               while P_{13} remains unstable.
     - Test C: Compute individual eigenvector overlaps |<u_j^{(N)}, u_j^{(N+4)}>|
               for j in {10, ..., 14} after zero-padding, to determine whether subspace
               convergence reflects individual mode stability or collective cluster rotation.
               Also audit the 2x2 cross-overlap block between {u_{11}, u_{12}}.
     - Test D: Fine-grained low-energy mode counting N(E; N) for thresholds near E_{11},
               and audit the deceleration trajectory of E_{11}(N).

THE FOUR INVESTIGATIVE TESTS OF CELL 91:
----------------------------------------
- Test A: Boundary Spectral Gaps & Energy Trajectories (j in {10, ..., 14})
- Test B: Boundary Projector Cauchy Differences ||Delta P_K||_{op} & ||Delta P_K||_F (K in {10, 11, 12, 13})
- Test C: Individual Eigenvector Overlaps & Cluster Rotation Analysis (j in {10, ..., 14})
- Test D: Fine-Grained Low-Energy Counting & Mode 11 Limit Audit

OUTPUT CONSTRAINTS:
-------------------
- Dispassionate, dry, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 91 EXECUTION COMPLETE.
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

# Sweep dimensions up to N = 44
N_LIST = [20, 24, 28, 32, 36, 40, 44]
BOUNDARY_MODES = [10, 11, 12, 13, 14]
PROJECTOR_K_LIST = [10, 11, 12, 13]

FINE_THRESHOLDS = [
    mp.mpf('0.001'),
    mp.mpf('0.005'),
    mp.mpf('0.010'),
    mp.mpf('0.015'),
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
# Overlap & Projector Metrics
# -----------------------------------------------------------------------------

def compute_projector_metrics(
    V_low: mp.matrix, N_low: int, V_high: mp.matrix, N_high: int, K: int
) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    """
    Compute operator norm and Frobenius norm of projector difference
    Delta P_K = P_K^{(N_high)} - tilde{P}_K^{(N_low)} via overlap Gram matrix.
    """
    dim_k = K + 1
    S = mp.matrix(dim_k, dim_k)

    for i in range(dim_k):
        for j in range(dim_k):
            val = mp.mpf('0')
            for m in range(N_low + 1):
                val += V_low[m, i] * V_high[m, j]
            S[i, j] = val

    G = S.T * S
    G = mp.mpf('0.5') * (G + G.T)

    evals_g, _ = mp.eigsy(G)
    evals_g_sorted = sorted(evals_g)

    sigma_min_sq = max(mp.mpf('0'), min(mp.mpf('1'), evals_g_sorted[0]))
    sin_theta_max = mp.sqrt(max(mp.mpf('0'), mp.mpf('1') - sigma_min_sq))

    trace_overlap = sum(evals_g_sorted)
    frob_sq = max(mp.mpf('0'), mp.mpf('2') * (mp.mpf(dim_k) - trace_overlap))
    frob_norm = mp.sqrt(frob_sq)

    return sin_theta_max, frob_norm, evals_g_sorted[0]


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

def run_cell91() -> None:
    print("=" * 135)
    print("CELL 91 — BOUNDARY CLUSTER DYNAMICS, INDIVIDUAL EIGENVECTOR OVERLAPS,")
    print("          BOUNDARY GAPS, AND PROJECTOR ENCLOSURE")
    print("=" * 135)
    print(f"Configuration: c = {C_PARAM}, L = log(c) = {mp.nstr(L_PARAM, 12)}, T = {T_PARAM}")
    print(f"Working Precision: {mp.mp.dps} decimal digits")
    print(f"Dimensions Sweep: N in {N_LIST}")
    print(f"Boundary Modes Focus: j in {BOUNDARY_MODES}")
    print(f"Projector Cutoffs: K in {PROJECTOR_K_LIST}")
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
        e11_str = mp.nstr(evals_e[11], 6) if 11 <= N else "N/A"
        e12_str = mp.nstr(evals_e[12], 6) if 12 <= N else "N/A"
        print(f"   Completed N = {N} in {elapsed:.2f}s | E_11 = {e11_str} | E_12 = {e12_str}")

        all_data[N] = {
            'N': N,
            'evals_e': evals_e,
            'V_e': V_e,
        }

    total_elapsed = time.time() - start_total_time
    print(f"\nAll dimensions processed in {total_elapsed:.2f}s.")

    # =========================================================================
    # TEST A: Boundary Spectral Gaps & Energy Trajectories
    # =========================================================================
    print("\n" + "=" * 135)
    print("TEST A: Boundary Spectral Gaps & Energy Trajectories for Modes j in {10, 11, 12, 13, 14}")
    print("Auditing the Stabilization of E_{12} and the Massive Barrier Gap g_{11} = E_{12} - E_{11}")
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

    step_pairs = [(N_LIST[i], N_LIST[i + 1]) for i in range(len(N_LIST) - 1)]

    for n_low, n_high in step_pairs:
        evals_low = all_data[n_low]['evals_e']
        evals_high = all_data[n_high]['evals_e']
        diffs = [mp.nstr(evals_high[j] - evals_low[j], 7) for j in BOUNDARY_MODES]
        diffs_str = " | ".join(f"{d:>12}" for d in diffs)
        print(f"{f'{n_low}->{n_high}':>7} | {diffs_str}")

    print("-" * 135)
    print("Key Diagnostic Summary for Test A:")
    print("1. Gap Isolation: Confirms whether the massive gap g_{11} remains stable around ~ 0.56.")
    print("2. E_{12} Stabilization: Measures whether the decrement of E_{12} continues to contract.")

    # =========================================================================
    # TEST B: Boundary Projector Cauchy Differences
    # =========================================================================
    print("\n" + "=" * 135)
    print("TEST B: Boundary Projector Cauchy Differences: Operator Norm ||Delta P_K||_{op} & Frobenius Norm ||Delta P_K||_F")
    print("Testing the Boundary Dichotomy: Rapid Convergence for K in {10, 11, 12} vs Instability for K = 13")
    print("=" * 135)

    proj_headers = []
    for K in PROJECTOR_K_LIST:
        proj_headers.append(f"K={K} ||Delta P||_op")
        proj_headers.append(f"K={K} ||Delta P||_F")
    proj_hdr_str = " | ".join(f"{ph:>16}" for ph in proj_headers)
    print(f"{'Step (N -> N+4)':>16} | {proj_hdr_str}")
    print("-" * 135)

    for n_low, n_high in step_pairs:
        v_low = all_data[n_low]['V_e']
        v_high = all_data[n_high]['V_e']

        row_metrics = []
        for K in PROJECTOR_K_LIST:
            op_norm, frob_norm, _ = compute_projector_metrics(v_low, n_low, v_high, n_high, K)
            row_metrics.append(mp.nstr(op_norm, 7))
            row_metrics.append(mp.nstr(frob_norm, 7))

        row_str = " | ".join(f"{m:>16}" for m in row_metrics)
        print(f"{f'{n_low} -> {n_high}':>16} | {row_str}")

    print("-" * 135)
    print("Key Diagnostic Summary for Test B:")
    print("1. Boundary Enclosure: Verifies whether ||Delta P_{11}|| and ||Delta P_{12}|| continue to decay.")
    print("2. Boundary Dichotomy: Confirms whether ||Delta P_{13}|| remains stubbornly unstable near 1.0.")

    # =========================================================================
    # TEST C: Individual Eigenvector Overlaps & Cluster Rotation Analysis
    # =========================================================================
    print("\n" + "=" * 135)
    print("TEST C: Individual Eigenvector Overlaps |<u_j^{(N)}, u_j^{(N+4)}>| for Modes j in {10, 11, 12, 13, 14}")
    print("Auditing Mode-by-Mode Directional Stabilization vs Collective Subspace Rotation")
    print("=" * 135)

    overlap_headers = [f"|<u_{j}, u_{j}>|" for j in BOUNDARY_MODES]
    overlap_hdr_str = " | ".join(f"{oh:>14}" for oh in overlap_headers)
    print(f"{'Step (N -> N+4)':>16} | {overlap_hdr_str}")
    print("-" * 135)

    for n_low, n_high in step_pairs:
        v_low = all_data[n_low]['V_e']
        v_high = all_data[n_high]['V_e']

        overlaps = []
        for j in BOUNDARY_MODES:
            s_jj = compute_individual_overlap(v_low, n_low, j, v_high, j)
            overlaps.append(mp.nstr(abs(s_jj), 8))

        row_str = " | ".join(f"{o:>14}" for o in overlaps)
        print(f"{f'{n_low} -> {n_high}':>16} | {row_str}")

    print("-" * 135)
    print("Cross-Mode Boundary Overlaps between Modes 11 and 12 across Steps:")
    print("-" * 135)
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

    print("-" * 135)
    print("Key Diagnostic Summary for Test C:")
    print("1. Individual Convergence: Distinguishes whether u_{11} and u_{12} converge individually (|s_{jj}| -> 1)")
    print("   or rotate into each other (|s_{11, 12}| > 0).")
    print("2. Cluster Autonomy: Measures the independence of mode 12 from transition mode 11.")

    # =========================================================================
    # TEST D: Fine-Grained Low-Energy Counting & Mode 11 Limit Audit
    # =========================================================================
    print("\n" + "=" * 135)
    print("TEST D: Fine-Grained Low-Energy Counting Function N(E; N) = #{l : E_l < E} & Mode 11 Limit Audit")
    print("Tracking the Deceleration Trajectory and Potential Floor of Transition Mode E_{11}")
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
    print("CELL 91 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell91()
