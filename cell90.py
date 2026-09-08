#!/usr/bin/env python3
"""
================================================================================
CELL 90 — NESTED-N SPECTRAL SUBSPACE OVERLAP, PROJECTOR CONVERGENCE,
          AND TRANSITION CLUSTER DYNAMICS
================================================================================

Strategic Roadmap Milestone M36 (Paper 4B Section 8.25 & Remark 8.32):
----------------------------------------------------------------------
Following the execution of Cell 89 and the reviewer's diagnostic:
  1. Deceleration of E_{13} Drift:
     At N = 36, the downward decrement of E_{13} drops by over an order of magnitude:
         E_{13}^{(32)} = 0.813235 -> E_{13}^{(36)} = 0.781610  (decrement Delta = 0.031625),
     compared with ~ 0.43, 0.26, 0.49, 0.50 at previous steps.
     This provides strong empirical indication of a stabilizing Ritz floor E_{13}^{(infty)} > 0.
  2. Persistent Coordinate Delocalization Invariant:
     The low-energy spectral projector P_{12} = sum_{j=0}^{12} u_j u_j^T has ~ 30.69% of its mass
     in modes m > 12 across N in {28, 32, 36}:
         30.694% (N=28), 30.657% (N=32), 30.695% (N=36).
     Coordinate delocalization is an asymptotic property of the continuous eigenspace.
  3. Severe Subspace Tilt in Coordinate Bases:
     cos theta_{min}(M=12) ~ 8.23e-5 and cos theta_{min}(M=20) ~ 0.0433 at N = 36.
     The coordinate basis is fundamentally ill-suited for separating the low-energy sector.
  4. Three-Zone Spectral Hierarchy:
     Deep tunneling (E < 0.01) stabilizes at 11 modes, while the transition sector (E < 1.0)
     grows slowly (10 -> 15).
  5. Strategic Shift to Milestone M36:
     Instead of projecting onto coordinate subspaces, investigate the convergence of the
     INTRINSIC SPECTRAL SUBSPACE across nested Galerkin embeddings V_N subset V_{N+4}:
       - Embed U_K^{(N)} into R^{(N+5) x (K+1)} via zero-padding.
       - Form the overlap matrix S_K^{(N, N+4)} = (U_K^{(N)})^T U_K^{(N+4)} in R^{(K+1) x (K+1)}.
       - Diagonalize S_K^T S_K to obtain principal cosines sigma_j and minimum cosine sigma_{min}.
       - Compute projector difference operator norm ||Delta P_K||_{op} = sin theta_{max} = sqrt(1 - sigma_{min}^2)
         and Frobenius norm ||Delta P_K||_F = sqrt(2 * ((K+1) - Tr(S_K^T S_K))).
       - Sweep K in {8, 9, 10, 11, 12, 13, 14, 15} to test whether low-energy modes converge faster.
       - Track transition cluster energies E_j^{(N)} for j in {7, ..., 16} across N in {16, ..., 40}.
       - Audit extended mode counting function N(E; N) up to N = 40.

THE FOUR INVESTIGATIVE TESTS OF CELL 90:
----------------------------------------
- Test A: Nested-N Spectral Subspace Overlap & Principal Angles
          Audit sigma_{min}((U_K^{(N)})^T U_K^{(N+4)}) and subspace tilt sin theta_{max}
          for K in {8, 10, 11, 12, 13, 14, 15} across dimension steps N -> N+4.

- Test B: Spectral Projector Difference Norms ||P_K^{(N+4)} - P_K^{(N)}||
          Compute Frobenius norm ||Delta P_K||_F and operator norm ||Delta P_K||_{op}
          for K in {10, 12, 14} across dimension steps N -> N+4.

- Test C: Transition Cluster Energy Tracking E_j^{(N)} for j in {7, ..., 16}
          Track individual Ritz energies across the transition barrier to audit whether
          E_{13} establishes a stable plateau and how adjacent modes evolve.

- Test D: Extended Low-Energy Mode Counting Function N(E; N) up to N = 40
          Audit N(E; N) for E in {0.001, 0.01, 0.10, 0.25, 0.50, 0.75, 1.00, 1.50}
          across N in {16, 20, 24, 28, 32, 36, 40}.

OUTPUT CONSTRAINTS:
-------------------
- Dispassionate, dry, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 90 EXECUTION COMPLETE.
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

# Sweep dimensions up to N = 40
N_LIST = [16, 20, 24, 28, 32, 36, 40]
K_SWEEP = [8, 9, 10, 11, 12, 13, 14, 15]
FOCUS_K_LIST = [10, 12, 14]
CLUSTER_MODES = list(range(7, 17))  # Modes 7 ... 16

ENERGY_THRESHOLDS = [
    mp.mpf('0.001'),
    mp.mpf('0.01'),
    mp.mpf('0.10'),
    mp.mpf('0.25'),
    mp.mpf('0.50'),
    mp.mpf('0.75'),
    mp.mpf('1.00'),
    mp.mpf('1.50'),
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
# Subspace Overlap & Projector Metrics
# -----------------------------------------------------------------------------

def compute_nested_subspace_metrics(
    V_low: mp.matrix, N_low: int, V_high: mp.matrix, N_high: int, K: int
) -> tuple[mp.mpf, mp.mpf, mp.mpf, mp.mpf]:
    """
    Compute overlap metrics between U_K^{(N_low)} (dim K+1) and U_K^{(N_high)} (dim K+1).
    Since U_K^{(N_low)} is embedded into R^{N_high+1} by zero-padding rows m in {N_low+1, ..., N_high},
    the overlap matrix S = (U_K^{(N_low)})^T U_K^{(N_high)} is:
        S[i, j] = sum_{m=0}^{N_low} V_low[m, i] * V_high[m, j].

    Returns:
        sigma_min: smallest singular value of S (minimum principal cosine)
        sin_theta_max: sqrt(1 - sigma_min^2) (worst-case subspace tilt / operator norm)
        frob_norm: sqrt(2 * ((K+1) - Tr(S^T S))) (projector difference Frobenius norm)
        trace_overlap: Tr(S^T S)
    """
    dim_k = K + 1
    S = mp.matrix(dim_k, dim_k)

    for i in range(dim_k):
        for j in range(dim_k):
            val = mp.mpf('0')
            for m in range(N_low + 1):
                val += V_low[m, i] * V_high[m, j]
            S[i, j] = val

    # Form Gram matrix G = S^T S
    G = S.T * S
    G = mp.mpf('0.5') * (G + G.T)

    evals_g, _ = mp.eigsy(G)
    evals_g_sorted = sorted(evals_g)

    # Singular values are sqrt(evals)
    sigma_min_sq = max(mp.mpf('0'), min(mp.mpf('1'), evals_g_sorted[0]))
    sigma_min = mp.sqrt(sigma_min_sq)
    sin_theta_max = mp.sqrt(max(mp.mpf('0'), mp.mpf('1') - sigma_min_sq))

    trace_overlap = sum(evals_g_sorted)
    frob_sq = max(mp.mpf('0'), mp.mpf('2') * (mp.mpf(dim_k) - trace_overlap))
    frob_norm = mp.sqrt(frob_sq)

    return sigma_min, sin_theta_max, frob_norm, trace_overlap


# -----------------------------------------------------------------------------
# Main Execution Function
# -----------------------------------------------------------------------------

def run_cell90() -> None:
    print("=" * 135)
    print("CELL 90 — NESTED-N SPECTRAL SUBSPACE OVERLAP, PROJECTOR CONVERGENCE,")
    print("          AND TRANSITION CLUSTER DYNAMICS")
    print("=" * 135)
    print(f"Configuration: c = {C_PARAM}, L = log(c) = {mp.nstr(L_PARAM, 12)}, T = {T_PARAM}")
    print(f"Working Precision: {mp.mp.dps} decimal digits")
    print(f"Dimensions Sweep: N in {N_LIST}")
    print(f"Subspace Cutoffs K in {K_SWEEP} (Subspace Dimensions {K_SWEEP[0]+1} ... {K_SWEEP[-1]+1})")
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
        e13_str = mp.nstr(evals_e[13], 6) if 13 <= N else "N/A"
        print(f"   Completed N = {N} in {elapsed:.2f}s | E_13 = {e13_str}")

        all_data[N] = {
            'N': N,
            'evals_e': evals_e,
            'V_e': V_e,
        }

    total_elapsed = time.time() - start_total_time
    print(f"\nAll dimensions processed in {total_elapsed:.2f}s.")

    # =========================================================================
    # TEST A: Nested-N Spectral Subspace Overlap & Principal Angles
    # =========================================================================
    print("\n" + "=" * 135)
    print("TEST A: Nested-N Spectral Subspace Overlap: Minimum Principal Cosine sigma_{min}((U_K^{(N)})^T U_K^{(N+4)})")
    print("Zero-Padded Embedding U_K^{(N)} into R^{(N+5) x (K+1)} across Dimension Steps N -> N+4")
    print("=" * 135)

    k_header = " | ".join(f"K={k:>2} (d={k+1:>2})" for k in K_SWEEP)
    print(f"{'Step (N -> N+4)':>16} | {k_header}")
    print("-" * 135)

    step_pairs = [(N_LIST[i], N_LIST[i + 1]) for i in range(len(N_LIST) - 1)]

    for n_low, n_high in step_pairs:
        v_low = all_data[n_low]['V_e']
        v_high = all_data[n_high]['V_e']

        sigmas = []
        for K in K_SWEEP:
            if K <= n_low:
                s_min, _, _, _ = compute_nested_subspace_metrics(v_low, n_low, v_high, n_high, K)
                sigmas.append(mp.nstr(s_min, 7))
            else:
                sigmas.append("---")

        row_str = " | ".join(f"{s:>11}" for s in sigmas)
        print(f"{f'{n_low} -> {n_high}':>16} | {row_str}")

    print("-" * 135)
    print("Subspace Tilt / Gap: sin(theta_{max}) = sqrt(1 - sigma_{min}^2) across Dimension Steps:")
    print("-" * 135)

    for n_low, n_high in step_pairs:
        v_low = all_data[n_low]['V_e']
        v_high = all_data[n_high]['V_e']

        tilts = []
        for K in K_SWEEP:
            if K <= n_low:
                _, sin_tilt, _, _ = compute_nested_subspace_metrics(v_low, n_low, v_high, n_high, K)
                tilts.append(mp.nstr(sin_tilt, 7))
            else:
                tilts.append("---")

        row_str = " | ".join(f"{t:>11}" for t in tilts)
        print(f"{f'{n_low} -> {n_high}':>16} | {row_str}")

    print("-" * 135)
    print("Key Diagnostic Summary for Test A:")
    print("1. Subspace Convergence: Tracks whether the low-energy spectral subspace stabilizes as N increases.")
    print("2. Mode Boundary Sensitivity: Compares whether K <= 11/12 converges faster than K >= 13.")

    # =========================================================================
    # TEST B: Spectral Projector Difference Norms ||P_K^{(N+4)} - P_K^{(N)}||
    # =========================================================================
    print("\n" + "=" * 135)
    print("TEST B: Spectral Projector Difference Norms: Operator Norm ||Delta P_K||_{op} and Frobenius Norm ||Delta P_K||_F")
    print("Delta P_K = P_K^{(N+4)} - tilde{P}_K^{(N)} in R^{(N+5) x (N+5)} across Canonical Cutoffs K in {10, 12, 14}")
    print("=" * 135)
    print(f"{'Step (N -> N+4)':>16} | {'K=10 ||Delta P||_op':>18} | {'K=10 ||Delta P||_F':>18} | {'K=12 ||Delta P||_op':>18} | {'K=12 ||Delta P||_F':>18} | {'K=14 ||Delta P||_op':>18} | {'K=14 ||Delta P||_F':>18}")
    print("-" * 135)

    for n_low, n_high in step_pairs:
        v_low = all_data[n_low]['V_e']
        v_high = all_data[n_high]['V_e']

        metrics_cols = []
        for K in FOCUS_K_LIST:
            if K <= n_low:
                _, op_norm, frob_norm, _ = compute_nested_subspace_metrics(v_low, n_low, v_high, n_high, K)
                metrics_cols.append(mp.nstr(op_norm, 8))
                metrics_cols.append(mp.nstr(frob_norm, 8))
            else:
                metrics_cols.append("---")
                metrics_cols.append("---")

        row_str = " | ".join(f"{m:>18}" for m in metrics_cols)
        print(f"{f'{n_low} -> {n_high}':>16} | {row_str}")

    print("-" * 135)
    print("Key Diagnostic Summary for Test B:")
    print("1. Cauchy Convergence Rate: ||Delta P_K|| measures the absolute distance between embedded spectral projectors.")
    print("2. Geometric Stabilization: Rapid decay of ||Delta P_K|| demonstrates that the spectral subspace converges in Hilbert space.")

    # =========================================================================
    # TEST C: Transition Cluster Energy Tracking E_j^{(N)} for j in {7, ..., 16}
    # =========================================================================
    print("\n" + "=" * 135)
    print("TEST C: Transition Cluster Energy Tracking: Ritz Energies E_j^{(N)} for Modes j in {7, ..., 16}")
    print("Tracking Downward Monotonic Drift and Investigating the Emergence of a Stabilizing Plateau")
    print("=" * 135)

    mode_headers = [f"E_{j}" for j in CLUSTER_MODES]
    header_str = " | ".join(f"{h:>11}" for h in mode_headers)
    print(f"{'N':>3} | {header_str}")
    print("-" * 135)

    for N in N_LIST:
        evals = all_data[N]['evals_e']
        row_vals = []
        for j in CLUSTER_MODES:
            if j <= N:
                row_vals.append(mp.nstr(evals[j], 6))
            else:
                row_vals.append("---")
        row_str = " | ".join(f"{v:>11}" for v in row_vals)
        print(f"{N:>3} | {row_str}")

    print("-" * 135)
    print("Step Decrements Delta E_j = E_j^{(N)} - E_j^{(N-4)} across the Dimension Steps:")
    print("-" * 135)

    for i in range(1, len(N_LIST)):
        n_prev = N_LIST[i - 1]
        n_curr = N_LIST[i]
        evals_prev = all_data[n_prev]['evals_e']
        evals_curr = all_data[n_curr]['evals_e']

        diffs = []
        for j in CLUSTER_MODES:
            if j <= n_prev and j <= n_curr:
                diff = evals_curr[j] - evals_prev[j]
                diffs.append(mp.nstr(diff, 6))
            else:
                diffs.append("---")
        row_str = " | ".join(f"{d:>11}" for d in diffs)
        print(f"{f'{n_prev}->{n_curr}':>7} | {row_str}")

    print("-" * 135)
    print("Key Diagnostic Summary for Test C:")
    print("1. Downward Deceleration Audit: Confirms whether the deceleration observed at N = 36 continues into N = 40.")
    print("2. Mode 13 Plateau: Quantifies the stability of the continuum base candidate E_{13}.")

    # =========================================================================
    # TEST D: Extended Low-Energy Mode Counting Function N(E; N) up to N = 40
    # =========================================================================
    print("\n" + "=" * 135)
    print("TEST D: Extended Low-Energy Mode Counting Function N(E; N) = #{l in {0, ..., N} : E_l(N) < E}")
    print("Auditing Whether Deep Tunneling (E < 0.01) Remains Saturated at 11 Modes up to N = 40")
    print("=" * 135)

    header_cols = [f"E < {mp.nstr(thresh, 3)}" for thresh in ENERGY_THRESHOLDS]
    header_str = " | ".join(f"{col:>12}" for col in header_cols)
    print(f"{'N':>3} | {header_str}")
    print("-" * 135)

    for N in N_LIST:
        evals = all_data[N]['evals_e']
        counts = [sum(1 for ev in evals if ev < thresh) for thresh in ENERGY_THRESHOLDS]
        counts_str = " | ".join(f"{cnt:>12}" for cnt in counts)
        print(f"{N:>3} | {counts_str}")

    print("-" * 135)
    print("Key Diagnostic Summary for Test D:")
    print("1. Deep Tunneling Sector: Tests whether N(E < 0.001) and N(E < 0.01) remain strictly locked at 11 modes at N = 40.")
    print("2. Transition Sector Growth: Tracks the growth of N(E < 1.0) and N(E < 1.5) as N increases to 40.")

    print("\n" + "=" * 80)
    print("CELL 90 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell90()
