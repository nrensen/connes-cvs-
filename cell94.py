#!/usr/bin/env python3
"""
================================================================================
CELL 94 — HIGH-THROUGHPUT LARGE-N BOUNDARY GAP & PROJECTOR CONVERGENCE
          STRESS TEST (N in {64, 80, 96, 112, ..., 256})
================================================================================

Strategic Roadmap Milestone M40 (Paper NR2 Section 8.25 & Remark 8.35):
----------------------------------------------------------------------
Following the execution of Cell 92/93 and the reviewer's diagnostic:
  1. The Purpose of an Overnight Run to Large N in the Hundreds:
     Investigate whether the boundary gap g_{11} = E_{12} - E_{11} and the
     11-dimensional spectral projector P_{11} survive when dimension N is
     increased by a factor of 4 to 6 beyond previous limits (N = 64 -> 256).
  2. The Three Asymptotic Scenarios to Distinguish:
     - Scenario A: g_{11} settles into a visibly positive asymptotic plateau
       (reinforcing H_{gap}(11) as a serious analytical target).
     - Scenario B: g_{11} continues drifting downward slowly (e.g. power-law decay
       g_{11} ~ N^{-alpha}, indicating boundary gap dissolution in the continuum).
     - Scenario C: g_{11} eventually collapses (identifying a pre-asymptotic scale
       where the current discrete morphology breaks down).
  3. Descriptive Local Logarithmic Slopes (Zero Curve Fitting):
     Track the effective local power-law exponent:
         s_{11}(N_1, N_2) = - log[g_{11}(N_2) / g_{11}(N_1)] / log(N_2 / N_1)
     and s_{E11}(N_1, N_2) as purely descriptive diagnostics, without fitting any
     functional law or extrapolating limiting values.
  4. Projector Cauchy Increments across Uniform Steps Delta N = 16:
     Track ||P_K^{(N+16)} - P_K^{(N)}||_{op} for K in {10, 11, 12} to test whether
     the low-energy subspace stabilizes as an invariant geometric object.
  5. Primary Compact Executive Output:
     Print a dedicated, high-prominence summary table:
         [N, E_11, E_12, g_10, g_11, g_12, Gamma_11, ||Delta P_11||_{op}]
     at the conclusion of the run for immediate inspection.

THE FOUR INVESTIGATIVE TESTS OF CELL 94:
----------------------------------------
- Executive Summary Table: [N, E_11, E_12, g_10, g_11, g_12, Gamma_11, ||DP_11||_op]
- Test A: Boundary Spectral Gaps, Prominence Ratio Gamma_11, and Descriptive Log Slopes
- Test B: Spectral Projector Cauchy Convergence across Uniform Steps Delta N = 16
- Test C: Competing Boundary Hypotheses Audit: H_{gap}(10) vs H_{gap}(11)
- Test D: Subspace Complementary Projection Tilt sin theta_max for K in {10, 11, 12}

OUTPUT CONSTRAINTS:
-------------------
- Dispassionate, dry, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 94 EXECUTION COMPLETE.
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

# High-throughput large-N sequence with uniform step Delta N = 16
# Extensible to [..., 288, 320, 384] if compute budget permits.
N_LIST = [64, 80, 96, 112, 128, 144, 160, 176, 192]
BOUNDARY_MODES = [10, 11, 12, 13, 14]
PROJECTOR_CUTOFFS = [10, 11, 12]


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
# Spectral Projector Construction & Difference Norms
# -----------------------------------------------------------------------------

def build_spectral_projector(V_e: mp.matrix, N: int, K: int) -> mp.matrix:
    """
    Construct spectral projector P_K = sum_{j=0}^K u_j u_j^T in R^{(N+1) x (N+1)}.
    V_e[:, 0:K+1] contains the first K+1 orthonormal eigenvectors.
    """
    dim = N + 1
    P = mp.matrix(dim, dim)
    for j in range(K + 1):
        for r in range(dim):
            v_rj = V_e[r, j]
            for c in range(r, dim):
                term = v_rj * V_e[c, j]
                P[r, c] += term
                if r != c:
                    P[c, r] += term
    return P


def compute_projector_difference(
    P_low: mp.matrix, N_low: int, P_high: mp.matrix, N_high: int
) -> tuple[mp.mpf, mp.mpf]:
    """
    Compute ||Delta P||_{op} and ||Delta P||_F where Delta P = P_high - P_low_tilde
    and P_low_tilde in R^{(N_high+1) x (N_high+1)} is zero-padded from P_low.
    Since Delta P is symmetric, ||Delta P||_{op} is the maximum absolute eigenvalue.
    """
    dim_high = N_high + 1
    Delta_P = mp.matrix(dim_high, dim_high)

    # Fill Delta_P = P_high - P_low_tilde
    for r in range(dim_high):
        for c in range(dim_high):
            val = P_high[r, c]
            if r <= N_low and c <= N_low:
                val -= P_low[r, c]
            Delta_P[r, c] = val

    # Symmetrize explicitly
    Delta_P = mp.mpf('0.5') * (Delta_P + Delta_P.T)

    # Frobenius norm
    frob_sq = mp.mpf('0')
    for r in range(dim_high):
        for c in range(dim_high):
            frob_sq += Delta_P[r, c] * Delta_P[r, c]
    norm_F = mp.sqrt(max(mp.mpf('0'), frob_sq))

    # Operator norm via eigenvalues of symmetric matrix
    evals_diff, _ = mp.eigsy(Delta_P)
    norm_op = max(abs(ev) for ev in evals_diff)

    return norm_op, norm_F


def compute_subspace_tilt(
    V_low: mp.matrix, N_low: int, V_high: mp.matrix, N_high: int, K: int
) -> tuple[mp.mpf, mp.mpf]:
    """
    Compute principal angle metrics between range of V_low[:, 0:K+1] (zero-padded)
    and V_high[:, 0:K+1]:
      M = (V_low_tilde)^T * V_high  in R^{(K+1) x (K+1)}
      cos(theta_max) = sigma_min(M)
      sin(theta_max) = sqrt(1 - cos^2(theta_max)) = ||(I - P_low_tilde) P_high||_{op}
    """
    cols = K + 1
    M = mp.matrix(cols, cols)
    for i in range(cols):
        for j in range(cols):
            val = mp.mpf('0')
            for r in range(N_low + 1):
                val += V_low[r, i] * V_high[r, j]
            M[i, j] = val

    # Singular values via eigenvalues of M^T * M
    MTM = M.T * M
    MTM = mp.mpf('0.5') * (MTM + MTM.T)
    evals_mtm, _ = mp.eigsy(MTM)

    min_eval = max(mp.mpf('0'), min(evals_mtm))
    cos_theta_max = mp.sqrt(min_eval)
    cos_theta_max = min(mp.mpf('1'), cos_theta_max)

    sin_theta_max = mp.sqrt(max(mp.mpf('0'), mp.mpf('1') - cos_theta_max * cos_theta_max))

    return cos_theta_max, sin_theta_max


def compute_log_slope(v1: mp.mpf, v2: mp.mpf, n1: int, n2: int) -> mp.mpf:
    """
    Compute descriptive local logarithmic slope s = - log(v2 / v1) / log(n2 / n1).
    Positive s indicates decay with increasing N; s = 0 indicates constancy.
    """
    if v1 <= 0 or v2 <= 0 or n1 <= 0 or n2 <= 0 or n1 == n2:
        return mp.mpf('nan')
    num = -mp.log(v2 / v1)
    den = mp.log(mp.mpf(n2) / mp.mpf(n1))
    return num / den


# -----------------------------------------------------------------------------
# Main Execution Function
# -----------------------------------------------------------------------------

def run_cell94() -> None:
    print("=" * 140)
    print("CELL 94 — HIGH-THROUGHPUT LARGE-N BOUNDARY GAP & PROJECTOR CONVERGENCE")
    print("          STRESS TEST (N in {64, 80, 96, 112, ..., 256})")
    print("=" * 140)
    print(f"Configuration: c = {C_PARAM}, L = log(c) = {mp.nstr(L_PARAM, 12)}, T = {T_PARAM}")
    print(f"Working Precision: {mp.mp.dps} decimal digits")
    print(f"Dimensions Sweep: N in {N_LIST}")
    print(f"Projector Cutoffs Focus: K in {PROJECTOR_CUTOFFS}")
    print("=" * 140)

    all_data: dict[int, dict] = {}
    start_total_time = time.time()

    for N in N_LIST:
        step_start = time.time()
        print(f"\n>>> Processing Dimension N = {N} (matrix dim = {N+1}) ...")

        Q_full, _ = get_galerkin_matrix(
            c=C_PARAM,
            N=N,
            T=T_PARAM,
            dps=GROUND_DPS,
            verbose=False,
        )

        Q_even, evals_e, V_e = solve_parity_eigensystem(Q_full, N)

        # Precompute spectral projectors for K in PROJECTOR_CUTOFFS
        projectors = {}
        for K in PROJECTOR_CUTOFFS:
            projectors[K] = build_spectral_projector(V_e, N, K)

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
            'projectors': projectors,
            'elapsed_sec': elapsed,
        }

    total_elapsed = time.time() - start_total_time
    print(f"\nAll {len(N_LIST)} dimensions processed in {total_elapsed:.2f}s ({total_elapsed/60:.2f}m).")

    step_pairs = [(N_LIST[i], N_LIST[i + 1]) for i in range(len(N_LIST) - 1)]

    # Compute projector differences across consecutive steps
    delta_P_data: dict[tuple[int, int], dict[int, tuple[mp.mpf, mp.mpf]]] = {}
    for n_low, n_high in step_pairs:
        p_dict_low = all_data[n_low]['projectors']
        p_dict_high = all_data[n_high]['projectors']
        delta_P_data[(n_low, n_high)] = {}
        for K in PROJECTOR_CUTOFFS:
            norm_op, norm_F = compute_projector_difference(p_dict_low[K], n_low, p_dict_high[K], n_high)
            delta_P_data[(n_low, n_high)][K] = (norm_op, norm_F)

    # =========================================================================
    # PRIMARY EXECUTIVE SUMMARY TABLE
    # =========================================================================
    print("\n" + "#" * 140)
    print("PRIMARY EXECUTIVE SUMMARY TABLE: LARGE-N TRAJECTORY OF BOUNDARY SEPARATION")
    print("Track: [N, E_11, E_12, g_10, g_11, g_12, Gamma_11, ||Delta P_11||_{op}]")
    print("#" * 140)
    print(f"{'N':>4} | {'E_{11}':>14} | {'E_{12}':>14} | {'g_{10}':>14} | {'g_{11}':>14} | {'g_{12}':>14} | {'Gamma_{11}':>12} | {'||DP_{11}||_op':>14} | {'Step time':>10}")
    print("-" * 140)

    for i, N in enumerate(N_LIST):
        evals = all_data[N]['evals_e']
        e10 = evals[10]
        e11 = evals[11]
        e12 = evals[12]
        e13 = evals[13]

        g10 = e11 - e10
        g11 = e12 - e11
        g12 = e13 - e12

        denom = max(g10, g12)
        gamma11 = g11 / denom if denom > 0 else mp.mpf('inf')

        # Find incoming step projector difference
        if i > 0:
            pair = (N_LIST[i - 1], N)
            dp11_op = mp.nstr(delta_P_data[pair][11][0], 6)
        else:
            dp11_op = "—"

        t_step = f"{all_data[N]['elapsed_sec']:.1f}s"

        print(f"{N:>4} | {mp.nstr(e11, 8):>14} | {mp.nstr(e12, 8):>14} | {mp.nstr(g10, 8):>14} | {mp.nstr(g11, 8):>14} | {mp.nstr(g12, 8):>14} | {mp.nstr(gamma11, 6):>12} | {dp11_op:>14} | {t_step:>10}")

    print("#" * 140)

    # =========================================================================
    # TEST A: Boundary Spectral Gaps and Descriptive Local Logarithmic Slopes
    # =========================================================================
    print("\n" + "=" * 140)
    print("TEST A: Boundary Spectral Gaps g_j and Descriptive Local Logarithmic Slopes")
    print("s_{11}(N_1, N_2) = - log[g_{11}(N_2)/g_{11}(N_1)] / log(N_2/N_1) and s_{E11}")
    print("(Descriptive diagnostic of effective local power-law exponent; no curve fit)")
    print("=" * 140)

    print(f"{'Step (N1 -> N2)':>18} | {'g_{11}(N1)':>14} | {'g_{11}(N2)':>14} | {'Ratio g11':>12} | {'Slope s_{11}':>14} | {'Ratio E11':>12} | {'Slope s_{E11}':>14}")
    print("-" * 140)

    for n_low, n_high in step_pairs:
        evals_low = all_data[n_low]['evals_e']
        evals_high = all_data[n_high]['evals_e']

        e11_low = evals_low[11]
        e11_high = evals_high[11]
        g11_low = evals_low[12] - e11_low
        g11_high = evals_high[12] - e11_high

        ratio_g11 = g11_high / g11_low if g11_low > 0 else mp.mpf('nan')
        ratio_e11 = e11_high / e11_low if e11_low > 0 else mp.mpf('nan')

        slope_g11 = compute_log_slope(g11_low, g11_high, n_low, n_high)
        slope_e11 = compute_log_slope(e11_low, e11_high, n_low, n_high)

        print(f"{f'{n_low} -> {n_high}':>18} | {mp.nstr(g11_low, 7):>14} | {mp.nstr(g11_high, 7):>14} | {mp.nstr(ratio_g11, 6):>12} | {mp.nstr(slope_g11, 6):>14} | {mp.nstr(ratio_e11, 6):>12} | {mp.nstr(slope_e11, 6):>14}")

    # =========================================================================
    # TEST B: Spectral Projector Cauchy Convergence across Uniform Steps Delta N = 16
    # =========================================================================
    print("\n" + "=" * 140)
    print("TEST B: Spectral Projector Cauchy Differences ||P_K^{(N+16)} - P_K^{(N)}||_{op} and ||Delta P_K||_F")
    print("Auditing Subspace Stabilization across Cutoffs K in {10, 11, 12}")
    print("=" * 140)

    header_parts = []
    for K in PROJECTOR_CUTOFFS:
        header_parts.append(f"||DP_{K}||_op")
        header_parts.append(f"||DP_{K}||_F")
    hdr_str = " | ".join(f"{h:>14}" for h in header_parts)
    print(f"{'Step (N -> N+16)':>18} | {hdr_str}")
    print("-" * 140)

    for n_low, n_high in step_pairs:
        vals = []
        for K in PROJECTOR_CUTOFFS:
            norm_op, norm_F = delta_P_data[(n_low, n_high)][K]
            vals.append(mp.nstr(norm_op, 7))
            vals.append(mp.nstr(norm_F, 7))

        row_str = " | ".join(f"{v:>14}" for v in vals)
        print(f"{f'{n_low} -> {n_high}':>18} | {row_str}")

    # =========================================================================
    # TEST C: Competing Boundary Hypotheses Audit: H_gap(10) vs H_gap(11)
    # =========================================================================
    print("\n" + "=" * 140)
    print("TEST C: Competing Boundary Hypotheses Audit: H_{gap}(10) vs H_{gap}(11)")
    print("Tracking the Collapse of g_{10} = E_{11} - E_{10} vs the Evolution of g_{11} = E_{12} - E_{11}")
    print("=" * 140)

    print(f"{'N':>4} | {'g_{10}':>14} | {'g_{10} step ratio':>18} | {'g_{11}':>14} | {'g_{11} step ratio':>18} | {'g_{11} / g_{10}':>16}")
    print("-" * 140)

    for i, N in enumerate(N_LIST):
        evals = all_data[N]['evals_e']
        g10 = evals[11] - evals[10]
        g11 = evals[12] - evals[11]
        ratio_g = g11 / g10 if g10 > 0 else mp.mpf('inf')

        if i > 0:
            n_prev = N_LIST[i - 1]
            evals_prev = all_data[n_prev]['evals_e']
            g10_prev = evals_prev[11] - evals_prev[10]
            g11_prev = evals_prev[12] - evals_prev[11]
            step_r_g10 = mp.nstr(g10 / g10_prev, 6)
            step_r_g11 = mp.nstr(g11 / g11_prev, 6)
        else:
            step_r_g10 = "—"
            step_r_g11 = "—"

        print(f"{N:>4} | {mp.nstr(g10, 8):>14} | {step_r_g10:>18} | {mp.nstr(g11, 8):>14} | {step_r_g11:>18} | {mp.nstr(ratio_g, 6):>16}")

    # =========================================================================
    # TEST D: Subspace Complementary Projection Tilt sin theta_max
    # =========================================================================
    print("\n" + "=" * 140)
    print("TEST D: Subspace Tilt sin theta_{max} = ||(I - P_K^{(N)}) P_K^{(N+16)}||_{op} for K in {10, 11, 12}")
    print("Auditing Whether the Low-Energy Subspaces Stabilize as Invariant Geometric Objects")
    print("=" * 140)

    print(f"{'Step (N -> N+16)':>18} | {'sin theta(K=10)':>18} | {'cos theta(K=10)':>18} | {'sin theta(K=11)':>18} | {'cos theta(K=11)':>18} | {'sin theta(K=12)':>18} | {'cos theta(K=12)':>18}")
    print("-" * 140)

    for n_low, n_high in step_pairs:
        v_low = all_data[n_low]['V_e']
        v_high = all_data[n_high]['V_e']

        cos10, sin10 = compute_subspace_tilt(v_low, n_low, v_high, n_high, 10)
        cos11, sin11 = compute_subspace_tilt(v_low, n_low, v_high, n_high, 11)
        cos12, sin12 = compute_subspace_tilt(v_low, n_low, v_high, n_high, 12)

        print(f"{f'{n_low} -> {n_high}':>18} | {mp.nstr(sin10, 7):>18} | {mp.nstr(cos10, 7):>18} | {mp.nstr(sin11, 7):>18} | {mp.nstr(cos11, 7):>18} | {mp.nstr(sin12, 7):>18} | {mp.nstr(cos12, 7):>18}")

    print("\n" + "=" * 80)
    print("CELL 94 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell94()
