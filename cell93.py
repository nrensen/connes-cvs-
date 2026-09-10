#!/usr/bin/env python3
"""
================================================================================
CELL 93 — SPECTRAL PROJECTOR CONVERGENCE ||P_K^{(N+4)} - P_K^{(N)}||,
          BOUNDARY GAP PROMINENCE GAMMA_11, AND EXTENDED SWEEPS (N <= 64)
================================================================================

Strategic Roadmap Milestone M39 (Paper NR2 Section 8.25 & Remark 8.35):
----------------------------------------------------------------------
Following the execution of Cell 92 and the reviewer's diagnostic:
  1. Downward Drift of the 11|12 Gap:
     While individual eigenvalues drifted downward, g_{11} = E_{12} - E_{11}
     descended from 0.542 (N=44) to 0.4273 (N=56), with the latest values
     (0.4523, 0.4294, 0.4273) showing signs of flattening without locking.
  2. Persistence of Structural Gap Prominence:
     The qualitative distinction remains conspicuous:
         g_{10} = E_{11} - E_{10} collapses to 6.26 x 10^{-3},
         g_{11} = E_{12} - E_{11} ~ 0.4273 remains macroscopic,
         g_{12} = E_{13} - E_{12} ~ 0.1541, g_{13} ~ 0.1139.
     The dimensionless prominence ratio:
         Gamma_{11}(N) = g_{11} / max(g_{10}, g_{12}) ~ 2.77
     confirms that 11|12 is by far the most distinguished spectral gap.
  3. Reorientation from Individual Modes to Spectral Subspace Projectors:
     In Cell 92, the Ritz residual ratio r_{11} / delta_{11} was confounded by
     the collapsing internal gap delta_{11} = E_{11} - E_{10} ~ 10^{-2}.
     The true invariant object is the spectral projector:
         P_K^{(N)} = sum_{j=0}^K u_j^{(N)} (u_j^{(N)})^T
     which is completely invariant to internal cluster rotations.
  4. Strategic Shift to Milestone M39:
     Push dimensions to N in {44, 48, 52, 56, 60, 64}.
     Audit Cauchy convergence ||Delta P_K||_{op} for K in {10, 11, 12, 13},
     track Gamma_{11}(N), contrast competing boundary hypotheses
     H_{gap}(10) vs H_{gap}(11), and evaluate complementary subspace tilt.

THE FOUR INVESTIGATIVE TESTS OF CELL 93:
----------------------------------------
- Test A: Boundary Gaps & Dimensionless Prominence Ratio Gamma_{11}(N)
- Test B: Spectral Projector Cauchy Convergence ||Delta P_K||_{op} and ||Delta P_K||_F
- Test C: Competing Boundary Hypotheses Audit: H_{gap}(10) vs H_{gap}(11)
- Test D: Subspace Complementary Projection ||(I - P_K) P_K^{(N+4)}|| and Subspace Tilt

OUTPUT CONSTRAINTS:
-------------------
- Dispassionate, dry, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 93 EXECUTION COMPLETE.
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

# Extended dimensions sweep up to N = 64
N_LIST = [44, 48, 52, 56, 60, 64]
BOUNDARY_MODES = [10, 11, 12, 13, 14]
PROJECTOR_CUTOFFS = [10, 11, 12, 13]


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

    # sigma_min is sqrt(min eigenvalue)
    min_eval = max(mp.mpf('0'), min(evals_mtm))
    cos_theta_max = mp.sqrt(min_eval)
    cos_theta_max = min(mp.mpf('1'), cos_theta_max)

    sin_theta_max = mp.sqrt(max(mp.mpf('0'), mp.mpf('1') - cos_theta_max * cos_theta_max))

    return cos_theta_max, sin_theta_max


# -----------------------------------------------------------------------------
# Main Execution Function
# -----------------------------------------------------------------------------

def run_cell93() -> None:
    print("=" * 135)
    print("CELL 93 — SPECTRAL PROJECTOR CONVERGENCE ||P_K^{(N+4)} - P_K^{(N)}||,")
    print("          BOUNDARY GAP PROMINENCE GAMMA_11, AND EXTENDED SWEEPS (N <= 64)")
    print("=" * 135)
    print(f"Configuration: c = {C_PARAM}, L = log(c) = {mp.nstr(L_PARAM, 12)}, T = {T_PARAM}")
    print(f"Working Precision: {mp.mp.dps} decimal digits")
    print(f"Dimensions Sweep: N in {N_LIST}")
    print(f"Projector Cutoffs Focus: K in {PROJECTOR_CUTOFFS}")
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
        }

    total_elapsed = time.time() - start_total_time
    print(f"\nAll dimensions processed in {total_elapsed:.2f}s.")

    step_pairs = [(N_LIST[i], N_LIST[i + 1]) for i in range(len(N_LIST) - 1)]

    # =========================================================================
    # TEST A: Boundary Gaps & Dimensionless Prominence Ratio Gamma_11
    # =========================================================================
    print("\n" + "=" * 135)
    print("TEST A: Boundary Gaps g_j and Dimensionless Prominence Ratio Gamma_{11}(N) = g_{11} / max(g_{10}, g_{12})")
    print("Auditing Whether the 11|12 Boundary Remains the Structurally Distinguished Interface")
    print("=" * 135)

    print(f"{'N':>3} | {'E_{10}':>12} | {'E_{11}':>12} | {'E_{12}':>12} | {'E_{13}':>12} | {'g_{10}':>12} | {'g_{11}':>12} | {'g_{12}':>12} | {'Gamma_{11}':>12}")
    print("-" * 135)

    for N in N_LIST:
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

        print(f"{N:>3} | {mp.nstr(e10, 7):>12} | {mp.nstr(e11, 7):>12} | {mp.nstr(e12, 7):>12} | {mp.nstr(e13, 7):>12} | {mp.nstr(g10, 7):>12} | {mp.nstr(g11, 7):>12} | {mp.nstr(g12, 7):>12} | {mp.nstr(gamma11, 6):>12}")

    # =========================================================================
    # TEST B: Spectral Projector Cauchy Convergence
    # =========================================================================
    print("\n" + "=" * 135)
    print("TEST B: Spectral Projector Cauchy Differences ||P_K^{(N+4)} - P_K^{(N)}||_{op} and ||Delta P_K||_F")
    print("Auditing Subspace Stabilization across Cutoffs K in {10, 11, 12, 13}")
    print("=" * 135)

    header_parts = []
    for K in PROJECTOR_CUTOFFS:
        header_parts.append(f"||DP_{K}||_op")
        header_parts.append(f"||DP_{K}||_F")
    hdr_str = " | ".join(f"{h:>12}" for h in header_parts)
    print(f"{'Step (N -> N+4)':>16} | {hdr_str}")
    print("-" * 135)

    for n_low, n_high in step_pairs:
        p_dict_low = all_data[n_low]['projectors']
        p_dict_high = all_data[n_high]['projectors']

        vals = []
        for K in PROJECTOR_CUTOFFS:
            norm_op, norm_F = compute_projector_difference(p_dict_low[K], n_low, p_dict_high[K], n_high)
            vals.append(mp.nstr(norm_op, 6))
            vals.append(mp.nstr(norm_F, 6))

        row_str = " | ".join(f"{v:>12}" for v in vals)
        print(f"{f'{n_low} -> {n_high}':>16} | {row_str}")

    # =========================================================================
    # TEST C: Competing Boundary Hypotheses Audit: H_gap(10) vs H_gap(11)
    # =========================================================================
    print("\n" + "=" * 135)
    print("TEST C: Competing Boundary Hypotheses Audit: H_{gap}(10) vs H_{gap}(11)")
    print("Tracking the Collapse of g_{10} = E_{11} - E_{10} vs the Survival of g_{11} = E_{12} - E_{11}")
    print("=" * 135)

    print(f"{'N':>3} | {'g_{10}':>14} | {'g_{10} step ratio':>18} | {'g_{11}':>14} | {'g_{11} step ratio':>18} | {'g_{11} / g_{10}':>16}")
    print("-" * 135)

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

        print(f"{N:>3} | {mp.nstr(g10, 7):>14} | {step_r_g10:>18} | {mp.nstr(g11, 7):>14} | {step_r_g11:>18} | {mp.nstr(ratio_g, 6):>16}")

    # =========================================================================
    # TEST D: Subspace Complementary Projection & Principal Angles
    # =========================================================================
    print("\n" + "=" * 135)
    print("TEST D: Subspace Tilt sin theta_{max} = ||(I - P_K^{(N)}) P_K^{(N+4)}||_{op} and cos theta_{max} = sigma_{min}(U_K^T U_K)")
    print("Auditing Whether the Low-Energy Subspaces Stabilize as Invariant Geometric Objects")
    print("=" * 135)

    print(f"{'Step (N -> N+4)':>16} | {'sin theta(K=10)':>16} | {'cos theta(K=10)':>16} | {'sin theta(K=11)':>16} | {'cos theta(K=11)':>16} | {'sin theta(K=12)':>16} | {'cos theta(K=12)':>16}")
    print("-" * 135)

    for n_low, n_high in step_pairs:
        v_low = all_data[n_low]['V_e']
        v_high = all_data[n_high]['V_e']

        cos10, sin10 = compute_subspace_tilt(v_low, n_low, v_high, n_high, 10)
        cos11, sin11 = compute_subspace_tilt(v_low, n_low, v_high, n_high, 11)
        cos12, sin12 = compute_subspace_tilt(v_low, n_low, v_high, n_high, 12)

        print(f"{f'{n_low} -> {n_high}':>16} | {mp.nstr(sin10, 6):>16} | {mp.nstr(cos10, 6):>16} | {mp.nstr(sin11, 6):>16} | {mp.nstr(cos11, 6):>16} | {mp.nstr(sin12, 6):>16} | {mp.nstr(cos12, 6):>16}")

    print("\n" + "=" * 80)
    print("CELL 93 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell93()
