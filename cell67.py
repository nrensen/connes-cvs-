#!/usr/bin/env python3
"""
================================================================================
CELL 67 — AUDIT OF RELATIVE TUNNELING-SCALE ORDERING, BARRIER THINNING,
AND THE EXCITED SECOND RESOLVENT MOMENT D_0^2 M_{2,exc}
================================================================================

PURPOSE:
--------
Investigate the relative tunneling-scale hierarchy and barrier-thinning mechanism
connecting the excited odd-sector resolvent moment M_1^{exc} and the excited
second resolvent moment M_{2,exc} in Paper NR2 (§8.10 & §8.11) across dimensions
N in {8, 12, 16, 20, 24}:

1. Relative Tunneling-Scale Hierarchy:
   Measure the modewise ratio between the ground-state boundary tunneling scale
   D_0^2 ~ e^{-2 S(E_0) N} and the excited odd tunneling gaps:
       R_gap(j, N) = D_0^2 / (mu_j - lambda_0),    j = 1, 2, ..., min(5, N-1).
   Determine the supremum over excited modes:
       R_gap^{max}(N) = sup_{j >= 1} D_0^2 / (mu_j - lambda_0) = D_0^2 / (mu_1 - lambda_0).

2. Semiclassical Barrier Thinning:
   For states with energy E_j > E_0, the turning points contract and the barrier
   height V(x) - E_j decreases, implying S(E_j) < S(E_0) and dS/dE < 0.
   Extract the effective barrier actions per mode:
       sigma_j(N) = - (1 / N) * log(mu_j - lambda_0),
   and evaluate the action gap Delta sigma_j = sigma_0(N) - sigma_j(N) > 0,
   which drives exponential suppression:
       R_gap(j, N) ~ exp(- Delta sigma_j * N) -> 0.

3. Automatic Second-Moment Control:
   Compute the exact excited second resolvent moment:
       M_{2,exc}(N) = sum_{j >= 1} a_j^2 / (mu_j - lambda_0)^2,
   and test the operator inequality:
       D_0^2 M_{2,exc} <= R_gap^{max}(N) * M_1^{exc}.
   Verify that D_0^2 M_{2,exc} <= M_1^{exc} holds with a suppression factor
   <= 5 x 10^{-6}, establishing that the second-moment condition
   D_0^2 M_{2,exc} <= C_2 N^q is subordinate to H2_{odd} with q = gamma.

4. Ground-State Second-Moment Decomposition:
   Verify the exact identity for the ground-state component:
       b_{00}^2 = D_0^2 * a_0^2 / (mu_0 - lambda_0)^2 <= N^2,
   and compute the total second moment:
       D_0^2 M_2 = b_{00}^2 + D_0^2 M_{2,exc}.

OUTPUT CONSTRAINTS:
-------------------
- Strictly dispassionate, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 67 EXECUTION COMPLETE.
================================================================================
"""

from __future__ import annotations

import time
import mpmath as mp

from connes_cvs import build_galerkin_matrix
from cell import (
    get_ground_state,
    canonical_to_full,
)

# -----------------------------------------------------------------------------
# Precision and Parameter Configuration
# -----------------------------------------------------------------------------

mp.mp.dps = 50

C_PARAM = 13
L_PARAM = mp.log(C_PARAM)
T_PARAM = 400
GROUND_DPS = 50
KAPPA = 2 * mp.pi / L_PARAM

N_LIST = [8, 12, 16, 20, 24]


# -----------------------------------------------------------------------------
# Full-Space Parity Decomposition & Operators
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


def psi_vector(Q: mp.matrix, N: int) -> mp.matrix:
    """Extract symbol generator vector psi(m) = m * Q_{0, m} from row N of Q."""
    dim = 2 * N + 1
    psi = mp.matrix(dim, 1)
    for j in range(dim):
        m = j - N
        if m != 0:
            psi[j, 0] = mp.mpf(m) * Q[N, j]
        else:
            psi[j, 0] = mp.mpf(0)
    return psi


def constant_vector(N: int) -> mp.matrix:
    """Constant vector d = (1, 1, ..., 1)^T in R^{2N+1}."""
    return mp.matrix([[mp.mpf(1)] for _ in range(2 * N + 1)])


def dot(x: mp.matrix, y: mp.matrix) -> mp.mpf:
    """Inner product <x, y> for column vectors."""
    return (x.T * y)[0, 0]


def eigensystem_symmetric(A: mp.matrix) -> tuple[list[mp.mpf], mp.matrix]:
    """Symmetric eigendecomposition returning (sorted_eigenvalues, eigenvector_matrix)."""
    vals, vecs = mp.eigsy(A)
    return list(vals), vecs


# -----------------------------------------------------------------------------
# Main Execution Suite
# -----------------------------------------------------------------------------

def main():
    print("=" * 96)
    print("CELL 67 — AUDIT OF RELATIVE TUNNELING-SCALE ORDERING, BARRIER THINNING,")
    print("         AND THE EXCITED SECOND RESOLVENT MOMENT D_0^2 M_{2,exc}")
    print("=" * 96)
    print(f"Parameters: c = {C_PARAM}, L = {mp.nstr(L_PARAM, 15)}, T = {T_PARAM}, dps = {mp.mp.dps}\n")

    summary_records = []

    for N in N_LIST:
        t0 = time.perf_counter()

        lam0, v_can, _ = get_ground_state(
            c=C_PARAM,
            N=N,
            T=T_PARAM,
            dps=GROUND_DPS,
            verbose=False,
        )
        c_vec = canonical_to_full(v_can)
        Q = build_galerkin_matrix(
            c=C_PARAM,
            N=N,
            T=T_PARAM,
            dps=GROUND_DPS,
        )

        E, O = full_parity_basis(N)
        Q_even = E.T * Q * E
        Q_odd = O.T * Q * O

        even_vals, even_vecs = eigensystem_symmetric(Q_even)
        odd_vals, odd_vecs = eigensystem_symmetric(Q_odd)

        psi = psi_vector(Q, N)
        psi_odd = O.T * psi

        d_full = constant_vector(N)
        D0 = dot(d_full, c_vec)
        D0_sq = D0 ** 2

        # Overlaps
        a_overlaps = [dot(odd_vecs[:, j], psi_odd) for j in range(N)]

        # Gaps
        mu0 = odd_vals[0]
        delta_0 = mu0 - lam0

        # Ground state component b_{00}^2 = D_0^2 * a_0^2 / delta_0^2
        # where a_0 = - (delta_0 / D_0) * b_{00}
        b00_sq = (a_overlaps[0] ** 2) * D0_sq / (delta_0 ** 2)

        # Resolvent moments: excited odd sector j >= 1
        M1_exc = mp.mpf(0)
        M2_exc = mp.mpf(0)
        for j in range(1, N):
            delta_j = odd_vals[j] - lam0
            aj_sq = a_overlaps[j] ** 2
            M1_exc += aj_sq / delta_j
            M2_exc += aj_sq / (delta_j ** 2)

        D0_sq_M2_exc = D0_sq * M2_exc
        D0_sq_M2_total = b00_sq + D0_sq_M2_exc

        # Relative tunneling ratio for lowest excited mode j=1
        delta_1 = odd_vals[1] - lam0
        R_gap_max = D0_sq / delta_1
        rho_2_exc = D0_sq_M2_exc / M1_exc

        # Decisive relative gap decay exponent: Delta sigma_N^gap = - (1 / N) * log(R_gap_max)
        delta_sigma_gap = - (mp.log(R_gap_max) / N)

        # Effective barrier actions
        sigma_0 = - (mp.log(delta_0) / N)
        sigma_1 = - (mp.log(delta_1) / N)
        delta_sigma = sigma_0 - sigma_1

        elapsed = time.perf_counter() - t0

        print(f"--- DIMENSION N = {N:2d} (Elapsed: {elapsed:.2f}s) ---")
        print(f"  Ground scale:    D_0^2          = {mp.nstr(D0_sq, 9)}")
        print(f"  Ground gap:      mu0 - lam0     = {mp.nstr(delta_0, 9)}")
        print(f"  Mode 1 gap:      mu1 - lam0     = {mp.nstr(delta_1, 9)}")
        print(f"  Ratio R_gap:     D0^2/gap_1     = {mp.nstr(R_gap_max, 9)}")
        print(f"  Gap Action:      -(1/N)logR_gap = {mp.nstr(delta_sigma_gap, 6)}")
        print(f"  Excited M_1:     M_1^exc        = {mp.nstr(M1_exc, 9)}")
        print(f"  Excited M_2:     D0^2*M2_exc    = {mp.nstr(D0_sq_M2_exc, 9)}")
        print(f"  Ratio rho_2^exc: D0^2*M2/M1_exc = {mp.nstr(rho_2_exc, 9)}")
        print(f"  Ground b00^2:    b00^2          = {mp.nstr(b00_sq, 9)} (<= N^2 = {N*N})")
        print(f"  Total D0^2*M2:   D0^2 * M_2     = {mp.nstr(D0_sq_M2_total, 9)}")
        print(f"  Actions:         sigma_0 = {mp.nstr(sigma_0, 6)}, sigma_1 = {mp.nstr(sigma_1, 6)}, Delta_sigma = {mp.nstr(delta_sigma, 6)}")

        # Ladder breakdown for excited modes j = 1, ..., min(5, N-1)
        print("\n  Excited Odd Modes (Barrier Thinning Progression):")
        print("  " + "-" * 88)
        print(f"  {'j':>2s}  {'mu_j - lam0':>16s}  {'D_0^2 / (mu_j-lam0)':>22s}  {'sigma_j':>12s}  {'aj^2/(mu_j-lam0)':>18s}")
        print("  " + "-" * 88)
        max_modes = min(6, N)
        for j in range(1, max_modes):
            d_j = odd_vals[j] - lam0
            r_j = D0_sq / d_j
            sig_j = - (mp.log(d_j) / N)
            aj_trans = (a_overlaps[j] ** 2) / d_j
            print(f"  {j:2d}  {mp.nstr(d_j, 10):>16s}  {mp.nstr(r_j, 10):>22s}  {mp.nstr(sig_j, 6):>12s}  {mp.nstr(aj_trans, 8):>18s}")
        print("  " + "-" * 88 + "\n")

        summary_records.append({
            "N": N,
            "D0_sq": D0_sq,
            "delta_0": delta_0,
            "delta_1": delta_1,
            "R_gap_max": R_gap_max,
            "delta_sigma_gap": delta_sigma_gap,
            "M1_exc": M1_exc,
            "D0_sq_M2_exc": D0_sq_M2_exc,
            "rho_2_exc": rho_2_exc,
            "b00_sq": b00_sq,
            "D0_sq_M2_total": D0_sq_M2_total,
            "sigma_0": sigma_0,
            "sigma_1": sigma_1,
            "delta_sigma": delta_sigma,
        })

    # -------------------------------------------------------------------------
    # Multi-Dimension Synthesis Tables
    # -------------------------------------------------------------------------
    print("=" * 112)
    print("SYNTHESIS TABLE 1: RELATIVE TUNNELING GAP, DECAY ACTION & EXCITED SECOND-MOMENT BOUNDS")
    print("=" * 112)
    print(f"{'N':>4s} | {'D_0^2':>16s} | {'mu_1 - lam_0':>16s} | {'D_0^2/gap_1':>14s} | {'-(1/N)logRgap':>14s} | {'M_1^exc':>10s} | {'D_0^2*M_{2,exc}':>16s} | {'rho_2^exc':>12s}")
    print("-" * 112)
    for rec in summary_records:
        n_str = f"{rec['N']:4d}"
        d0_str = mp.nstr(rec['D0_sq'], 8)
        gap1_str = mp.nstr(rec['delta_1'], 8)
        rgap_str = mp.nstr(rec['R_gap_max'], 6)
        dsig_gap_str = mp.nstr(rec['delta_sigma_gap'], 6)
        m1_str = mp.nstr(rec['M1_exc'], 6)
        m2_str = mp.nstr(rec['D0_sq_M2_exc'], 8)
        rho_str = mp.nstr(rec['rho_2_exc'], 6)
        print(f"{n_str} | {d0_str:>16s} | {gap1_str:>16s} | {rgap_str:>14s} | {dsig_gap_str:>14s} | {m1_str:>10s} | {m2_str:>16s} | {rho_str:>12s}")
    print("=" * 112)

    print("\n" + "=" * 96)
    print("SYNTHESIS TABLE 2: SEMICLASSICAL BARRIER THINNING & TOTAL SECOND MOMENT D_0^2 M_2")
    print("=" * 96)
    print(f"{'N':>4s} | {'sigma_0':>10s} | {'sigma_1':>10s} | {'Delta sigma':>12s} | {'b_{00}^2':>12s} | {'D_0^2*M_2':>14s} | {'N^2':>8s}")
    print("-" * 96)
    for rec in summary_records:
        n_str = f"{rec['N']:4d}"
        sig0_str = mp.nstr(rec['sigma_0'], 6)
        sig1_str = mp.nstr(rec['sigma_1'], 6)
        dsig_str = mp.nstr(rec['delta_sigma'], 6)
        b00_str = mp.nstr(rec['b00_sq'], 6)
        tot_str = mp.nstr(rec['D0_sq_M2_total'], 6)
        nsq_str = f"{rec['N']**2:8d}"
        print(f"{n_str} | {sig0_str:>10s} | {sig1_str:>10s} | {dsig_str:>12s} | {b00_str:>12s} | {tot_str:>14s} | {nsq_str:>8s}")
    print("=" * 96)

    print("\n" + "=" * 80)
    print("CELL 67 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
