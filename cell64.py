#!/usr/bin/env python3
"""
================================================================================
CELL 64 — HIGH-PRECISION SCHUR COMPLEMENT BLOCK POSITIVITY VERIFICATION
================================================================================

PURPOSE:
--------
Execute the Stage III high-precision Schur positivity protocol specified in ROADMAP.md:
Following the Cell 63 diagnosis that generalized eigenvalue whitening
(Q_-^{-1/2} Q_pos Q_-^{-1/2}) becomes numerically ill-posed for N >= 12,
Cell 64 investigates finite-rank Weil positivity directly via the symmetric
block decomposition of the canonical Weil matrix:

    Q_Weil = [ A    B   ]
             [ B^T  C   ]

By the Schur complement theorem, Q_Weil is strictly positive definite if and
only if:
    1. The high-mode block is positive definite: C > 0.
    2. The low-mode Schur complement is positive definite:
           S_low = A - B * C^{-1} * B^T > 0.

This formulation avoids inverting the singular Gram matrix Q_-^{(N)} entirely
and operates directly on the native symmetric Galerkin operator at high precision.

PROTOCOL:
---------
1. Canonical Basis: Project the full Galerkin matrix Q_full to the canonical
   even v-basis of size (N+1) using the orthonormal parity matrix E:
       Q_Weil = E^T * Q_full * E.
2. Diagnostic 1 (High-Mode Block C Definiteness):
   Partition at m_cut = 3: A is 3x3, B is 3x(N-2), and C is (N-2)x(N-2).
   Compute the full spectrum, min eig(C), condition number kappa(C), and an
   exact LDL^T factorization of C, auditing all diagonal pivots D_{ii}(C) > 0
   and computing the relative backward error ||C - L D L^T||_inf / ||C||_inf.
3. Diagnostic 2 (3x3 Low-Mode Schur Complement S_low):
   Compute S_low = A - B * C^{-1} * B^T at 80-digit precision.
   Evaluate its spectrum (sigma_0 <= sigma_1 <= sigma_2), condition number,
   and LDL^T factorization pivots.
4. Diagnostic 3 (Cutoff Sensitivity Sweep at N = 24):
   Sweep the partition cutoff m_cut in {1, 2, 3, 4, 6, 8, 12 = N/2}.
   Test how min eig(C) and min eig(S) scale as the high-mode boundary shifts.
   At m_cut = 12 (barrier top), test whether the scattering continuum block
   C_scatt possesses an O(1) gap (lambda_min >= 1.30), connecting to
   Proposition 8.6 of Paper NR2.
5. Synthesis Table & Clean Termination Sentinel.

DIMENSIONS:
-----------
Sweep N in {4, 8, 12, 16, 20, 24} for cutoff c = 13 (L = log 13, T = 400).
================================================================================
"""

from __future__ import annotations

import time
import mpmath as mp

from connes_cvs import build_galerkin_matrix


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

mp.mp.dps = 80

C_PARAM = 13
L_PARAM = mp.log(C_PARAM)
T_PARAM = 400
GROUND_DPS = 80

N_LIST = [4, 8, 12, 16, 20, 24]
N_SWEEP_CUTOFF = 24
CUTOFF_LIST = [1, 2, 3, 4, 6, 8, 12]


# -----------------------------------------------------------------------------
# Canonical Parity Basis Transformation (c in R^{2N+1} -> v in R^{N+1})
# -----------------------------------------------------------------------------

def canonical_even_basis(N: int) -> mp.matrix:
    """
    Construct the (2N+1) x (N+1) orthonormal matrix E mapping canonical v-basis
    v in R^{N+1} to the full exponential basis c in R^{2N+1}: c = E v.
    """
    dim_full = 2 * N + 1
    dim_can = N + 1
    E = mp.matrix(dim_full, dim_can)

    E[N, 0] = mp.mpf(1)
    inv_sqrt2 = 1 / mp.sqrt(2)
    for m in range(1, dim_can):
        E[N + m, m] = inv_sqrt2
        E[N - m, m] = inv_sqrt2

    return E


# -----------------------------------------------------------------------------
# Exact LDL^T Symmetric Factorization
# -----------------------------------------------------------------------------

def ldl_factorization(A: mp.matrix) -> tuple[mp.matrix, list[mp.mpf], mp.mpf]:
    """
    Compute the LDL^T factorization of a real symmetric matrix A:
        A = L * diag(D) * L^T,
    where L is unit lower triangular and D is a vector of diagonal pivots.
    Returns:
        L: unit lower triangular matrix
        D: list of diagonal pivots
        rel_residual: ||A - L * D * L^T||_inf / ||A||_inf
    """
    n = A.rows
    L = mp.matrix(n, n)
    D = [mp.mpf('0')] * n

    for i in range(n):
        L[i, i] = mp.mpf('1')

    for j in range(n):
        # Compute D[j]
        d_val = A[j, j]
        for k in range(j):
            d_val -= (L[j, k] ** 2) * D[k]
        D[j] = d_val

        # Abort if pivot lacks sufficient separation above numerical floor
        pivot_tol = mp.mpf('1e-75')
        if abs(D[j]) < pivot_tol:
            raise ZeroDivisionError(
                f"Pivoting breakdown at index {j}: |D[{j}]| = {abs(D[j])} < {pivot_tol}. "
                "Insufficient spectral separation for stable LDL^T factorization."
            )
        inv_d = 1 / D[j]

        # Compute L[i, j] for i > j
        for i in range(j + 1, n):
            l_val = A[i, j]
            for k in range(j):
                l_val -= L[i, k] * L[j, k] * D[k]
            L[i, j] = l_val * inv_d

    # Compute reconstruction and relative backward error
    norm_A = max(sum(abs(A[r, c]) for c in range(n)) for r in range(n))
    if norm_A == 0:
        norm_A = mp.mpf('1')

    # Reconstruct R = L * D * L^T
    max_diff = mp.mpf('0')
    for r in range(n):
        row_sum = mp.mpf('0')
        for c in range(n):
            recon = sum(L[r, k] * D[k] * L[c, k] for k in range(min(r, c) + 1))
            diff = abs(A[r, c] - recon)
            row_sum += diff
        if row_sum > max_diff:
            max_diff = row_sum

    rel_residual = max_diff / norm_A
    return L, D, rel_residual


# -----------------------------------------------------------------------------
# Main Execution Protocol
# -----------------------------------------------------------------------------

def main():
    print("=" * 80)
    print("CELL 64 — CERTIFIED POSITIVITY VIA SCHUR COMPLEMENT BLOCK DECOUPLING")
    print("=" * 80)
    print(f"Parameters: c = {C_PARAM}, L = {mp.nstr(L_PARAM, 15)}, T = {T_PARAM}, dps = {mp.mp.dps}")
    print("Primary Partition: m_cut = 3 (Low modes: {e_0, e_1, e_2}, High modes: {e_3, ..., e_N})\n")

    summary_records = []

    for N in N_LIST:
        t0_N = time.perf_counter()
        dim = N + 1
        print("-" * 80)
        print(f"DIMENSION N = {N} (Canonical Basis Dimension = {dim})")
        print("-" * 80)

        # 1. Build full Galerkin matrix and project to canonical v-basis
        t0_Q = time.perf_counter()
        Q_full = build_galerkin_matrix(c=C_PARAM, N=N, T=T_PARAM, dps=GROUND_DPS)
        E = canonical_even_basis(N)
        Q_weil = E.T * Q_full * E
        Q_weil = mp.mpf('0.5') * (Q_weil + Q_weil.T)
        t1_Q = time.perf_counter()
        print(f"  Matrix build elapsed time: {t1_Q - t0_Q:.3f} s")

        # Full matrix eigenvalues
        full_evals, _ = mp.eigsy(Q_weil)
        full_evals_sorted = sorted(full_evals)
        lam_0_full = full_evals_sorted[0]
        lam_1_full = full_evals_sorted[1]
        print(f"  Full matrix Q_Weil smallest eigenvalue lambda_0 = {mp.nstr(lam_0_full, 25)}")
        print(f"  Full matrix Q_Weil second eigenvalue   lambda_1 = {mp.nstr(lam_1_full, 15)}")

        # ---------------------------------------------------------------------
        # Diagnostic 1: High-Mode Block C Definiteness (m_cut = 3)
        # ---------------------------------------------------------------------
        if dim > 3:
            m_cut = 3
            A = Q_weil[:m_cut, :m_cut]
            B = Q_weil[:m_cut, m_cut:]
            C = Q_weil[m_cut:, m_cut:]

            C_sym = mp.mpf('0.5') * (C + C.T)
            c_evals, _ = mp.eigsy(C_sym)
            c_evals_sorted = sorted(c_evals)
            lam_min_C = c_evals_sorted[0]
            lam_max_C = c_evals_sorted[-1]
            cond_C = lam_max_C / lam_min_C if lam_min_C > 0 else mp.inf

            # LDL^T Factorization of C
            L_c, D_c, res_C = ldl_factorization(C_sym)
            min_pivot_C = min(D_c)
            pivots_C_positive = all(p > 0 for p in D_c)

            print("\nDIAGNOSTIC 1: High-Mode Block C Audit (dimension: {0}x{0}):".format(C.rows))
            print(f"  Smallest eigenvalue min eig(C) = {mp.nstr(lam_min_C, 25)}")
            print(f"  Largest eigenvalue  max eig(C) = {mp.nstr(lam_max_C, 15)}")
            print(f"  Condition number kappa(C)      = {mp.nstr(cond_C, 12)}")
            print(f"  LDL^T factorization relative residual = {mp.nstr(res_C, 10)}")
            print(f"  Minimum diagonal pivot min(D_C)       = {mp.nstr(min_pivot_C, 25)}")
            print(f"  All pivots D_{{ii}}(C) strictly positive: {pivots_C_positive}")
            print(f"  Bottom 3 eigenvalues of C : {[mp.nstr(e, 8) for e in c_evals_sorted[:min(3, len(c_evals_sorted))]]}")

            # -----------------------------------------------------------------
            # Diagnostic 2: 3x3 Schur Complement S_low = A - B C^{-1} B^T
            # -----------------------------------------------------------------
            if lam_min_C > 0:
                C_inv = mp.inverse(C_sym)
                S_low = A - (B * C_inv * B.T)
                S_low_sym = mp.mpf('0.5') * (S_low + S_low.T)

                s_evals, _ = mp.eigsy(S_low_sym)
                s_evals_sorted = sorted(s_evals)
                lam_min_S = s_evals_sorted[0]
                lam_max_S = s_evals_sorted[-1]
                cond_S = lam_max_S / lam_min_S if lam_min_S > 0 else mp.inf

                # LDL^T Factorization of 3x3 S_low
                L_s, D_s, res_S = ldl_factorization(S_low_sym)
                min_pivot_S = min(D_s)
                pivots_S_positive = all(p > 0 for p in D_s)

                print("\nDIAGNOSTIC 2: 3x3 Schur Complement S_low Audit:")
                print(f"  Eigenvalues of S_low: {[mp.nstr(e, 25) for e in s_evals_sorted]}")
                print(f"  Smallest eigenvalue min eig(S_low) = {mp.nstr(lam_min_S, 25)}")
                print(f"  Condition number kappa(S_low)      = {mp.nstr(cond_S, 12)}")
                print(f"  LDL^T factorization relative residual = {mp.nstr(res_S, 10)}")
                print(f"  Pivots of S_low [D_00, D_11, D_22] = {[mp.nstr(p, 20) for p in D_s]}")
                print(f"  All pivots D_{{ii}}(S_low) strictly positive: {pivots_S_positive}")
                print(f"  Schur complement condition satisfied (C > 0 and S_low > 0): {pivots_C_positive and pivots_S_positive}")
            else:
                lam_min_S = mp.mpf('nan')
                cond_S = mp.mpf('inf')
                min_pivot_S = mp.mpf('nan')
                pivots_S_positive = False
                res_S = mp.mpf('nan')
                print("  High-mode block C is NOT positive definite.")
        else:
            # Dimension <= 3 (N = 2, 4)
            lam_min_C = mp.mpf('nan')
            cond_C = mp.mpf('nan')
            min_pivot_C = mp.mpf('nan')
            pivots_C_positive = True
            res_C = mp.mpf('0')
            lam_min_S = lam_0_full
            cond_S = mp.mpf('1')
            min_pivot_S = lam_0_full
            pivots_S_positive = True
            res_S = mp.mpf('0')
            print("  N <= 2: Dimension <= 3 (pure low-mode system).")

        t1_N = time.perf_counter()
        print(f"\nTotal elapsed time for N = {N}: {t1_N - t0_N:.3f} s\n")

        summary_records.append({
            'N': N,
            'dim': dim,
            'lam_0_full': lam_0_full,
            'lam_min_C': lam_min_C,
            'cond_C': cond_C,
            'min_pivot_C': min_pivot_C,
            'pivots_C_pos': pivots_C_positive,
            'lam_min_S': lam_min_S,
            'cond_S': cond_S,
            'min_pivot_S': min_pivot_S,
            'pivots_S_pos': pivots_S_positive,
        })

    # -------------------------------------------------------------------------
    # Diagnostic 3: Cutoff Sensitivity Sweep at N = 24
    # -------------------------------------------------------------------------
    print("=" * 80)
    print(f"DIAGNOSTIC 3: Cutoff Sensitivity Sweep at N = {N_SWEEP_CUTOFF}")
    print("=" * 80)
    print(f"{'m_cut':>5} | {'dim(A)':>6} | {'dim(C)':>6} | {'min_eig(C)':>22} | {'min_eig(S)':>22} | {'All Pivots > 0':>14}")
    print("-" * 86)

    dim_24 = N_SWEEP_CUTOFF + 1
    Q_full_24 = build_galerkin_matrix(c=C_PARAM, N=N_SWEEP_CUTOFF, T=T_PARAM, dps=GROUND_DPS)
    E_24 = canonical_even_basis(N_SWEEP_CUTOFF)
    Q_weil_24 = E_24.T * Q_full_24 * E_24
    Q_weil_24 = mp.mpf('0.5') * (Q_weil_24 + Q_weil_24.T)

    for m_cut in CUTOFF_LIST:
        if m_cut >= dim_24:
            continue
        A_cut = Q_weil_24[:m_cut, :m_cut]
        B_cut = Q_weil_24[:m_cut, m_cut:]
        C_cut = Q_weil_24[m_cut:, m_cut:]

        C_cut_sym = mp.mpf('0.5') * (C_cut + C_cut.T)
        c_vals, _ = mp.eigsy(C_cut_sym)
        c_min = min(c_vals)

        _, D_c_cut, _ = ldl_factorization(C_cut_sym)
        c_pos = all(p > 0 for p in D_c_cut)

        if c_min > 0:
            C_inv_cut = mp.inverse(C_cut_sym)
            S_cut = A_cut - (B_cut * C_inv_cut * B_cut.T)
            S_cut_sym = mp.mpf('0.5') * (S_cut + S_cut.T)
            s_vals, _ = mp.eigsy(S_cut_sym)
            s_min = min(s_vals)

            _, D_s_cut, _ = ldl_factorization(S_cut_sym)
            s_pos = all(p > 0 for p in D_s_cut)
            all_pos = c_pos and s_pos
        else:
            s_min = mp.mpf('nan')
            all_pos = False

        c_str = mp.nstr(c_min, 12)
        s_str = mp.nstr(s_min, 12)
        print(f"{m_cut:5d} | {m_cut:6d} | {C_cut.rows:6d} | {c_str:>22} | {s_str:>22} | {str(all_pos):>14}")

    # -------------------------------------------------------------------------
    # Multi-Dimension Synthesis Table
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("MULTI-DIMENSION SYNTHESIS TABLE (c = 13, Primary m_cut = 3)")
    print("=" * 80)
    header = (
        f"{'N':>3} | {'lambda_0(Q_Weil)':>22} | {'min_eig(C)':>20} | "
        f"{'min_eig(S_low)':>22} | {'C > 0':>6} | {'S > 0':>6}"
    )
    print(header)
    print("-" * len(header))
    for rec in summary_records:
        n_val = rec['N']
        l0_val = mp.nstr(rec['lam_0_full'], 12)
        c_val = mp.nstr(rec['lam_min_C'], 12)
        s_val = mp.nstr(rec['lam_min_S'], 12)
        c_pos = str(rec['pivots_C_pos'])
        s_pos = str(rec['pivots_S_pos'])
        print(f"{n_val:3d} | {l0_val:>22} | {c_val:>20} | {s_val:>22} | {c_pos:>6} | {s_pos:>6}")

    print("=" * 80)
    print("CELL 64 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
