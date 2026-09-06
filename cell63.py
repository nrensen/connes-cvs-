#!/usr/bin/env python3
"""
================================================================================
CELL 63 — FINITE-RANK OPERATOR DOMINANCE & SPECTRUM OF THE NEGATIVE
ARCHIMEDEAN OPERATOR Q_arch^{(-)}
================================================================================

PURPOSE:
--------
Execute the canonical Stage III dominance reconnaissance specified in ROADMAP.md:
1. Decompose the Weil quadratic form into positive and negative operators:
       Q_Weil = Q_positive - Q_arch^{(-)},
   where:
       [Q_arch^{(-)}]_{mn} = (1 / pi) * int_0^{r_*} |h_+(r)| * phi_m(r) * phi_n(r) dr >= 0
   is the compact negative Archimedean Gram matrix supported on [0, r_*], with
   r_* approx 6.289835 denoting the first positive root of the Archimedean weight
   h_+(r) = Re psi(1/4 + i r / 2) - log pi.
2. Form the positive operator:
       Q_positive = Q_Weil + Q_arch^{(-)} = Q_pole + Q_prime + Q_arch^{(+)}.
3. Audit the Gram matrix definiteness, condition number, and effective rank of
   Q_- = Q_arch^{(-)} on the finite-dimensional canonical v-basis.
4. Solve the generalized eigenvalue problem:
       Q_positive x = lambda * Q_arch^{(-)} x.
   Evaluate the dominance margin mu_min = lambda_min - 1 to test whether
   Q_positive dominates Q_arch^{(-)} (lambda_min > 1 <===> Q_Weil > 0).
5. Execute the three core diagnostics:
   - Diagnostic 1: Spectrum, conditioning, and effective rank of Q_arch^{(-)}.
   - Diagnostic 2: Modal energy projection of the dangerous eigenvector x_min
     onto the 3-mode sector span{e_0, e_1, e_2} vs high modes span{e_3, ..., e_N}.
   - Diagnostic 3: Schur complement test of Q_Weil partitioned into 3x3 low-mode
     block A and high-mode block C: S_low = A - B * C^{-1} * B^T.

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

mp.mp.dps = 50

C_PARAM = 13
L_PARAM = mp.log(C_PARAM)
T_PARAM = 400
GROUND_DPS = 50

N_LIST = [4, 8, 12, 16, 20, 24]
QUAD_ORDER = 80


# -----------------------------------------------------------------------------
# Archimedean Multiplier and Positive Zero r_*
# -----------------------------------------------------------------------------

def h_plus(r: mp.mpf) -> mp.mpf:
    """Archimedean weight h_+(r) = Re psi(1/4 + i r / 2) - log pi."""
    z = mp.mpc(mp.mpf('0.25'), mp.mpf('0.5') * r)
    return mp.re(mp.digamma(z)) - mp.log(mp.pi)


def compute_r_star() -> mp.mpf:
    """Compute the unique positive root r_* of h_+(r) = 0 via root finding."""
    return mp.findroot(h_plus, mp.mpf('6.289835499'))


# -----------------------------------------------------------------------------
# Canonical Fourier Basis Amplitudes phi_m(r)
# -----------------------------------------------------------------------------

def phi_basis(m: int, r: mp.mpf, L: mp.mpf) -> mp.mpf:
    """
    Canonical Fourier basis amplitude phi_m(r) on R.
    phi_0(r) = (2 / sqrt(L)) * (sin(r L / 2) / r)
    phi_m(r) = (2 * sqrt(2) / sqrt(L)) * (r * sin(r L / 2) / (r^2 - a_m^2))  for m >= 1.
    """
    half_L = L / 2
    if m == 0:
        if abs(r) < mp.mpf('1e-25'):
            return mp.sqrt(L)
        return (2 / mp.sqrt(L)) * (mp.sin(r * half_L) / r)

    a_m = 2 * mp.pi * m / L
    diff_pos = r - a_m
    diff_neg = r + a_m

    # Removable singularity at r = a_m
    if abs(diff_pos) < mp.mpf('1e-12'):
        u = diff_pos * half_L
        # Taylor expansion of sin(u) / u
        sinc_u = 1 - (u**2) / 6 + (u**4) / 120 - (u**6) / 5040
        geom = (a_m + diff_pos) / (2 * a_m + diff_pos)
        sgn = (-1) ** m
        val = sgn * half_L * geom * sinc_u
        return (2 * mp.sqrt(2) / mp.sqrt(L)) * val

    # Removable singularity at r = -a_m
    if abs(diff_neg) < mp.mpf('1e-12'):
        u = diff_neg * half_L
        sinc_u = 1 - (u**2) / 6 + (u**4) / 120 - (u**6) / 5040
        geom = (-a_m + diff_neg) / (-2 * a_m + diff_neg)
        sgn = (-1) ** m
        val = sgn * half_L * geom * sinc_u
        return (2 * mp.sqrt(2) / mp.sqrt(L)) * val

    denom = r**2 - a_m**2
    return (2 * mp.sqrt(2) / mp.sqrt(L)) * (r * mp.sin(r * half_L) / denom)


# -----------------------------------------------------------------------------
# Gauss-Legendre Quadrature Generator (Golub-Welsch)
# -----------------------------------------------------------------------------

def gauss_legendre_nodes_weights(order: int, a: mp.mpf, b: mp.mpf) -> tuple[list[mp.mpf], list[mp.mpf]]:
    """
    Compute Gauss-Legendre quadrature nodes and weights on [a, b] using
    the Golub-Welsch tridiagonal eigenvalue method at current mpmath dps.
    """
    J = mp.matrix(order, order)
    for i in range(order - 1):
        k = i + 1
        b_k = mp.mpf(k) / mp.sqrt(4 * k * k - 1)
        J[i, i + 1] = b_k
        J[i + 1, i] = b_k

    nodes_std, V = mp.eigsy(J)

    # Map [-1, 1] to [a, b]
    mid = (b + a) / 2
    half_width = (b - a) / 2

    nodes = []
    weights = []
    for i in range(order):
        x_i = nodes_std[i]
        w_i = 2 * (V[0, i] ** 2)
        nodes.append(mid + half_width * x_i)
        weights.append(half_width * w_i)

    return nodes, weights


# -----------------------------------------------------------------------------
# Negative Archimedean Gram Matrix Construction
# -----------------------------------------------------------------------------

def build_Q_arch_minus(N: int, r_star: mp.mpf, L: mp.mpf, order: int = QUAD_ORDER) -> mp.matrix:
    """
    Construct the (N+1) x (N+1) Gram matrix Q_arch^{(-)}:
        [Q_arch^{(-)}]_{mn} = (1 / pi) * int_0^{r_*} |h_+(r)| * phi_m(r) * phi_n(r) dr.
    """
    nodes, weights = gauss_legendre_nodes_weights(order, mp.mpf('0'), r_star)
    dim = N + 1
    Q_minus = mp.matrix(dim, dim)

    # Precompute basis amplitudes at quadrature nodes
    phi_table = []
    weights_h = []
    for k in range(order):
        rk = nodes[k]
        hk = abs(h_plus(rk))
        weights_h.append((weights[k] * hk) / mp.pi)
        phi_row = [phi_basis(m, rk, L) for m in range(dim)]
        phi_table.append(phi_row)

    # Accumulate Gram matrix
    for k in range(order):
        w_eff = weights_h[k]
        phi_k = phi_table[k]
        for m in range(dim):
            phi_km = phi_k[m]
            w_eff_phi_m = w_eff * phi_km
            for n in range(m, dim):
                Q_minus[m, n] += w_eff_phi_m * phi_k[n]

    # Symmetrize
    for m in range(dim):
        for n in range(m + 1, dim):
            val = Q_minus[m, n]
            Q_minus[n, m] = val

    return Q_minus


# -----------------------------------------------------------------------------
# Parity Basis Transformation
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
# Generalized Eigensolver with Nullspace Audit
# -----------------------------------------------------------------------------

def solve_generalized_eigenproblem(
    A: mp.matrix,
    B: mp.matrix,
    tol_rel: mp.mpf = mp.mpf('1e-35'),
) -> tuple[list[mp.mpf], mp.matrix, list[mp.mpf], int]:
    """
    Solve the generalized eigenvalue problem A x = lambda B x for symmetric A, B with B >= 0.
    Returns:
        evals: generalized eigenvalues sorted ascending
        evecs: generalized eigenvectors normalized to unit Euclidean length
        B_evals: eigenvalues of B sorted descending
        null_dim: dimension of numerical nullspace of B
    """
    dim = A.rows
    b_vals, V_b = mp.eigsy(B)

    # Sort B eigenvalues descending
    idx_b = sorted(range(dim), key=lambda i: b_vals[i], reverse=True)
    B_evals = [b_vals[i] for i in idx_b]
    V_sorted = mp.matrix(dim, dim)
    for j, i in enumerate(idx_b):
        for r in range(dim):
            V_sorted[r, j] = V_b[r, i]

    sigma_max = B_evals[0]
    rank = 0
    for s in B_evals:
        if s > tol_rel * sigma_max and s > 0:
            rank += 1
        else:
            break

    null_dim = dim - rank

    # Form transformation matrix S = V_r * diag(1 / sqrt(sigma_i)) of shape (dim, rank)
    S = mp.matrix(dim, rank)
    for j in range(rank):
        inv_sqrt = 1 / mp.sqrt(B_evals[j])
        for r in range(dim):
            S[r, j] = V_sorted[r, j] * inv_sqrt

    # Projected standard symmetric problem M = S^T * A * S
    M = S.T * A * S
    M = mp.mpf('0.5') * (M + M.T)

    m_vals, Y = mp.eigsy(M)

    # Sort generalized eigenvalues ascending
    idx_m = sorted(range(rank), key=lambda i: m_vals[i])
    evals = [m_vals[i] for i in idx_m]

    # Transform back to original coordinates x = S * y and normalize ||x||_2 = 1
    evecs = mp.matrix(dim, rank)
    for j, i in enumerate(idx_m):
        y_col = Y[:, i]
        x_col = S * y_col
        norm_x = mp.sqrt(sum(x_col[r, 0] ** 2 for r in range(dim)))
        for r in range(dim):
            evecs[r, j] = x_col[r, 0] / norm_x

    return evals, evecs, B_evals, null_dim


# -----------------------------------------------------------------------------
# Main Execution Protocol
# -----------------------------------------------------------------------------

def main():
    print("=" * 80)
    print("CELL 63 — FINITE-RANK OPERATOR DOMINANCE & SPECTRUM OF Q_arch^{(-)}")
    print("=" * 80)
    print(f"Parameters: c = {C_PARAM}, L = {mp.nstr(L_PARAM, 15)}, T = {T_PARAM}, dps = {mp.mp.dps}")

    t0_rstar = time.perf_counter()
    r_star = compute_r_star()
    t1_rstar = time.perf_counter()
    print(f"Computed Archimedean zero r_* = {mp.nstr(r_star, 25)} ({t1_rstar - t0_rstar:.3f} s)\n")

    summary_records = []

    for N in N_LIST:
        t0_N = time.perf_counter()
        dim = N + 1
        print("-" * 80)
        print(f"DIMENSION N = {N} (Canonical Basis Dimension = {dim})")
        print("-" * 80)

        # 1. Build negative Archimedean Gram matrix Q_arch^{(-)}
        t0_minus = time.perf_counter()
        Q_minus = build_Q_arch_minus(N, r_star, L_PARAM, order=QUAD_ORDER)
        t1_minus = time.perf_counter()

        # 2. Build full-space Galerkin matrix Q and project to canonical v-basis
        t0_Q = time.perf_counter()
        Q_full = build_galerkin_matrix(c=C_PARAM, N=N, T=T_PARAM, dps=GROUND_DPS)
        E = canonical_even_basis(N)
        Q_weil = E.T * Q_full * E
        Q_weil = mp.mpf('0.5') * (Q_weil + Q_weil.T)
        t1_Q = time.perf_counter()

        # 3. Form positive operator Q_positive = Q_weil + Q_minus
        Q_pos = Q_weil + Q_minus
        Q_pos = mp.mpf('0.5') * (Q_pos + Q_pos.T)

        # ---------------------------------------------------------------------
        # Diagnostic 1: Spectrum, Conditioning, and Effective Rank of Q_-
        # ---------------------------------------------------------------------
        evals_gen, evecs_gen, sigma_minus, null_dim = solve_generalized_eigenproblem(
            Q_pos, Q_minus, tol_rel=mp.mpf('1e-35')
        )

        sigma_max = sigma_minus[0]
        sigma_min = sigma_minus[-1]
        cond_minus = sigma_max / sigma_min if sigma_min > 0 else mp.inf

        print("DIAGNOSTIC 1: Gram Matrix Audit of Q_arch^{(-)}:")
        print(f"  Build time: Q_minus = {t1_minus - t0_minus:.3f} s, Q_weil = {t1_Q - t0_Q:.3f} s")
        print(f"  sigma_max(Q_-) = {mp.nstr(sigma_max, 15)}")
        print(f"  sigma_min(Q_-) = {mp.nstr(sigma_min, 15)}")
        print(f"  Condition number kappa(Q_-) = {mp.nstr(cond_minus, 10)}")
        print(f"  Numerical nullspace dimension = {null_dim} (Rank = {dim - null_dim} / {dim})")
        print(f"  Top 5 eigenvalues of Q_- : {[mp.nstr(s, 8) for s in sigma_minus[:min(5, dim)]]}")
        print(f"  Bottom 3 eigenvalues of Q_- : {[mp.nstr(s, 8) for s in sigma_minus[-min(3, dim):]]}\n")

        # ---------------------------------------------------------------------
        # Generalized Eigenvalues & Dominance Margin
        # ---------------------------------------------------------------------
        lam_min = evals_gen[0]
        margin_min = lam_min - 1
        x_min = evecs_gen[:, 0]

        # Residual verification: || Q_pos x_min - lambda_min Q_minus x_min ||_2
        res_vec = (Q_pos * x_min) - (lam_min * (Q_minus * x_min))
        res_norm = mp.sqrt(sum(res_vec[r, 0] ** 2 for r in range(dim)))

        print("GENERALIZED EIGENSYSTEM AUDIT (Q_positive x = lambda * Q_arch^{(-)} x):")
        print(f"  Smallest generalized eigenvalue lambda_min = {mp.nstr(lam_min, 20)}")
        print(f"  Operator dominance margin mu_min = lambda_min - 1 = {mp.nstr(margin_min, 20)}")
        print(f"  Eigenpair residual ||Q_pos x_min - lam_min Q_- x_min||_2 = {mp.nstr(res_norm, 10)}")
        print(f"  Next 4 generalized eigenvalues : {[mp.nstr(l, 10) for l in evals_gen[1:min(5, len(evals_gen))]]}\n")

        # ---------------------------------------------------------------------
        # Diagnostic 2: Modal Energy Projection of Dangerous State x_min
        # ---------------------------------------------------------------------
        e_low = sum(x_min[r, 0] ** 2 for r in range(min(3, dim)))
        e_high = sum(x_min[r, 0] ** 2 for r in range(3, dim)) if dim > 3 else mp.mpf('0')
        e_ratio = e_low / (e_low + e_high)

        print("DIAGNOSTIC 2: Dangerous State x_min Modal Energy Projection:")
        print(f"  Energy in 3-mode sector {e_0, e_1, e_2} : E_low  = {mp.nstr(e_low, 15)}")
        print(f"  Energy in high modes {e_3, ..., e_N}   : E_high = {mp.nstr(e_high, 15)}")
        print(f"  Fractional 3-mode energy E_low / E_total = {mp.nstr(e_ratio, 12)} ({mp.nstr(e_ratio * 100, 6)} %)")
        print(f"  Leading coordinates of x_min : {[mp.nstr(x_min[r, 0], 8) for r in range(min(6, dim))]}\n")

        # ---------------------------------------------------------------------
        # Diagnostic 3: Schur Complement Test for High-Mode Elimination
        # ---------------------------------------------------------------------
        print("DIAGNOSTIC 3: Schur Complement Partition of Q_Weil:")
        if dim > 3:
            A = Q_weil[:3, :3]
            B = Q_weil[:3, 3:]
            C = Q_weil[3:, 3:]

            C_sym = mp.mpf('0.5') * (C + C.T)
            c_evals, _ = mp.eigsy(C_sym)
            min_eig_C = min(c_evals)

            if min_eig_C > 0:
                C_inv = mp.inverse(C_sym)
                S_low = A - (B * C_inv * B.T)
                S_low_sym = mp.mpf('0.5') * (S_low + S_low.T)
                s_evals, _ = mp.eigsy(S_low_sym)
                s_evals_sorted = sorted(s_evals)
                min_eig_S = s_evals_sorted[0]
                print(f"  High-mode block C dimension: {C.rows} x {C.cols}")
                print(f"  Minimum eigenvalue of high-mode block min eig(C) = {mp.nstr(min_eig_C, 15)}")
                print(f"  Eigenvalues of 3x3 Schur complement S_low = {[mp.nstr(v, 12) for v in s_evals_sorted]}")
                print(f"  Minimum eigenvalue of S_low = {mp.nstr(min_eig_S, 15)}")
                print(f"  High-mode elimination condition satisfied (C > 0 and S_low > 0): {min_eig_C > 0 and min_eig_S > 0}")
            else:
                min_eig_S = mp.mpf('nan')
                print(f"  High-mode block C is NOT strictly positive definite (min eig = {mp.nstr(min_eig_C, 10)})")
        else:
            c_evals, _ = mp.eigsy(Q_weil)
            min_eig_C = mp.mpf('0')
            min_eig_S = min(c_evals)
            print(f"  N = {N} <= 2: Dimension {dim} <= 3 (pure low-mode system, min eig = {mp.nstr(min_eig_S, 15)})")

        t1_N = time.perf_counter()
        print(f"\nTotal elapsed time for N = {N}: {t1_N - t0_N:.3f} s\n")

        summary_records.append({
            'N': N,
            'dim': dim,
            'sigma_min_minus': sigma_min,
            'cond_minus': cond_minus,
            'lambda_min': lam_min,
            'margin_min': margin_min,
            'e_ratio': e_ratio,
            'min_eig_C': min_eig_C,
            'min_eig_S': min_eig_S,
        })

    # -------------------------------------------------------------------------
    # Comprehensive Multi-Dimension Synthesis Table
    # -------------------------------------------------------------------------
    print("=" * 80)
    print("MULTI-DIMENSION SYNTHESIS TABLE (c = 13, Sweep N in {4, 8, 12, 16, 20, 24})")
    print("=" * 80)
    header = (
        f"{'N':>3} | {'sigma_min(Q_-)':>15} | {'kappa(Q_-)':>11} | "
        f"{'lambda_min':>15} | {'lambda_min - 1':>15} | {'E_low/E_tot':>11} | {'min_eig(S)':>15}"
    )
    print(header)
    print("-" * len(header))
    for rec in summary_records:
        n_val = rec['N']
        s_min = mp.nstr(rec['sigma_min_minus'], 8)
        k_val = mp.nstr(rec['cond_minus'], 6)
        l_min = mp.nstr(rec['lambda_min'], 9)
        m_min = mp.nstr(rec['margin_min'], 9)
        e_rat = mp.nstr(rec['e_ratio'] * 100, 5) + "%"
        s_low = mp.nstr(rec['min_eig_S'], 8)
        print(f"{n_val:3d} | {s_min:>15} | {k_val:>11} | {l_min:>15} | {m_min:>15} | {e_rat:>11} | {s_low:>15}")

    print("=" * 80)
    print("CELL 63 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
