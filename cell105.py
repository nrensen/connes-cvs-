# ============================================================
# CELL 105 — GROUND-STATE QUANTITATIVE LOCALIZATION & DECAY PROFILE AUDIT
# ============================================================
#
# Gate:       Gate 1 (Finite-N Spectral Mechanism & Asymptotic Tail Extinction)
#
# Target Propositions:
#
#   1. Route A (Bernstein-Paley-Wiener Continuum Embedding):
#      The ground-state Fourier coefficients v_m decay rapidly due to
#      infinite-order boundary contact and complex strip analyticity of
#      the solitary wave profile T_infty(t).
#
#   2. Route D (Sobolev Tail-Mass Enclosure - Theorem 3 in cell105.md):
#      For any Sobolev index s >= 1 and cutoff M >= 1:
#
#          ||v^(Q)(M)||^2 <= M^(-2s) * ||v_N||_{H^s}^2,
#
#      providing unconditional, deterministic polynomial tail bounds.
#
#   3. Combes-Thomas Locality Audit:
#      The Galerkin coupling matrix H_{j, k} exhibits off-diagonal decay
#      in distance d = |j - k|, enabling spatial resolvent localization.
#
#   4. Model Selection for Asymptotic Decay:
#      Determine whether the empirical decay of |v_m| is best described by:
#      (a) Pure exponential:        ln |v_m| ~ c_0 - sigma * m
#      (b) Modulated exponential:   ln |v_m| ~ c_0 - sigma * m - gamma * ln m
#      (c) Stretched exponential:   ln |v_m| ~ c_0 - sigma * m^beta
#      (d) Pure power law:          ln |v_m| ~ c_0 - p * ln m
#
# Falsification Criteria:
#
#   - If ||v^(Q)(M)||^2 > M^(-2s) * ||v||_{H^s}^2 for any M, s, Theorem 3 is refuted.
#   - If |v_m| fails to decay faster than m^(-4), super-polynomial decay is refuted.
#   - If the Sobolev norms ||v||_{H^s} diverge rapidly at small s, Route D is refuted.
#
# Design:
#
#   Evaluates at N=192, c=13, T=600 with mpmath dps = 70.
#   Extracts complete eigenvector v_N of canonical even matrix H.
#   Computes exact tail masses across dense grid M in {8, 12, 16, 20, 24, 32, 40, 48, 56, 64, 80, 96}.
#   Audits kinetic moments K_s = sum m^(2s) v_m^2 for s in {1, 2, 3, 4, 6}.
#   Measures off-diagonal decay profile of H_{j, k} as a function of separation d = |j - k|.
#
# ============================================================

import time
import mpmath as mp

from cell import get_galerkin_matrix

# ============================================================
# PARAMETERS & PRECISION
# ============================================================

mp.mp.dps = 70

C_PARAM = 13
T_PARAM = 600
N_PARAM = 192
GROUND_DPS = 50

DENSE_M_GRID = [8, 12, 16, 20, 24, 32, 40, 48, 56, 64, 80, 96]
SOBOLEV_INDICES = [1, 2, 3, 4, 6]

# ============================================================
# MATRIX UTILITIES
# ============================================================

def full_to_canonical_matrix(Q, N):
    """
    Convert full (2N+1) x (2N+1) Galerkin matrix Q into canonical
    even-basis (N+1) x (N+1) matrix H.
    """
    H = mp.matrix(N + 1, N + 1)
    centre = N

    H[0, 0] = Q[centre, centre]

    for k in range(1, N + 1):
        H[0, k] = mp.sqrt(2) * Q[centre, centre + k]
        H[k, 0] = H[0, k]

    for j in range(1, N + 1):
        for k in range(j, N + 1):
            value = (
                Q[centre + j, centre + k]
                + Q[centre + j, centre - k]
            )
            H[j, k] = value
            H[k, j] = value

    return H


def eigsys_sym(A):
    """
    Compute sorted eigenvalues and orthonormal eigenvectors of symmetric matrix A.
    Returns (evals, evecs) where evals is a list and evecs is an mp.matrix
    whose columns are the eigenvectors.
    """
    vals, V = mp.eigsy(A)
    dim = A.rows
    idx = sorted(range(dim), key=lambda i: vals[i])

    evals = [vals[i] for i in idx]
    evecs = mp.matrix(dim, dim)

    for col_out, col_in in enumerate(idx):
        for row in range(dim):
            evecs[row, col_out] = V[row, col_in]

    return evals, evecs


# ============================================================
# STATISTICAL REGRESSION UTILITIES (HIGH PRECISION)
# ============================================================

def fit_linear(x_vals, y_vals):
    """
    Fit y = a0 + a1 * x using standard linear least squares.
    Returns (a0, a1, r_squared, max_abs_res).
    """
    n = len(x_vals)
    sx = sum(x_vals)
    sy = sum(y_vals)
    sxx = sum(x * x for x in x_vals)
    syy = sum(y * y for y in y_vals)
    sxy = sum(x * y for x, y in zip(x_vals, y_vals))

    denom = n * sxx - sx * sx
    if abs(denom) < mp.mpf("1e-40"):
        return mp.mpf(0), mp.mpf(0), mp.mpf(0), mp.mpf(0)

    a1 = (n * sxy - sx * sy) / denom
    a0 = (sy - a1 * sx) / n

    # R^2 and residuals
    y_mean = sy / n
    ss_tot = sum((y - y_mean) ** 2 for y in y_vals)
    ss_res = sum((y - (a0 + a1 * x)) ** 2 for x, y in zip(x_vals, y_vals))
    r2 = mp.mpf(1) - (ss_res / ss_tot) if ss_tot > mp.mpf("1e-40") else mp.mpf(0)
    max_res = max(abs(y - (a0 + a1 * x)) for x, y in zip(x_vals, y_vals))

    return a0, a1, r2, max_res


def fit_bivariate_linear(x1_vals, x2_vals, y_vals):
    """
    Fit y = a0 + a1 * x1 + a2 * x2 using 2D linear least squares.
    Returns (a0, a1, a2, r_squared).
    """
    n = len(y_vals)
    # Normal equations: (X^T X) a = X^T y
    X = mp.matrix(n, 3)
    Y = mp.matrix(n, 1)
    for i in range(n):
        X[i, 0] = mp.mpf(1)
        X[i, 1] = x1_vals[i]
        X[i, 2] = x2_vals[i]
        Y[i, 0] = y_vals[i]

    XtX = X.T * X
    XtY = X.T * Y
    # Solve 3x3 system via LU
    a = mp.lu_solve(XtX, XtY)
    a0, a1, a2 = a[0, 0], a[1, 0], a[2, 0]

    y_mean = sum(y_vals) / n
    ss_tot = sum((y - y_mean) ** 2 for y in y_vals)
    ss_res = sum((y_vals[i] - (a0 + a1 * x1_vals[i] + a2 * x2_vals[i])) ** 2 for i in range(n))
    r2 = mp.mpf(1) - (ss_res / ss_tot) if ss_tot > mp.mpf("1e-40") else mp.mpf(0)

    return a0, a1, a2, r2


# ============================================================
# MAIN EXECUTION
# ============================================================

def main():
    t_start = time.time()

    print("=" * 80)
    print("CELL 105 — GROUND-STATE QUANTITATIVE LOCALIZATION & DECAY PROFILE AUDIT")
    print("  Functional Form Selection, Sobolev Enclosures, & Combes-Thomas Locality")
    print(f"  Configuration: N = {N_PARAM}, c = {C_PARAM}, T = {T_PARAM}")
    print(f"  mpmath dps: {mp.mp.dps}")
    print("=" * 80)
    print()

    # 1. Retrieve Galerkin matrix
    print("--- STEP 1: RETRIEVING GALERKIN MATRIX ---")
    t0 = time.time()
    Q_full, _ = get_galerkin_matrix(
        c=C_PARAM,
        N=N_PARAM,
        T=T_PARAM,
        dps=GROUND_DPS,
        verbose=False,
    )
    H = full_to_canonical_matrix(Q_full, N_PARAM)
    print(f"Canonical even matrix H: shape ({H.rows}, {H.cols}), time: {time.time() - t0:.2f} s")
    print()

    # 2. Ground-State Eigensolve on Full H
    print("--- STEP 2: GROUND-STATE EIGENSYSTEM OF FULL H ---")
    t0 = time.time()
    evals_H, evecs_H = eigsys_sym(H)
    E11 = evals_H[0]
    H_norm_op = evals_H[-1]
    v_ground = [evecs_H[i, 0] for i in range(H.rows)]
    norm_v = mp.sqrt(sum(x * x for x in v_ground))

    # Fix sign convention: ground state solitary wave is positive at midpoint / v_0 > 0
    if v_ground[0] < 0:
        v_ground = [-x for x in v_ground]

    print(f"Ground-state eigenvalue E11:      {mp.nstr(E11, 25)}")
    print(f"Spectral ceiling ||H||_op:        {mp.nstr(H_norm_op, 15)}")
    print(f"Ground-state norm:                {mp.nstr(norm_v, 20)}")
    print(f"Leading component v_0:            {mp.nstr(v_ground[0], 15)}")
    print(f"Eigensolve runtime: {time.time() - t0:.2f} s")
    print()

    # 3. Dense Tail Mass Sweep & Instantaneous Decay Rates
    print("--- STEP 3: DENSE TAIL MASS SWEEP & RATE PROFILE ---")
    tail_records = []
    for M in DENSE_M_GRID:
        tail_sq = sum(v_ground[m] ** 2 for m in range(M + 1, N_PARAM + 1))
        # Instantaneous rate sigma_M = -log(tail_sq) / (2*M)
        sigma_M = -mp.log(tail_sq) / (2 * M) if tail_sq > mp.mpf("1e-120") else mp.mpf(0)
        # Ratio of successive tails if available
        tail_records.append((M, tail_sq, sigma_M))
        print(f"  M = {M:>2} | ||v^(Q)||^2 = {mp.nstr(tail_sq, 12):>18} | sigma_M = {mp.nstr(sigma_M, 6):>8}")
    print()

    # 4. Sobolev Kinetic Moments and Polynomial Tail Enclosures (Route D)
    print("--- STEP 4: SOBOLEV KINETIC MOMENTS & THEOREM 3 ENCLOSURES ---")
    sobolev_norms = {}
    for s in SOBOLEV_INDICES:
        # ||v||_{H^s}^2 = sum_{m=0}^N (1 + m^2)^s * v_m^2
        norm_Hs_sq = sum(((mp.mpf(1) + mp.mpf(m) ** 2) ** s) * (v_ground[m] ** 2) for m in range(N_PARAM + 1))
        sobolev_norms[s] = norm_Hs_sq
        print(f"  Sobolev index s = {s:>1} | ||v||_{{H^{s}}}^2 = {mp.nstr(norm_Hs_sq, 12)}")

    print()
    print("  Audit of Theorem 3: ||v^(Q)(M)||^2 <= M^(-2s) * ||v||_{H^s}^2")
    all_thm3_verified = True
    for s in [1, 2, 4]:
        print(f"  [Sobolev Index s = {s} (Rate M^(-{2*s}))]")
        for M in [24, 32, 48, 64]:
            tail_actual = sum(v_ground[m] ** 2 for m in range(M + 1, N_PARAM + 1))
            bound_val = sobolev_norms[s] / (mp.mpf(M) ** (2 * s))
            holds = (tail_actual <= bound_val)
            slack = bound_val / tail_actual if tail_actual > mp.mpf("1e-120") else mp.inf
            if not holds:
                all_thm3_verified = False
            print(f"    M = {M:>2}: Actual = {mp.nstr(tail_actual, 8):>12} | Bound = {mp.nstr(bound_val, 8):>12} | "
                  f"Holds: {holds} (Slack: {mp.nstr(slack, 4)}x)")
    print()

    # 5. Functional Form Selection on the Eigenvector Tail |v_m|
    print("--- STEP 5: FUNCTIONAL FORM REGRESSION ON TAIL |v_m| ---")
    # Select mode window where |v_m| is above the 70 dps precision floor (~ 1e-35 to 1e-30)
    # Mode index 10 to 65
    FIT_MIN_M = 12
    FIT_MAX_M = 60
    fit_modes = list(range(FIT_MIN_M, FIT_MAX_M + 1))
    v_abs_fit = [abs(v_ground[m]) for m in fit_modes]
    ln_v_fit = [mp.log(x) for x in v_abs_fit]

    m_float = [mp.mpf(m) for m in fit_modes]
    ln_m_float = [mp.log(m) for m in fit_modes]

    # Model 1: Pure Exponential: ln |v_m| = c0 - sigma * m
    c0_exp, neg_sigma, r2_exp, max_res_exp = fit_linear(m_float, ln_v_fit)
    sigma_pure = -neg_sigma

    # Model 2: Modulated Exponential: ln |v_m| = c0 - sigma * m - gamma * ln(m)
    c0_mod, neg_sigma_mod, neg_gamma_mod, r2_mod = fit_bivariate_linear(m_float, ln_m_float, ln_v_fit)
    sigma_mod = -neg_sigma_mod
    gamma_mod = -neg_gamma_mod

    # Model 3: Pure Power Law: ln |v_m| = c0 - p * ln(m)
    c0_pow, neg_p, r2_pow, max_res_pow = fit_linear(ln_m_float, ln_v_fit)
    p_pure = -neg_p

    # Model 4: Stretched Exponential: ln |v_m| = c0 - sigma * m^beta
    # Sweep beta in [0.7, 1.3]
    best_beta = mp.mpf(1)
    best_r2_beta = mp.mpf(0)
    best_sigma_beta = mp.mpf(0)
    for b_step in range(70, 131, 5):
        b_val = mp.mpf(b_step) / mp.mpf(100)
        m_beta = [m ** b_val for m in m_float]
        _, neg_sig_b, r2_b, _ = fit_linear(m_beta, ln_v_fit)
        if r2_b > best_r2_beta:
            best_r2_beta = r2_b
            best_beta = b_val
            best_sigma_beta = -neg_sig_b

    print(f"  Fitting Window: modes m in [{FIT_MIN_M}, {FIT_MAX_M}] ({len(fit_modes)} points)")
    print(f"  -------------------------------------------------------------")
    print(f"  Model 1 (Pure Exponential):        R^2 = {mp.nstr(r2_exp, 6)}")
    print(f"    |v_m| ~ C * exp(-sigma * m):     sigma = {mp.nstr(sigma_pure, 6)}")
    print(f"  -------------------------------------------------------------")
    print(f"  Model 2 (Modulated Exponential):   R^2 = {mp.nstr(r2_mod, 6)}")
    print(f"    |v_m| ~ C * m^(-gamma) * e^(-sigma*m): sigma = {mp.nstr(sigma_mod, 6)}, gamma = {mp.nstr(gamma_mod, 6)}")
    print(f"  -------------------------------------------------------------")
    print(f"  Model 3 (Stretched Exponential):   R^2 = {mp.nstr(best_r2_beta, 6)}")
    print(f"    |v_m| ~ C * exp(-sigma * m^beta): beta = {mp.nstr(best_beta, 4)}, sigma = {mp.nstr(best_sigma_beta, 6)}")
    print(f"  -------------------------------------------------------------")
    print(f"  Model 4 (Pure Power Law):          R^2 = {mp.nstr(r2_pow, 6)}")
    print(f"    |v_m| ~ C * m^(-p):              p = {mp.nstr(p_pure, 6)}")
    print()

    # 6. Combes-Thomas Off-Diagonal Decay Profile in H
    print("--- STEP 6: COMBES-THOMAS OFF-DIAGONAL DECAY IN H ---")
    distances = [1, 2, 3, 4, 8, 12, 16, 24, 32, 48, 64]
    off_diag_records = []
    for d in distances:
        entries = []
        for i in range(H.rows - d):
            entries.append(abs(H[i, i + d]))
        mean_val = sum(entries) / len(entries)
        max_val = max(entries)
        off_diag_records.append((d, mean_val, max_val))
        print(f"  Distance d = {d:>2} | Mean |H_{{i, i+d}}| = {mp.nstr(mean_val, 8):>14} | Max |H_{{i, i+d}}| = {mp.nstr(max_val, 8):>14}")
    print()

    # 7. Diagonal Dominance and Diagonal Spectral Distances
    print("--- STEP 7: DIAGONAL VALUES AND SPECTRAL DISTANCES ---")
    diag_samples = [0, 1, 2, 4, 8, 16, 24, 32, 48, 64, 96, 128, 160, 192]
    for k in diag_samples:
        H_kk = H[k, k]
        row_sum_excl = sum(abs(H[k, j]) for j in range(H.cols) if j != k)
        ratio_diag = (H_kk - E11) / row_sum_excl if row_sum_excl > mp.mpf("1e-40") else mp.mpf(0)
        print(f"  k = {k:>3} | H_{{kk}} = {mp.nstr(H_kk, 8):>10} | (H_{{kk}}-E11) = {mp.nstr(H_kk - E11, 8):>10} | "
              f"Row sum excl: {mp.nstr(row_sum_excl, 8):>10} | Diag ratio: {mp.nstr(ratio_diag, 4):>6}")
    print()

    # ============================================================
    # SYNTHESIS TABLE
    # ============================================================
    print("=" * 110)
    print("CELL 105 SYNTHESIS: QUANTITATIVE LOCALIZATION SCORECARD")
    print("=" * 110)
    print(f"{'M':>4} | {'||v^(Q)||^2':>16} | {'sigma_M':>8} | {'s=1 Bound':>14} | {'s=2 Bound':>14} | {'s=4 Bound':>14} | {'Best Model':>12}")
    print("-" * 110)

    for M, tail_sq, sig in tail_records:
        b1 = sobolev_norms[1] / (mp.mpf(M) ** 2)
        b2 = sobolev_norms[2] / (mp.mpf(M) ** 4)
        b4 = sobolev_norms[4] / (mp.mpf(M) ** 8)
        model_name = "Modulated Exp" if r2_mod > r2_exp else "Pure Exp"
        print(
            f"{M:>4} | "
            f"{mp.nstr(tail_sq, 6):>16} | "
            f"{mp.nstr(sig, 4):>8} | "
            f"{mp.nstr(b1, 4):>14} | "
            f"{mp.nstr(b2, 4):>14} | "
            f"{mp.nstr(b4, 4):>14} | "
            f"{model_name:>12}"
        )
    print("=" * 110)
    print()

    # Final Conclusions
    print("--- PRECISE THEORETICAL CONCLUSIONS ---")
    print(f"1. Theorem 3 (Sobolev Tail Enclosures):          {'CERTIFIED' if all_thm3_verified else 'VIOLATED'}")
    best_model_name = "Modulated Exp" if r2_mod >= max(r2_exp, best_r2_beta, r2_pow) else ("Pure Exp" if r2_exp >= max(best_r2_beta, r2_pow) else ("Stretched Exp" if best_r2_beta >= r2_pow else "Pure Power Law"))
    best_r2 = max(r2_mod, r2_exp, best_r2_beta, r2_pow)
    print(f"2. Best-fit Functional Form:                     {best_model_name} (R^2 = {mp.nstr(best_r2, 6)})")
    print(f"   Model parameters (Modulated Exp):             sigma = {mp.nstr(sigma_mod, 5)}, gamma = {mp.nstr(gamma_mod, 5)}")
    print(f"   Model parameters (Pure Exp):                  sigma = {mp.nstr(sigma_pure, 5)}")
    print(f"3. Pure Power Law Fit (R^2):                     {mp.nstr(r2_pow, 6)} (Exponent p = {mp.nstr(p_pure, 5)})")
    print(f"Total script runtime: {time.time() - t_start:.2f} s")
    print()

    print("=" * 80)
    print("CELL 105 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
