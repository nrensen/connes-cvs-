# ============================================================
# CELL 115 — LOCALIZED BRANCH BOUNDARY DEFECT ASYMPTOTIC POWER-LAW REGRESSION
# ============================================================
#
# Gate:       Gate 1 (Finite-N Spectral Mechanism & Asymptotic Tail Extinction)
# Route:      Route D (Regularity Bridge / Boundary Defect Extinction)
# Milestone:  M12 (Boundary Defect Decoupling & Regularity Bridge)
#
# Target Propositions:
#
#   1. Module 1 — Multi-Dimension Spectral Extraction on Dense Grid:
#      Sweep N in {48, 64, 80, 96, 112, 128, 144, 160, 192} at c = 13, T = 600, dps = 70.
#      Extract canonical even Hamiltonian H from cached Q(192, 600).
#      Compute the lowest two eigenpairs (E_0, v^(0)) and (E_1, v^(1)).
#      For each N and each branch k in {0, 1}, evaluate:
#        - Eigenvalue E_k
#        - Central amplitude v_{k, 0}
#        - High-mode tail mass (1 - L_24)
#        - Kinetic Sobolev moment K_2
#        - Physical boundary contact defect T_{v_k}(0)
#        - Boundary coupling scalar alpha_N(k)
#        - Extinction product P_alpha(k) = |alpha_N(k)| * sqrt(N)
#        - Contact product P_T(k) = |T_{v_k}(0)| * N^(3/2)
#
#   2. Module 2 — Rescaled Boundary Observables & Branch Ratios:
#      On the localized branch (k = 1 for N >= 64, k = 0 for N = 48), evaluate:
#        |alpha_N|, N^(1/2) |alpha_N|, N |alpha_N|,
#        |T_v(0)|, N^(3/2) |T_v(0)|, N |T_v(0)|,
#      and the branch coupling ratio |alpha_N^(1)| / |alpha_N^(0)|.
#
#   3. Module 3 — Logarithmic Power-Law Regressions (N in [64, 192], n = 8):
#      Fit empirical power laws Y(N) ~ C * N^beta via ordinary least squares:
#        log10(Y) = beta * log10(N) + log10(C).
#      Compute beta, standard error SE(beta), log10(C), C, and R^2 for:
#        - Boundary coupling |alpha_N| on k = 1 (Target: beta < -0.5)
#        - Contact defect |T_v(0)| on k = 1 (Target: beta < -1.5)
#        - Extinction product P_alpha = |alpha_N| * sqrt(N) on k = 1 (Target: beta < 0)
#        - Contact product P_T = |T_v(0)| * N^(3/2) on k = 1 (Target: beta < 0)
#        - Branch coupling ratio |alpha_N^(1)| / |alpha_N^(0)|
#        - Corresponding regressions on edge branch k = 0 for comparison.
#
#   4. Module 4 — Scaling Exponents Summary & Quantitative Comparison:
#      Tabulate fitted exponents, empirical confidence intervals, and
#      distance from Gate 1 extinction criteria.
#
# Output Standard:
#   Strictly objective, dispassionate numerical reporting of computed values,
#   regression coefficients, standard errors, and residuals without editorial flourish.
#
# Configuration:
#   c = 13, Primary T = 600, N_max = 192, working dps = 70, matrix dps = 70.
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
T_PRIMARY = 600
N_MAX = 192
GROUND_DPS = 70

# Dense evaluation grid (pre-transition N=48, post-transition N in [64, 192])
N_GRID = [48, 64, 80, 96, 112, 128, 144, 160, 192]
N_POST_TRANSITION = [64, 80, 96, 112, 128, 144, 160, 192]

# ============================================================
# MATRIX & ALGEBRA UTILITIES
# ============================================================

def extract_canonical_H(Q_full, N_full, N_sub):
    """
    Extract canonical even-basis (N_sub + 1) x (N_sub + 1) matrix H
    from full (2*N_full + 1) x (2*N_full + 1) Galerkin matrix Q.
    """
    H = mp.matrix(N_sub + 1, N_sub + 1)
    centre = N_full

    H[0, 0] = Q_full[centre, centre]

    for k in range(1, N_sub + 1):
        H[0, k] = mp.sqrt(2) * Q_full[centre, centre + k]
        H[k, 0] = H[0, k]

    for j in range(1, N_sub + 1):
        for k in range(j, N_sub + 1):
            val = (
                Q_full[centre + j, centre + k]
                + Q_full[centre + j, centre - k]
            )
            H[j, k] = val
            H[k, j] = val

    return H


def eigsys_sym(A):
    """
    Compute sorted eigenvalues and orthonormal eigenvectors of symmetric matrix A.
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


def compute_invariants(v_vec, N_sub):
    """
    Compute localization observables on eigenvector v:
      v_0: central amplitude
      tail_24: high-mode mass sum_{m=25}^N v_m^2 = 1 - L_24
      K_2: kinetic Sobolev moment sum_{m=1}^N m^4 v_m^2
    """
    v0 = v_vec[0]
    dim = N_sub + 1

    L24 = sum(v_vec[m] ** 2 for m in range(min(25, dim)))
    tail_24 = max(mp.mpf("0"), mp.mpf("1") - L24)
    K2 = sum((mp.mpf(m) ** 4) * (v_vec[m] ** 2) for m in range(1, dim))

    return v0, tail_24, K2


def log_linear_regression(n_list, y_list):
    """
    Fit log10(Y) = beta * log10(N) + log10(C) via ordinary least squares.
    Returns: beta, se_beta, log10_C, C_val, R2
    """
    x_vals = [mp.log10(mp.mpf(n)) for n in n_list]
    y_vals = [mp.log10(abs(y)) for y in y_list]
    n_pts = len(x_vals)

    x_bar = sum(x_vals) / mp.mpf(n_pts)
    y_bar = sum(y_vals) / mp.mpf(n_pts)

    S_xx = sum((x - x_bar) ** 2 for x in x_vals)
    S_yy = sum((y - y_bar) ** 2 for y in y_vals)
    S_xy = sum((x - x_bar) * (y - y_bar) for x, y in zip(x_vals, y_vals))

    beta = S_xy / S_xx
    log10_C = y_bar - beta * x_bar
    C_val = mp.mpf("10") ** log10_C

    if S_xx > 0 and S_yy > 0:
        R2 = (S_xy ** 2) / (S_xx * S_yy)
    else:
        R2 = mp.mpf("0")

    if n_pts > 2 and S_xx > 0:
        res_var = max(mp.mpf("0"), (S_yy - beta * S_xy) / mp.mpf(n_pts - 2))
        se_beta = mp.sqrt(res_var / S_xx)
    else:
        se_beta = mp.mpf("0")

    return beta, se_beta, log10_C, C_val, R2


# ============================================================
# MAIN EXECUTION
# ============================================================

def main():
    t_start = time.time()

    print("=" * 80)
    print("CELL 115 — LOCALIZED BRANCH BOUNDARY DEFECT ASYMPTOTIC POWER-LAW REGRESSION")
    print(f"  Configuration: c = {C_PARAM}, Primary T = {T_PRIMARY}, N_max = {N_MAX}")
    print(f"  mpmath dps = {mp.mp.dps}, matrix dps = {GROUND_DPS}")
    print(f"  Evaluation Grid: N in {N_GRID}")
    print("=" * 80)
    print()

    # Step 1: Retrieve cached Galerkin matrix at N = 192, T = 600
    print("--- STEP 1: RETRIEVING CACHED GALERKIN MATRIX (T = 600, N = 192) ---")
    t0 = time.time()
    Q_full_600, _ = get_galerkin_matrix(
        c=C_PARAM,
        N=N_MAX,
        T=T_PRIMARY,
        dps=GROUND_DPS,
        verbose=False,
    )
    H_192 = extract_canonical_H(Q_full_600, N_MAX, N_MAX)
    print(f"Canonical even matrix H: shape ({H_192.rows}, {H_192.cols}), retrieval time: {time.time() - t0:.2f} s")
    print()

    # Modal multiplier vector a on N_MAX: a_k = sqrt(2) * k^2 * H_{0, k}
    dim_192 = N_MAX + 1
    a_vec = [mp.mpf("0")] * dim_192
    for k in range(1, dim_192):
        a_vec[k] = mp.sqrt(2) * (mp.mpf(k) ** 2) * H_192[0, k]

    # ============================================================
    # MODULE 1: MULTI-DIMENSION SPECTRAL EXTRACTION ON GRID
    # ============================================================
    print("--- MODULE 1: MULTI-DIMENSION SPECTRAL EXTRACTION (k = 0 and k = 1) ---")
    print(f"{'N':>4} | {'k':>2} | {'Energy E':>18} | {'v_0':>8} | {'1 - L_24':>12} | {'K_2':>10} | {'|T_v(0)|':>16} | {'|alpha_N|':>16} | {'P_alpha':>16}")
    print("-" * 118)

    records = {}

    for N_sub in N_GRID:
        H_sub = extract_canonical_H(Q_full_600, N_MAX, N_sub)
        evals_sub, evecs_sub = eigsys_sym(H_sub)
        dim_sub = N_sub + 1

        records[N_sub] = {}

        for k in [0, 1]:
            vec = [evecs_sub[r, k] for r in range(dim_sub)]
            if vec[0] < 0:
                vec = [-x for x in vec]

            E_val = evals_sub[k]
            v0_val, tail_24, K2_val = compute_invariants(vec, N_sub)

            T_zero = vec[0] + mp.sqrt(2) * sum(vec[m] for m in range(1, dim_sub))
            alpha_N = sum(a_vec[m] * vec[m] for m in range(1, dim_sub))

            abs_T = abs(T_zero)
            abs_alpha = abs(alpha_N)
            sqrt_N = mp.sqrt(mp.mpf(N_sub))
            P_alpha = abs_alpha * sqrt_N
            P_T = abs_T * (mp.mpf(N_sub) ** mp.mpf("1.5"))

            records[N_sub][k] = {
                "E": E_val,
                "v0": v0_val,
                "tail_24": tail_24,
                "K2": K2_val,
                "abs_T": abs_T,
                "abs_alpha": abs_alpha,
                "P_alpha": P_alpha,
                "P_T": P_T,
            }

            print(
                f"{N_sub:>4} | "
                f"{k:>2} | "
                f"{mp.nstr(E_val, 10):>18} | "
                f"{mp.nstr(v0_val, 5):>8} | "
                f"{mp.nstr(tail_24, 6):>12} | "
                f"{mp.nstr(K2_val, 5):>10} | "
                f"{mp.nstr(abs_T, 8):>16} | "
                f"{mp.nstr(abs_alpha, 8):>16} | "
                f"{mp.nstr(P_alpha, 8):>16}"
            )

    print("-" * 118)
    print()

    # ============================================================
    # MODULE 2: RESCALED BOUNDARY OBSERVABLES ON LOCALIZED CANDIDATE
    # ============================================================
    print("--- MODULE 2: RESCALED BOUNDARY OBSERVABLES ON LOCALIZED BRANCH (k = 1 for N >= 64) ---")
    print(f"{'N':>4} | {'|alpha_N|':>16} | {'N^(1/2)|alpha|':>16} | {'N|alpha|':>16} | {'|T_v(0)|':>16} | {'N^(3/2)|T|':>16} | {'N|T|':>16} | {'|a1|/|a0|':>10}")
    print("-" * 122)

    for N_sub in N_POST_TRANSITION:
        r1 = records[N_sub][1]
        r0 = records[N_sub][0]

        a_val = r1["abs_alpha"]
        T_val = r1["abs_T"]

        N_mp = mp.mpf(N_sub)
        p_half = a_val * mp.sqrt(N_mp)
        p_one = a_val * N_mp

        t_three_half = T_val * (N_mp ** mp.mpf("1.5"))
        t_one = T_val * N_mp

        ratio_a = a_val / r0["abs_alpha"] if r0["abs_alpha"] > mp.mpf("1e-65") else mp.mpf("nan")

        print(
            f"{N_sub:>4} | "
            f"{mp.nstr(a_val, 8):>16} | "
            f"{mp.nstr(p_half, 8):>16} | "
            f"{mp.nstr(p_one, 8):>16} | "
            f"{mp.nstr(T_val, 8):>16} | "
            f"{mp.nstr(t_three_half, 8):>16} | "
            f"{mp.nstr(t_one, 8):>16} | "
            f"{mp.nstr(ratio_a, 6):>10}"
        )

    print("-" * 122)
    print()

    # ============================================================
    # MODULE 3: LOGARITHMIC POWER-LAW REGRESSIONS (N in [64, 192])
    # ============================================================
    print("--- MODULE 3: LOGARITHMIC POWER-LAW REGRESSIONS (N in [64, 192], n = 8) ---")
    print("  Model: log10(Y) = beta * log10(N) + log10(C)  <=>  Y(N) ~ C * N^beta")
    print()

    # Data vectors for k = 1
    alpha_list_1 = [records[n][1]["abs_alpha"] for n in N_POST_TRANSITION]
    T_list_1 = [records[n][1]["abs_T"] for n in N_POST_TRANSITION]
    P_alpha_list_1 = [records[n][1]["P_alpha"] for n in N_POST_TRANSITION]
    P_T_list_1 = [records[n][1]["P_T"] for n in N_POST_TRANSITION]
    ratio_alpha_list = [records[n][1]["abs_alpha"] / records[n][0]["abs_alpha"] for n in N_POST_TRANSITION]

    # Data vectors for k = 0
    alpha_list_0 = [records[n][0]["abs_alpha"] for n in N_POST_TRANSITION]
    T_list_0 = [records[n][0]["abs_T"] for n in N_POST_TRANSITION]
    P_alpha_list_0 = [records[n][0]["P_alpha"] for n in N_POST_TRANSITION]
    P_T_list_0 = [records[n][0]["P_T"] for n in N_POST_TRANSITION]

    # Regressions on k = 1 (Localized Candidate)
    beta_a1, se_a1, logC_a1, C_a1, R2_a1 = log_linear_regression(N_POST_TRANSITION, alpha_list_1)
    beta_T1, se_T1, logC_T1, C_T1, R2_T1 = log_linear_regression(N_POST_TRANSITION, T_list_1)
    beta_Pa1, se_Pa1, logC_Pa1, C_Pa1, R2_Pa1 = log_linear_regression(N_POST_TRANSITION, P_alpha_list_1)
    beta_PT1, se_PT1, logC_PT1, C_PT1, R2_PT1 = log_linear_regression(N_POST_TRANSITION, P_T_list_1)
    beta_rat, se_rat, logC_rat, C_rat, R2_rat = log_linear_regression(N_POST_TRANSITION, ratio_alpha_list)

    # Regressions on k = 0 (Edge Candidate)
    beta_a0, se_a0, logC_a0, C_a0, R2_a0 = log_linear_regression(N_POST_TRANSITION, alpha_list_0)
    beta_T0, se_T0, logC_T0, C_T0, R2_T0 = log_linear_regression(N_POST_TRANSITION, T_list_0)
    beta_Pa0, se_Pa0, logC_Pa0, C_Pa0, R2_Pa0 = log_linear_regression(N_POST_TRANSITION, P_alpha_list_0)
    beta_PT0, se_PT0, logC_PT0, C_PT0, R2_PT0 = log_linear_regression(N_POST_TRANSITION, P_T_list_0)

    print(f"{'Quantity Y(N)':<28} | {'Branch':<11} | {'Exponent beta':>14} | {'SE(beta)':>11} | {'R^2':>10} | {'Pre-factor C':>16}")
    print("-" * 105)
    print(f"{'|alpha_N|':<28} | {'k = 1 (loc)':<11} | {mp.nstr(beta_a1, 7):>14} | {mp.nstr(se_a1, 5):>11} | {mp.nstr(R2_a1, 6):>10} | {mp.nstr(C_a1, 8):>16}")
    print(f"{'|alpha_N|':<28} | {'k = 0 (edge)':<11} | {mp.nstr(beta_a0, 7):>14} | {mp.nstr(se_a0, 5):>11} | {mp.nstr(R2_a0, 6):>10} | {mp.nstr(C_a0, 8):>16}")
    print(f"{'|T_v(0)|':<28} | {'k = 1 (loc)':<11} | {mp.nstr(beta_T1, 7):>14} | {mp.nstr(se_T1, 5):>11} | {mp.nstr(R2_T1, 6):>10} | {mp.nstr(C_T1, 8):>16}")
    print(f"{'|T_v(0)|':<28} | {'k = 0 (edge)':<11} | {mp.nstr(beta_T0, 7):>14} | {mp.nstr(se_T0, 5):>11} | {mp.nstr(R2_T0, 6):>10} | {mp.nstr(C_T0, 8):>16}")
    print(f"{'P_alpha = |alpha_N|*sqrt(N)':<28} | {'k = 1 (loc)':<11} | {mp.nstr(beta_Pa1, 7):>14} | {mp.nstr(se_Pa1, 5):>11} | {mp.nstr(R2_Pa1, 6):>10} | {mp.nstr(C_Pa1, 8):>16}")
    print(f"{'P_alpha = |alpha_N|*sqrt(N)':<28} | {'k = 0 (edge)':<11} | {mp.nstr(beta_Pa0, 7):>14} | {mp.nstr(se_Pa0, 5):>11} | {mp.nstr(R2_Pa0, 6):>10} | {mp.nstr(C_Pa0, 8):>16}")
    print(f"{'P_T = |T_v(0)|*N^(3/2)':<28} | {'k = 1 (loc)':<11} | {mp.nstr(beta_PT1, 7):>14} | {mp.nstr(se_PT1, 5):>11} | {mp.nstr(R2_PT1, 6):>10} | {mp.nstr(C_PT1, 8):>16}")
    print(f"{'P_T = |T_v(0)|*N^(3/2)':<28} | {'k = 0 (edge)':<11} | {mp.nstr(beta_PT0, 7):>14} | {mp.nstr(se_PT0, 5):>11} | {mp.nstr(R2_PT0, 6):>10} | {mp.nstr(C_PT0, 8):>16}")
    print(f"{'Ratio |alpha_1| / |alpha_0|':<28} | {'k=1 vs k=0':<11} | {mp.nstr(beta_rat, 7):>14} | {mp.nstr(se_rat, 5):>11} | {mp.nstr(R2_rat, 6):>10} | {mp.nstr(C_rat, 8):>16}")
    print("-" * 105)
    print()

    # ============================================================
    # MODULE 4: EXTINCTION CRITERIA COMPARISON
    # ============================================================
    print("--- MODULE 4: QUANTITATIVE DISTANCE TO GATE 1 EXTINCTION CRITERIA ---")
    print(f"{'Observable':<24} | {'Fitted Exponent':>16} | {'Target Extinction Exponent':>28} | {'Exponent Gap (Fitted - Target)':>32}")
    print("-" * 108)

    gap_alpha = beta_a1 - mp.mpf("-0.5")
    gap_T = beta_T1 - mp.mpf("-1.5")
    gap_Pa = beta_Pa1 - mp.mpf("0.0")
    gap_PT = beta_PT1 - mp.mpf("0.0")

    print(f"{'|alpha_N| (k=1)':<24} | {mp.nstr(beta_a1, 6):>16} | {'< -0.500000':>28} | {mp.nstr(gap_alpha, 6):>32}")
    print(f"{'|T_v(0)| (k=1)':<24} | {mp.nstr(beta_T1, 6):>16} | {'< -1.500000':>28} | {mp.nstr(gap_T, 6):>32}")
    print(f"{'P_alpha(N) (k=1)':<24} | {mp.nstr(beta_Pa1, 6):>16} | {'<  0.000000':>28} | {mp.nstr(gap_Pa, 6):>32}")
    print(f"{'P_T(N) (k=1)':<24} | {mp.nstr(beta_PT1, 6):>16} | {'<  0.000000':>28} | {mp.nstr(gap_PT, 6):>32}")
    print("-" * 108)
    print()

    # Step 5: Summary of Execution Metrics
    print("=" * 80)
    print("CELL 115 NUMERICAL SUMMARY OF COMPUTED METRICS")
    print("=" * 80)
    print(f"Dimension grid span:               N in {N_GRID}")
    print(f"Regression span (post-transition):  N in {N_POST_TRANSITION} (n = {len(N_POST_TRANSITION)})")
    print(f"Total script runtime:               {time.time() - t_start:.2f} s")
    print()

    # Unambiguous Completion Sentinel
    print("=" * 80)
    print("CELL 115 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
