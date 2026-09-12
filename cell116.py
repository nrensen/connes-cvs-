# ============================================================
# CELL 116 — BOUNDARY-FLUX KERNEL ASYMPTOTICS & AIRY BOUNDARY LAYER MECHANISM
# ============================================================
#
# Gate:       Gate 1 (Finite-N Spectral Mechanism & Asymptotic Tail Extinction)
# Route:      Route D (Regularity Bridge / Boundary Defect Extinction)
# Milestone:  M12 (Boundary Defect Decoupling & Regularity Bridge)
#
# Target Propositions:
#
#   1. Module 1 — Boundary Proportionality Lock & Constant Ratio:
#      Audit the exact ratio kappa_alpha(N) = |alpha_N| / |T_v(0)| across
#      N in {64, 80, 96, 112, 128, 144, 160, 192}. Verify that alpha_N
#      and T_v(0) are locked together with kappa_alpha(N) = 283 +/- 6,
#      proving that both observables share the exact same asymptotic scaling.
#
#   2. Module 2 — Semiclassical Airy Boundary Layer Model Testing:
#      Test the universal Airy Dirichlet boundary layer scaling:
#        |T_v(0)| ~ C_T * N^(-1/3),   |alpha_N| ~ C_alpha * N^(-1/3),
#        P_alpha(N) ~ C_alpha * N^(1/6),   P_T(N) ~ C_T * N^(7/6).
#      Compute rescaled invariant sequences:
#        N^(1/3) * |alpha_N|,   N^(1/3) * |T_v(0)|,
#        N^(-1/6) * P_alpha(N),  N^(-7/6) * P_T(N),
#      and report their means, standard deviations, and maximum relative variations.
#
#   3. Module 3 — Boundary-Row Modal Flux Decomposition (N = 192):
#      On the cached canonical Hamiltonian at N = 192, solve the localized
#      eigenvector v_192 and decompose the boundary row sum:
#        (H u_N)_N = sum_{k=1}^N H_{Nk} k^2 v_{N, k}.
#      Partition into:
#        - Bulk sector: k <= M (with M in {24, 48, 96})
#        - Boundary layer sector: M < k < N
#        - Endpoint diagonal term: k = N.
#      Measure the fractional contribution of each sector to (H u_N)_N.
#
#   4. Module 4 — Equipartition Identity Residuals:
#      Verify the exact Theorem 1 identity:
#        (H u_N)_N = alpha_N - (a_N / sqrt(2)) * T_v(0) + E_11 * N^2 * v_{N, N}.
#      Compute the numerical equipartition ratio:
#        R_equip = |(H u_N)_N| / |(a_N / sqrt(2)) * T_v(0)|.
#
# Output Standard:
#   Strictly objective, dispassionate numerical reporting of computed values,
#   residuals, and ratios without editorial flourish.
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

N_GRID = [64, 80, 96, 112, 128, 144, 160, 192]

# Exact raw data recorded from Cell 115 (localized branch k = 1)
CELL115_RAW_DATA = {
    64:  {"alpha": mp.mpf("2.4218845e-22"), "T": mp.mpf("8.7211935e-25")},
    80:  {"alpha": mp.mpf("1.9705339e-22"), "T": mp.mpf("7.0170728e-25")},
    96:  {"alpha": mp.mpf("1.8162911e-22"), "T": mp.mpf("6.4268511e-25")},
    112: {"alpha": mp.mpf("1.7951172e-22"), "T": mp.mpf("6.3200914e-25")},
    128: {"alpha": mp.mpf("1.7286674e-22"), "T": mp.mpf("6.0646725e-25")},
    144: {"alpha": mp.mpf("1.7476383e-22"), "T": mp.mpf("6.1060114e-25")},
    160: {"alpha": mp.mpf("1.7116185e-22"), "T": mp.mpf("5.9438973e-25")},
    192: {"alpha": mp.mpf("1.6798055e-22"), "T": mp.mpf("5.8116828e-25")},
}

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


# ============================================================
# MAIN EXECUTION
# ============================================================

def main():
    t_start = time.time()

    print("=" * 80)
    print("CELL 116 — BOUNDARY-FLUX KERNEL ASYMPTOTICS & AIRY BOUNDARY LAYER MECHANISM")
    print(f"  Configuration: c = {C_PARAM}, Primary T = {T_PRIMARY}, N_max = {N_MAX}")
    print(f"  mpmath dps = {mp.mp.dps}, matrix dps = {GROUND_DPS}")
    print("=" * 80)
    print()

    # ============================================================
    # MODULE 1: BOUNDARY PROPORTIONALITY LOCK
    # ============================================================
    print("--- MODULE 1: BOUNDARY PROPORTIONALITY LOCK kappa_alpha(N) = |alpha_N| / |T_v(0)| ---")
    print(f"{'N':>4} | {'|alpha_N|':>16} | {'|T_v(0)|':>16} | {'kappa_alpha = |a|/|T|':>22} | {'a_N / sqrt(2)':>16} | {'kappa / (a_N/sqrt(2))':>24}")
    print("-" * 105)

    # Step 1a: Retrieve cached Galerkin matrix at N = 192, T = 600
    t0 = time.time()
    Q_full_600, _ = get_galerkin_matrix(
        c=C_PARAM,
        N=N_MAX,
        T=T_PRIMARY,
        dps=GROUND_DPS,
        verbose=False,
    )
    H_192 = extract_canonical_H(Q_full_600, N_MAX, N_MAX)
    print(f"  [Retrieved cached H_192 in {time.time() - t0:.2f} s]")
    print()

    dim_192 = N_MAX + 1
    a_vec = [mp.mpf("0")] * dim_192
    for k in range(1, dim_192):
        a_vec[k] = mp.sqrt(2) * (mp.mpf(k) ** 2) * H_192[0, k]

    kappa_list = []
    ratio_kappa_a = []

    for N_sub in N_GRID:
        d = CELL115_RAW_DATA[N_sub]
        a_val = d["alpha"]
        T_val = d["T"]

        kappa = a_val / T_val
        kappa_list.append(kappa)

        a_N_term = a_vec[N_sub] / mp.sqrt(mp.mpf("2"))
        ratio_term = kappa / a_N_term
        ratio_kappa_a.append(ratio_term)

        print(
            f"{N_sub:>4} | "
            f"{mp.nstr(a_val, 8):>16} | "
            f"{mp.nstr(T_val, 8):>16} | "
            f"{mp.nstr(kappa, 8):>22} | "
            f"{mp.nstr(a_N_term, 8):>16} | "
            f"{mp.nstr(ratio_term, 8):>24}"
        )

    print("-" * 105)

    kappa_mean = sum(kappa_list) / len(kappa_list)
    kappa_std = mp.sqrt(sum((x - kappa_mean) ** 2 for x in kappa_list) / len(kappa_list))
    kappa_rel_var = (max(kappa_list) - min(kappa_list)) / kappa_mean

    print(f"  Mean kappa_alpha:                {mp.nstr(kappa_mean, 8)}")
    print(f"  Std Dev kappa_alpha:             {mp.nstr(kappa_std, 6)}")
    print(f"  Relative variation (max - min):  {mp.nstr(kappa_rel_var * 100, 4)} %")
    print()

    # ============================================================
    # MODULE 2: AIRY BOUNDARY LAYER MODEL RESCALING (beta = -1/3)
    # ============================================================
    print("--- MODULE 2: AIRY BOUNDARY LAYER MODEL RESCALING (Exponent beta = -1/3) ---")
    print("  Evaluating invariant products: S_alpha = N^(1/3)*|alpha_N|,  S_T = N^(1/3)*|T_v(0)|")
    print(f"{'N':>4} | {'N^(1/3)*|alpha_N|':>20} | {'N^(1/3)*|T_v(0)|':>20} | {'N^(-1/6)*P_alpha':>20} | {'N^(-7/6)*P_T':>20}")
    print("-" * 92)

    s_alpha_list = []
    s_T_list = []
    s_Pa_list = []
    s_PT_list = []

    for N_sub in N_GRID:
        d = CELL115_RAW_DATA[N_sub]
        a_val = d["alpha"]
        T_val = d["T"]

        N_mp = mp.mpf(N_sub)
        n_one_third = N_mp ** (mp.mpf("1") / mp.mpf("3"))
        n_one_sixth = N_mp ** (mp.mpf("1") / mp.mpf("6"))
        n_seven_sixth = N_mp ** (mp.mpf("7") / mp.mpf("6"))

        s_a = a_val * n_one_third
        s_T = T_val * n_one_third

        P_a = a_val * mp.sqrt(N_mp)
        P_T = T_val * (N_mp ** mp.mpf("1.5"))

        res_Pa = P_a / n_one_sixth
        res_PT = P_T / n_seven_sixth

        s_alpha_list.append(s_a)
        s_T_list.append(s_T)
        s_Pa_list.append(res_Pa)
        s_PT_list.append(res_PT)

        print(
            f"{N_sub:>4} | "
            f"{mp.nstr(s_a, 8):>20} | "
            f"{mp.nstr(s_T, 8):>20} | "
            f"{mp.nstr(res_Pa, 8):>20} | "
            f"{mp.nstr(res_PT, 8):>20}"
        )

    print("-" * 92)

    mean_s_a = sum(s_alpha_list) / len(s_alpha_list)
    var_s_a = (max(s_alpha_list) - min(s_alpha_list)) / mean_s_a
    mean_s_T = sum(s_T_list) / len(s_T_list)
    var_s_T = (max(s_T_list) - min(s_T_list)) / mean_s_T

    print(f"  Invariant C_alpha = N^(1/3)*|alpha_N|:  mean = {mp.nstr(mean_s_a, 8)},  max spread = {mp.nstr(var_s_a * 100, 4)} %")
    print(f"  Invariant C_T     = N^(1/3)*|T_v(0)|:   mean = {mp.nstr(mean_s_T, 8)},  max spread = {mp.nstr(var_s_T * 100, 4)} %")
    print()

    # ============================================================
    # MODULE 3: BOUNDARY-ROW MODAL FLUX DECOMPOSITION AT N = 192
    # ============================================================
    print("--- MODULE 3: BOUNDARY-ROW MODAL FLUX DECOMPOSITION AT N = 192 ---")
    print("  Solving localized eigenvector v_192 and decomposing (H u_N)_N = sum_{k=1}^N H_{Nk} k^2 v_k")

    t_diag = time.time()
    evals_192, evecs_192 = eigsys_sym(H_192)
    print(f"  [Eigenpair extraction at N = 192 completed in {time.time() - t_diag:.2f} s]")

    # State k = 1 is the localized solitary wave
    v_loc = [evecs_192[r, 1] for r in range(dim_192)]
    if v_loc[0] < 0:
        v_loc = [-x for x in v_loc]

    E_192 = evals_192[1]
    v0_192 = v_loc[0]

    # Full boundary row sum
    N = N_MAX
    H_row_N = [H_192[N, k] for k in range(dim_192)]

    flux_terms = [H_row_N[k] * (mp.mpf(k) ** 2) * v_loc[k] for k in range(1, dim_192)]
    total_flux = sum(flux_terms)

    # Partition across sectors
    M_CUTS = [24, 48, 96, 144, 180, 191]

    print()
    print(f"{'Sector / Partition':<32} | {'Index Range':<14} | {'Partial Sum':>18} | {'% of Total Flux':>18}")
    print("-" * 90)

    for M_val in M_CUTS:
        part_bulk = sum(flux_terms[:M_val])
        pct_bulk = (part_bulk / total_flux) * 100
        print(
            f"{f'Bulk Sector k <= {M_val}':<32} | "
            f"{f'k in [1, {M_val}]':<14} | "
            f"{mp.nstr(part_bulk, 8):>18} | "
            f"{mp.nstr(pct_bulk, 6):>17} %"
        )

    # Diagonal endpoint term k = N
    diag_term = flux_terms[-1]
    pct_diag = (diag_term / total_flux) * 100
    print(
        f"{'Endpoint Diagonal k = N':<32} | "
        f"{f'k = {N}':<14} | "
        f"{mp.nstr(diag_term, 8):>18} | "
        f"{mp.nstr(pct_diag, 6):>17} %"
    )

    print(
        f"{'Total Boundary Flux (H u_N)_N':<32} | "
        f"{f'k in [1, {N}]':<14} | "
        f"{mp.nstr(total_flux, 8):>18} | "
        f"{'100.0000 %':>18}"
    )
    print("-" * 90)
    print()

    # ============================================================
    # MODULE 4: EQUIPARTITION IDENTITY RESIDUALS & BALANCE
    # ============================================================
    print("--- MODULE 4: THEOREM 1 EQUIPARTITION IDENTITY RESIDUALS ---")
    print("  Identity: (H u_N)_N - alpha_N + (a_N / sqrt(2)) * T_v(0) - E_11 * N^2 * v_N = 0")

    T_zero_192 = v_loc[0] + mp.sqrt(2) * sum(v_loc[m] for m in range(1, dim_192))
    alpha_192 = sum(a_vec[m] * v_loc[m] for m in range(1, dim_192))

    contact_term = (a_vec[N] / mp.sqrt(2)) * T_zero_192
    eig_term = E_192 * (mp.mpf(N) ** 2) * v_loc[N]

    identity_rhs = alpha_192 - contact_term + eig_term
    residual = abs(total_flux - identity_rhs)

    ratio_flux_contact = abs(total_flux) / abs(contact_term)
    ratio_flux_alpha = abs(total_flux) / abs(alpha_192)

    print(f"  Computed boundary flux (H u_N)_N:      {mp.nstr(total_flux, 10)}")
    print(f"  Boundary coupling alpha_N:              {mp.nstr(alpha_192, 10)}")
    print(f"  Contact flux (a_N / sqrt(2)) * T_v(0):  {mp.nstr(contact_term, 10)}")
    print(f"  Eigenvalue term E_11 * N^2 * v_{{N, N}}:  {mp.nstr(eig_term, 10)}")
    print(f"  Exact Identity Residual:                {mp.nstr(residual, 8)}")
    print()
    print(f"  Flux / Contact Ratio:                   {mp.nstr(ratio_flux_contact, 8)}")
    print(f"  Flux / Alpha Ratio:                     {mp.nstr(ratio_flux_alpha, 8)}")
    print("-" * 90)
    print()

    # Step 5: Summary of Execution Metrics
    print("=" * 80)
    print("CELL 116 NUMERICAL SUMMARY OF COMPUTED METRICS")
    print("=" * 80)
    print(f"Total script runtime:               {time.time() - t_start:.2f} s")
    print()

    # Unambiguous Completion Sentinel
    print("=" * 80)
    print("CELL 116 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
