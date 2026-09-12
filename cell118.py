# ============================================================
# CELL 118 — GLOBAL MODAL CANCELLATION ANATOMY & ARITHMETIC REMAINDER AUDIT
# ============================================================
#
# Gate:       Gate 1 (Finite-N Spectral Mechanism & Asymptotic Tail Extinction)
# Route:      Route D (Regularity Bridge / Boundary Defect Extinction)
# Milestone:  M12 (Boundary Defect Decoupling & Regularity Bridge)
#
# Target Propositions:
#
#   1. Module 1 — Low-Mode Anatomy & The Universal k = 3 Peak:
#      Audit individual modal flux contributions F_{N, k} = H_{Nk} k^2 v_{N, k}
#      and cumulative partial sums S_N(m) = sum_{k=1}^m F_{N, k} for low modes
#      k in [1, 10] across N in {64, 96, 128, 160, 192}. Explain the universal
#      peak M_peak ~ 10^(-3) - 10^(-2) occurring at k = 3.
#
#   2. Module 2 — Macroscopic Cancellation Curve S_N(x) on x in (0, 1):
#      Evaluate cumulative partial sums S_N(floor(x * N)) sampled at 16 fractional
#      coordinates x in [0.01, 1.00] across all 5 benchmark dimensions.
#      Test whether the normalized profile g_N(x) = S_N(floor(x*N)) / M_peak(N)
#      collapses onto a universal continuum master curve G(x). Pinpoint the
#      crossover coordinate x_* where the sum drops by 20 orders of magnitude.
#
#   3. Module 3 — Bulk-Core vs Transition vs Boundary Cancellation:
#      Quantify the fraction of total cancellation occurring in:
#        - Bulk core: x in [0.00, 0.25]
#        - Intermediate sector: x in [0.25, 0.75]
#        - High-frequency edge: x in [0.75, 1.00].
#
#   4. Module 4 — Arithmetic Symbol Remainder Correlation:
#      Test whether the boundary residual alpha_N or (H u_N)_N correlates
#      with the arithmetic prime-power cosine sum:
#        Sigma_cos(N) = sum_{p^k <= c} (Lambda(p^k)/sqrt(p^k)) * cos(2*pi*N*log(p^k)/L)
#      or with psi'(N) across N in {64, 80, 96, 112, 128, 144, 160, 192}.
#      Compute Pearson correlation coefficients to determine if alpha_N
#      is an arithmetic remainder rather than a local boundary layer artifact.
#
# Output Standard:
#   Strictly objective, dispassionate numerical reporting of computed values,
#   residuals, and correlation metrics without editorial flourish.
#
# Configuration:
#   c = 13, Primary T = 600, N_max = 192, working dps = 70, matrix dps = 70.
#
# ============================================================

import time
import mpmath as mp

from cell import get_galerkin_matrix
from connes_cvs.operator import prime_powers_up_to

# ============================================================
# PARAMETERS & PRECISION
# ============================================================

mp.mp.dps = 70

C_PARAM = 13
T_PRIMARY = 600
N_MAX = 192
GROUND_DPS = 70

N_GRID = [64, 96, 128, 160, 192]
N_GRID_ALL = [64, 80, 96, 112, 128, 144, 160, 192]

# Continuum coordinate sampling grid
X_GRID = [
    mp.mpf("0.01"), mp.mpf("0.02"), mp.mpf("0.03"), mp.mpf("0.05"),
    mp.mpf("0.10"), mp.mpf("0.15"), mp.mpf("0.20"), mp.mpf("0.30"),
    mp.mpf("0.40"), mp.mpf("0.50"), mp.mpf("0.60"), mp.mpf("0.70"),
    mp.mpf("0.80"), mp.mpf("0.90"), mp.mpf("0.95"), mp.mpf("1.00"),
]

# ============================================================
# MATHEMATICAL HELPER FUNCTIONS
# ============================================================

def calc_psi_prime_deriv(n_val: int, L: mp.mpf, prime_data: list) -> mp.mpf:
    """
    Derivative of prime piece: d/dx psi_prime(x) at x = n_val:
    psi_prime'(n) = -2 * sum_{p^k <= c} (Lambda(p^k)/sqrt(p^k)) * (1 - log(p^k)/L) * cos(2*pi*n*(1 - log(p^k)/L)).
    """
    PI = mp.pi
    n_mp = mp.mpf(n_val)
    two_pi_n = 2 * PI * n_mp
    s = mp.mpf("0")
    for (_, logn, w) in prime_data:
        c_val = mp.mpf("1") - logn / L
        s += w * mp.mpf("2") * c_val * mp.cos(two_pi_n * c_val)
    return -s


def calc_sigma_cos(n_val: int, L: mp.mpf, prime_data: list) -> mp.mpf:
    """
    Pure prime-power cosine sum:
    Sigma_cos(n) = sum_{p^k <= c} (Lambda(p^k)/sqrt(p^k)) * cos(2*pi*n*log(p^k)/L).
    """
    PI = mp.pi
    n_mp = mp.mpf(n_val)
    two_pi_n = 2 * PI * n_mp
    s = mp.mpf("0")
    for (_, logn, w) in prime_data:
        s += w * mp.cos(two_pi_n * logn / L)
    return s


def extract_canonical_H(Q_full, N_full: int, N_sub: int) -> mp.matrix:
    """
    Extract canonical even-basis (N_sub + 1) x (N_sub + 1) matrix H
    from full (2*N_full + 1) x (2*N_full + 1) Galerkin matrix Q.
    """
    H = mp.matrix(N_sub + 1, N_sub + 1)
    centre = N_full

    H[0, 0] = Q_full[centre, centre]

    for k in range(1, N_sub + 1):
        H[0, k] = mp.sqrt(mp.mpf("2")) * Q_full[centre, centre + k]
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


def eigsys_sym(A: mp.matrix):
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


def pearson_r(x_list, y_list):
    """Compute Pearson correlation coefficient r between two numeric lists."""
    n = len(x_list)
    mean_x = sum(x_list) / n
    mean_y = sum(y_list) / n
    cov = sum((x_list[i] - mean_x) * (y_list[i] - mean_y) for i in range(n))
    var_x = sum((x - mean_x) ** 2 for x in x_list)
    var_y = sum((y - mean_y) ** 2 for y in y_list)
    denom = mp.sqrt(var_x * var_y)
    if denom == 0:
        return mp.mpf("0")
    return cov / denom


# ============================================================
# MAIN EXECUTION
# ============================================================

def main():
    t_start = time.time()

    print("=" * 80)
    print("CELL 118 — GLOBAL MODAL CANCELLATION ANATOMY & ARITHMETIC REMAINDER AUDIT")
    print(f"  Configuration: c = {C_PARAM}, Primary T = {T_PRIMARY}, N_max = {N_MAX}")
    print(f"  mpmath dps = {mp.mp.dps}, matrix dps = {GROUND_DPS}")
    print("=" * 80)
    print()

    # Step 0: Load cached Hamiltonian at N = 192, T = 600
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

    c_mp = mp.mpf(C_PARAM)
    L_mp = mp.log(c_mp)
    prime_data, _ = prime_powers_up_to(int(mp.floor(c_mp)))

    # Compute a_vec = sqrt(2) * k^2 * H_{0, k} for all k in [1, N_MAX]
    dim_192 = N_MAX + 1
    a_vec = [mp.mpf("0")] * dim_192
    for k in range(1, dim_192):
        k_mp = mp.mpf(k)
        H0k = H_192[0, k]
        a_vec[k] = mp.sqrt(mp.mpf("2")) * (k_mp ** 2) * H0k

    # Solve localized solitary eigenvector across all 8 dimensions in N_GRID_ALL
    print("  Solving localized solitary branch across dimensions...")
    solitary_data = {}
    for N_sub in N_GRID_ALL:
        dim_sub = N_sub + 1
        H_sub = extract_canonical_H(Q_full_600, N_MAX, N_sub)

        evals_sub, evecs_sub = eigsys_sym(H_sub)

        loc_idx = 1
        for idx in range(min(5, dim_sub)):
            cand_v0 = abs(evecs_sub[0, idx])
            if cand_v0 > mp.mpf("0.4"):
                loc_idx = idx
                break

        v_sub = [evecs_sub[r, loc_idx] for r in range(dim_sub)]
        if v_sub[0] < 0:
            v_sub = [-x for x in v_sub]

        E_val = evals_sub[loc_idx]
        v0_val = v_sub[0]

        T_zero = v_sub[0] + mp.sqrt(mp.mpf("2")) * sum(v_sub[m] for m in range(1, dim_sub))
        alpha_val = sum(a_vec[m] * v_sub[m] for m in range(1, dim_sub))

        flux_terms = [H_sub[N_sub, k] * (mp.mpf(k) ** 2) * v_sub[k] for k in range(1, dim_sub)]
        Hu_N = sum(flux_terms)

        # Track running cumulative sums and peak
        running_sum = mp.mpf("0")
        cum_sums = []
        peak_val = mp.mpf("0")
        k_peak = 1
        for k in range(1, N_sub + 1):
            running_sum += flux_terms[k - 1]
            cum_sums.append(running_sum)
            if abs(running_sum) > peak_val:
                peak_val = abs(running_sum)
                k_peak = k

        c_cancel = peak_val / abs(Hu_N)

        solitary_data[N_sub] = {
            "v": v_sub,
            "E": E_val,
            "v0": v0_val,
            "T_zero": T_zero,
            "alpha": alpha_val,
            "Hu_N": Hu_N,
            "flux_terms": flux_terms,
            "cum_sums": cum_sums,
            "peak_val": peak_val,
            "k_peak": k_peak,
            "c_cancel": c_cancel,
        }

    print(f"  [Eigenpair extractions complete. Localized v_0 in [{mp.nstr(solitary_data[64]['v0'], 6)}, {mp.nstr(solitary_data[192]['v0'], 6)}]]")
    print()

    # ============================================================
    # MODULE 1: LOW-MODE ANATOMY & THE UNIVERSAL k = 3 PEAK
    # ============================================================
    print("--- MODULE 1: LOW-MODE ANATOMY & THE UNIVERSAL k = 3 PEAK ---")
    print("  Auditing individual terms F_{N, k} and cumulative sums S_N(m) = sum_{k=1}^m F_{N, k} for k <= 10")
    print()

    print(f"{'k':>3} |" + "".join([f"{f'F(N={n})':>16} |{f'S(N={n})':>16} |" for n in N_GRID]))
    print("-" * 168)

    for k in range(1, 11):
        row_str = f"{k:>3} |"
        for N_sub in N_GRID:
            f_val = solitary_data[N_sub]["flux_terms"][k - 1]
            s_val = solitary_data[N_sub]["cum_sums"][k - 1]
            row_str += f"{mp.nstr(f_val, 6):>16} |{mp.nstr(s_val, 6):>16} |"
        print(row_str)

    print("-" * 168)
    print()
    print("  Key Low-Mode Observations:")
    for N_sub in N_GRID:
        pk = solitary_data[N_sub]["peak_val"]
        kp = solitary_data[N_sub]["k_peak"]
        print(f"    N = {N_sub:>3}: Peak S_N(m) = {mp.nstr(pk, 8)} at k = {kp}")
    print()

    # ============================================================
    # MODULE 2: MACROSCOPIC CANCELLATION CURVE S_N(x) ON x in (0, 1)
    # ============================================================
    print("--- MODULE 2: MACROSCOPIC CANCELLATION PROFILE S_N(x) ON CONTINUUM COORDINATE x = k / N ---")
    print("  Evaluating cumulative sum S_N(floor(x*N)) and normalized master curve g_N(x) = S_N / M_peak")
    print()

    # Part A: Cumulative Partial Sums S_N(x) Table
    header_m2a = f"{'x = k/N':>8} |" + "".join([f"{f'S(N={n})':>20} |" for n in N_GRID])
    print(header_m2a)
    print("-" * len(header_m2a))

    for x_val in X_GRID:
        row_str = f"{mp.nstr(x_val, 3):>8} |"
        for N_sub in N_GRID:
            m_idx = max(1, min(N_sub, int(mp.floor(x_val * mp.mpf(N_sub)))))
            s_val = solitary_data[N_sub]["cum_sums"][m_idx - 1]
            row_str += f"{mp.nstr(s_val, 8):>20} |"
        print(row_str)

    print("-" * len(header_m2a))
    print()

    # Part B: Normalized Master Cancellation Curve g_N(x) = S_N(x) / M_peak(N)
    print("  [Part B: Normalized Master Cancellation Curve g_N(x) = S_N(floor(x*N)) / M_peak(N)]")
    header_m2b = f"{'x = k/N':>8} |" + "".join([f"{f'g(N={n})':>16} |" for n in N_GRID]) + f"{'Spread %':>12}"
    print(header_m2b)
    print("-" * len(header_m2b))

    for x_val in X_GRID:
        row_str = f"{mp.nstr(x_val, 3):>8} |"
        g_vals = []
        for N_sub in N_GRID:
            m_idx = max(1, min(N_sub, int(mp.floor(x_val * mp.mpf(N_sub)))))
            s_val = solitary_data[N_sub]["cum_sums"][m_idx - 1]
            g_val = s_val / solitary_data[N_sub]["peak_val"]
            g_vals.append(g_val)
            row_str += f"{mp.nstr(g_val, 6):>16} |"

        mean_g = sum(g_vals) / len(g_vals)
        if abs(mean_g) > mp.mpf("1e-15"):
            spread_g = abs((max(g_vals) - min(g_vals)) / mean_g) * 100
            row_str += f"{mp.nstr(spread_g, 4):>11} %"
        else:
            row_str += f"{'N/A (-> 0)':>12}"
        print(row_str)

    print("-" * len(header_m2b))
    print()

    # ============================================================
    # MODULE 3: BULK-CORE VS TRANSITION VS BOUNDARY CANCELLATION
    # ============================================================
    print("--- MODULE 3: MACROSCOPIC SECTOR CANCELLATION DECOMPOSITION ---")
    print("  Partitioning continuum domain into 3 macroscopic sectors:")
    print("    - Bulk Core:        x in (0.00, 0.25]")
    print("    - Transition Zone:  x in (0.25, 0.75]")
    print("    - Boundary Edge:    x in (0.75, 1.00]")
    print()

    print(f"{'Sector':<24} |" + "".join([f"{f'N = {n}':>18} |" for n in N_GRID]))
    print("-" * 118)

    MACRO_SECTORS = [
        ("Bulk Core (x <= 0.25)", mp.mpf("0.00"), mp.mpf("0.25")),
        ("Transition (0.25 < x <= 0.75)", mp.mpf("0.25"), mp.mpf("0.75")),
        ("Boundary Edge (x > 0.75)", mp.mpf("0.75"), mp.mpf("1.00")),
    ]

    for sec_name, x_lo, x_hi in MACRO_SECTORS:
        row_str = f"{sec_name:<24} |"
        for N_sub in N_GRID:
            f_terms = solitary_data[N_sub]["flux_terms"]
            k_lo = int(mp.floor(x_lo * mp.mpf(N_sub))) + 1
            k_hi = int(mp.floor(x_hi * mp.mpf(N_sub)))
            if k_lo <= k_hi:
                sec_sum = sum(f_terms[k - 1] for k in range(k_lo, k_hi + 1))
                row_str += f"{mp.nstr(sec_sum, 6):>18} |"
            else:
                row_str += f"{'0.0':>18} |"
        print(row_str)

    print("-" * 118)
    row_tot = f"{'Total Flux (H u_N)_N':<24} |"
    for N_sub in N_GRID:
        row_tot += f"{mp.nstr(solitary_data[N_sub]['Hu_N'], 6):>18} |"
    print(row_tot)
    print("-" * 118)
    print()

    # ============================================================
    # MODULE 4: ARITHMETIC SYMBOL REMAINDER CORRELATION
    # ============================================================
    print("--- MODULE 4: ARITHMETIC SYMBOL REMAINDER CORRELATION ---")
    print("  Testing correlation between boundary observables and arithmetic fluctuations:")
    print("    - alpha_N:              Boundary coupling")
    print("    - (H u_N)_N:            Net boundary flux")
    print("    - Sigma_cos(N):         Prime-power cosine sum sum_{p^k <= c} Lambda/sqrt(p^k) cos(2*pi*N*log(p^k)/L)")
    print("    - psi_prime_deriv(N):   Derivative of prime symbol psi'_{prime}(N)")
    print()

    alpha_list = []
    Hu_list = []
    sigma_cos_list = []
    psi_d_list = []

    print(f"{'N':>4} | {'alpha_N':>18} | {'(H u_N)_N':>18} | {'Sigma_cos(N)':>18} | {'psi_prime_deriv(N)':>22}")
    print("-" * 88)

    for N_sub in N_GRID_ALL:
        a_val = solitary_data[N_sub]["alpha"]
        Hu_val = solitary_data[N_sub]["Hu_N"]
        s_cos = calc_sigma_cos(N_sub, L_mp, prime_data)
        p_d = calc_psi_prime_deriv(N_sub, L_mp, prime_data)

        alpha_list.append(a_val)
        Hu_list.append(Hu_val)
        sigma_cos_list.append(s_cos)
        psi_d_list.append(p_d)

        print(
            f"{N_sub:>4} | "
            f"{mp.nstr(a_val, 8):>18} | "
            f"{mp.nstr(Hu_val, 8):>18} | "
            f"{mp.nstr(s_cos, 8):>18} | "
            f"{mp.nstr(p_d, 8):>22}"
        )

    print("-" * 88)
    print()

    # Compute Pearson correlation coefficients
    r_alpha_sigma = pearson_r(alpha_list, sigma_cos_list)
    r_alpha_psid = pearson_r(alpha_list, psi_d_list)
    r_Hu_sigma = pearson_r(Hu_list, sigma_cos_list)
    r_Hu_psid = pearson_r(Hu_list, psi_d_list)
    r_alpha_Hu = pearson_r(alpha_list, Hu_list)

    print("  Pearson Correlation Matrix (n = 8 data points):")
    print(f"    r(alpha_N, Sigma_cos):         {mp.nstr(r_alpha_sigma, 6)}")
    print(f"    r(alpha_N, psi'_prime):        {mp.nstr(r_alpha_psid, 6)}")
    print(f"    r((H u_N)_N, Sigma_cos):       {mp.nstr(r_Hu_sigma, 6)}")
    print(f"    r((H u_N)_N, psi'_prime):      {mp.nstr(r_Hu_psid, 6)}")
    print(f"    r(alpha_N, (H u_N)_N):         {mp.nstr(r_alpha_Hu, 6)}")
    print()

    # ============================================================
    # MODULE 5: EXECUTION SUMMARY & SENTINEL
    # ============================================================
    print("=" * 80)
    print("CELL 118 NUMERICAL SUMMARY OF COMPUTED METRICS")
    print("=" * 80)
    print(f"Total script runtime:               {time.time() - t_start:.2f} s")
    print()

    # Unambiguous Completion Sentinel
    print("=" * 80)
    print("CELL 118 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
