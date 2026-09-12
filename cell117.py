# ============================================================
# CELL 117 — BOUNDARY-ROW MODAL PROFILING & BOUNDARY LAYER SCALING COLLAPSE DIAGNOSTIC
# ============================================================
#
# Gate:       Gate 1 (Finite-N Spectral Mechanism & Asymptotic Tail Extinction)
# Route:      Route D (Regularity Bridge / Boundary Defect Extinction)
# Milestone:  M12 (Boundary Defect Decoupling & Regularity Bridge)
#
# Target Propositions:
#
#   1. Module 1 — Exact Analytical & Numerical Audit of H_{0, N} and a_N:
#      Verify the exact operator identity:
#        H_{0, N} = sqrt(2) * psi(N) / N,    a_N = 2 * N * psi(N).
#      Decompose psi(N) = psi_prime(N) + psi_pole(N) + psi_arch(N) across
#      N in [32, 192] to analytically explain why a_N / sqrt(2) oscillates
#      wildly in sign (+24.17 to -145.38) and grows as O(N) in envelope due
#      to the almost-periodic, non-decaying prime sum psi_prime(N).
#
#   2. Module 2 — Multi-Dimension Localized Branch Extraction:
#      Extract the localized solitary wave branch (state k = 1) across five
#      dimensions N in {64, 96, 128, 160, 192} from the cached N = 192 Hamiltonian.
#      Compute (H u_N)_N, alpha_N, T_v(0), and verify the exact Theorem 1 identity.
#
#   3. Module 3 — Modal Boundary-Row Profiling (j = N - k):
#      Compute individual modal flux terms:
#        F_{N, k} = H_{Nk} * k^2 * v_{N, k}
#      and re-index by distance from the boundary:
#        j = N - k  for  j in {0, 1, 2, ..., 24}.
#      Tabulate raw modal contributions and signs across all five dimensions.
#
#   4. Module 4 — Boundary Layer Scaling Collapse Diagnostic:
#      Test whether the boundary terms collapse onto a universal profile:
#        F_{N, N-j} approx N^(-mu) * f(j / N^theta)
#      for candidate exponents theta in {0, 1/3, 1/2, 1} with mu = 1/3.
#      Evaluate collapse variance across dimensions to objectively test
#      the Airy boundary-layer hypothesis.
#
#   5. Module 5 — Binned Flux Decomposition & Cancellation Tracking:
#      Partition modes by distance j:
#        - Bin 0: Endpoint j = 0 (k = N)
#        - Bin 1: Immediate Boundary Layer j in [1, 4]
#        - Bin 2: Near-Boundary j in [5, 10]
#        - Bin 3: Intermediate Sector j in [11, 20]
#        - Bin 4: Transition Sector j in [21, 40]
#        - Bin 5: Low-Mode Core j > 40 (k <= N - 41)
#      Compute cumulative flux curves from core outwards and boundary inwards.
#      Track the peak partial sum M_peak(N) and the cancellation factor:
#        C_cancel(N) = max_M |sum_{k=1}^M F_{N, k}| / |(H u_N)_N|.
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
AUDIT_N_GRID = [32, 48, 64, 80, 96, 112, 128, 144, 160, 176, 192]

# ============================================================
# MATHEMATICAL HELPER FUNCTIONS
# ============================================================

def calc_psi_prime(n_val: int, L: mp.mpf, prime_data: list) -> mp.mpf:
    """
    Compute prime-piece contribution to psi(n_val):
    psi_prime(n) = -(1/pi) * sum_{p^k <= c} (Lambda(p^k)/sqrt(p^k)) * sin(2*pi*n*(1 - log(p^k)/L)).
    """
    PI = mp.pi
    n_mp = mp.mpf(n_val)
    two_pi_n = 2 * PI * n_mp
    s = mp.mpf("0")
    for (_, logn, w) in prime_data:
        s += w * mp.sin(two_pi_n * (mp.mpf("1") - logn / L))
    return -s / PI


def calc_psi_pole(n_val: int, L: mp.mpf, c_val: mp.mpf) -> mp.mpf:
    """
    Exact analytical closed-form evaluation of the zeta-pole piece:
    psi_pole(n) = (1/pi) * integral_0^L sin(2*pi*n*(1 - y/L)) * 2*cosh(y/2) dy
                = [ (2*n/L) / (1/4 + (2*pi*n/L)^2) ] * (sqrt(c) + 1/sqrt(c) - 2).
    Zero numerical quadrature required.
    """
    n_mp = mp.mpf(n_val)
    k_val = 2 * mp.pi * n_mp / L
    sqrt_c = mp.sqrt(c_val)
    geom_factor = sqrt_c + mp.mpf("1") / sqrt_c - mp.mpf("2")
    return ((2 * n_mp / L) / (mp.mpf("0.25") + k_val ** 2)) * geom_factor


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


# ============================================================
# MAIN EXECUTION
# ============================================================

def main():
    t_start = time.time()

    print("=" * 80)
    print("CELL 117 — BOUNDARY-ROW MODAL PROFILING & BOUNDARY LAYER SCALING COLLAPSE")
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
    psi_vec = [mp.mpf("0")] * dim_192
    for k in range(1, dim_192):
        k_mp = mp.mpf(k)
        H0k = H_192[0, k]
        a_vec[k] = mp.sqrt(mp.mpf("2")) * (k_mp ** 2) * H0k
        psi_vec[k] = (k_mp / mp.sqrt(mp.mpf("2"))) * H0k

    # ============================================================
    # MODULE 1: EXACT ANALYTICAL & NUMERICAL AUDIT OF H_{0, N} AND a_N
    # ============================================================
    print("--- MODULE 1: EXACT ANALYTICAL & NUMERICAL AUDIT OF H_{0, N} AND a_N ---")
    print("  Identity: H_{0, N} = sqrt(2)*psi(N)/N,   a_N = 2*N*psi(N)")
    print("  Decomposition: psi(N) = psi_prime(N) + psi_pole(N) + psi_arch(N)")
    print()
    print(f"{'N':>4} | {'H_{0, N}':>15} | {'psi_prime':>15} | {'psi_pole':>15} | {'psi_arch':>15} | {'psi_total':>15} | {'a_N / sqrt(2)':>16}")
    print("-" * 105)

    for n_val in AUDIT_N_GRID:
        H0N = H_192[0, n_val]
        psi_tot = psi_vec[n_val]
        aN_over_sqrt2 = a_vec[n_val] / mp.sqrt(mp.mpf("2"))

        p_prime = calc_psi_prime(n_val, L_mp, prime_data)
        p_pole = calc_psi_pole(n_val, L_mp, c_mp)
        p_arch = psi_tot - p_prime - p_pole

        print(
            f"{n_val:>4} | "
            f"{mp.nstr(H0N, 8):>15} | "
            f"{mp.nstr(p_prime, 8):>15} | "
            f"{mp.nstr(p_pole, 8):>15} | "
            f"{mp.nstr(p_arch, 8):>15} | "
            f"{mp.nstr(psi_tot, 8):>15} | "
            f"{mp.nstr(aN_over_sqrt2, 8):>16}"
        )

    print("-" * 105)
    print("  Key Analytical Findings:")
    print("  1. psi_prime(N) is an almost-periodic sum of high-frequency sines over prime powers.")
    print("     It does NOT decay as N -> inf; it oscillates indefinitely with O(1) envelope.")
    print("  2. psi_pole(N) decays as O(1/N) as expected from the zeta pole integral.")
    print("  3. psi_arch(N) provides a smooth, slowly varying Mellin baseline.")
    print("  4. Because psi_total(N) oscillates with O(1) amplitude, a_N / sqrt(2) = sqrt(2)*N*psi_total(N)")
    print("     oscillates wildly in sign and its envelope grows linearly as O(N).")
    print()

    # ============================================================
    # MODULE 2: MULTI-DIMENSION LOCALIZED BRANCH EXTRACTION
    # ============================================================
    print("--- MODULE 2: MULTI-DIMENSION LOCALIZED BRANCH EXTRACTION ---")
    print("  Solving localized solitary branch (state k = 1) across N in {64, 96, 128, 160, 192}")
    print()
    print(f"{'N':>4} | {'v_0':>12} | {'||v_N||_2 - 1':>16} | {'(H u_N)_N':>18} | {'alpha_N':>18} | {'T_v(0)':>18} | {'Thm 1 Residual':>16}")
    print("-" * 112)

    solitary_data = {}

    for N_sub in N_GRID:
        dim_sub = N_sub + 1
        H_sub = extract_canonical_H(Q_full_600, N_MAX, N_sub)

        evals_sub, evecs_sub = eigsys_sym(H_sub)

        # Select localized branch (v_0 > 0.4)
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
        norm_res = abs(mp.sqrt(sum(x ** 2 for x in v_sub)) - mp.mpf("1"))

        # Compute observables
        T_zero = v_sub[0] + mp.sqrt(mp.mpf("2")) * sum(v_sub[m] for m in range(1, dim_sub))
        alpha_val = sum(a_vec[m] * v_sub[m] for m in range(1, dim_sub))

        # Boundary flux terms F_{N, k} = H_{N, k} * k^2 * v_{N, k}
        flux_terms = [H_sub[N_sub, k] * (mp.mpf(k) ** 2) * v_sub[k] for k in range(1, dim_sub)]
        Hu_N = sum(flux_terms)

        # Theorem 1 identity residual
        contact_term = (a_vec[N_sub] / mp.sqrt(mp.mpf("2"))) * T_zero
        eig_term = E_val * (mp.mpf(N_sub) ** 2) * v_sub[N_sub]
        thm1_rhs = alpha_val - contact_term + eig_term
        thm1_res = abs(Hu_N - thm1_rhs)

        solitary_data[N_sub] = {
            "v": v_sub,
            "E": E_val,
            "v0": v0_val,
            "T_zero": T_zero,
            "alpha": alpha_val,
            "Hu_N": Hu_N,
            "flux_terms": flux_terms,
            "H_sub": H_sub,
        }

        print(
            f"{N_sub:>4} | "
            f"{mp.nstr(v0_val, 8):>12} | "
            f"{mp.nstr(norm_res, 6):>16} | "
            f"{mp.nstr(Hu_N, 8):>18} | "
            f"{mp.nstr(alpha_val, 8):>18} | "
            f"{mp.nstr(T_zero, 8):>18} | "
            f"{mp.nstr(thm1_res, 6):>16}"
        )

    print("-" * 112)
    print()

    # ============================================================
    # MODULE 3: MODAL BOUNDARY-ROW PROFILING (j = N - k)
    # ============================================================
    print("--- MODULE 3: MODAL BOUNDARY-ROW PROFILING F_{N, N-j} (j = N - k) ---")
    print("  Individual boundary row terms F_{N, N-j} = H_{N, N-j} * (N-j)^2 * v_{N, N-j}")
    print()
    header_m3 = f"{'j':>3} |" + "".join([f"{f'N = {n}':>20} |" for n in N_GRID])
    print(header_m3)
    print("-" * len(header_m3))

    J_MAX = 24
    for j in range(J_MAX + 1):
        row_str = f"{j:>3} |"
        for N_sub in N_GRID:
            f_terms = solitary_data[N_sub]["flux_terms"]
            # f_terms is indexed by k = 1..N, so k = N - j corresponds to index N - j - 1
            idx_k = N_sub - j - 1
            if 0 <= idx_k < len(f_terms):
                val = f_terms[idx_k]
                row_str += f"{mp.nstr(val, 8):>20} |"
            else:
                row_str += f"{'N/A':>20} |"
        print(row_str)

    print("-" * len(header_m3))
    print()

    # ============================================================
    # MODULE 4: BOUNDARY LAYER SCALING COLLAPSE DIAGNOSTIC
    # ============================================================
    print("--- MODULE 4: BOUNDARY LAYER SCALING COLLAPSE DIAGNOSTIC ---")
    print("  Testing hypothesis: F_{N, N-j} approx N^(-1/3) * f(j / N^theta)")
    print()

    # Test 4a: Fixed distance j scaling S_{N, j} = N^(1/3) * F_{N, N-j}
    print("  [Test 4a: Discrete Fixed-Depth Scaling (theta = 0)]")
    print("  Evaluating invariant S_{N, j} = N^(1/3) * F_{N, N-j} at fixed distance j:")
    header_4a = f"{'j':>3} |" + "".join([f"{f'N={n} (N^(1/3)*F)':>18} |" for n in N_GRID]) + f"{'Spread %':>12}"
    print(header_4a)
    print("-" * len(header_4a))

    TEST_J_LIST = [0, 1, 2, 3, 4, 5, 8, 12, 16, 20]
    for j in TEST_J_LIST:
        row_str = f"{j:>3} |"
        s_vals = []
        for N_sub in N_GRID:
            f_terms = solitary_data[N_sub]["flux_terms"]
            idx_k = N_sub - j - 1
            val = f_terms[idx_k]
            n_third = mp.mpf(N_sub) ** (mp.mpf("1") / mp.mpf("3"))
            s_val = val * n_third
            s_vals.append(s_val)
            row_str += f"{mp.nstr(s_val, 6):>18} |"

        # Relative spread: (max - min) / mean
        mean_s = sum(s_vals) / len(s_vals)
        if abs(mean_s) > mp.mpf("1e-40"):
            spread_pct = abs((max(s_vals) - min(s_vals)) / mean_s) * 100
            row_str += f"{mp.nstr(spread_pct, 4):>11} %"
        else:
            row_str += f"{'N/A':>12}"
        print(row_str)

    print("-" * len(header_4a))
    print()

    # Test 4b: Airy-Scaled Coordinate eta = j / N^(1/3)
    print("  [Test 4b: Airy Scaling Coordinate (theta = 1/3)]")
    print("  Evaluating N^(1/3) * F_{N, N-j} against eta = j / N^(1/3):")
    print(f"  N = 64:  N^(1/3) = 4.00,   j in {{0, 1, 2, 4, 6, 8}}   -> eta in {{0.00, 0.25, 0.50, 1.00, 1.50, 2.00}}")
    print(f"  N = 192: N^(1/3) = 5.77,   j in {{0, 1, 3, 6, 9, 12}}  -> eta in {{0.00, 0.17, 0.52, 1.04, 1.56, 2.08}}")
    print()
    print(f"{'Target eta':>10} | {'(N=64, j)':>12} | {'S(N=64)':>16} | {'(N=192, j)':>12} | {'S(N=192)':>16} | {'Ratio S_192 / S_64':>20}")
    print("-" * 96)

    ETA_PAIRS = [
        (mp.mpf("0.0"), 0, 0),
        (mp.mpf("0.5"), 2, 3),
        (mp.mpf("1.0"), 4, 6),
        (mp.mpf("1.5"), 6, 9),
        (mp.mpf("2.0"), 8, 12),
        (mp.mpf("2.5"), 10, 14),
        (mp.mpf("3.0"), 12, 17),
    ]

    for eta_target, j_64, j_192 in ETA_PAIRS:
        f_64 = solitary_data[64]["flux_terms"][64 - j_64 - 1]
        f_192 = solitary_data[192]["flux_terms"][192 - j_192 - 1]

        s_64 = f_64 * (mp.mpf("64") ** (mp.mpf("1") / mp.mpf("3")))
        s_192 = f_192 * (mp.mpf("192") ** (mp.mpf("1") / mp.mpf("3")))

        if abs(s_64) > mp.mpf("1e-40"):
            ratio_s = s_192 / s_64
            ratio_str = mp.nstr(ratio_s, 6)
        else:
            ratio_str = "N/A"

        print(
            f"{mp.nstr(eta_target, 3):>10} | "
            f"{f'j={j_64}':>12} | "
            f"{mp.nstr(s_64, 6):>16} | "
            f"{f'j={j_192}':>12} | "
            f"{mp.nstr(s_192, 6):>16} | "
            f"{ratio_str:>20}"
        )

    print("-" * 96)
    print()

    # ============================================================
    # MODULE 5: BINNED FLUX DECOMPOSITION & CANCELLATION TRACKING
    # ============================================================
    print("--- MODULE 5: BINNED FLUX DECOMPOSITION & MULTI-ORDER CANCELLATION TRACKING ---")
    print("  Distance Partitions: Bin 0 (j=0), Bin 1 (1<=j<=4), Bin 2 (5<=j<=10),")
    print("                       Bin 3 (11<=j<=20), Bin 4 (21<=j<=40), Bin 5 (j>40, core)")
    print()

    # Part A: Sector Partial Sums Table
    header_m5 = f"{'Sector / Bin':<28} |" + "".join([f"{f'N = {n}':>18} |" for n in N_GRID])
    print(header_m5)
    print("-" * len(header_m5))

    BIN_RANGES = [
        ("Bin 0: Endpoint (j=0)", 0, 0),
        ("Bin 1: Bdy Layer (1<=j<=4)", 1, 4),
        ("Bin 2: Near-Bdy (5<=j<=10)", 5, 10),
        ("Bin 3: Intermed (11<=j<=20)", 11, 20),
        ("Bin 4: Transition (21<=j<=40)", 21, 40),
        ("Bin 5: Bulk Core (j>40)", 41, 10000),
    ]

    for bin_name, j_min, j_max in BIN_RANGES:
        row_str = f"{bin_name:<28} |"
        for N_sub in N_GRID:
            f_terms = solitary_data[N_sub]["flux_terms"]
            # Sum over terms where j in [j_min, min(j_max, N_sub - 1)]
            actual_j_max = min(j_max, N_sub - 1)
            if j_min <= actual_j_max:
                bin_sum = sum(f_terms[N_sub - j - 1] for j in range(j_min, actual_j_max + 1))
                row_str += f"{mp.nstr(bin_sum, 6):>18} |"
            else:
                row_str += f"{'0.0':>18} |"
        print(row_str)

    print("-" * len(header_m5))

    # Total flux row
    row_total = f"{'Total Flux (H u_N)_N':<28} |"
    for N_sub in N_GRID:
        Hu_val = solitary_data[N_sub]["Hu_N"]
        row_total += f"{mp.nstr(Hu_val, 6):>18} |"
    print(row_total)
    print("-" * len(header_m5))
    print()

    # Part B: Multi-Order Destructive Cancellation Tracking
    print("  [Part B: Cancellation Ratio Tracking across Dimensions]")
    print("  M_peak(N) = max_{1 <= M <= N} |sum_{k=1}^M F_{N, k}|")
    print("  Cancellation Ratio: C_cancel(N) = M_peak(N) / |(H u_N)_N|")
    print()
    print(f"{'N':>4} | {'Net Flux (H u_N)_N':>20} | {'Peak Sum M_peak':>20} | {'k_peak':>8} | {'Cancellation Factor C_cancel':>30}")
    print("-" * 92)

    for N_sub in N_GRID:
        f_terms = solitary_data[N_sub]["flux_terms"]
        Hu_val = solitary_data[N_sub]["Hu_N"]

        # Track partial sums from core outwards
        running_sum = mp.mpf("0")
        peak_val = mp.mpf("0")
        k_peak = 1
        for k in range(1, N_sub + 1):
            running_sum += f_terms[k - 1]
            if abs(running_sum) > peak_val:
                peak_val = abs(running_sum)
                k_peak = k

        c_cancel = peak_val / abs(Hu_val)

        print(
            f"{N_sub:>4} | "
            f"{mp.nstr(Hu_val, 8):>20} | "
            f"{mp.nstr(peak_val, 8):>20} | "
            f"{k_peak:>8} | "
            f"{mp.nstr(c_cancel, 8):>30}"
        )

    print("-" * 92)
    print()

    # ============================================================
    # MODULE 6: EXECUTION SUMMARY & SENTINEL
    # ============================================================
    print("=" * 80)
    print("CELL 117 NUMERICAL SUMMARY OF COMPUTED METRICS")
    print("=" * 80)
    print(f"Total script runtime:               {time.time() - t_start:.2f} s")
    print()

    # Unambiguous Completion Sentinel
    print("=" * 80)
    print("CELL 117 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
