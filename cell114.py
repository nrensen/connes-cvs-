# ============================================================
# CELL 114 — TWO-STATE SPECTRAL REORDERING ANATOMY & LOCALIZATION INVARIANTS
# ============================================================
#
# Gate:       Gate 1 (Finite-N Spectral Mechanism & Asymptotic Tail Extinction)
# Route:      Route D (Regularity Bridge / Boundary Defect Extinction)
# Milestone:  M12 (Boundary Defect Decoupling & Regularity Bridge)
#
# Target Propositions:
#
#   1. Module 1 — Two-Dimensional Low-Energy Subspace Tracking:
#      Sweep N in {40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60} (step Delta N = 2).
#      Compute the lowest eigenpairs (E_0, E_1) and the spectral gap Delta E = E_1 - E_0.
#      Measure the 2D invariant subspace overlap via SVD of the 2x2 inter-dimensional
#      overlap matrix M_{ij} = <v_{N_prev}^(i), v_N^(j)>, extracting principal singular
#      values (sigma_1, sigma_2) and the internal rotation angle phi(N) across the transition.
#
#   2. Module 2 — Operator-Independent Physical Localization Invariants:
#      For each eigenstate k in {0, 1} across the transition, evaluate:
#        - Central amplitude: v_0
#        - Core mass concentration: L_M(v) = sum_{m=0}^M v_m^2 for M in {8, 16, 24}
#        - Kinetic Sobolev moment: K_2(v) = sum_{m=1}^N m^4 v_m^2
#        - Curvature moment: S_2(v) = sum_{m=1}^N m^2 v_m
#      identifying the localized solitary wave independently of eigenvalue ordering.
#
#   3. Module 3 — High-N Branch Comparison (N in {64, 80, 96, 128, 192}):
#      Concurrently evaluate both states k = 0 and k = 1 at high discrete dimensions,
#      measuring the boundary defect T_v(0), coupling scalar alpha_N, and extinction
#      products P_T(N) = |T(0)| * N^(3/2) and P_alpha(N) = |alpha_N| * sqrt(N),
#      testing whether the localized state exhibits extinction decay or a plateau.
#
#   4. Module 4 — Finite-T Spectral Structure at Fixed N = 48:
#      At N = 48, evaluate eigenpairs, core masses L_24, and kinetic moments K_2
#      across cached cutoffs T in {400, 500, 600}, testing the response of both
#      states to the finite-T Archimedean truncation.
#
# Output Standard:
#   Objective, dispassionate output reporting computed numerical values,
#   singular values, gaps, moments, and extinction metrics without qualitative labels.
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

N_FINE_SWEEP = [40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60]
N_HIGH_SWEEP = [64, 80, 96, 128, 192]
T_SWEEP_LIST = [400, 500, 600]
N_T_FIXED = 48
NUM_EIGENPAIRS = 5

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


def compute_2x2_subspace_svd(v_prev_0, v_prev_1, N_prev, v_curr_0, v_curr_1, N_curr):
    """
    Compute 2x2 inter-dimensional overlap matrix M and its singular values:
      M_{ij} = <v_prev^(i), v_curr^(j)> for i, j in {0, 1}.
    v_prev is embedded into R^(N_curr + 1) by zero padding.
    """
    # Dot products over common indices 0 .. N_prev
    m00 = sum(v_prev_0[m] * v_curr_0[m] for m in range(N_prev + 1))
    m01 = sum(v_prev_0[m] * v_curr_1[m] for m in range(N_prev + 1))
    m10 = sum(v_prev_1[m] * v_curr_0[m] for m in range(N_prev + 1))
    m11 = sum(v_prev_1[m] * v_curr_1[m] for m in range(N_prev + 1))

    # M^T M components
    a = m00 ** 2 + m10 ** 2
    b = m00 * m01 + m10 * m11
    c = m01 ** 2 + m11 ** 2

    # Eigenvalues of M^T M
    tr = a + c
    det_val = a * c - b ** 2
    disc = mp.sqrt(max(mp.mpf("0"), (a - c) ** 2 + mp.mpf("4") * b ** 2))

    lam1 = (tr + disc) / mp.mpf("2")
    lam2 = (tr - disc) / mp.mpf("2")

    sig1 = mp.sqrt(max(mp.mpf("0"), lam1))
    sig2 = mp.sqrt(max(mp.mpf("0"), lam2))

    # Rotation angle within the subspace: phi = arctan(|m01| / |m00|)
    phi_rad = mp.atan2(abs(m01), abs(m00))
    phi_deg = phi_rad * mp.mpf("180") / mp.pi

    return (m00, m01, m10, m11), (sig1, sig2), phi_deg


def compute_invariants(v_vec, N_sub):
    """
    Compute physical localization invariants on eigenvector v:
      v_0: central amplitude
      L_8, L_16, L_24: core masses
      K_2: kinetic Sobolev moment sum m^4 v_m^2
      S_2: curvature moment sum m^2 v_m
    """
    v0 = v_vec[0]
    dim = N_sub + 1

    L8 = sum(v_vec[m] ** 2 for m in range(min(9, dim)))
    L16 = sum(v_vec[m] ** 2 for m in range(min(17, dim)))
    L24 = sum(v_vec[m] ** 2 for m in range(min(25, dim)))

    K2 = sum((mp.mpf(m) ** 4) * (v_vec[m] ** 2) for m in range(1, dim))
    S2 = sum((mp.mpf(m) ** 2) * v_vec[m] for m in range(1, dim))

    return v0, L8, L16, L24, K2, S2


# ============================================================
# MAIN EXECUTION
# ============================================================

def main():
    t_start = time.time()

    print("=" * 80)
    print("CELL 114 — TWO-STATE SPECTRAL REORDERING ANATOMY & LOCALIZATION INVARIANTS")
    print(f"  Configuration: c = {C_PARAM}, Primary T = {T_PRIMARY}, N_max = {N_MAX}")
    print(f"  mpmath dps = {mp.mp.dps}, matrix dps = {GROUND_DPS}")
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
    # MODULE 1: TWO-DIMENSIONAL LOW-ENERGY SUBSPACE TRACKING (Delta N = 2)
    # ============================================================
    print("--- MODULE 1: TWO-DIMENSIONAL SUBSPACE TRACKING (N in [40, 60], Delta N = 2) ---")
    print("  Subspace SVD: M_{ij} = <v_prev^(i), v_curr^(j)>, singular values sigma_1, sigma_2, rotation angle phi")
    print(f"{'N':>4} | {'E_0':>18} | {'E_1':>18} | {'Delta E = E_1 - E_0':>20} | {'sigma_1':>9} | {'sigma_2':>9} | {'phi (deg)':>9} | {'|M_00|':>8} | {'|M_01|':>8}")
    print("-" * 122)

    fine_data = {}
    v_prev_0 = None
    v_prev_1 = None
    N_prev = None

    for N_sub in N_FINE_SWEEP:
        H_sub = extract_canonical_H(Q_full_600, N_MAX, N_sub)
        evals_sub, evecs_sub = eigsys_sym(H_sub)
        dim_sub = N_sub + 1

        v0_curr = [evecs_sub[r, 0] for r in range(dim_sub)]
        if v0_curr[0] < 0:
            v0_curr = [-x for x in v0_curr]

        v1_curr = [evecs_sub[r, 1] for r in range(dim_sub)]
        if v1_curr[0] < 0:
            v1_curr = [-x for x in v1_curr]

        E0 = evals_sub[0]
        E1 = evals_sub[1]
        delta_E = E1 - E0

        if v_prev_0 is None:
            sig1_str = "1.000000"
            sig2_str = "1.000000"
            phi_str = "0.00"
            m00_str = "1.0000"
            m01_str = "0.0000"
        else:
            (m00, m01, m10, m11), (sig1, sig2), phi_deg = compute_2x2_subspace_svd(
                v_prev_0, v_prev_1, N_prev, v0_curr, v1_curr, N_sub
            )
            sig1_str = mp.nstr(sig1, 7)
            sig2_str = mp.nstr(sig2, 7)
            phi_str = f"{float(phi_deg):.2f}"
            m00_str = mp.nstr(abs(m00), 5)
            m01_str = mp.nstr(abs(m01), 5)

        fine_data[N_sub] = {
            "evals": evals_sub[:2],
            "v0": v0_curr,
            "v1": v1_curr,
            "delta_E": delta_E,
        }

        v_prev_0 = v0_curr
        v_prev_1 = v1_curr
        N_prev = N_sub

        print(
            f"{N_sub:>4} | "
            f"{mp.nstr(E0, 10):>18} | "
            f"{mp.nstr(E1, 10):>18} | "
            f"{mp.nstr(delta_E, 10):>20} | "
            f"{sig1_str:>9} | "
            f"{sig2_str:>9} | "
            f"{phi_str:>9} | "
            f"{m00_str:>8} | "
            f"{m01_str:>8}"
        )

    print("-" * 122)
    print()

    # ============================================================
    # MODULE 2: PHYSICAL LOCALIZATION INVARIANTS ACROSS THE TRANSITION
    # ============================================================
    print("--- MODULE 2: PHYSICAL LOCALIZATION INVARIANTS ACROSS THE TRANSITION ---")
    print("  State k = 0 vs State k = 1: Central Amplitude v_0, Core Mass L_24, Kinetic Moment K_2")
    print(f"{'N':>4} | {'v_0(k=0)':>10} | {'v_0(k=1)':>10} | {'L_24(k=0)':>10} | {'L_24(k=1)':>10} | {'K_2(k=0)':>14} | {'K_2(k=1)':>14} | {'Localized':>10}")
    print("-" * 96)

    for N_sub in N_FINE_SWEEP:
        d = fine_data[N_sub]
        v0_0, _, _, L24_0, K2_0, _ = compute_invariants(d["v0"], N_sub)
        v0_1, _, _, L24_1, K2_1, _ = compute_invariants(d["v1"], N_sub)

        # Localized state identified by core mass L_24
        loc_label = "k = 0" if L24_0 > L24_1 else "k = 1"

        print(
            f"{N_sub:>4} | "
            f"{mp.nstr(v0_0, 6):>10} | "
            f"{mp.nstr(v0_1, 6):>10} | "
            f"{mp.nstr(L24_0, 6):>10} | "
            f"{mp.nstr(L24_1, 6):>10} | "
            f"{mp.nstr(K2_0, 6):>14} | "
            f"{mp.nstr(K2_1, 6):>14} | "
            f"{loc_label:>10}"
        )

    print("-" * 96)
    print()

    # ============================================================
    # MODULE 3: HIGH-N BRANCH COMPARISON & EXTINCTION METRICS
    # ============================================================
    print("--- MODULE 3: HIGH-N BRANCH COMPARISON & EXTINCTION METRICS (T = 600) ---")
    print("  Concurrent evaluation of State k = 0 (Edge Candidate) and State k = 1 (Localized Candidate)")
    print(f"{'N':>4} | {'k':>2} | {'Energy E':>18} | {'v_0':>8} | {'L_24':>9} | {'K_2':>10} | {'|T_v(0)|':>16} | {'|alpha_N|':>16} | {'P_alpha = |a|*sqrt(N)':>22}")
    print("-" * 122)

    high_records = {}

    for N_sub in N_HIGH_SWEEP:
        H_sub = extract_canonical_H(Q_full_600, N_MAX, N_sub)
        evals_sub, evecs_sub = eigsys_sym(H_sub)
        dim_sub = N_sub + 1
        N_mp = mp.mpf(N_sub)

        high_records[N_sub] = []

        for k in [0, 1]:
            vec = [evecs_sub[r, k] for r in range(dim_sub)]
            if vec[0] < 0:
                vec = [-x for x in vec]

            E_val = evals_sub[k]
            v0_val, _, _, L24_val, K2_val, _ = compute_invariants(vec, N_sub)

            T_zero = vec[0] + mp.sqrt(2) * sum(vec[m] for m in range(1, dim_sub))
            alpha_N = sum(a_vec[m] * vec[m] for m in range(1, dim_sub))

            abs_T = abs(T_zero)
            abs_alpha = abs(alpha_N)
            P_alpha = abs_alpha * mp.sqrt(N_mp)

            high_records[N_sub].append({
                "k": k,
                "E": E_val,
                "v0": v0_val,
                "L24": L24_val,
                "K2": K2_val,
                "abs_T": abs_T,
                "abs_alpha": abs_alpha,
                "P_alpha": P_alpha,
            })

            print(
                f"{N_sub:>4} | "
                f"{k:>2} | "
                f"{mp.nstr(E_val, 10):>18} | "
                f"{mp.nstr(v0_val, 5):>8} | "
                f"{mp.nstr(L24_val, 5):>9} | "
                f"{mp.nstr(K2_val, 5):>10} | "
                f"{mp.nstr(abs_T, 8):>16} | "
                f"{mp.nstr(abs_alpha, 8):>16} | "
                f"{mp.nstr(P_alpha, 8):>22}"
            )

    print("-" * 122)
    print()

    # Step 3b: Extinction Ratio Table between k = 1 and k = 0
    print("--- MODULE 3b: EXTINCTION COMPARISON RATIOS (k = 1 vs k = 0) ---")
    print(f"{'N':>4} | {'|T(k=1)| / |T(k=0)|':>22} | {'|alpha(k=1)| / |alpha(k=0)|':>26} | {'P_alpha(k=1) / P_alpha(k=0)':>28}")
    print("-" * 88)

    for N_sub in N_HIGH_SWEEP:
        r0 = high_records[N_sub][0]
        r1 = high_records[N_sub][1]

        ratio_T = r1["abs_T"] / r0["abs_T"] if r0["abs_T"] > mp.mpf("1e-65") else mp.mpf("nan")
        ratio_alpha = r1["abs_alpha"] / r0["abs_alpha"] if r0["abs_alpha"] > mp.mpf("1e-65") else mp.mpf("nan")
        ratio_P = r1["P_alpha"] / r0["P_alpha"] if r0["P_alpha"] > mp.mpf("1e-65") else mp.mpf("nan")

        print(
            f"{N_sub:>4} | "
            f"{mp.nstr(ratio_T, 8):>22} | "
            f"{mp.nstr(ratio_alpha, 8):>26} | "
            f"{mp.nstr(ratio_P, 8):>28}"
        )

    print("-" * 88)
    print()

    # ============================================================
    # MODULE 4: FINITE-T SPECTRAL STRUCTURE AT FIXED N = 48
    # ============================================================
    print(f"--- MODULE 4: FINITE-T SPECTRAL STRUCTURE AT FIXED N = {N_T_FIXED} ---")
    print("  Evaluation across cached cutoffs T in {400, 500, 600}:")
    print(f"{'T':>4} | {'k':>2} | {'Energy E':>18} | {'v_0':>8} | {'L_24':>9} | {'K_2':>10} | {'|T_v(0)|':>16} | {'|alpha_48|':>16} | {'Solve Time':>10}")
    print("-" * 107)

    t_mod4 = time.time()

    for T_val in T_SWEEP_LIST:
        t_sub = time.time()
        Q_T, _ = get_galerkin_matrix(
            c=C_PARAM,
            N=N_T_FIXED,
            T=T_val,
            dps=GROUND_DPS,
            verbose=False,
        )
        H_T = extract_canonical_H(Q_T, N_T_FIXED, N_T_FIXED)
        evals_T, evecs_T = eigsys_sym(H_T)
        dim_T = N_T_FIXED + 1

        a_T = [mp.mpf("0")] * dim_T
        for k in range(1, dim_T):
            a_T[k] = mp.sqrt(2) * (mp.mpf(k) ** 2) * H_T[0, k]

        dt_str = f"{time.time() - t_sub:.2f} s"

        for k in [0, 1]:
            vec_T = [evecs_T[r, k] for r in range(dim_T)]
            if vec_T[0] < 0:
                vec_T = [-x for x in vec_T]

            E_val = evals_T[k]
            v0_val, _, _, L24_val, K2_val, _ = compute_invariants(vec_T, N_T_FIXED)

            T_zero = vec_T[0] + mp.sqrt(2) * sum(vec_T[m] for m in range(1, dim_T))
            alpha_T = sum(a_T[m] * vec_T[m] for m in range(1, dim_T))

            print(
                f"{T_val:>4} | "
                f"{k:>2} | "
                f"{mp.nstr(E_val, 10):>18} | "
                f"{mp.nstr(v0_val, 5):>8} | "
                f"{mp.nstr(L24_val, 5):>9} | "
                f"{mp.nstr(K2_val, 5):>10} | "
                f"{mp.nstr(abs(T_zero), 8):>16} | "
                f"{mp.nstr(abs(alpha_T), 8):>16} | "
                f"{dt_str:>10}"
            )
            # Only print time on first row for given T
            dt_str = ' ' * 10

    print("-" * 107)
    print(f"Total Module 4 runtime: {time.time() - t_mod4:.2f} s")
    print()

    # Step 5: Summary Table of Computed Metrics
    print("=" * 80)
    print("CELL 114 NUMERICAL SUMMARY OF COMPUTED METRICS")
    print("=" * 80)
    print(f"Fine transition span:               N in [{N_FINE_SWEEP[0]}, {N_FINE_SWEEP[-1]}] (Delta N = 2)")
    print(f"High-N evaluation span:             N in {N_HIGH_SWEEP}")
    print(f"Finite-T evaluation points:         T in {T_SWEEP_LIST} at N = {N_T_FIXED}")
    print(f"Total script runtime:               {time.time() - t_start:.2f} s")
    print()

    # Unambiguous Completion Sentinel
    print("=" * 80)
    print("CELL 114 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
