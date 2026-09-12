# ============================================================
# CELL 112A — MULTI-EIGENVALUE BRANCH TRACKING & FAST 70-DPS EXTINCTION AUDIT
# ============================================================
#
# Gate:       Gate 1 (Finite-N Spectral Mechanism & Asymptotic Tail Extinction)
# Route:      Route D (Regularity Bridge / Boundary Defect Extinction)
# Milestone:  M12 (Boundary Defect Decoupling & Regularity Bridge)
#
# Target Propositions:
#
#   1. Spectrum & Branch Tracking (Solitary Wave vs Edge Mode Crossing):
#      In Cell 111, for N >= 64, the lowest eigenvalue E_0 turned negative
#      (-1.063e-51) matching the finite-T=600 Archimedean tail leakage floor
#      (-delta_T^tail), accompanied by a sudden drop of v_0 from 0.456 to 0.064.
#      Track the lowest K = 5 eigenpairs (E_k, v^(k)) to detect whether an
#      eigenvalue crossing occurred and isolate both:
#        (a) The literal ground state (Branch k = 0)
#        (b) The solitary wave branch (Branch k = k_sol with max v_0 ~ 0.5)
#
#   2. Multi-Route Evaluation of alpha_N on Both Branches:
#      Verify that alpha_N computed via:
#        - Route 1: Direct modal sum sum_{k=1}^N a_k v_{N, k}
#        - Route 2: Boundary row (H u_N)_N + (T(0) / sqrt(2)) * a_N - E * N^2 * v_N
#        - Route 3: High-sector average (1 / (N - M)) * sum_{m=M+1}^N [(H u)_m + (T(0)/sqrt(2))*a_m - E*m^2*v_m]
#      coincide to numerical precision (< 1e-60 at 70 dps) across both branches.
#
#   3. High-Precision Eigensolver Residual Certification:
#      Certify that computed states v satisfy:
#
#          ||(H_N - E * I) v||_2 < 1e-60,
#
#      verifying that all observed properties are genuine spectral features
#      rather than numerical noise or eigensolver ill-conditioning.
#
#   4. Extinction Progression P_alpha(N) Comparison:
#      Compare the extinction products P_T(N) = |T(0)| * N^(3/2) and
#      P_alpha(N) = |alpha_N| * sqrt(N) for both the ground state and the
#      solitary wave branch, testing whether the plateau at 1.92e-22 is
#      confined to the edge mode while the solitary wave continues decaying.
#
# Falsification Criteria:
#
#   - If multi-route alpha_N discrepancy exceeds 1e-60, algebraic consistency is refuted.
#   - If ||(H - E * I) v||_2 > 1e-55, the eigensolution is numerically uncertified.
#   - If no distinct solitary wave branch (v_0 ~ 0.5) exists in the low-lying spectrum
#     for N >= 64, the branch crossing hypothesis is refuted.
#
# Configuration:
#
#   c = 13, T = 600, N_max = 192
#   mpmath dps = 70 (utilizes cached 70-dps Galerkin matrix, runtime < 1 minute)
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
N_MAX = 192
GROUND_DPS = 70

DIMENSION_SWEEP = [16, 24, 32, 48, 64, 80, 96, 128, 192]
CUTOFF_GRID = [24, 32, 48, 64]
NUM_EIGENPAIRS = 5

# ============================================================
# MATRIX UTILITIES
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
    print("CELL 112A — MULTI-EIGENVALUE BRANCH TRACKING & FAST 70-DPS EXTINCTION AUDIT")
    print("  Low-Lying Spectrum Tracking, Solitary Wave vs Edge Mode Crossing & Multi-Route Forensics")
    print(f"  Configuration: c = {C_PARAM}, T = {T_PARAM}, N_max = {N_MAX}")
    print(f"  mpmath dps: {mp.mp.dps}, matrix cache dps: {GROUND_DPS}")
    print("=" * 80)
    print()

    # 1. Retrieve cached Galerkin matrix at N=192
    print("--- STEP 1: RETRIEVING CACHED GALERKIN MATRIX (N = 192, DPS = 70) ---")
    t0 = time.time()
    Q_full, _ = get_galerkin_matrix(
        c=C_PARAM,
        N=N_MAX,
        T=T_PARAM,
        dps=GROUND_DPS,
        verbose=False,
    )
    H_192 = extract_canonical_H(Q_full, N_MAX, N_MAX)
    print(f"Canonical even matrix H: shape ({H_192.rows}, {H_192.cols}), retrieval time: {time.time() - t0:.2f} s")
    print()

    # 2. Modal multiplier vector a on N_MAX: a_k = sqrt(2) * k^2 * H_{0, k}
    dim_192 = N_MAX + 1
    a_vec = [mp.mpf("0")] * dim_192
    for k in range(1, dim_192):
        a_vec[k] = mp.sqrt(2) * (mp.mpf(k) ** 2) * H_192[0, k]

    # Pre-solve eigensystems and track low-lying spectrum
    print("--- STEP 2: LOW-LYING SPECTRUM & BRANCH IDENTIFICATION (K = 5) ---")
    print(f"{'N':>4} | {'E_0':>18} | {'v_0(0)':>10} | {'E_1':>18} | {'v_0(1)':>10} | {'E_2':>18} | {'v_0(2)':>10} | {'k_sol':>5} | {'Solve Time':>10}")
    print("-" * 116)

    t0 = time.time()
    spectrum_data = {}

    for N_sub in DIMENSION_SWEEP:
        t_sub = time.time()
        H_sub = extract_canonical_H(Q_full, N_MAX, N_sub)
        evals_sub, evecs_sub = eigsys_sym(H_sub)

        dim_sub = N_sub + 1
        num_k = min(NUM_EIGENPAIRS, dim_sub)

        vectors = []
        for k in range(num_k):
            vec = [evecs_sub[r, k] for r in range(dim_sub)]
            # Fix phase: ensure v[0] >= 0
            if vec[0] < 0:
                vec = [-x for x in vec]
            vectors.append(vec)

        # Identify solitary wave branch: candidate with highest v_0 component
        k_sol = 0
        max_v0 = vectors[0][0]
        for k in range(1, num_k):
            if vectors[k][0] > max_v0:
                max_v0 = vectors[k][0]
                k_sol = k

        spectrum_data[N_sub] = {
            "evals": evals_sub[:num_k],
            "vectors": vectors,
            "H": H_sub,
            "k_sol": k_sol,
        }

        e0_str = mp.nstr(evals_sub[0], 10)
        v0_0_str = mp.nstr(vectors[0][0], 6)
        e1_str = mp.nstr(evals_sub[1], 10) if num_k > 1 else "N/A"
        v0_1_str = mp.nstr(vectors[1][0], 6) if num_k > 1 else "N/A"
        e2_str = mp.nstr(evals_sub[2], 10) if num_k > 2 else "N/A"
        v0_2_str = mp.nstr(vectors[2][0], 6) if num_k > 2 else "N/A"
        dt_str = f"{time.time() - t_sub:.2f} s"

        print(
            f"{N_sub:>4} | "
            f"{e0_str:>18} | "
            f"{v0_0_str:>10} | "
            f"{e1_str:>18} | "
            f"{v0_1_str:>10} | "
            f"{e2_str:>18} | "
            f"{v0_2_str:>10} | "
            f"{k_sol:>5} | "
            f"{dt_str:>10}"
        )

    print("-" * 116)
    print(f"Total eigensolves runtime: {time.time() - t0:.2f} s")
    print()

    # Step 3: Independent Eigensolver Residual Certification
    print("--- STEP 3: EIGENSOLVER RESIDUAL & RAYLEIGH QUOTIENT CERTIFICATION ---")
    print("  Residual: ||(H - E * I) v||_2, Rayleigh Error: |E - <v, Hv>/||v||^2|")
    print(f"{'N':>4} | {'Branch':>14} | {'Energy E':>20} | {'||(H - E I) v||_2':>22} | {'Rayleigh Error':>20} | {'<v, xi> Orthog':>20} | {'Certified':>9}")
    print("-" * 120)

    for N_sub in DIMENSION_SWEEP:
        spec = spectrum_data[N_sub]
        H_sub = spec["H"]
        dim_sub = N_sub + 1

        branches_to_check = [("Ground (k=0)", 0)]
        if spec["k_sol"] != 0:
            branches_to_check.append((f"Solitary (k={spec['k_sol']})", spec["k_sol"]))

        for b_name, k_idx in branches_to_check:
            E_val = spec["evals"][k_idx]
            v_vec = spec["vectors"][k_idx]

            # Operator residual: r = H v - E v
            Hv = [mp.mpf("0")] * dim_sub
            for r in range(dim_sub):
                Hv[r] = sum(H_sub[r, c] * v_vec[c] for c in range(dim_sub))

            res_vec = [Hv[r] - E_val * v_vec[r] for r in range(dim_sub)]
            norm_res = mp.sqrt(sum(x ** 2 for x in res_vec))

            # Rayleigh quotient error
            norm_v_sq = sum(x ** 2 for x in v_vec)
            vHv = sum(v_vec[r] * Hv[r] for r in range(dim_sub))
            rayleigh_err = abs(E_val - vHv / norm_v_sq)

            # Commutator orthogonality: <v, xi> = 0
            T_zero = v_vec[0] + mp.sqrt(2) * sum(v_vec[m] for m in range(1, dim_sub))
            alpha_N = sum(a_vec[k] * v_vec[k] for k in range(1, dim_sub))
            xi_0 = alpha_N / mp.sqrt(2)
            xi_vec = [xi_0] + [alpha_N - (T_zero / mp.sqrt(2)) * a_vec[m] for m in range(1, dim_sub)]
            orthog = abs(sum(v_vec[r] * xi_vec[r] for r in range(dim_sub)))

            certified = bool(norm_res < mp.mpf("1e-58") and rayleigh_err < mp.mpf("1e-58"))

            print(
                f"{N_sub:>4} | "
                f"{b_name:>14} | "
                f"{mp.nstr(E_val, 10):>20} | "
                f"{mp.nstr(norm_res, 10):>22} | "
                f"{mp.nstr(rayleigh_err, 10):>20} | "
                f"{mp.nstr(orthog, 10):>20} | "
                f"{str(certified):>9}"
            )

    print("-" * 120)
    print()

    # Step 4: Multi-Route Evaluation of alpha_N Across Branches
    print("--- STEP 4: MULTI-ROUTE EVALUATION OF alpha_N (AGREEMENT AUDIT) ---")
    print("  Route 1: Direct sum sum a_k v_k")
    print("  Route 2: Boundary mode (H u)_N + (T(0)/sqrt(2))*a_N - E*N^2*v_N")
    print("  Route 3: High-sector average over m in {N/2+1, ..., N}")
    print(f"{'N':>4} | {'Branch':>14} | {'alpha_N (Route 1)':>22} | {'|R1 - R2|':>18} | {'|R1 - R3|':>18} | {'Max Discrepancy':>18} | {'Agreement':>9}")
    print("-" * 115)

    for N_sub in DIMENSION_SWEEP:
        spec = spectrum_data[N_sub]
        H_sub = spec["H"]
        dim_sub = N_sub + 1

        branches_to_check = [("Ground (k=0)", 0)]
        if spec["k_sol"] != 0:
            branches_to_check.append((f"Solitary (k={spec['k_sol']})", spec["k_sol"]))

        for b_name, k_idx in branches_to_check:
            E_val = spec["evals"][k_idx]
            v_vec = spec["vectors"][k_idx]

            # Route 1: direct sum
            alpha_1 = sum(a_vec[k] * v_vec[k] for k in range(1, dim_sub))

            # Kinetic vector u
            u_vec = [mp.mpf(m) ** 2 * v_vec[m] for m in range(dim_sub)]
            Hu_vec = [mp.mpf("0")] * dim_sub
            for r in range(dim_sub):
                Hu_vec[r] = sum(H_sub[r, c] * u_vec[c] for c in range(dim_sub))

            T_zero = v_vec[0] + mp.sqrt(2) * sum(v_vec[m] for m in range(1, dim_sub))

            # Route 2: boundary mode m = N
            alpha_2 = Hu_vec[N_sub] + (T_zero / mp.sqrt(2)) * a_vec[N_sub] - E_val * (mp.mpf(N_sub) ** 2) * v_vec[N_sub]

            # Route 3: high-sector average (M = N_sub // 2)
            M_sub = N_sub // 2
            high_rows = [
                Hu_vec[m] + (T_zero / mp.sqrt(2)) * a_vec[m] - E_val * (mp.mpf(m) ** 2) * v_vec[m]
                for m in range(M_sub + 1, dim_sub)
            ]
            alpha_3 = sum(high_rows) / mp.mpf(len(high_rows))

            diff_12 = abs(alpha_1 - alpha_2)
            diff_13 = abs(alpha_1 - alpha_3)
            max_diff = max(diff_12, diff_13)
            agrees = bool(max_diff < mp.mpf("1e-58"))

            print(
                f"{N_sub:>4} | "
                f"{b_name:>14} | "
                f"{mp.nstr(alpha_1, 10):>22} | "
                f"{mp.nstr(diff_12, 8):>18} | "
                f"{mp.nstr(diff_13, 8):>18} | "
                f"{mp.nstr(max_diff, 8):>18} | "
                f"{str(agrees):>9}"
            )

    print("-" * 115)
    print()

    # Step 5: Boundary Kinetic Flux (H u_N)_N Forensics
    print("--- STEP 5: BOUNDARY KINETIC FLUX (H u_N)_N FORENSICS ---")
    print("  Decomposition: alpha_N = (H u_N)_N + (T(0)/sqrt(2))*a_N - E*N^2*v_N")
    print(f"{'N':>4} | {'Branch':>14} | {'(H u_N)_N':>22} | {'(T(0)/sqrt(2))*a_N':>22} | {'alpha_N':>22} | {'Flux/Contact':>14}")
    print("-" * 108)

    for N_sub in DIMENSION_SWEEP:
        spec = spectrum_data[N_sub]
        H_sub = spec["H"]
        dim_sub = N_sub + 1

        branches_to_check = [("Ground (k=0)", 0)]
        if spec["k_sol"] != 0:
            branches_to_check.append((f"Solitary (k={spec['k_sol']})", spec["k_sol"]))

        for b_name, k_idx in branches_to_check:
            E_val = spec["evals"][k_idx]
            v_vec = spec["vectors"][k_idx]

            u_vec = [mp.mpf(m) ** 2 * v_vec[m] for m in range(dim_sub)]
            Hu_N = sum(H_sub[N_sub, c] * u_vec[c] for c in range(dim_sub))

            T_zero = v_vec[0] + mp.sqrt(2) * sum(v_vec[m] for m in range(1, dim_sub))
            contact = (T_zero / mp.sqrt(2)) * a_vec[N_sub]
            alpha_N = sum(a_vec[k] * v_vec[k] for k in range(1, dim_sub))

            ratio = abs(Hu_N) / abs(contact) if abs(contact) > mp.mpf("1e-65") else mp.mpf("nan")

            print(
                f"{N_sub:>4} | "
                f"{b_name:>14} | "
                f"{mp.nstr(Hu_N, 10):>22} | "
                f"{mp.nstr(contact, 10):>22} | "
                f"{mp.nstr(alpha_N, 10):>22} | "
                f"{mp.nstr(ratio, 5):>14}"
            )

    print("-" * 108)
    print()

    # Step 6: Extinction Products & Branch Comparison
    print("--- STEP 6: EXTINCTION PRODUCTS & BRANCH COMPARISON ---")
    print(f"{'N':>4} | {'Branch':>14} | {'v_0':>10} | {'|T_v(0)|':>18} | {'|alpha_N|':>18} | {'kappa_alpha':>12} | {'P_T(N)':>18} | {'P_alpha(N)':>18}")
    print("-" * 122)

    for N_sub in DIMENSION_SWEEP:
        spec = spectrum_data[N_sub]
        dim_sub = N_sub + 1

        branches_to_check = [("Ground (k=0)", 0)]
        if spec["k_sol"] != 0:
            branches_to_check.append((f"Solitary (k={spec['k_sol']})", spec["k_sol"]))

        for b_name, k_idx in branches_to_check:
            v_vec = spec["vectors"][k_idx]

            T_zero = v_vec[0] + mp.sqrt(2) * sum(v_vec[m] for m in range(1, dim_sub))
            alpha_N = sum(a_vec[k] * v_vec[k] for k in range(1, dim_sub))

            abs_T = abs(T_zero)
            abs_alpha = abs(alpha_N)

            kappa = abs_alpha / abs_T if abs_T > mp.mpf("1e-65") else mp.mpf("nan")

            N_mp = mp.mpf(N_sub)
            P_T = abs_T * (N_mp ** mp.mpf("1.5"))
            P_alpha = abs_alpha * mp.sqrt(N_mp)

            print(
                f"{N_sub:>4} | "
                f"{b_name:>14} | "
                f"{mp.nstr(v_vec[0], 6):>10} | "
                f"{mp.nstr(abs_T, 8):>18} | "
                f"{mp.nstr(abs_alpha, 8):>18} | "
                f"{mp.nstr(kappa, 5):>12} | "
                f"{mp.nstr(P_T, 8):>18} | "
                f"{mp.nstr(P_alpha, 8):>18}"
            )

    print("-" * 122)
    print()

    # Step 7: Exact High-Sector Boundary Defect Norm ||xi_N^(Q)||_2 Grid
    print("--- STEP 7: EXACT HIGH-SECTOR DEFECT NORM ||xi_N^(Q)||_2 GRID ---")
    print(f"{'N':>4} | {'Branch':>14} | {'M':>4} | {'||xi_N^(Q)||_2 (Exact)':>24} | {'Triangle Bound':>22} | {'Ratio':>10} | {'Holds':>7}")
    print("-" * 97)

    for N_sub in DIMENSION_SWEEP:
        if N_sub <= 24:
            continue
        spec = spectrum_data[N_sub]
        dim_sub = N_sub + 1

        branches_to_check = [("Ground (k=0)", 0)]
        if spec["k_sol"] != 0:
            branches_to_check.append((f"Solitary (k={spec['k_sol']})", spec["k_sol"]))

        for b_name, k_idx in branches_to_check:
            v_vec = spec["vectors"][k_idx]

            T_zero = v_vec[0] + mp.sqrt(2) * sum(v_vec[m] for m in range(1, dim_sub))
            alpha_N = sum(a_vec[k] * v_vec[k] for k in range(1, dim_sub))

            for M in CUTOFF_GRID:
                if M >= N_sub:
                    continue

                xi_sq = sum(
                    (alpha_N - (T_zero / mp.sqrt(2)) * a_vec[m]) ** 2
                    for m in range(M + 1, dim_sub)
                )
                norm_xi_Q = mp.sqrt(xi_sq)

                term1 = abs(alpha_N) * mp.sqrt(mp.mpf(N_sub - M))
                norm_a_Q = mp.sqrt(sum(a_vec[k] ** 2 for k in range(M + 1, dim_sub)))
                term2 = (abs(T_zero) / mp.sqrt(2)) * norm_a_Q
                tri_bound = term1 + term2

                ratio = tri_bound / norm_xi_Q if norm_xi_Q > mp.mpf("1e-65") else mp.mpf("1.0")
                holds = bool(norm_xi_Q <= tri_bound + mp.mpf("1e-65"))

                print(
                    f"{N_sub:>4} | "
                    f"{b_name:>14} | "
                    f"{M:>4} | "
                    f"{mp.nstr(norm_xi_Q, 10):>24} | "
                    f"{mp.nstr(tri_bound, 10):>22} | "
                    f"{mp.nstr(ratio, 4):>10} | "
                    f"{str(holds):>7}"
                )

    print("-" * 97)
    print()

    # Step 8: Synthesis Scorecard
    print("=" * 80)
    print("CELL 112A SYNTHESIS: MULTI-EIGENVALUE BRANCH TRACKING SCORECARD")
    print("=" * 80)
    print("1. Spectrum Tracking (Lowest 5 Eigenpairs Tracked):                 COMPLETE")
    print("2. Branch Identification (Ground State vs Solitary Wave):          TRACKED")
    print("3. Eigensolver Residual Certification (||(H - E I) v|| < 1e-58):   EVALUATED")
    print("4. Multi-Route alpha_N Consistency (|R1 - R2| < 1e-58):             EVALUATED")
    print("5. Boundary Flux/Contact Coupling Forensics:                       TRACKED")
    print("6. Extinction Progression P_alpha(N) Across Branches:               AUDITED")
    print(f"Total script runtime: {time.time() - t_start:.2f} s")
    print()

    # Unambiguous Completion Sentinel
    print("=" * 80)
    print("CELL 112A EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
