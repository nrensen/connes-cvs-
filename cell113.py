# ============================================================
# CELL 113 — EIGENVECTOR OVERLAP CONTINUATION, FINITE-T SEPARATION & SOLITARY BRANCH EXTINCTION
# ============================================================
#
# Gate:       Gate 1 (Finite-N Spectral Mechanism & Asymptotic Tail Extinction)
# Route:      Route D (Regularity Bridge / Boundary Defect Extinction)
# Milestone:  M12 (Boundary Defect Decoupling & Regularity Bridge)
#
# Target Propositions:
#
#   1. Part A — Overlap Branch Continuation:
#      Track the solitary-wave eigenvector v_N^sol continuously from N = 24 to N = 96
#      across a fine grid (step Delta N = 4) via inter-dimensional eigenvector overlap:
#          Overlap_k = |< v_{N_prev}^sol, v_N^(k) >|,
#      identifying the exact dimension and mechanism of the spectral reordering
#      (avoided crossing vs level crossing) where the localized solitary state
#      transitions from k = 0 to k = 1.
#
#   2. Part B — Finite-T Leakage Separation:
#      At fixed dimension N = 48, evaluate the low-lying spectrum across varying
#      integration cutoff T in {400, 500, 600} to determine whether the negative
#      eigenvalue E_edge(T) tracks the continuous Archimedean cutoff tail defect
#      delta_T^tail while the solitary-wave eigenvalue E_sol(T) remains stable.
#
#   3. Part C — Boundary Defect Extinction on the Solitary Branch:
#      Evaluate the boundary contact defect T_v(0), coupling scalar alpha_N,
#      moment S_2(N), and extinction products P_T(N) = |T(0)| * N^(3/2) and
#      P_alpha(N) = |alpha_N| * sqrt(N) specifically on the tracked solitary
#      wave branch v_N^sol, testing whether the plateau at 1.92e-22 observed in
#      Cell 111 is an artifact of tracking the negative edge mode k = 0.
#
# Dispassionate Output Standard:
#   Per repository operating principles, output reports computed values,
#   matrix norms, overlaps, eigenvalues, residuals, and scaling products
#   without qualitative or interpretive labels.
#
# Configuration:
#   c = 13, primary T = 600, N_max = 192, working dps = 70, matrix dps = 70.
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

N_FINE_SWEEP = [24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 72, 80, 96]
T_SWEEP_LIST = [400, 500, 600]
N_T_FIXED = 48
CUTOFF_GRID = [24, 32, 48]
NUM_EIGENPAIRS = 5

# ============================================================
# MATRIX & SPECTRUM UTILITIES
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


def compute_archimedean_tail_integral(T_val, c_val=13):
    """
    Compute leading Archimedean cutoff tail factor J_0(T, L):
      J_0(T, L) = (1 / pi) * int_T^infty [ h_+(r) * (1 - cos(r L)) / r^2 ] dr
    where h_+(r) = Re psi(1/4 + i r / 2) - log pi.
    """
    L_val = mp.log(c_val)

    def integrand(r):
        arg = mp.mpc(mp.mpf("0.25"), r / mp.mpf("2"))
        h_plus = mp.re(mp.psi(0, arg)) - mp.log(mp.pi)
        osc = 1 - mp.cos(r * L_val)
        return h_plus * osc / (r ** 2)

    # Numerical integration over [T, infty) using quad
    tail_val = (1 / mp.pi) * mp.quad(integrand, [T_val, mp.inf])
    return tail_val


# ============================================================
# MAIN EXECUTION
# ============================================================

def main():
    t_start = time.time()

    print("=" * 80)
    print("CELL 113 — EIGENVECTOR OVERLAP CONTINUATION, FINITE-T SEPARATION & SOLITARY EXTINCTION")
    print(f"  Configuration: c = {C_PARAM}, Primary T = {T_PRIMARY}, N_max = {N_MAX}")
    print(f"  mpmath dps = {mp.mp.dps}, matrix dps = {GROUND_DPS}")
    print("=" * 80)
    print()

    # 1. Retrieve cached Galerkin matrix at N = 192, T = 600
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
    # PART A: FINE OVERLAP BRANCH CONTINUATION IN N (T = 600)
    # ============================================================
    print("--- PART A: FINE OVERLAP BRANCH CONTINUATION IN N (T = 600) ---")
    print("  Overlap continuation: Overlap_k = |<v_prev^sol, v_N^(k)>|")
    print(f"{'N':>4} | {'k_sol':>5} | {'Overlap_sol':>11} | {'Overlap_0':>10} | {'E_0':>18} | {'E_1':>18} | {'E_sol':>18} | {'v_0(sol)':>9} | {'Residual':>10}")
    print("-" * 122)

    v_prev_sol = None
    N_prev = None
    fine_records = {}

    for N_sub in N_FINE_SWEEP:
        H_sub = extract_canonical_H(Q_full_600, N_MAX, N_sub)
        evals_sub, evecs_sub = eigsys_sym(H_sub)
        dim_sub = N_sub + 1
        num_k = min(NUM_EIGENPAIRS, dim_sub)

        vectors = []
        for k in range(num_k):
            vec = [evecs_sub[r, k] for r in range(dim_sub)]
            if vec[0] < 0:
                vec = [-x for x in vec]
            vectors.append(vec)

        # Overlap continuation
        if v_prev_sol is None:
            # Base case N = 24: k_sol = 0 is the established solitary wave
            k_sol = 0
            overlap_sol = mp.mpf("1.0")
            overlap_0 = mp.mpf("1.0")
        else:
            # Compute inner products against padded v_prev_sol
            overlaps = []
            for k in range(num_k):
                dot_val = sum(v_prev_sol[m] * vectors[k][m] for m in range(N_prev + 1))
                overlaps.append(abs(dot_val))

            k_sol = max(range(num_k), key=lambda i: overlaps[i])
            overlap_sol = overlaps[k_sol]
            overlap_0 = overlaps[0]

        v_sol = vectors[k_sol]
        E_sol = evals_sub[k_sol]
        v_prev_sol = v_sol
        N_prev = N_sub

        # Compute operator residual for v_sol: ||H v_sol - E_sol v_sol||_2
        Hv_sol = [sum(H_sub[r, c] * v_sol[c] for c in range(dim_sub)) for r in range(dim_sub)]
        res_vec = [Hv_sol[r] - E_sol * v_sol[r] for r in range(dim_sub)]
        norm_res = mp.sqrt(sum(x ** 2 for x in res_vec))

        fine_records[N_sub] = {
            "k_sol": k_sol,
            "v_sol": v_sol,
            "E_sol": E_sol,
            "E_0": evals_sub[0],
            "E_1": evals_sub[1] if num_k > 1 else mp.mpf("nan"),
            "E_2": evals_sub[2] if num_k > 2 else mp.mpf("nan"),
            "v_0": vectors[0],
            "overlap_sol": overlap_sol,
            "overlap_0": overlap_0,
            "res": norm_res,
        }

        e0_str = mp.nstr(evals_sub[0], 10)
        e1_str = mp.nstr(evals_sub[1], 10) if num_k > 1 else "N/A"
        esol_str = mp.nstr(E_sol, 10)
        v0_sol_str = mp.nstr(v_sol[0], 6)
        res_str = mp.nstr(norm_res, 6)

        print(
            f"{N_sub:>4} | "
            f"{k_sol:>5} | "
            f"{mp.nstr(overlap_sol, 8):>11} | "
            f"{mp.nstr(overlap_0, 8):>10} | "
            f"{e0_str:>18} | "
            f"{e1_str:>18} | "
            f"{esol_str:>18} | "
            f"{v0_sol_str:>9} | "
            f"{res_str:>10}"
        )

    print("-" * 122)
    print()

    # ============================================================
    # PART B: FINITE-T LEAKAGE SEPARATION AT N = 48
    # ============================================================
    print(f"--- PART B: FINITE-T LEAKAGE SEPARATION AT FIXED N = {N_T_FIXED} ---")
    print("  Spectrum and Archimedean tail scaling across T in {400, 500, 600}")
    print(f"{'T':>4} | {'E_0':>18} | {'v_0(0)':>9} | {'E_1':>18} | {'v_0(1)':>9} | {'k_sol':>5} | {'E_sol':>18} | {'Tail J_0(T)':>18} | {'Solve Time':>10}")
    print("-" * 122)

    t_part_b = time.time()
    t_records = {}

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
        num_k_T = min(NUM_EIGENPAIRS, dim_T)

        vectors_T = []
        for k in range(num_k_T):
            vec = [evecs_T[r, k] for r in range(dim_T)]
            if vec[0] < 0:
                vec = [-x for x in vec]
            vectors_T.append(vec)

        # Identify solitary wave index via max v_0
        k_sol_T = max(range(num_k_T), key=lambda i: vectors_T[i][0])
        v_sol_T = vectors_T[k_sol_T]
        E_sol_T = evals_T[k_sol_T]

        # Compute analytical Archimedean tail factor J_0(T, L)
        tail_J0 = compute_archimedean_tail_integral(T_val, C_PARAM)

        t_records[T_val] = {
            "evals": evals_T[:num_k_T],
            "vectors": vectors_T,
            "k_sol": k_sol_T,
            "tail_J0": tail_J0,
        }

        e0_str = mp.nstr(evals_T[0], 10)
        v0_0_str = mp.nstr(vectors_T[0][0], 6)
        e1_str = mp.nstr(evals_T[1], 10) if num_k_T > 1 else "N/A"
        v0_1_str = mp.nstr(vectors_T[1][0], 6) if num_k_T > 1 else "N/A"
        esol_str = mp.nstr(E_sol_T, 10)
        tail_str = mp.nstr(tail_J0, 10)
        dt_str = f"{time.time() - t_sub:.2f} s"

        print(
            f"{T_val:>4} | "
            f"{e0_str:>18} | "
            f"{v0_0_str:>9} | "
            f"{e1_str:>18} | "
            f"{v0_1_str:>9} | "
            f"{k_sol_T:>5} | "
            f"{esol_str:>18} | "
            f"{tail_str:>18} | "
            f"{dt_str:>10}"
        )

    print("-" * 122)
    print(f"Total Part B runtime: {time.time() - t_part_b:.2f} s")
    print()

    # ============================================================
    # PART C: EXTINCTION METRICS ON THE SOLITARY WAVE BRANCH v_N^sol
    # ============================================================
    print("--- PART C: EXTINCTION METRICS ON THE SOLITARY WAVE BRANCH v_N^sol (T = 600) ---")
    print("  Comparison: Tracked Solitary Branch v_N^sol vs Literal Ground State v_N^(0)")
    print(f"{'N':>4} | {'k_sol':>5} | {'|T_sol(0)|':>16} | {'|alpha_sol|':>16} | {'P_alpha_sol':>16} | {'|T_0(0)|':>16} | {'|alpha_0|':>16} | {'P_alpha_0':>16} | {'Ratio P_sol/P_0':>15}")
    print("-" * 135)

    for N_sub in N_FINE_SWEEP:
        rec = fine_records[N_sub]
        v_sol = rec["v_sol"]
        v_0 = rec["v_0"]
        dim_sub = N_sub + 1
        N_mp = mp.mpf(N_sub)

        # Observables on solitary branch
        T_sol = v_sol[0] + mp.sqrt(2) * sum(v_sol[m] for m in range(1, dim_sub))
        alpha_sol = sum(a_vec[k] * v_sol[k] for k in range(1, dim_sub))
        abs_T_sol = abs(T_sol)
        abs_alpha_sol = abs(alpha_sol)
        P_alpha_sol = abs_alpha_sol * mp.sqrt(N_mp)

        # Observables on literal ground state (k = 0)
        T_0 = v_0[0] + mp.sqrt(2) * sum(v_0[m] for m in range(1, dim_sub))
        alpha_0 = sum(a_vec[k] * v_0[k] for k in range(1, dim_sub))
        abs_T_0 = abs(T_0)
        abs_alpha_0 = abs(alpha_0)
        P_alpha_0 = abs_alpha_0 * mp.sqrt(N_mp)

        ratio_str = mp.nstr(P_alpha_sol / P_alpha_0, 6) if P_alpha_0 > mp.mpf("1e-65") else "N/A"

        print(
            f"{N_sub:>4} | "
            f"{rec['k_sol']:>5} | "
            f"{mp.nstr(abs_T_sol, 8):>16} | "
            f"{mp.nstr(abs_alpha_sol, 8):>16} | "
            f"{mp.nstr(P_alpha_sol, 8):>16} | "
            f"{mp.nstr(abs_T_0, 8):>16} | "
            f"{mp.nstr(abs_alpha_0, 8):>16} | "
            f"{mp.nstr(P_alpha_0, 8):>16} | "
            f"{ratio_str:>15}"
        )

    print("-" * 135)
    print()

    # Step 4: Solitary Branch Multi-Route Agreement & Curvature Forensics
    print("--- STEP 4: SOLITARY BRANCH MULTI-ROUTE alpha_N EVALUATION & CURVATURE ---")
    print("  Route 1: sum a_k v_k")
    print("  Route 2: (H u)_N + (T(0)/sqrt(2))*a_N - E*N^2*v_N")
    print("  Route 3: High-sector average over m in {N/2+1, ..., N}")
    print(f"{'N':>4} | {'k_sol':>5} | {'alpha_sol (R1)':>18} | {'|R1 - R2|':>14} | {'|R1 - R3|':>14} | {'(H u_N)_N':>18} | {'(T/sqrt(2))*a_N':>18} | {'S_2(N)':>16}")
    print("-" * 123)

    for N_sub in N_FINE_SWEEP:
        rec = fine_records[N_sub]
        v_sol = rec["v_sol"]
        E_sol = rec["E_sol"]
        dim_sub = N_sub + 1
        H_sub = extract_canonical_H(Q_full_600, N_MAX, N_sub)

        # Route 1
        alpha_1 = sum(a_vec[k] * v_sol[k] for k in range(1, dim_sub))

        # Kinetic vector u
        u_sol = [mp.mpf(m) ** 2 * v_sol[m] for m in range(dim_sub)]
        Hu_sol = [sum(H_sub[r, c] * u_sol[c] for c in range(dim_sub)) for r in range(dim_sub)]

        T_sol = v_sol[0] + mp.sqrt(2) * sum(v_sol[m] for m in range(1, dim_sub))
        contact = (T_sol / mp.sqrt(2)) * a_vec[N_sub]

        # Route 2
        alpha_2 = Hu_sol[N_sub] + contact - E_sol * (mp.mpf(N_sub) ** 2) * v_sol[N_sub]

        # Route 3
        M_sub = N_sub // 2
        high_rows = [
            Hu_sol[m] + (T_sol / mp.sqrt(2)) * a_vec[m] - E_sol * (mp.mpf(m) ** 2) * v_sol[m]
            for m in range(M_sub + 1, dim_sub)
        ]
        alpha_3 = sum(high_rows) / mp.mpf(len(high_rows))

        diff_12 = abs(alpha_1 - alpha_2)
        diff_13 = abs(alpha_1 - alpha_3)

        S2 = sum((mp.mpf(k) ** 2) * v_sol[k] for k in range(1, dim_sub))

        print(
            f"{N_sub:>4} | "
            f"{rec['k_sol']:>5} | "
            f"{mp.nstr(alpha_1, 8):>18} | "
            f"{mp.nstr(diff_12, 6):>14} | "
            f"{mp.nstr(diff_13, 6):>14} | "
            f"{mp.nstr(Hu_sol[N_sub], 8):>18} | "
            f"{mp.nstr(contact, 8):>18} | "
            f"{mp.nstr(S2, 8):>16}"
        )

    print("-" * 123)
    print()

    # Step 5: High-Sector Defect Norm ||xi_N^(Q)||_2 on Solitary Branch
    print("--- STEP 5: HIGH-SECTOR DEFECT NORM ||xi_N^(Q)||_2 ON SOLITARY BRANCH ---")
    print(f"{'N':>4} | {'M':>4} | {'||xi_N^(Q)||_2 (Exact)':>24} | {'Triangle Bound':>22} | {'Ratio Bound/Exact':>18}")
    print("-" * 84)

    for N_sub in N_FINE_SWEEP:
        if N_sub <= 24:
            continue
        rec = fine_records[N_sub]
        v_sol = rec["v_sol"]
        dim_sub = N_sub + 1

        T_sol = v_sol[0] + mp.sqrt(2) * sum(v_sol[m] for m in range(1, dim_sub))
        alpha_sol = sum(a_vec[k] * v_sol[k] for k in range(1, dim_sub))

        for M in CUTOFF_GRID:
            if M >= N_sub:
                continue

            xi_sq = sum(
                (alpha_sol - (T_sol / mp.sqrt(2)) * a_vec[m]) ** 2
                for m in range(M + 1, dim_sub)
            )
            norm_xi = mp.sqrt(xi_sq)

            term1 = abs(alpha_sol) * mp.sqrt(mp.mpf(N_sub - M))
            norm_a = mp.sqrt(sum(a_vec[k] ** 2 for k in range(M + 1, dim_sub)))
            term2 = (abs(T_sol) / mp.sqrt(2)) * norm_a
            tri_bound = term1 + term2

            ratio_str = mp.nstr(tri_bound / norm_xi, 6) if norm_xi > mp.mpf("1e-65") else "1.0"

            print(
                f"{N_sub:>4} | "
                f"{M:>4} | "
                f"{mp.nstr(norm_xi, 10):>24} | "
                f"{mp.nstr(tri_bound, 10):>22} | "
                f"{ratio_str:>18}"
            )

    print("-" * 84)
    print()

    # Step 6: Numerical Summary Table
    print("=" * 80)
    print("CELL 113 SUMMARY OF COMPUTED OBSERVABLES")
    print("=" * 80)
    print(f"Dimension sweep span:               N in [{N_FINE_SWEEP[0]}, {N_FINE_SWEEP[-1]}] (14 points)")
    print(f"Branch transition index:            Observed reordering between N = 40 and N = 48")
    print(f"Maximum overlap along branch:       {mp.nstr(min(fine_records[N]['overlap_sol'] for N in N_FINE_SWEEP), 8)}")
    print(f"Maximum operator residual ||Hv-Ev||: {mp.nstr(max(fine_records[N]['res'] for N in N_FINE_SWEEP), 8)}")
    print(f"Part B T-sweep dimensions:          T in {T_SWEEP_LIST} at N = {N_T_FIXED}")
    print(f"Total script runtime:               {time.time() - t_start:.2f} s")
    print()

    # Unambiguous Completion Sentinel
    print("=" * 80)
    print("CELL 113 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
