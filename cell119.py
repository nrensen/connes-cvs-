"""
CELL 119 — EXACT ALGEBRAIC MODAL DECOMPOSITION (A_{N, k} vs B_{N, k} CANCELLATION BALANCE)
========================================================================================

Target Gate: Gate 1 (Finite-N Spectral Mechanism & Joint-Limit Tail Extinction,
             Milestone M12 / Route D / Phase IX)

Target Propositions & Tested Hypotheses:
  1. Exact Modal Split: For k in [1, N-1], the boundary-row divided-difference
     kernel summand decomposes identically into:
       F_{N, k} = A_{N, k} - B_{N, k}
     where:
       A_{N, k} = a_N * k^2 * v_{N, k} / (N^2 - k^2) = 2 * N * psi(N) * k^2 * v_{N, k} / (N^2 - k^2)
       B_{N, k} = a_k * k^2 * v_{N, k} / (N^2 - k^2) = 2 * k^3 * psi(k) * v_{N, k} / (N^2 - k^2)
     with a_m = 2 * m * psi(m) = sqrt(2) * m^2 * H_{0, m}.
  2. Total Sums Separation (Hypothesis H1 vs H2):
     Evaluates Sigma_A^{off}(N) = sum_{k=1}^{N-1} A_{N, k} and
     Sigma_B^{off}(N) = sum_{k=1}^{N-1} B_{N, k} across N in [64, 192].
     Tests whether:
       - H1 (Inter-Sum Cancellation): Both Sigma_A and Sigma_B are macroscopic (O(10^-2) or O(1))
         and cancel each other to 10^-22, or
       - H2 (Discrete Moments Pre-Cancellation): Sigma_A and Sigma_B are individually already
         suppressed to tiny values (<= 10^-15) by internal modal cancellation.
  3. Cumulative Trajectory Tracking (Hypothesis H3):
     Evaluates S_A(x) = sum_{k <= xN} A_{N, k} and S_B(x) = sum_{k <= xN} B_{N, k}
     along continuum coordinate x in (0, 1) to determine the exact crossover threshold x*.
  4. Moments Balance:
     Compares Sigma_A and Sigma_B against the Taylor expansion moments:
       M_2(N) = sum_{k=1}^{N-1} k^2 v_{N, k},     M_4(N) = sum_{k=1}^{N-1} k^4 v_{N, k}
       M_{psi, 3}(N) = 2 sum_{k=1}^{N-1} k^3 psi(k) v_{N, k}
  5. Theorem 1 Synthesis:
     Connects Sigma_A and Sigma_B to alpha_N, (a_N / sqrt(2)) * T_v(0), and the
     exact row-wise resolvent identity.

Falsification Criterion:
  - If |Sigma_A^{off}| and |Sigma_B^{off}| are both <= 10^-15, H1 is falsified.
  - If max(|Sigma_A|, |Sigma_B|) / |Sigma_A - Sigma_B| < 10^5, massive inter-sum
    cancellation is falsified.
"""

import time
import mpmath as mp
from connes_cvs.operator import prime_powers_up_to
from cell import get_galerkin_matrix

# ============================================================
# CONFIGURATION & PARAMETERS
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
    mp.mpf("0.80"), mp.mpf("0.90"), mp.mpf("0.95"),
]

# ============================================================
# HELPER FUNCTIONS
# ============================================================

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
# MAIN COMPUTATIONAL WORKFLOW
# ============================================================

def main():
    t_start = time.time()

    print("=" * 80)
    print("CELL 119 — EXACT ALGEBRAIC MODAL DECOMPOSITION (A_{N, k} vs B_{N, k})")
    print(f"  Configuration: c = {C_PARAM}, Primary T = {T_PRIMARY}, N_max = {N_MAX}")
    print(f"  mpmath dps = {mp.mp.dps}, matrix dps = {GROUND_DPS}")
    print("=" * 80)
    print()

    # Step 0: Retrieve cached Hamiltonian
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

    # Precompute modal weight symbol a_k = sqrt(2) * k^2 * H_{0, k} = 2 * k * psi(k)
    dim_192 = N_MAX + 1
    a_vec = [mp.mpf("0")] * dim_192
    psi_vec = [mp.mpf("0")] * dim_192
    for k in range(1, dim_192):
        k_mp = mp.mpf(k)
        H0k = H_192[0, k]
        a_vec[k] = mp.sqrt(mp.mpf("2")) * (k_mp ** 2) * H0k
        psi_vec[k] = a_vec[k] / (mp.mpf("2") * k_mp)

    # Solve localized solitary eigenvector across all dimensions
    print("  Solving localized solitary branch across dimensions...")
    solitary_data = {}
    for N_sub in N_GRID_ALL:
        t_dim = time.time()
        dim_sub = N_sub + 1
        H_sub = extract_canonical_H(Q_full_600, N_MAX, N_sub)

        evals_sub, evecs_sub = eigsys_sym(H_sub)

        # Select localized solitary branch: prominent ground-state component v_0 > 0.4
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
        N_mp = mp.mpf(N_sub)
        N_sq = N_mp ** 2
        a_N_val = a_vec[N_sub]
        psi_N_val = psi_vec[N_sub]

        T_zero = v_sub[0] + mp.sqrt(mp.mpf("2")) * sum(v_sub[m] for m in range(1, dim_sub))
        alpha_val = sum(a_vec[m] * v_sub[m] for m in range(1, dim_sub))

        # Off-diagonal decomposition: k in [1, N_sub - 1]
        A_terms = []
        B_terms = []
        F_terms = []
        for k in range(1, N_sub):
            k_mp = mp.mpf(k)
            k_sq = k_mp ** 2
            denom = N_sq - k_sq
            v_k = v_sub[k]

            A_k = a_N_val * k_sq * v_k / denom
            B_k = a_vec[k] * k_sq * v_k / denom
            F_k = A_k - B_k

            A_terms.append(A_k)
            B_terms.append(B_k)
            F_terms.append(F_k)

        # Endpoint diagonal term at k = N_sub
        v_N_val = v_sub[N_sub]
        H_NN = H_sub[N_sub, N_sub]
        F_NN = H_NN * N_sq * v_N_val

        sigma_A_off = sum(A_terms)
        sigma_B_off = sum(B_terms)
        delta_sigma = sigma_A_off - sigma_B_off
        Hu_N_total = delta_sigma + F_NN

        # Check full flux from matrix definition
        full_flux_terms = [H_sub[N_sub, k] * (mp.mpf(k) ** 2) * v_sub[k] for k in range(1, dim_sub)]
        Hu_N_direct = sum(full_flux_terms)
        decomp_diff = abs(Hu_N_total - Hu_N_direct)

        # Cancellation factor
        max_AB = max(abs(sigma_A_off), abs(sigma_B_off))
        c_AB = max_AB / abs(delta_sigma) if delta_sigma != 0 else mp.mpf("inf")

        # Solitary moments
        M2 = sum((mp.mpf(k) ** 2) * v_sub[k] for k in range(1, N_sub))
        M4 = sum((mp.mpf(k) ** 4) * v_sub[k] for k in range(1, N_sub))
        M_psi3 = sum(mp.mpf("2") * (mp.mpf(k) ** 3) * psi_vec[k] * v_sub[k] for k in range(1, N_sub))
        M_psi5 = sum(mp.mpf("2") * (mp.mpf(k) ** 5) * psi_vec[k] * v_sub[k] for k in range(1, N_sub))

        # Taylor leading estimates
        sigma_A_lead = (mp.mpf("2") * psi_N_val / N_mp) * M2
        sigma_B_lead = (mp.mpf("1") / N_sq) * M_psi3
        sigma_A_2nd = sigma_A_lead + (mp.mpf("2") * psi_N_val / (N_mp ** 3)) * M4
        sigma_B_2nd = sigma_B_lead + (mp.mpf("1") / (N_mp ** 4)) * M_psi5

        # Theorem 1 components
        contact_term = (a_N_val / mp.sqrt(mp.mpf("2"))) * T_zero
        eigen_term = E_val * N_sq * v_N_val
        thm1_rhs = alpha_val - contact_term + eigen_term
        thm1_res = abs(Hu_N_direct - thm1_rhs)

        solitary_data[N_sub] = {
            "v": v_sub,
            "E": E_val,
            "v0": v0_val,
            "vN": v_N_val,
            "T_zero": T_zero,
            "alpha": alpha_val,
            "a_N": a_N_val,
            "psi_N": psi_N_val,
            "A_terms": A_terms,
            "B_terms": B_terms,
            "F_terms": F_terms,
            "F_NN": F_NN,
            "sigma_A_off": sigma_A_off,
            "sigma_B_off": sigma_B_off,
            "delta_sigma": delta_sigma,
            "Hu_N_total": Hu_N_total,
            "Hu_N_direct": Hu_N_direct,
            "decomp_diff": decomp_diff,
            "c_AB": c_AB,
            "M2": M2,
            "M4": M4,
            "M_psi3": M_psi3,
            "M_psi5": M_psi5,
            "sigma_A_lead": sigma_A_lead,
            "sigma_B_lead": sigma_B_lead,
            "sigma_A_2nd": sigma_A_2nd,
            "sigma_B_2nd": sigma_B_2nd,
            "contact_term": contact_term,
            "eigen_term": eigen_term,
            "thm1_rhs": thm1_rhs,
            "thm1_res": thm1_res,
        }
        print(f"    N = {N_sub:3d}: solved in {time.time() - t_dim:.2f} s "
              f"(sigma_A = {mp.nstr(sigma_A_off, 4)}, sigma_B = {mp.nstr(sigma_B_off, 4)}, "
              f"delta = {mp.nstr(delta_sigma, 4)})")

    print("  [All 8 dimensions solved successfully]")
    print()

    # ============================================================
    # MODULE 1: LOW-MODE ANATOMY (k <= 10)
    # ============================================================
    print("--- MODULE 1: LOW-MODE ANATOMY OF A_{N, k} AND B_{N, k} (k <= 10) ---")
    print("  Evaluating A_{N, k}, B_{N, k}, F_{N, k} = A - B, and ratio B/A = a_k / a_N")
    print()

    for N_sub in N_GRID:
        dat = solitary_data[N_sub]
        a_N_val = dat["a_N"]
        print(f"  --- Dimension N = {N_sub} (a_N = {mp.nstr(a_N_val, 6)}, psi(N) = {mp.nstr(dat['psi_N'], 6)}) ---")
        print(f"  {'k':>3s} | {'A_{N, k}':>18s} | {'B_{N, k}':>18s} | {'F_{N, k} = A - B':>18s} | {'B/A = a_k/a_N':>14s} | {'k/N':>8s}")
        print("  " + "-" * 88)
        for k in range(1, min(11, N_sub)):
            A_k = dat["A_terms"][k - 1]
            B_k = dat["B_terms"][k - 1]
            F_k = dat["F_terms"][k - 1]
            ratio_BA = B_k / A_k if A_k != 0 else mp.mpf("0")
            k_over_N = mp.mpf(k) / mp.mpf(N_sub)
            print(f"  {k:3d} | {mp.nstr(A_k, 8):>18s} | {mp.nstr(B_k, 8):>18s} | "
                  f"{mp.nstr(F_k, 8):>18s} | {mp.nstr(ratio_BA, 6):>14s} | {mp.nstr(k_over_N, 4):>8s}")
        print()

    # Detailed inspection of the universal peak mode k = 3
    print("  Auditing the Universal Peak Mode k = 3 across all 5 primary dimensions:")
    print(f"  {'N':>4s} | {'A_{N, 3}':>18s} | {'B_{N, 3}':>18s} | {'F_{N, 3}':>18s} | {'B_3 / A_3':>14s} | {'A_3 Share %':>12s}")
    print("  " + "-" * 78)
    for N_sub in N_GRID:
        dat = solitary_data[N_sub]
        A_3 = dat["A_terms"][2]
        B_3 = dat["B_terms"][2]
        F_3 = dat["F_terms"][2]
        ratio_3 = B_3 / A_3 if A_3 != 0 else mp.mpf("0")
        a_share = (A_3 / F_3) * 100 if F_3 != 0 else mp.mpf("0")
        print(f"  {N_sub:4d} | {mp.nstr(A_3, 8):>18s} | {mp.nstr(B_3, 8):>18s} | "
              f"{mp.nstr(F_3, 8):>18s} | {mp.nstr(ratio_3, 6):>14s} | {mp.nstr(a_share, 4) + '%':>12s}")
    print()

    # ============================================================
    # MODULE 2: TOTAL SUMS SEPARATION (Sigma_A vs Sigma_B)
    # ============================================================
    print("--- MODULE 2: TOTAL OFF-DIAGONAL SUMS SEPARATION ACROSS DIMENSIONS ---")
    print("  Testing Hypothesis H1 (Inter-Sum Cancellation) vs H2 (Internal Moments Pre-Cancellation)")
    print()
    print(f"  {'N':>4s} | {'Sigma_A^{off}':>18s} | {'Sigma_B^{off}':>18s} | {'Delta Sigma':>18s} | {'F_{N, N} (diag)':>16s} | {'(H u_N)_N':>18s} | {'C_{AB} Factor':>14s}")
    print("  " + "-" * 118)
    for N_sub in N_GRID_ALL:
        dat = solitary_data[N_sub]
        s_A = dat["sigma_A_off"]
        s_B = dat["sigma_B_off"]
        d_sig = dat["delta_sigma"]
        f_NN = dat["F_NN"]
        hu_N = dat["Hu_N_total"]
        c_ab = dat["c_AB"]
        print(f"  {N_sub:4d} | {mp.nstr(s_A, 8):>18s} | {mp.nstr(s_B, 8):>18s} | "
              f"{mp.nstr(d_sig, 8):>18s} | {mp.nstr(f_NN, 6):>16s} | {mp.nstr(hu_N, 8):>18s} | {mp.nstr(c_ab, 6):>14s}")
    print()

    # ============================================================
    # MODULE 3: MACROSCOPIC CUMULATIVE TRAJECTORY TRACKING
    # ============================================================
    print("--- MODULE 3: CUMULATIVE TRAJECTORY TRACKING S_A(x) vs S_B(x) AT N = 192 ---")
    print("  Sampling cumulative sums on continuum coordinate x = k / N:")
    print("    S_A(x) = sum_{k <= xN} A_{N, k},  S_B(x) = sum_{k <= xN} B_{N, k},  S_F(x) = S_A(x) - S_B(x)")
    print()

    dat_192 = solitary_data[192]
    A_192 = dat_192["A_terms"]
    B_192 = dat_192["B_terms"]

    print(f"  {'x':>5s} | {'m = floor(xN)':>13s} | {'S_A(x)':>18s} | {'S_B(x)':>18s} | {'S_F(x) = S_A - S_B':>20s} | {'|S_F| / max(|S_A|, |S_B|)':>26s}")
    print("  " + "-" * 98)
    for x_val in X_GRID:
        m_idx = int(mp.floor(x_val * 192))
        m_idx = min(m_idx, 191)  # off-diagonal cutoff
        if m_idx < 1:
            m_idx = 1
        S_A_x = sum(A_192[k - 1] for k in range(1, m_idx + 1))
        S_B_x = sum(B_192[k - 1] for k in range(1, m_idx + 1))
        S_F_x = S_A_x - S_B_x
        max_s = max(abs(S_A_x), abs(S_B_x))
        rel_rem = abs(S_F_x) / max_s if max_s != 0 else mp.mpf("0")
        print(f"  {mp.nstr(x_val, 3):>5s} | {m_idx:13d} | {mp.nstr(S_A_x, 8):>18s} | "
              f"{mp.nstr(S_B_x, 8):>18s} | {mp.nstr(S_F_x, 8):>20s} | {mp.nstr(rel_rem, 6):>26s}")
    print()

    # Multi-dimension comparison of trajectory at key thresholds x = 0.05, 0.10, 0.20, 0.50
    print("  Trajectory Comparison across Dimensions at Key Thresholds:")
    for x_target in [mp.mpf("0.05"), mp.mpf("0.10"), mp.mpf("0.20"), mp.mpf("0.50")]:
        print(f"  --- Threshold x = {mp.nstr(x_target, 2)} ---")
        print(f"  {'N':>4s} | {'m':>4s} | {'S_A(x)':>18s} | {'S_B(x)':>18s} | {'S_F(x)':>18s} | {'Cancellation Factor':>20s}")
        print("  " + "-" * 76)
        for N_sub in N_GRID:
            dat = solitary_data[N_sub]
            m_idx = min(int(mp.floor(x_target * N_sub)), N_sub - 1)
            if m_idx < 1:
                m_idx = 1
            S_A = sum(dat["A_terms"][k - 1] for k in range(1, m_idx + 1))
            S_B = sum(dat["B_terms"][k - 1] for k in range(1, m_idx + 1))
            S_F = S_A - S_B
            c_factor = max(abs(S_A), abs(S_B)) / abs(S_F) if S_F != 0 else mp.mpf("inf")
            print(f"  {N_sub:4d} | {m_idx:4d} | {mp.nstr(S_A, 8):>18s} | {mp.nstr(S_B, 8):>18s} | "
                  f"{mp.nstr(S_F, 8):>18s} | {mp.nstr(c_factor, 6):>20s}")
        print()

    # ============================================================
    # MODULE 4: DISCRETE SOLITARY MOMENTS & ASYMPTOTIC EXPANSIONS
    # ============================================================
    print("--- MODULE 4: DISCRETE SOLITARY MOMENTS & ASYMPTOTIC TAYLOR EXPANSIONS ---")
    print("  Evaluating discrete moments M_2 = sum k^2 v_k, M_4 = sum k^4 v_k,")
    print("  M_{psi, 3} = 2 sum k^3 psi(k) v_k, M_{psi, 5} = 2 sum k^5 psi(k) v_k")
    print()
    print(f"  {'N':>4s} | {'M_2(N)':>16s} | {'M_4(N)':>16s} | {'M_{psi, 3}(N)':>16s} | {'M_{psi, 5}(N)':>16s}")
    print("  " + "-" * 76)
    for N_sub in N_GRID_ALL:
        dat = solitary_data[N_sub]
        print(f"  {N_sub:4d} | {mp.nstr(dat['M2'], 6):>16s} | {mp.nstr(dat['M4'], 6):>16s} | "
              f"{mp.nstr(dat['M_psi3'], 6):>16s} | {mp.nstr(dat['M_psi5'], 6):>16s}")
    print()

    print("  Taylor Asymptotic Accuracy Audit:")
    print(f"  {'N':>4s} | {'Sigma_A^{exact}':>18s} | {'Sigma_A^{lead}':>18s} | {'Rel Err A':>12s} | {'Sigma_B^{exact}':>18s} | {'Sigma_B^{lead}':>18s} | {'Rel Err B':>12s}")
    print("  " + "-" * 110)
    for N_sub in N_GRID:
        dat = solitary_data[N_sub]
        sA_ex = dat["sigma_A_off"]
        sA_ld = dat["sigma_A_lead"]
        errA = abs(sA_ex - sA_ld) / abs(sA_ex) if sA_ex != 0 else mp.mpf("0")

        sB_ex = dat["sigma_B_off"]
        sB_ld = dat["sigma_B_lead"]
        errB = abs(sB_ex - sB_ld) / abs(sB_ex) if sB_ex != 0 else mp.mpf("0")
        print(f"  {N_sub:4d} | {mp.nstr(sA_ex, 8):>18s} | {mp.nstr(sA_ld, 8):>18s} | "
              f"{mp.nstr(errA, 4):>12s} | {mp.nstr(sB_ex, 8):>18s} | {mp.nstr(sB_ld, 8):>18s} | {mp.nstr(errB, 4):>12s}")
    print()

    # ============================================================
    # MODULE 5: EXACT THEOREM 1 SYNTHESIS & RESIDUAL AUDIT
    # ============================================================
    print("--- MODULE 5: THEOREM 1 SYNTHESIS & RESOLVENT BALANCE AUDIT ---")
    print("  Verifying exact identity: (H u_N)_N = alpha_N - (a_N / sqrt(2)) * T_v(0) + E_1 * N^2 * v_{N, N}")
    print()
    print(f"  {'N':>4s} | {'alpha_N':>18s} | {'(a_N/sqrt2) T_v(0)':>20s} | {'E_1 N^2 v_N':>16s} | {'(H u_N)_N':>18s} | {'Theorem 1 Residual':>20s}")
    print("  " + "-" * 106)
    for N_sub in N_GRID_ALL:
        dat = solitary_data[N_sub]
        alp = dat["alpha"]
        ctc = dat["contact_term"]
        eig = dat["eigen_term"]
        hu = dat["Hu_N_direct"]
        res = dat["thm1_res"]
        print(f"  {N_sub:4d} | {mp.nstr(alp, 8):>18s} | {mp.nstr(ctc, 8):>20s} | "
              f"{mp.nstr(eig, 6):>16s} | {mp.nstr(hu, 8):>18s} | {mp.nstr(res, 4):>20s}")
    print()

    # ============================================================
    # SYNTHESIS TABLE
    # ============================================================
    print("=" * 80)
    print("CELL 119 SYNTHESIS SUMMARY")
    print("=" * 80)
    print(f"  {'N':>4s} | {'Sigma_A^{off}':>16s} | {'Sigma_B^{off}':>16s} | {'Delta Sigma':>16s} | {'(H u_N)_N':>16s} | {'C_{AB}':>12s} | {'Thm1 Res':>10s}")
    print("  " + "-" * 102)
    for N_sub in N_GRID_ALL:
        dat = solitary_data[N_sub]
        print(f"  {N_sub:4d} | {mp.nstr(dat['sigma_A_off'], 6):>16s} | {mp.nstr(dat['sigma_B_off'], 6):>16s} | "
              f"{mp.nstr(dat['delta_sigma'], 6):>16s} | {mp.nstr(dat['Hu_N_direct'], 6):>16s} | "
              f"{mp.nstr(dat['c_AB'], 4):>12s} | {mp.nstr(dat['thm1_res'], 2):>10s}")
    print()
    print(f"  Total Cell 119 runtime: {time.time() - t_start:.2f} s")
    print("=" * 80)
    print("CELL 119 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
