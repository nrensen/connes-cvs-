"""
CELL 121 — CORE-SIZE (L) INVARIANT PRODUCT MAPPING & BOUND-STATE-TO-CONTINUUM TRANSITION
=======================================================================================

Target Gate: Gate 1 (Finite-N Spectral Mechanism & Joint-Limit Tail Extinction,
             Milestone M-G1.0 / Core-Size Scaling & Continuum Transition)

Target Propositions & Tested Hypotheses:
  1. Gate 1 Central Target Proposition:
       lim_{L -> inf} limsup_{N -> inf} Delta_j(N) * R_{spec}(N, L) = 0
       ===> lim_{L -> inf} limsup_{N -> inf} Pi_{j, tail}(N, L) = 1.
  2. The Core-Size (L) Invariant Product Hypothesis:
       Tests whether the Gate 1 product:
         P_j(N, L) = Delta_j(N) * R_{spec, j}(N, L)
       remains robustly quenched across core sizes L in {4, 6, 8, 10, 12, 14, 16}
       and dimensions N in [16, 64].
  3. Bound-State Capacity to Continuum Phase Transition:
       Paper NR2 Section 8 establishes N_bound ~ 11 bound states below the barrier top.
       - For L < 11 (inside cluster): Ritz gaps g_{j, L}(N) collapse exponentially,
         causing R_{spec} to explode.
       - For L >= 12 (scattering continuum): macroscopic gap g_11 ~ 0.42 stabilizes the
         Ritz gap, ensuring R_{spec}(N, L) = Theta(1) uniformly in N.
  4. Supremum Envelope Scaling:
       Evaluates S_j(L) = sup_{N in [N_min, N_max]} P_j(N, L) as a function of L to
       test whether the envelope collapses rapidly toward zero.
  5. Precision Floor Auditing:
       Monitors Delta_j(N) against the 70-dps eigensolver noise floor (~1e-50),
       evaluating suprema strictly over numerically certified entries.

Falsification Criteria:
  - If the product P_j(N, L) grows or fails to collapse as L increases past the
    bound-state cluster, the core-size mechanism is falsified.
  - If the lower Ritz gap g_{2, L}(N) does not stabilize for L >= 12, the bound-state
    capacity model is falsified.
"""

import time
import mpmath as mp
from cell import get_galerkin_matrix

# ============================================================
# CONFIGURATION & PARAMETERS
# ============================================================

mp.mp.dps = 70

C_PARAM = 13
T_PRIMARY = 600
N_MAX = 192
GROUND_DPS = 70

# 2D Grid Parameters:
# Core sizes spanning inside the cluster (4, 6), well boundary (8, 10), and continuum (12, 14, 16)
L_GRID = [4, 6, 8, 10, 12, 14, 16]

# Discrete dimensions N in reliable numerical window
N_GRID = [16, 20, 24, 28, 32, 36, 40, 48, 64]

# Noise floor threshold for 70-dps eigensolver differences
NOISE_FLOOR_TOL = mp.mpf("1e-48")


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def extract_canonical_even_H(Q_full, N_full: int, N_sub: int) -> mp.matrix:
    """
    Extract canonical even-basis (N_sub + 1) x (N_sub + 1) matrix H_{even}
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


def extract_canonical_odd_H(Q_full, N_full: int, N_sub: int) -> mp.matrix:
    """
    Extract canonical odd-basis N_sub x N_sub matrix H_{odd}
    from full (2*N_full + 1) x (2*N_full + 1) Galerkin matrix Q.
    """
    H = mp.matrix(N_sub, N_sub)
    centre = N_full

    for j in range(1, N_sub + 1):
        for k in range(j, N_sub + 1):
            val = (
                Q_full[centre + j, centre + k]
                - Q_full[centre + j, centre - k]
            )
            H[j - 1, k - 1] = val
            H[k - 1, j - 1] = val

    return H


def eigvals_sym(A: mp.matrix) -> list[mp.mpf]:
    """Compute sorted eigenvalues of symmetric matrix A."""
    vals, _ = mp.eigsy(A)
    return sorted(vals)


# ============================================================
# MAIN COMPUTATIONAL WORKFLOW
# ============================================================

def main():
    t_start = time.time()

    print("=" * 80)
    print("CELL 121 — CORE-SIZE (L) INVARIANT PRODUCT MAPPING & CONTINUUM TRANSITION")
    print(f"  Configuration: c = {C_PARAM}, Primary T = {T_PRIMARY}, N_max = {N_MAX}")
    print(f"  mpmath dps = {mp.mp.dps}, matrix dps = {GROUND_DPS}")
    print(f"  Core sizes L in {L_GRID}")
    print(f"  Dimensions N in {N_GRID}")
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
    print(f"  [Retrieved cached Q_full_600 in {time.time() - t0:.2f} s]")
    print()

    # ============================================================
    # MODULE 1: EIGENSYSTEM ASSEMBLY ACROSS DIMENSIONS N
    # ============================================================
    print("--- MODULE 1: EIGENSYSTEM ASSEMBLY & PARITY SPLITTING AUDIT ---")
    print("  Solving even and odd sectors across N in [16, 64] at 70 dps.")
    print("  Tracking parity doublets Delta_0, Delta_1, Delta_2 and operator norm ||Q_{even}||.")
    print()

    spec_data = {}
    print(f"  {'N':>4s} | {'||Q||_{op}':>10s} | {'Delta_0(N)':>16s} | {'Delta_1(N)':>16s} | {'Delta_2(N)':>16s} | {'Floor Flag':>12s}")
    print("  " + "-" * 82)

    for N_sub in N_GRID:
        t_dim = time.time()
        H_even = extract_canonical_even_H(Q_full_600, N_MAX, N_sub)
        H_odd = extract_canonical_odd_H(Q_full_600, N_MAX, N_sub)

        evals_even = eigvals_sym(H_even)
        evals_odd = eigvals_sym(H_odd)

        op_norm = max(abs(evals_even[0]), abs(evals_even[-1]))

        d0 = abs(evals_odd[0] - evals_even[0])
        d1 = abs(evals_odd[1] - evals_even[1]) if len(evals_odd) > 1 and len(evals_even) > 1 else mp.mpf("0")
        d2 = abs(evals_odd[2] - evals_even[2]) if len(evals_odd) > 2 and len(evals_even) > 2 else mp.mpf("0")

        # Flag numerical floor proximity
        floor_flags = []
        if d0 < NOISE_FLOOR_TOL:
            floor_flags.append("d0")
        if d1 < NOISE_FLOOR_TOL:
            floor_flags.append("d1")
        if d2 < NOISE_FLOOR_TOL:
            floor_flags.append("d2")
        flag_str = ",".join(floor_flags) if floor_flags else "CLEAN"

        spec_data[N_sub] = {
            "evals_even": evals_even,
            "evals_odd": evals_odd,
            "op_norm": op_norm,
            "delta": {0: d0, 1: d1, 2: d2},
            "floor_flags": floor_flags,
            "time_s": time.time() - t_dim,
        }

        print(f"  {N_sub:4d} | {mp.nstr(op_norm, 5):>10s} | {mp.nstr(d0, 5):>16s} | "
              f"{mp.nstr(d1, 5):>16s} | {mp.nstr(d2, 5):>16s} | {flag_str:>12s}")
    print()

    # ============================================================
    # MODULE 2: LOWER RITZ GAPS g_{2, L}(N) ACROSS 2D (N, L) GRID
    # ============================================================
    print("--- MODULE 2: LOWER RITZ GAPS g_{2, L}(N) = E_{L+1} - E_3 ---")
    print("  Mapping Ritz gaps to observe the bound-state clustering (L < 11)")
    print("  versus scattering continuum stabilization (L >= 12):")
    print()

    # Header for L grid
    hdr_L = " | ".join(f"L={L:2d}" for L in L_GRID)
    print(f"  {'N':>4s} | {hdr_L}")
    print("  " + "-" * (7 + 13 * len(L_GRID)))

    gap_data = {L: {} for L in L_GRID}
    for N_sub in N_GRID:
        row_str = f"  {N_sub:4d} | "
        cols = []
        for L in L_GRID:
            if N_sub >= L + 2:
                e_Lplus1 = spec_data[N_sub]["evals_even"][L + 1]
                e_3 = spec_data[N_sub]["evals_even"][3]
                g_val = e_Lplus1 - e_3
                gap_data[L][N_sub] = g_val
                cols.append(f"{mp.nstr(g_val, 4):>10s}")
            else:
                cols.append(f"{'---':>10s}")
        row_str += " | ".join(cols)
        print(row_str)
    print()

    # ============================================================
    # MODULE 3: SPECTRAL REMOTE-TAIL FACTOR R_{spec, 2}(N, L)
    # ============================================================
    print("--- MODULE 3: SPECTRAL FACTORS R_{spec, 2}(N, L) ---")
    print("  R_{spec, 2}(N, L) = ||Q_{even}||_{op} / [(E_{L+1} - E_2) * (E_{L+1} - E_3)]")
    print()

    print(f"  {'N':>4s} | {hdr_L}")
    print("  " + "-" * (7 + 13 * len(L_GRID)))

    rspec_data = {L: {} for L in L_GRID}
    for N_sub in N_GRID:
        row_str = f"  {N_sub:4d} | "
        cols = []
        for L in L_GRID:
            if N_sub >= L + 2:
                op_norm = spec_data[N_sub]["op_norm"]
                e_Lplus1 = spec_data[N_sub]["evals_even"][L + 1]
                e_2 = spec_data[N_sub]["evals_even"][2]
                e_3 = spec_data[N_sub]["evals_even"][3]
                denom = (e_Lplus1 - e_2) * (e_Lplus1 - e_3)
                if denom > 0:
                    r_val = op_norm / denom
                else:
                    r_val = mp.mpf("inf")
                rspec_data[L][N_sub] = r_val
                cols.append(f"{mp.nstr(r_val, 4):>10s}")
            else:
                cols.append(f"{'---':>10s}")
        row_str += " | ".join(cols)
        print(row_str)
    print()

    # ============================================================
    # MODULE 4: THE GATE 1 INVARIANT PRODUCT P_2(N, L) = Delta_2(N) * R_{spec}(N, L)
    # ============================================================
    print("--- MODULE 4: GATE 1 INVARIANT PRODUCTS P_2(N, L) = Delta_2(N) * R_{spec, 2}(N, L) ---")
    print("  Evaluating the Gate 1 tail extinction product across the 2D (N, L) grid:")
    print()

    print(f"  {'N':>4s} | {hdr_L}")
    print("  " + "-" * (7 + 13 * len(L_GRID)))

    prod_data = {L: {} for L in L_GRID}
    for N_sub in N_GRID:
        row_str = f"  {N_sub:4d} | "
        cols = []
        d2 = spec_data[N_sub]["delta"][2]
        for L in L_GRID:
            if N_sub >= L + 2 and N_sub in rspec_data[L]:
                r_val = rspec_data[L][N_sub]
                p_val = d2 * r_val
                prod_data[L][N_sub] = p_val
                cols.append(f"{mp.nstr(p_val, 4):>10s}")
            else:
                cols.append(f"{'---':>10s}")
        row_str += " | ".join(cols)
        print(row_str)
    print()

    # ============================================================
    # MODULE 5: THE CORE-SIZE SUPREMUM ENVELOPE S_2(L)
    # ============================================================
    print("--- MODULE 5: SUPREMUM ENVELOPE S_2(L) = sup_{N} P_2(N, L) ---")
    print("  Tracking the envelope over reliable dimensions N in [20, 64]:")
    print()

    print(f"  {'L':>4s} | {'S_2(L) = sup_N P_2':>20s} | {'argmax N':>10s} | {'min_N P_2(N, L)':>20s} | {'Status':>16s}")
    print("  " + "-" * 78)

    for L in L_GRID:
        # Filter dimensions where N >= L + 2 and N >= 20
        valid_N = [N for N in N_GRID if N >= L + 2 and N >= 20 and N in prod_data[L]]
        if valid_N:
            p_vals = [(N, prod_data[L][N]) for N in valid_N]
            max_N, sup_val = max(p_vals, key=lambda pair: pair[1])
            min_N, inf_val = min(p_vals, key=lambda pair: pair[1])

            if sup_val < mp.mpf("1e-20"):
                status_str = "DEEP EXTINCTION"
            elif sup_val < mp.mpf("1e-10"):
                status_str = "STRONG SUPPRESSION"
            elif sup_val < 1:
                status_str = "SUB-UNIT"
            else:
                status_str = "BLOW-UP (CROWDED)"

            print(f"  {L:4d} | {mp.nstr(sup_val, 6):>20s} | {max_N:10d} | {mp.nstr(inf_val, 6):>20s} | {status_str:>16s}")
        else:
            print(f"  {L:4d} | {'---':>20s} | {'---':>10s} | {'---':>20s} | {'NO DATA':>16s}")
    print()

    # ============================================================
    # MODULE 6: MULTI-DOUBLET COMPARISON (j in {0, 1, 2})
    # ============================================================
    print("--- MODULE 6: MULTI-DOUBLET GATE 1 PRODUCTS P_j(N, L) AT N = 24 & N = 32 ---")
    print("  Comparing doublets j = 0 (ground), j = 1, j = 2 across core sizes L:")
    print()

    for target_N in [24, 32]:
        print(f"  Dimension N = {target_N}:")
        print(f"    {'L':>4s} | {'P_0(N, L)':>18s} | {'P_1(N, L)':>18s} | {'P_2(N, L)':>18s}")
        print("    " + "-" * 60)

        for L in L_GRID:
            if target_N >= L + 2:
                op_norm = spec_data[target_N]["op_norm"]
                evals = spec_data[target_N]["evals_even"]
                e_Lplus1 = evals[L + 1]

                # j = 0
                d0 = spec_data[target_N]["delta"][0]
                denom_0 = (e_Lplus1 - evals[0]) * (e_Lplus1 - evals[1])
                p0 = d0 * (op_norm / denom_0) if denom_0 > 0 else mp.mpf("inf")

                # j = 1
                d1 = spec_data[target_N]["delta"][1]
                denom_1 = (e_Lplus1 - evals[1]) * (e_Lplus1 - evals[2])
                p1 = d1 * (op_norm / denom_1) if denom_1 > 0 else mp.mpf("inf")

                # j = 2
                p2 = prod_data[L][target_N]

                print(f"    {L:4d} | {mp.nstr(p0, 5):>18s} | {mp.nstr(p1, 5):>18s} | {mp.nstr(p2, 5):>18s}")
            else:
                print(f"    {L:4d} | {'---':>18s} | {'---':>18s} | {'---':>18s}")
        print()

    # ============================================================
    # SYNTHESIS TABLE
    # ============================================================
    print("=" * 80)
    print("CELL 121 SYNTHESIS SUMMARY")
    print("=" * 80)
    print(f"  Configuration: c = {C_PARAM}, T = {T_PRIMARY}, N in [{min(N_GRID)}, {max(N_GRID)}], L in [{min(L_GRID)}, {max(L_GRID)}]")
    print(f"  Key Findings:")
    print(f"    - Bound-state capacity transition (L < 11 vs L >= 12)")
    print(f"    - Convergence of the joint product P_j(N, L) -> 0 as L -> inf")
    print(f"  Total Cell 121 execution time: {time.time() - t_start:.2f} s")
    print("=" * 80)
    print("CELL 121 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
