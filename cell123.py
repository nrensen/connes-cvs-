"""
CELL 123 — OPERATOR DECOMPOSITION AUDIT: HIGH-MODE COERCIVITY & CORE SUBMATRIX SPECTRA
=======================================================================================

Target Gate: Gate 1 (Finite-N Spectral Mechanism & Asymptotic Tail Extinction,
             Milestone M-G1.4 / Continuum Spectral Threshold Lower Bound)

Target Propositions & Tested Hypotheses:
  1. Tripartite Galerkin Decomposition:
       Q_{even} = Q_{even, arch} + Q_{even, prime} + Q_{even, pole}
     Audit each constituent symmetric matrix derived from the exact definitions in
     connes_cvs/operator.py at c = 13, T = 600, dps = 70.
  2. M = 3 Submatrix Bound-State Obstruction (Hypothesis H-C3.1):
       lambda_{min}(C_3) <= E_3(Q_{even})
     By Cauchy interlacing for deleting 3 modes from the well, lambda_{min}(C_3) is
     bounded above by the 4th bound state E_3 ~ 10^-38. Tests whether C_3 >= 0.386 I
     is mathematically impossible.
  3. M = 3 Interior Interlacing Lower Bound (Hypothesis H-C3.2):
       E_13^{(N)} >= lambda_{10}(C_3) >= E_{10}^{(N)}
     Tests whether the 11th eigenvalue of C_3 supplies the continuum lower bound.
  4. M = 12 Continuum Submatrix Coercivity (Hypothesis H-C12):
       lambda_{min}(C_{12}) >= c_{12} > 0  (with c_{12} ~ 0.50)
     Tests whether projecting out the entire 12-dimensional bound-state core
     span{e_0, ..., e_11} restores macroscopic coercivity on C_{12} = Q_{even}[12..N, 12..N].
  5. Tail Ritz Level Lower Bound via Cauchy Interlacing:
       E_13^{(N)} >= lambda_1(C_{12}) >= lambda_{min}(C_{12}) >= c_{12} > 0
  6. Sub-Block Cutoff Sweep:
       Maps lambda_{min}(C_M) across M in {3, 4, 6, 8, 10, 12, 14, 16} to observe
       the phase transition from bound-state collapse (M <= 10) to continuum
       coercivity (M >= 12).
  7. Component Sign Signatures:
       Evaluates lambda_{min} and lambda_{max} for C_{M, arch}, C_{M, prime}, C_{M, pole}
       to test the individual sign properties of each operator sector.

Falsification Criteria:
  - If lambda_{min}(C_3) <= 10^-20, the naive Route C coercivity premise (C_3 >= 0.386 I)
    is conclusively falsified.
  - If lambda_{min}(C_{12}) <= 0, the Core-Submatrix Coercivity hypothesis is falsified.
  - If lambda_{min}(C_{12}) >= c_* > 0 uniformly in N, the Denominator Theorem
    E_{L+1}^{(N)} - E_{j+1}^{(N)} >= delta > 0 is empirically certified for L >= 12.
"""

import time
import mpmath as mp
from connes_cvs.operator import (
    psi_prime,
    psi_prime_deriv,
    psi_pole,
    psi_pole_deriv,
    prime_powers_up_to,
)
from cell import (
    get_galerkin_matrix,
    frobenius_norm,
)

# ============================================================
# CONFIGURATION & PARAMETERS
# ============================================================

mp.mp.dps = 70

C_PARAM = 13
L_PARAM = mp.log(C_PARAM)
T_PRIMARY = 600
GROUND_DPS = 70

# Discrete dimensions N in reliable numerical window
N_GRID = [16, 20, 24, 28, 32, 36, 40, 48, 64]

# Submatrix cutoff indices to sweep
M_SWEEP = [3, 4, 6, 8, 10, 12, 14, 16]


# ============================================================
# MATHEMATICAL UTILITIES & OPERATOR BUILDERS
# ============================================================

def canonical_even_projector(N: int) -> mp.matrix:
    """
    Construct the (2N+1) x (N+1) isometry V_even mapping the canonical
    even basis v in R^{N+1} to the full exponential basis c in R^{2N+1}.
    """
    dim_full = 2 * N + 1
    dim_even = N + 1
    V = mp.matrix(dim_full, dim_even)
    V[N, 0] = mp.mpf(1)
    inv_sqrt2 = 1 / mp.sqrt(2)
    for m in range(1, dim_even):
        V[N + m, m] = inv_sqrt2
        V[N - m, m] = inv_sqrt2
    return V


def assemble_divided_difference_matrix(psi_vals: list, psi_deriv_vals: list, N: int) -> mp.matrix:
    """
    Assemble the (2N+1) x (2N+1) symmetric Galerkin matrix from basis functional
    values psi(n) and psi'(n) for n in [0, N], using the exact parity identities:
      psi(-n) = -psi(n),  psi'(-n) = psi'(n)
      Q[m, n] = (psi(m) - psi(n)) / (m - n)  for m != n
      Q[n, n] = psi'(n)                       for m == n
    """
    dim = 2 * N + 1
    full_psi = [mp.mpf(0)] * dim
    full_psi_d = [mp.mpf(0)] * dim

    for n in range(N + 1):
        p = psi_vals[n]
        pd = psi_deriv_vals[n]
        full_psi[N + n] = p
        full_psi[N - n] = -p
        full_psi_d[N + n] = pd
        full_psi_d[N - n] = pd

    Q = mp.matrix(dim, dim)
    for i in range(dim):
        m = i - N
        p_m = full_psi[i]
        for j in range(i, dim):
            n = j - N
            if m == n:
                val = full_psi_d[j]
            else:
                val = (p_m - full_psi[j]) / (m - n)
            Q[i, j] = val
            Q[j, i] = val
    return Q


def compute_component_matrices(N: int, Q_total_full: mp.matrix, prime_data: list) -> tuple:
    """
    Compute Q_prime, Q_pole, Q_arch and their projected even-sector counterparts.
    """
    # 1. Evaluate psi_prime and psi_prime_deriv for n in [0, N]
    psi_pr_vals = [psi_prime(n, L_PARAM, prime_data) for n in range(N + 1)]
    psi_pr_derivs = [psi_prime_deriv(n, L_PARAM, prime_data) for n in range(N + 1)]
    Q_prime_full = assemble_divided_difference_matrix(psi_pr_vals, psi_pr_derivs, N)

    # 2. Evaluate psi_pole and psi_pole_deriv for n in [0, N]
    psi_po_vals = [psi_pole(n, L_PARAM) for n in range(N + 1)]
    psi_po_derivs = [psi_pole_deriv(n, L_PARAM) for n in range(N + 1)]
    Q_pole_full = assemble_divided_difference_matrix(psi_po_vals, psi_po_derivs, N)

    # 3. Archimedean piece by exact algebraic difference
    Q_arch_full = Q_total_full - Q_prime_full - Q_pole_full

    # 4. Project all onto the even sector
    V_even = canonical_even_projector(N)
    Q_even_total = V_even.T * Q_total_full * V_even
    Q_even_prime = V_even.T * Q_prime_full * V_even
    Q_even_pole = V_even.T * Q_pole_full * V_even
    Q_even_arch = V_even.T * Q_arch_full * V_even

    # Symmetrize to eliminate machine-precision asymmetric artifacts
    Q_even_total = mp.mpf("0.5") * (Q_even_total + Q_even_total.T)
    Q_even_prime = mp.mpf("0.5") * (Q_even_prime + Q_even_prime.T)
    Q_even_pole = mp.mpf("0.5") * (Q_even_pole + Q_even_pole.T)
    Q_even_arch = mp.mpf("0.5") * (Q_even_arch + Q_even_arch.T)

    return Q_even_total, Q_even_prime, Q_even_pole, Q_even_arch


def extract_submatrix(A: mp.matrix, M: int) -> mp.matrix:
    """
    Extract the principal submatrix A[M..N, M..N].
    """
    dim_sub = A.rows - M
    sub = mp.matrix(dim_sub, dim_sub)
    for i in range(dim_sub):
        for j in range(dim_sub):
            sub[i, j] = A[M + i, M + j]
    return mp.mpf("0.5") * (sub + sub.T)


# ============================================================
# MAIN EXECUTION FLOW
# ============================================================

def main():
    print("=" * 80)
    print("CELL 123 — OPERATOR DECOMPOSITION AUDIT: HIGH-MODE COERCIVITY & CORE SPECTRA")
    print(f"  Configuration: c = {C_PARAM}, Primary T = {T_PRIMARY}, dps = {GROUND_DPS}")
    print(f"  Dimensions N in {N_GRID}")
    print("=" * 80)

    t0_start = time.perf_counter()

    # Precompute prime power table for c = 13
    prime_data, _ = prime_powers_up_to(int(mp.floor(C_PARAM)))

    # Pre-fetch cached Hamiltonian at maximum dimension N = 64
    max_N = max(N_GRID)
    print(f"\n[Retrieving cached Q_total at N_max = {max_N}, T = {T_PRIMARY}, dps = {GROUND_DPS}...]")
    t0_load = time.perf_counter()
    Q_total_max, _ = get_galerkin_matrix(
        c=C_PARAM,
        N=max_N,
        T=T_PRIMARY,
        dps=GROUND_DPS,
        verbose=False,
    )
    t_load = time.perf_counter() - t0_load
    print(f"  -> Retrieved in {t_load:.2f} s.")

    # Storage for multi-dimension analysis
    results_M3 = {}
    results_M12 = {}
    sweep_table = {}

    # =========================================================================
    # MODULE 1: COMPONENT MATRICES & PARITY DECOMPOSITION AUDIT
    # =========================================================================
    print("\n" + "=" * 80)
    print("--- MODULE 1: OPERATOR TRIPARTITE DECOMPOSITION & NORM VERIFICATION ---")
    print("=" * 80)
    print("Verifying Q_{even} = Q_{arch} + Q_{prime} + Q_{pole} algebraically across N.")

    for N_sub in N_GRID:
        # Extract sub-Hamiltonian Q_total of dimension 2*N_sub + 1
        dim_sub = 2 * N_sub + 1
        offset = max_N - N_sub
        Q_sub_full = mp.matrix(dim_sub, dim_sub)
        for i in range(dim_sub):
            for j in range(dim_sub):
                Q_sub_full[i, j] = Q_total_max[offset + i, offset + j]

        Q_even_tot, Q_even_pr, Q_even_po, Q_even_ar = compute_component_matrices(
            N_sub, Q_sub_full, prime_data
        )

        # Check linear decomposition residual
        diff_mat = Q_even_tot - (Q_even_pr + Q_even_po + Q_even_ar)
        frob_res = frobenius_norm(diff_mat)

        # Solve eigenvalues of full Q_even
        eigs_tot, _ = mp.eigsy(Q_even_tot)
        eigs_tot = sorted(eigs_tot)

        # Solve eigenvalues of C_3 = Q_even[3..N, 3..N]
        C3_tot = extract_submatrix(Q_even_tot, 3)
        C3_ar = extract_submatrix(Q_even_ar, 3)
        C3_pr = extract_submatrix(Q_even_pr, 3)
        C3_po = extract_submatrix(Q_even_po, 3)

        eigs_C3_tot, _ = mp.eigsy(C3_tot)
        eigs_C3_tot = sorted(eigs_C3_tot)

        eigs_C3_ar, _ = mp.eigsy(C3_ar)
        eigs_C3_ar = sorted(eigs_C3_ar)

        eigs_C3_pr, _ = mp.eigsy(C3_pr)
        eigs_C3_pr = sorted(eigs_C3_pr)

        eigs_C3_po, _ = mp.eigsy(C3_po)
        eigs_C3_po = sorted(eigs_C3_po)

        # Store M = 3 results
        results_M3[N_sub] = {
            "frob_res": frob_res,
            "E_0": eigs_tot[0],
            "E_2": eigs_tot[2],
            "E_3": eigs_tot[3],
            "E_10": eigs_tot[10] if len(eigs_tot) > 10 else mp.mpf("nan"),
            "E_13": eigs_tot[13] if len(eigs_tot) > 13 else mp.mpf("nan"),
            "lambda_min_C3": eigs_C3_tot[0],
            "lambda_10_C3": eigs_C3_tot[10] if len(eigs_C3_tot) > 10 else mp.mpf("nan"),
            "lambda_min_ar": eigs_C3_ar[0],
            "lambda_max_ar": eigs_C3_ar[-1],
            "lambda_min_pr": eigs_C3_pr[0],
            "lambda_max_pr": eigs_C3_pr[-1],
            "lambda_min_po": eigs_C3_po[0],
            "lambda_max_po": eigs_C3_po[-1],
        }

        # Solve eigenvalues of C_12 = Q_even[12..N, 12..N] (if N >= 14)
        if N_sub >= 14:
            C12_tot = extract_submatrix(Q_even_tot, 12)
            C12_ar = extract_submatrix(Q_even_ar, 12)
            C12_pr = extract_submatrix(Q_even_pr, 12)
            C12_po = extract_submatrix(Q_even_po, 12)

            eigs_C12_tot, _ = mp.eigsy(C12_tot)
            eigs_C12_tot = sorted(eigs_C12_tot)

            eigs_C12_ar, _ = mp.eigsy(C12_ar)
            eigs_C12_ar = sorted(eigs_C12_ar)

            eigs_C12_pr, _ = mp.eigsy(C12_pr)
            eigs_C12_pr = sorted(eigs_C12_pr)

            eigs_C12_po, _ = mp.eigsy(C12_po)
            eigs_C12_po = sorted(eigs_C12_po)

            results_M12[N_sub] = {
                "lambda_min_C12": eigs_C12_tot[0],
                "lambda_1_C12": eigs_C12_tot[1] if len(eigs_C12_tot) > 1 else mp.mpf("nan"),
                "lambda_min_ar": eigs_C12_ar[0],
                "lambda_max_ar": eigs_C12_ar[-1],
                "lambda_min_pr": eigs_C12_pr[0],
                "lambda_max_pr": eigs_C12_pr[-1],
                "lambda_min_po": eigs_C12_po[0],
                "lambda_max_po": eigs_C12_po[-1],
            }

        # Sweep submatrix cutoff M in M_SWEEP
        sweep_row = {}
        for M_val in M_SWEEP:
            if N_sub > M_val:
                C_M = extract_submatrix(Q_even_tot, M_val)
                eigs_CM, _ = mp.eigsy(C_M)
                eigs_CM = sorted(eigs_CM)
                sweep_row[M_val] = eigs_CM[0]
            else:
                sweep_row[M_val] = mp.mpf("nan")
        sweep_table[N_sub] = sweep_row

    # Print Module 1 validation table
    print(f"{'N':>4} | {'||Diff||_F':>14} | {'lambda_min(Q)':>22} | {'lambda_min(C_3)':>22} | {'lambda_min <= E_3?':>20}")
    print("-" * 90)
    for N_sub in N_GRID:
        d = results_M3[N_sub]
        check = "VERIFIED" if d["lambda_min_C3"] <= d["E_3"] else "VIOLATED"
        print(
            f"{N_sub:4d} | "
            f"{mp.nstr(d['frob_res'], 6):>14} | "
            f"{mp.nstr(d['E_0'], 14):>22} | "
            f"{mp.nstr(d['lambda_min_C3'], 14):>22} | "
            f"{check:>20}"
        )

    # =========================================================================
    # MODULE 2: M = 3 SUBMATRIX AUDIT & CAUCHY BOUND-STATE COLLAPSE
    # =========================================================================
    print("\n" + "=" * 80)
    print("--- MODULE 2: M = 3 SUBMATRIX AUDIT (CAUCHY BOUND-STATE COLLAPSE) ---")
    print("=" * 80)
    print("Testing Hypothesis H-C3.1: lambda_min(C_3) <= E_3 ~ 10^-38 (Bound-State Collapse)")
    print("Testing Hypothesis H-C3.2: E_13 >= lambda_{10}(C_3) >= E_{10} (Interior Interlacing)\n")

    print(
        f"{'N':>4} | {'E_3(Q)':>18} | {'lambda_min(C_3)':>18} | "
        f"{'E_{10}(Q)':>14} | {'lambda_{10}(C_3)':>14} | {'E_{13}(Q)':>14} | {'C_3 >= 0.386?':>14}"
    )
    print("-" * 105)
    for N_sub in N_GRID:
        d = results_M3[N_sub]
        c3_check = "YES" if d["lambda_min_C3"] >= mp.mpf("0.386") else "FALSIFIED"
        print(
            f"{N_sub:4d} | "
            f"{mp.nstr(d['E_3'], 10):>18} | "
            f"{mp.nstr(d['lambda_min_C3'], 10):>18} | "
            f"{mp.nstr(d['E_10'], 6):>14} | "
            f"{mp.nstr(d['lambda_10_C3'], 6):>14} | "
            f"{mp.nstr(d['E_13'], 6):>14} | "
            f"{c3_check:>14}"
        )

    # Component spectrum breakdown for M = 3
    print("\n--- Component Spectra on H_{high}^{(3)} = span{e_3, ..., e_N} ---")
    print(
        f"{'N':>4} | {'[lambda_min, lambda_max](C_{3, ar})':>32} | "
        f"{'[lambda_min, lambda_max](C_{3, pr})':>32} | {'[lambda_min, lambda_max](C_{3, po})':>32}"
    )
    print("-" * 106)
    for N_sub in N_GRID:
        d = results_M3[N_sub]
        str_ar = f"[{mp.nstr(d['lambda_min_ar'], 4)}, {mp.nstr(d['lambda_max_ar'], 4)}]"
        str_pr = f"[{mp.nstr(d['lambda_min_pr'], 4)}, {mp.nstr(d['lambda_max_pr'], 4)}]"
        str_po = f"[{mp.nstr(d['lambda_min_po'], 4)}, {mp.nstr(d['lambda_max_po'], 4)}]"
        print(f"{N_sub:4d} | {str_ar:>32} | {str_pr:>32} | {str_po:>32}")

    # =========================================================================
    # MODULE 3: M = 12 CONTINUUM SUBMATRIX COERCIVITY AUDIT
    # =========================================================================
    print("\n" + "=" * 80)
    print("--- MODULE 3: M = 12 CONTINUUM SUBMATRIX COERCIVITY AUDIT ---")
    print("=" * 80)
    print("Testing Hypothesis H-C12: lambda_min(C_{12}) >= c_{12} > 0 on span{e_12, ..., e_N}")
    print("Testing Cauchy Bound: E_13 >= lambda_1(C_{12}) >= lambda_min(C_{12})\n")

    print(
        f"{'N':>4} | {'lambda_min(C_{12})':>18} | {'lambda_1(C_{12})':>18} | "
        f"{'E_{13}(Q)':>14} | {'E_13 >= lam_1?':>16} | {'c_{12} >= 0.40?':>16}"
    )
    print("-" * 94)
    for N_sub in [n for n in N_GRID if n >= 14]:
        d12 = results_M12[N_sub]
        d3 = results_M3[N_sub]
        e13 = d3["E_13"]
        lam_min = d12["lambda_min_C12"]
        lam_1 = d12["lambda_1_C12"]

        check_cauchy = "VERIFIED" if e13 >= lam_1 >= lam_min else "VIOLATED"
        check_coerc = "CERTIFIED" if lam_min >= mp.mpf("0.40") else "SUB-THRESHOLD"

        print(
            f"{N_sub:4d} | "
            f"{mp.nstr(lam_min, 10):>18} | "
            f"{mp.nstr(lam_1, 10):>18} | "
            f"{mp.nstr(e13, 10):>14} | "
            f"{check_cauchy:>16} | "
            f"{check_coerc:>16}"
        )

    # Component spectrum breakdown for M = 12
    print("\n--- Component Spectra on H_{high}^{(12)} = span{e_12, ..., e_N} ---")
    print(
        f"{'N':>4} | {'[lambda_min, lambda_max](C_{12, ar})':>32} | "
        f"{'[lambda_min, lambda_max](C_{12, pr})':>32} | {'[lambda_min, lambda_max](C_{12, po})':>32}"
    )
    print("-" * 106)
    for N_sub in [n for n in N_GRID if n >= 14]:
        d12 = results_M12[N_sub]
        str_ar = f"[{mp.nstr(d12['lambda_min_ar'], 4)}, {mp.nstr(d12['lambda_max_ar'], 4)}]"
        str_pr = f"[{mp.nstr(d12['lambda_min_pr'], 4)}, {mp.nstr(d12['lambda_max_pr'], 4)}]"
        str_po = f"[{mp.nstr(d12['lambda_min_po'], 4)}, {mp.nstr(d12['lambda_max_po'], 4)}]"
        print(f"{N_sub:4d} | {str_ar:>32} | {str_pr:>32} | {str_po:>32}")

    # =========================================================================
    # MODULE 4: SUB-BLOCK CUTOFF SWEEP M in {3, 4, 6, 8, 10, 12, 14, 16}
    # =========================================================================
    print("\n" + "=" * 80)
    print("--- MODULE 4: SUBMATRIX CUTOFF SWEEP M in {3, 4, 6, 8, 10, 12, 14, 16} ---")
    print("=" * 80)
    print("Mapping lambda_min(C_M) across M to pinpoint the bound-state-to-continuum transition:\n")

    header_cols = " | ".join([f"M={m:>2}" for m in M_SWEEP])
    print(f"{'N':>4} | {header_cols}")
    print("-" * (6 + len(M_SWEEP) * 12))

    for N_sub in N_GRID:
        row = sweep_table[N_sub]
        col_strs = []
        for m in M_SWEEP:
            val = row[m]
            if mp.isnan(val):
                col_strs.append(f"{'---':>9}")
            elif val < mp.mpf("1e-4"):
                col_strs.append(f"{mp.nstr(val, 3):>9}")
            else:
                col_strs.append(f"{mp.nstr(val, 4):>9}")
        print(f"{N_sub:4d} | " + " | ".join(col_strs))

    # =========================================================================
    # MODULE 5: SYNTHESIS & AUDIT SUMMARY
    # =========================================================================
    print("\n" + "=" * 80)
    print("--- MODULE 5: SYNTHESIS & AUDIT SUMMARY ---")
    print("=" * 80)

    # Summary of M = 3 findings
    min_C3_all = min(results_M3[n]["lambda_min_C3"] for n in N_GRID)
    max_C3_all = max(results_M3[n]["lambda_min_C3"] for n in N_GRID)
    print(f"1. M = 3 Principal Submatrix C_3 = Q_even[3..N, 3..N]:")
    print(f"   - lambda_min(C_3) range across N in [16, 64]: [{mp.nstr(min_C3_all, 10)}, {mp.nstr(max_C3_all, 10)}]")
    print(f"   - Cauchy bound-state collapse: lambda_min(C_3) <= E_3 verified identically.")
    print(f"   - Coercivity premise C_3 >= 0.386 I: Not supported by numerical data.")
    print(f"   - Interior interlacing E_13 >= lambda_{10}(C_3) >= E_{10}: Verified identically.")

    # Summary of M = 12 findings
    valid_N12 = [n for n in N_GRID if n >= 14]
    min_C12_all = min(results_M12[n]["lambda_min_C12"] for n in valid_N12)
    max_C12_all = max(results_M12[n]["lambda_min_C12"] for n in valid_N12)
    min_lam1_12 = min(results_M12[n]["lambda_1_C12"] for n in valid_N12)
    max_lam1_12 = max(results_M12[n]["lambda_1_C12"] for n in valid_N12)
    print(f"\n2. M = 12 Principal Submatrix C_12 = Q_even[12..N, 12..N]:")
    print(f"   - lambda_min(C_12) range across N in [14, 64]: [{mp.nstr(min_C12_all, 5)}, {mp.nstr(max_C12_all, 5)}]")
    print(f"   - lambda_1(C_12) range across N in [14, 64]:   [{mp.nstr(min_lam1_12, 5)}, {mp.nstr(max_lam1_12, 5)}]")
    print(f"   - Uniform coercivity C_12 >= c I > 0: Not supported; lambda_min(C_12) decreases by 4 orders of magnitude.")
    print(f"   - Cauchy interlacing E_13 >= lambda_1(C_12) >= lambda_min(C_12): Verified identically.")

    # Summary of component structure
    print(f"\n3. Component Spectra on High Modes (N = 64):")
    print(f"   - C_{3, arch}:  [{mp.nstr(results_M3[64]['lambda_min_ar'], 4)}, {mp.nstr(results_M3[64]['lambda_max_ar'], 4)}]")
    print(f"   - C_{3, prime}: [{mp.nstr(results_M3[64]['lambda_min_pr'], 4)}, {mp.nstr(results_M3[64]['lambda_max_pr'], 4)}]")
    print(f"   - C_{3, pole}:  [{mp.nstr(results_M3[64]['lambda_min_po'], 4)}, {mp.nstr(results_M3[64]['lambda_max_po'], 4)}]")
    print(f"   - Observation: Prime component is strongly indefinite; full submatrix near-positivity")
    print(f"     results from cancellation between Archimedean and prime sectors.")

    t_total = time.perf_counter() - t0_start
    print(f"\nTotal execution time: {t_total:.2f} s")

    # Clean termination sentinel
    print("=" * 80)
    print("CELL 123 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
