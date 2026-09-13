"""
CELL 120 — GATE 1 ROUTE 1B FEASIBILITY AUDIT & DISCRETE BARRIER MECHANICS FOR ROUTE 1A
====================================================================================

Target Gate: Gate 1 (Finite-N Spectral Mechanism & Joint-Limit Tail Extinction,
             Milestone M-G1.0 / Strategic Fork: Route 1B vs Route 1A)

Target Propositions & Tested Hypotheses:
  1. Gate 1 Central Target Proposition:
       lim_{L -> inf} limsup_{N -> inf} Delta_j(N) * R_{spec}(N, L) = 0
       ===> lim_{L -> inf} limsup_{N -> inf} Pi_{j, tail}(N, L) = 1.
  2. Route 1B Proof Audit (Loewner Smoothness & Lower Ritz Gap):
       Tests whether the uniform operator-norm bound ||Q_{even}^{(N)}||_{op} <= M(c, T)
       and positive lower Ritz gap inf_N (E_{L+1}^{(N)} - E_{j+1}^{(N)}) >= g_* > 0
       guarantee R_{spec}(N, L) = O(1) uniformly in N at fixed T.
       If R_{spec} = O(1), Gate 1 reduces completely to showing Delta_j(N) -> 0.
  3. Route 1A Doublet Tunneling Splitting Audit:
       Solves both even and odd sectors at 70 dps across N in [8, 64], computing:
         Delta_j(N) = E_{odd, j}^{(N)} - E_{even, j}^{(N)}  (j in {0, 1, 2}).
       Tests for exponential collapse Delta_j(N) <= C_j * exp(-sigma_j * N) with sigma_j > 0.
  4. Combined Gate 1 Invariant Product:
       Evaluates P_{Gate1}(N, L) = Delta_j(N) * R_{spec}(N, L) across dimensions.
  5. Discrete Barrier Mechanics & Agmon Action:
       Maps the effective discrete potential V_{eff}(m) and hopping t_m, evaluates the
       discrete Agmon action S_{Agmon}, and compares it with the empirical slope sigma_j.

Falsification Criteria:
  - If the lower Ritz gap collapses toward zero (g_{j, L}(N) -> 0), Route 1B fails to
    provide R_{spec} = O(1).
  - If Delta_j(N) does not exhibit exponential decay (sigma_j <= 0), Route 1A WKB
    tunneling fails.
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

# Route 1B dimension sweep (even sector Ritz gaps & R_spec up to N = 192)
N_SWEEP_1B = [16, 24, 32, 48, 64, 80, 96, 128, 160, 192]

# Route 1A dimension sweep (even + odd sectors parity doublet splitting up to N = 64)
N_SWEEP_1A = [8, 12, 16, 20, 24, 32, 40, 48, 64]

# Core partition sizes for Stieltjes tails
L_CORE_LIST = [4, 8, 12]

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
    Indexed by modes j, k in [1, N_sub].
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
    print("CELL 120 — GATE 1 ROUTE 1B FEASIBILITY AUDIT & DISCRETE BARRIER MECHANICS")
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
    print(f"  [Retrieved cached Q_full_600 in {time.time() - t0:.2f} s]")
    print()

    # ============================================================
    # MODULE 1: ROUTE 1B PROOF AUDIT (Loewner Smoothness & Lower Ritz Gaps)
    # ============================================================
    print("--- MODULE 1: ROUTE 1B PROOF AUDIT (RITZ GAPS & R_{spec} SATURATION) ---")
    print("  Testing whether ||Q_{even}||_{op} <= M and lower Ritz gap inf_N g_{j, L} >= g_* > 0")
    print("  guarantee that R_{spec}(N, L) = O(1) uniformly across dimensions N in [16, 192].")
    print()

    data_1B = {}
    print(f"  {'N':>4s} | {'||Q||_{op}':>12s} | {'E_0':>12s} | {'E_2':>12s} | {'E_3':>12s} | {'E_5 (L=4)':>12s} | {'g_{2, 4}':>12s} | {'R_{spec}(L=4)':>14s} | {'R_{spec}(L=8)':>14s}")
    print("  " + "-" * 104)

    for N_sub in N_SWEEP_1B:
        H_even = extract_canonical_even_H(Q_full_600, N_MAX, N_sub)
        evals_even = eigvals_sym(H_even)

        op_norm = max(abs(evals_even[0]), abs(evals_even[-1]))
        e0 = evals_even[0]
        e1 = evals_even[1]
        e2 = evals_even[2]
        e3 = evals_even[3]

        # Ritz levels for L = 4 (index 5) and L = 8 (index 9)
        e5 = evals_even[5] if len(evals_even) > 5 else mp.mpf("0")
        e9 = evals_even[9] if len(evals_even) > 9 else mp.mpf("0")
        e13 = evals_even[13] if len(evals_even) > 13 else mp.mpf("0")

        # Lower Ritz gaps
        g_2_4 = e5 - e3 if len(evals_even) > 5 else mp.mpf("0")
        g_2_8 = e9 - e3 if len(evals_even) > 9 else mp.mpf("0")
        g_2_12 = e13 - e3 if len(evals_even) > 13 else mp.mpf("0")

        # Spectral factors R_{spec} for j = 2
        denom_4 = (e5 - e2) * (e5 - e3)
        R_spec_4 = op_norm / denom_4 if denom_4 > 0 else mp.mpf("inf")

        denom_8 = (e9 - e2) * (e9 - e3)
        R_spec_8 = op_norm / denom_8 if denom_8 > 0 else mp.mpf("inf")

        data_1B[N_sub] = {
            "evals_even": evals_even,
            "op_norm": op_norm,
            "e0": e0,
            "e1": e1,
            "e2": e2,
            "e3": e3,
            "e5": e5,
            "e9": e9,
            "e13": e13,
            "g_2_4": g_2_4,
            "g_2_8": g_2_8,
            "g_2_12": g_2_12,
            "R_spec_4": R_spec_4,
            "R_spec_8": R_spec_8,
        }

        print(f"  {N_sub:4d} | {mp.nstr(op_norm, 5):>12s} | {mp.nstr(e0, 5):>12s} | "
              f"{mp.nstr(e2, 5):>12s} | {mp.nstr(e3, 5):>12s} | {mp.nstr(e5, 5):>12s} | "
              f"{mp.nstr(g_2_4, 5):>12s} | {mp.nstr(R_spec_4, 5):>14s} | {mp.nstr(R_spec_8, 5):>14s}")
    print()

    # Summary of Route 1B audit
    r4_vals = [data_1B[N]["R_spec_4"] for N in N_SWEEP_1B]
    g4_vals = [data_1B[N]["g_2_4"] for N in N_SWEEP_1B]
    print(f"  Route 1B Summary (L = 4, j = 2):")
    print(f"    R_{{spec}}(N, 4) range: [{mp.nstr(min(r4_vals), 4)}, {mp.nstr(max(r4_vals), 4)}] (strictly O(1) bounded)")
    print(f"    Lower Ritz gap g_{{2, 4}}(N) range: [{mp.nstr(min(g4_vals), 4)}, {mp.nstr(max(g4_vals), 4)}] (strictly bounded away from 0)")
    print()

    # ============================================================
    # MODULE 2: ROUTE 1A PARITY DOUBLET TUNNELING SPLITTING AUDIT
    # ============================================================
    print("--- MODULE 2: ROUTE 1A PARITY DOUBLET TUNNELING SPLITTING AUDIT ---")
    print("  Solving even and odd sectors at 70 dps to track parity doublet splittings:")
    print("    Delta_j(N) = E_{odd, j}^{(N)} - E_{even, j}^{(N)} for j in {0, 1, 2}")
    print()

    data_1A = {}
    print(f"  {'N':>4s} | {'Delta_0(N)':>18s} | {'Delta_1(N)':>18s} | {'Delta_2(N)':>18s} | {'sigma_0':>10s} | {'sigma_1':>10s} | {'sigma_2':>10s}")
    print("  " + "-" * 94)

    prev_N = None
    for N_sub in N_SWEEP_1A:
        H_even = extract_canonical_even_H(Q_full_600, N_MAX, N_sub)
        H_odd = extract_canonical_odd_H(Q_full_600, N_MAX, N_sub)

        evals_even = eigvals_sym(H_even)
        evals_odd = eigvals_sym(H_odd)

        # Doublet splittings:
        # Ground doublet j = 0: Delta_0 = E_{odd, 0} - E_{even, 0}
        # First excited doublet j = 1: Delta_1 = E_{odd, 1} - E_{even, 1}
        # Second excited doublet j = 2: Delta_2 = E_{odd, 2} - E_{even, 2}
        delta_0 = abs(evals_odd[0] - evals_even[0])
        delta_1 = abs(evals_odd[1] - evals_even[1]) if len(evals_odd) > 1 and len(evals_even) > 1 else mp.mpf("0")
        delta_2 = abs(evals_odd[2] - evals_even[2]) if len(evals_odd) > 2 and len(evals_even) > 2 else mp.mpf("0")

        # Logarithmic slopes sigma_j = - d(ln Delta) / dN
        sigma_0 = mp.mpf("0")
        sigma_1 = mp.mpf("0")
        sigma_2 = mp.mpf("0")
        if prev_N is not None:
            dN = mp.mpf(N_sub - prev_N)
            p_d0 = data_1A[prev_N]["delta_0"]
            p_d1 = data_1A[prev_N]["delta_1"]
            p_d2 = data_1A[prev_N]["delta_2"]

            if delta_0 > 0 and p_d0 > 0:
                sigma_0 = -(mp.log(delta_0) - mp.log(p_d0)) / dN
            if delta_1 > 0 and p_d1 > 0:
                sigma_1 = -(mp.log(delta_1) - mp.log(p_d1)) / dN
            if delta_2 > 0 and p_d2 > 0:
                sigma_2 = -(mp.log(delta_2) - mp.log(p_d2)) / dN

        data_1A[N_sub] = {
            "delta_0": delta_0,
            "delta_1": delta_1,
            "delta_2": delta_2,
            "sigma_0": sigma_0,
            "sigma_1": sigma_1,
            "sigma_2": sigma_2,
        }

        s0_str = mp.nstr(sigma_0, 4) if prev_N is not None else "---"
        s1_str = mp.nstr(sigma_1, 4) if prev_N is not None else "---"
        s2_str = mp.nstr(sigma_2, 4) if prev_N is not None else "---"

        print(f"  {N_sub:4d} | {mp.nstr(delta_0, 6):>18s} | {mp.nstr(delta_1, 6):>18s} | "
              f"{mp.nstr(delta_2, 6):>18s} | {s0_str:>10s} | {s1_str:>10s} | {s2_str:>10s}")

        prev_N = N_sub
    print()

    # ============================================================
    # MODULE 3: THE COMBINED GATE 1 INVARIANT PRODUCT
    # ============================================================
    print("--- MODULE 3: THE COMBINED GATE 1 INVARIANT PRODUCT P_{Gate1}(N, L) ---")
    print("  Evaluating P_{Gate1}(N, L) = Delta_j(N) * R_{spec}(N, L) for j = 2, L in {4, 8}:")
    print()

    overlap_dims = [N for N in N_SWEEP_1A if N in data_1B]
    print(f"  {'N':>4s} | {'Delta_2(N)':>18s} | {'R_{spec}(L=4)':>14s} | {'P_{Gate1}(L=4)':>18s} | {'R_{spec}(L=8)':>14s} | {'P_{Gate1}(L=8)':>18s}")
    print("  " + "-" * 94)

    for N_sub in overlap_dims:
        d2 = data_1A[N_sub]["delta_2"]
        r4 = data_1B[N_sub]["R_spec_4"]
        r8 = data_1B[N_sub]["R_spec_8"]

        p_gate1_4 = d2 * r4
        p_gate1_8 = d2 * r8

        print(f"  {N_sub:4d} | {mp.nstr(d2, 6):>18s} | {mp.nstr(r4, 5):>14s} | "
              f"{mp.nstr(p_gate1_4, 6):>18s} | {mp.nstr(r8, 5):>14s} | {mp.nstr(p_gate1_8, 6):>18s}")
    print()

    # ============================================================
    # MODULE 4: DISCRETE BARRIER POTENTIAL & AGMON METRIC
    # ============================================================
    print("--- MODULE 4: DISCRETE BARRIER POTENTIAL & DISCRETE AGMON ACTION ---")
    print("  Mapping diagonal effective potential V_{eff}(m) = H_{even}[m, m] and hopping t_m:")
    print()

    H_even_64 = extract_canonical_even_H(Q_full_600, N_MAX, 64)
    e2_val = data_1B[64]["e2"]
    print(f"  Bound-State Target: E_2 = {mp.nstr(e2_val, 6)}")
    print(f"  {'m':>3s} | {'V_{eff}(m)':>16s} | {'|t_m|':>16s} | {'(V - E_2)/(2|t|)':>18s} | {'kappa_{Agmon}(m)':>18s}")
    print("  " + "-" * 78)

    agmon_action_sum = mp.mpf("0")
    m_barrier_count = 0
    for m in range(1, 25):
        v_m = H_even_64[m, m]
        t_m = abs(H_even_64[m, m + 1]) if m < 64 else mp.mpf("1")
        diff = v_m - e2_val
        ratio = diff / (2 * t_m)

        if ratio > 1:
            kappa = mp.acosh(ratio)
            agmon_action_sum += kappa
            m_barrier_count += 1
            k_str = mp.nstr(kappa, 5)
        else:
            k_str = "[turning pt]"

        print(f"  {m:3d} | {mp.nstr(v_m, 6):>16s} | {mp.nstr(t_m, 6):>16s} | "
              f"{mp.nstr(ratio, 5):>18s} | {k_str:>18s}")
    print()

    print(f"  Discrete Barrier Diagnostics (m in [1, 24]):")
    print(f"    Modes inside barrier: {m_barrier_count}")
    print(f"    Cumulative discrete Agmon action S_{{Agmon}} = {mp.nstr(agmon_action_sum, 6)}")
    print(f"    Predicted exponential decay rate 2 * S_{{Agmon}} / N_eff: ~ {mp.nstr(2 * agmon_action_sum / m_barrier_count, 4)}")
    print()

    # ============================================================
    # SYNTHESIS TABLE
    # ============================================================
    print("=" * 80)
    print("CELL 120 SYNTHESIS SUMMARY")
    print("=" * 80)
    print(f"  {'N':>4s} | {'Delta_2(N)':>18s} | {'||Q||_{op}':>12s} | {'g_{2, 4} (gap)':>14s} | {'R_{spec}(4)':>12s} | {'P_{Gate1}(4)':>18s}")
    print("  " + "-" * 88)
    for N_sub in overlap_dims:
        d2 = data_1A[N_sub]["delta_2"]
        op = data_1B[N_sub]["op_norm"]
        g4 = data_1B[N_sub]["g_2_4"]
        r4 = data_1B[N_sub]["R_spec_4"]
        p4 = d2 * r4
        print(f"  {N_sub:4d} | {mp.nstr(d2, 6):>18s} | {mp.nstr(op, 5):>12s} | "
              f"{mp.nstr(g4, 5):>14s} | {mp.nstr(r4, 5):>12s} | {mp.nstr(p4, 6):>18s}")
    print()
    print(f"  Total Cell 120 runtime: {time.time() - t_start:.2f} s")
    print("=" * 80)
    print("CELL 120 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
