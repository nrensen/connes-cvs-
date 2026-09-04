"""
CELL 54 — ANALYTIC ANATOMY OF THE FIRST-JET CANCELLATION SCALE D_0 / D_1,
SOBOLEV TRACE BOUNDS, AND EXPONENTIAL FACTOR CANCELLATION

Following the findings of Cell 53:
1. Cell 53 found a surprisingly good cross-N collapse of the normalized heat profile:
       Theta_N(theta) = H_N(theta u_1) / D_0
   when heat time u is scaled by the first-jet scale u_1 = D_0 / D_1
   (matching within ~1.5% at theta = 1.0 across the sampled range N in {8, ..., 24}).
2. The dimensionless ratio:
       s_N = (kappa * N)^2 * (|D_0| / |D_1|) = u_1 / u_edge
   drifts slowly from 9.19e-3 (N=8) to 6.65e-3 (N=24). Since s_N is drifting,
   we do not yet know whether this prefactor has a non-zero continuum limit s_infty > 0
   or vanishes as N -> infinity.
3. Because D_0 = T_N(0) and D_1 = T_N''(0) both vanish rapidly, the central
   analytic question proposed by the reviewer is:
       "What determines u_1 = D_0 / D_1?"

Cell 54 executes a four-part mathematical dissection:

PART 1: MODE-BY-MODE ANATOMY & SIGNED/DESTRUCTIVE CANCELLATION MECHANICS
  Decomposes D_0 = v_0 + sqrt(2) sum v_m and D_1 = -sqrt(2) kappa^2 sum m^2 v_m
  into positive and negative sub-sums S^+, S^-.
  Measures the degree of signed cancellation:
      epsilon_0(N) = |D_0| / (v_0 + sqrt(2) S_0^+)
      epsilon_1(N) = |D_1| / (sqrt(2) kappa^2 S_1^+)
  and tracks low-frequency bulk (m <= N/2) vs. edge boundary layer (m > N/2).

PART 2: SOBOLEV NORMS, TRACE BOUNDS, AND QUADRATIC FORM
  Computes the L^2 norms of derivatives:
      ||T_v||_{L^2}^2 = v_0^2 + sum v_m^2 = 1.0
      ||T'_v||_{L^2}^2 = kappa^2 sum m^2 v_m^2 (Dirichlet kinetic energy)
      ||T''_v||_{L^2}^2 = kappa^4 sum m^4 v_m^2 (curvature energy)
  Compares D_1 against the Cauchy-Schwarz trace bound:
      |D_1| <= sqrt(2) * kappa * ||T'_v|| * sqrt(sum_{m=1}^N m^2).
  Conceptual note: Because ||T_N||_2 = 1 while T_N(0) = D_0 -> 0, the endpoint value
  is not controlled by the global L^2 norm in any useful asymptotic sense. Similarly,
  trace inequalities give upper bounds but cannot explain the extraordinary cancellation;
  endpoint suppression is encoded in the specific ground-state mode vector.

PART 3: ASYMPTOTIC EXPONENTIAL RATES & DIRECT DIFFERENCE Delta_N
  Tests whether D_0 and D_1 share a common leading exponential decay rate:
      alpha_0(N) = -log|D_0(N)| / N   vs.   alpha_1(N) = -log|D_1(N)| / N.
  If alpha_0 = alpha_1, the leading exponential factor cancels in D_0 / D_1,
  leaving subexponential algebraic factors (e.g. D_0 / D_1 ~ N^{p_0 - p_1}).
  At this stage, this tests an empirical exponential-rate hypothesis.
  Directly tracks the exponential difference:
      Delta_N = -log|D_0| + log|D_1|
  such that log(s_N) = 2*log(kappa*N) - Delta_N, isolating the subexponential remainder.

PART 4: ASYMPTOTIC SCALING DIAGNOSTICS FOR s_N = (kappa*N)^2 * (D_0 / D_1)
  Connects Part 3 and Part 4 via log(s_N) = 2*log(kappa*N) - Delta_N, and applies
  three exploratory scaling diagnostics:
      - Diagnostic A (Quadratic extrapolation diagnostic on N in {16, 20, 24}):
            s_N = s_infty + c_1 / N + c_2 / N^2 (3-point interpolation)
      - Diagnostic B (Local power-law diagnostic on N in {20, 24}):
            s_N ~ A * N^{-p}
      - Diagnostic C (Local logarithmic diagnostic on N in {20, 24}):
            s_N ~ A / (log N)^p
  These exploratory diagnostics provide qualitative indications to guide higher-N runs.
"""

from __future__ import annotations

import mpmath as mp

from cell import get_ground_state


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

mp.mp.dps = 50

C_PARAM = 13
L_PARAM = mp.log(C_PARAM)
KAPPA = 2 * mp.pi / L_PARAM
T_QUAD = 400
DPS_RUN = 50

N_LIST = [8, 12, 16, 20, 24]


# ---------------------------------------------------------------------------
# Evaluators and Derivatives
# ---------------------------------------------------------------------------

def compute_jets_and_subsums(v, kappa):
    """
    Computes D_0, D_1, positive/negative subsums, and Sobolev norms.
    """
    N = len(v) - 1

    # D_0 decomposition
    v0 = v[0]
    s0_pos = mp.mpf(0)
    s0_neg = mp.mpf(0)
    bulk_s0 = mp.mpf(0)
    edge_s0 = mp.mpf(0)

    # D_1 decomposition: D_1 = -sqrt(2)*kappa^2 sum_{m=1}^N m^2 v_m
    s1_pos = mp.mpf(0)
    s1_neg = mp.mpf(0)
    bulk_s1 = mp.mpf(0)
    edge_s1 = mp.mpf(0)

    # Sobolev norms
    l2_norm2 = v0 ** 2
    h1_semi2 = mp.mpf(0)
    h2_semi2 = mp.mpf(0)

    for m in range(1, len(v)):
        vm = v[m]
        l2_norm2 += vm ** 2
        am2 = (kappa * m) ** 2
        am4 = am2 ** 2

        h1_semi2 += am2 * (vm ** 2)
        h2_semi2 += am4 * (vm ** 2)

        # D_0 terms
        if vm >= 0:
            s0_pos += vm
        else:
            s0_neg += abs(vm)

        # D_1 terms: term is m^2 * vm
        m2_vm = (m ** 2) * vm
        if m2_vm >= 0:
            s1_pos += m2_vm
        else:
            s1_neg += abs(m2_vm)

        # Bulk vs edge (split at N/2)
        if m <= N // 2:
            bulk_s0 += vm
            bulk_s1 += m2_vm
        else:
            edge_s0 += vm
            edge_s1 += m2_vm

    d0 = v0 + mp.sqrt(2) * (s0_pos - s0_neg)
    sum_m2_vm = s1_pos - s1_neg
    d1 = -mp.sqrt(2) * (kappa ** 2) * sum_m2_vm

    # Cancellation factors
    denom_0 = v0 + mp.sqrt(2) * s0_pos
    eps_0 = abs(d0) / denom_0 if denom_0 > 0 else mp.mpf(0)

    denom_1 = mp.sqrt(2) * (kappa ** 2) * s1_pos
    eps_1 = abs(d1) / denom_1 if denom_1 > 0 else mp.mpf(0)

    return {
        "d0": d0,
        "d1": d1,
        "v0": v0,
        "s0_pos": s0_pos,
        "s0_neg": s0_neg,
        "bulk_s0": bulk_s0,
        "edge_s0": edge_s0,
        "eps_0": eps_0,
        "s1_pos": s1_pos,
        "s1_neg": s1_neg,
        "bulk_s1": bulk_s1,
        "edge_s1": edge_s1,
        "eps_1": eps_1,
        "h1_norm": mp.sqrt(h1_semi2),
        "h2_norm": mp.sqrt(h2_semi2),
        "l2_norm": mp.sqrt(l2_norm2),
    }


# ---------------------------------------------------------------------------
# MAIN INVESTIGATION
# ---------------------------------------------------------------------------

def main():
    print("=" * 80)
    print("CELL 54 — ANALYTIC ANATOMY OF D_0 / D_1 & THE CANCELLATION SCALE")
    print("=" * 80)
    print(f"Parameters: c = {C_PARAM}, L = ln(13) = {mp.nstr(L_PARAM, 10)}, kappa = {mp.nstr(KAPPA, 10)}")
    print(f"Working precision: {DPS_RUN} decimal digits")
    print()

    # Step 1: Load ground states via persistent cache (.cell_cache)
    data = {}
    print("Loading ground states via persistent cache (.cell_cache)...")
    for N in N_LIST:
        lam, v, meta = get_ground_state(c=C_PARAM, N=N, T=T_QUAD, dps=DPS_RUN, verbose=False)
        res = compute_jets_and_subsums(v, KAPPA)
        res["lam"] = lam
        data[N] = res
        hit_status = "cache hit" if meta.get("cache_hit") else "computed"
        sec_str = f" ({meta.get('total_seconds', 0):.2f}s)" if 'total_seconds' in meta else ""
        print(f"    N = {N:2d} [{hit_status}{sec_str}]: lambda_min = {mp.nstr(lam, 8)}, D_0 = {mp.nstr(res['d0'], 6)}, D_1 = {mp.nstr(res['d1'], 6)}")

    print()

    # -----------------------------------------------------------------------
    # Part 1: Mode-by-Mode Anatomy & Signed Cancellation Mechanics
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("1. MODE-BY-MODE ANATOMY & SIGNED/DESTRUCTIVE CANCELLATION MECHANICS")
    print("=" * 80)
    print("Decomposing D_0 and D_1 into positive/negative sums to quantify destructive interference:")
    print("    epsilon_0 = |D_0| / (v_0 + sqrt(2)*S_0^+),   epsilon_1 = |D_1| / (sqrt(2)*kappa^2*S_1^+)")
    print()
    print(f"{'N':>4}  {'D_0':>14}  {'Positive S_0^+':>16}  {'Negative S_0^-':>16}  {'Cancellation eps_0':>20}")
    print("-" * 76)
    for N in N_LIST:
        d = data[N]
        print(f"{N:>4d}  {mp.nstr(d['d0'], 6):>14}  {mp.nstr(d['s0_pos'], 6):>16}  {mp.nstr(d['s0_neg'], 6):>16}  {mp.nstr(d['eps_0'], 6):>20}")
    print("-" * 76)
    print()

    print(f"{'N':>4}  {'D_1':>14}  {'Positive S_1^+':>16}  {'Negative S_1^-':>16}  {'Cancellation eps_1':>20}")
    print("-" * 76)
    for N in N_LIST:
        d = data[N]
        print(f"{N:>4d}  {mp.nstr(d['d1'], 6):>14}  {mp.nstr(d['s1_pos'], 6):>16}  {mp.nstr(d['s1_neg'], 6):>16}  {mp.nstr(d['eps_1'], 6):>20}")
    print("-" * 76)
    print()

    print("Bulk (m <= N/2) vs. Edge (m > N/2) contributions to the Fourier sums:")
    print(f"{'N':>4}  {'Bulk Sum (D_0)':>16}  {'Edge Sum (D_0)':>16}  {'Bulk Sum (D_1)':>16}  {'Edge Sum (D_1)':>16}")
    print("-" * 74)
    for N in N_LIST:
        d = data[N]
        print(f"{N:>4d}  {mp.nstr(d['bulk_s0'], 6):>16}  {mp.nstr(d['edge_s0'], 6):>16}  {mp.nstr(d['bulk_s1'], 6):>16}  {mp.nstr(d['edge_s1'], 6):>16}")
    print("-" * 74)
    print()

    # -----------------------------------------------------------------------
    # Part 2: Sobolev Norms, Trace Bounds, and Kinetic Energy
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("2. SOBOLEV NORMS & CAUCHY-SCHWARZ TRACE BOUNDS")
    print("=" * 80)
    print("Dirichlet kinetic energy ||T'_v||, curvature ||T''_v||, and Cauchy-Schwarz trace bounds:")
    print("    |D_1| <= sqrt(2) * kappa * ||T'_v||_2 * sqrt(sum_{m=1}^N m^2)")
    print()
    print(f"{'N':>4}  {'||T_v||_{L^2}':>14}  {'||T_v\'||_{L^2}':>16}  {'||T_v\'\'||_{L^2}':>16}  {'CS Bound Ratio':>18}")
    print("-" * 74)
    for N in N_LIST:
        d = data[N]
        sum_m2 = mp.mpf(N * (N + 1) * (2 * N + 1)) / 6
        cs_bound = mp.sqrt(2) * KAPPA * d["h1_norm"] * mp.sqrt(sum_m2)
        cs_ratio = abs(d["d1"]) / cs_bound if cs_bound > 0 else mp.mpf(0)
        print(f"{N:>4d}  {mp.nstr(d['l2_norm'], 6):>14}  {mp.nstr(d['h1_norm'], 6):>16}  {mp.nstr(d['h2_norm'], 6):>16}  {mp.nstr(cs_ratio, 6):>18}")
    print("-" * 74)
    print("Note: The elementary Cauchy–Schwarz bound is extremely non-sharp for D_1.")
    print("This provides evidence that signed cancellation is active, but because the bound")
    print("itself grows with N (scaling as ~ N^{3/2}), its ratio to D_1 becoming tiny is an")
    print("expected consequence of D_1 -> 0 rather than a separate proof of mechanism.")
    print("Furthermore, because ||T_N||_{L^2} = 1 while T_N(0) = D_0 -> 0, the endpoint value")
    print("is not controlled by the global L^2 norm in any useful asymptotic sense.")
    print("Endpoint suppression is encoded in the specific ground-state mode vector rather")
    print("than being a generic Sobolev phenomenon.")
    print()

    # -----------------------------------------------------------------------
    # Part 3: Asymptotic Exponential Rates & Direct Difference Delta_N
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("3. ASYMPTOTIC EXPONENTIAL RATES & DIRECT DIFFERENCE Delta_N")
    print("=" * 80)
    print("Testing whether D_0 and D_1 share a common leading exponential decay rate:")
    print("    alpha_0(N) = -log|D_0(N)| / N   vs.   alpha_1(N) = -log|D_1(N)| / N")
    print("and tracking the direct difference Delta_N = -log|D_0| + log|D_1|:")
    print("    If alpha_0 = alpha_1, the leading exponential factor cancels in D_0 / D_1.")
    print()
    print(f"{'N':>4}  {'-log|D_0|':>14}  {'-log|D_1|':>14}  {'alpha_0(N)':>12}  {'alpha_1(N)':>12}  {'Delta_N':>14}  {'s_N':>14}")
    print("-" * 88)
    for N in N_LIST:
        d = data[N]
        log_d0 = -mp.log(abs(d["d0"]))
        log_d1 = -mp.log(abs(d["d1"]))
        a0 = log_d0 / N
        a1 = log_d1 / N
        delta_n = log_d0 - log_d1  # -log|D_0| - (-log|D_1|) = -log|D_0| + log|D_1|
        s_val = ((KAPPA * N) ** 2) * (abs(d["d0"]) / abs(d["d1"]))
        print(f"{N:>4d}  {mp.nstr(log_d0, 6):>14}  {mp.nstr(log_d1, 6):>14}  {mp.nstr(a0, 6):>12}  {mp.nstr(a1, 6):>12}  {mp.nstr(delta_n, 6):>14}  {mp.nstr(s_val, 6):>14}")
    print("-" * 88)
    print()

    print("Two-point consecutive logarithmic decay rates alpha_N = -[log|D(N)| - log|D(N-4)|] / 4:")
    print(f"{'Step N -> N+4':>16}  {'alpha_N(D_0)':>16}  {'alpha_N(D_1)':>16}  {'|alpha_0 - alpha_1|':>22}")
    print("-" * 60)
    for i in range(len(N_LIST) - 1):
        n_prev = N_LIST[i]
        n_curr = N_LIST[i + 1]
        dn = n_curr - n_prev
        rate_0 = (-mp.log(abs(data[n_curr]["d0"])) - (-mp.log(abs(data[n_prev]["d0"])))) / dn
        rate_1 = (-mp.log(abs(data[n_curr]["d1"])) - (-mp.log(abs(data[n_prev]["d1"])))) / dn
        diff_rate = abs(rate_0 - rate_1)
        print(f"{n_prev:>5d} -> {n_curr:<5d}  {mp.nstr(rate_0, 6):>16}  {mp.nstr(rate_1, 6):>16}  {mp.nstr(diff_rate, 6):>22}")
    print("-" * 60)
    print("Note: The finite-N effective exponential rates can be compared to test whether")
    print("D_0 and D_1 share a common leading exponential factor. The present N-range does not")
    print("establish asymptotic equality.")
    print()

    # -----------------------------------------------------------------------
    # Part 4: Asymptotic Scaling Diagnostics for s_N = (kappa*N)^2 * (D_0 / D_1)
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("4. SCALING DIAGNOSTICS FOR s_N = (kappa*N)^2 * (D_0 / D_1)")
    print("=" * 80)
    print("Connection between Part 3 and Part 4 via the subexponential decomposition:")
    print("    log(s_N) = 2*log(kappa*N) - [-log|D_1| + log|D_0|] = 2*log(kappa*N) - Delta_N")
    print("where Delta_N = -log|D_0| + log|D_1| removes the leading exponential order.")
    print()
    print(f"{'N':>4}  {'s_N':>14}  {'Delta_N':>14}  {'2*log(kappa*N)':>16}  {'log(s_N)':>14}  {'Identity Check':>16}")
    print("-" * 82)
    s_vals = {}
    for N in N_LIST:
        d = data[N]
        log_d0 = -mp.log(abs(d["d0"]))
        log_d1 = -mp.log(abs(d["d1"]))
        delta_n = log_d0 - log_d1
        s_val = ((KAPPA * N) ** 2) * (abs(d["d0"]) / abs(d["d1"]))
        s_vals[N] = s_val
        two_log_kn = 2 * mp.log(KAPPA * N)
        log_s = mp.log(s_val)
        decomp_diff = abs(log_s - (two_log_kn - delta_n))
        print(f"{N:>4d}  {mp.nstr(s_val, 8):>14}  {mp.nstr(delta_n, 8):>14}  {mp.nstr(two_log_kn, 8):>16}  {mp.nstr(log_s, 8):>14}  {mp.nstr(decomp_diff, 4):>16}")
    print("-" * 82)
    print()

    # Diagnostic A: Quadratic extrapolation diagnostic on N in {16, 20, 24}
    print("Diagnostic A: Quadratic extrapolation diagnostic on N in {16, 20, 24}:")
    print("    Assumes s_N = s_infty + c_1 / N + c_2 / N^2 (exploratory 3-point interpolation):")
    n_pts = [16, 20, 24]
    M = mp.matrix(3, 3)
    Y = mp.matrix(3, 1)
    for row, n_val in enumerate(n_pts):
        M[row, 0] = mp.mpf(1)
        M[row, 1] = mp.mpf(1) / n_val
        M[row, 2] = mp.mpf(1) / (n_val ** 2)
        Y[row, 0] = s_vals[n_val]

    sol = mp.lu_solve(M, Y)
    s_inf_poly = sol[0, 0]
    c1_poly = sol[1, 0]
    c2_poly = sol[2, 0]

    print(f"        Extrapolated s_infty = {mp.nstr(s_inf_poly, 6)}")
    print(f"        Fitted c_1           = {mp.nstr(c1_poly, 6)}")
    print(f"        Fitted c_2           = {mp.nstr(c2_poly, 6)}")
    print()

    # Diagnostic B: Local two-point power-law diagnostic on N in {20, 24}
    log_s24 = mp.log(s_vals[24])
    log_s20 = mp.log(s_vals[20])
    p_power = -(log_s24 - log_s20) / (mp.log(24) - mp.log(20))
    a_power = s_vals[24] * (mp.mpf(24) ** p_power)

    print("Diagnostic B: Local power-law diagnostic s_N ~ A * N^(-p) on N in {20, 24}:")
    print(f"        Local effective power p = {mp.nstr(p_power, 4)}")
    print(f"        Local prefactor A       = {mp.nstr(a_power, 6)}")
    print()

    # Diagnostic C: Local logarithmic diagnostic s_N ~ A / (log N)^p on N in {20, 24}
    ll_24 = mp.log(mp.log(24))
    ll_20 = mp.log(mp.log(20))
    p_log = -(log_s24 - log_s20) / (ll_24 - ll_20)
    a_log = s_vals[24] * (mp.log(24) ** p_log)

    print("Diagnostic C: Local logarithmic diagnostic s_N ~ A / (log N)^p on N in {20, 24}:")
    print(f"        Local logarithmic power p = {mp.nstr(p_log, 4)}")
    print(f"        Local prefactor A         = {mp.nstr(a_log, 6)}")
    print()

    # Comparison Table of Diagnostic Models
    print("Comparison of diagnostic models against observed s_N:")
    print(f"{'N':>4}  {'s_N observed':>14}  {'Diag A (Quad Extrap)':>22}  {'Diag B (Local Power)':>22}  {'Diag C (Local Log)':>20}")
    print("-" * 88)
    for N in N_LIST:
        obs = s_vals[N]
        pred_a = s_inf_poly + c1_poly / N + c2_poly / (N ** 2)
        pred_b = a_power * (mp.mpf(N) ** (-p_power))
        pred_c = a_log * (mp.log(N) ** (-p_log))
        print(f"{N:>4d}  {mp.nstr(obs, 6):>14}  {mp.nstr(pred_a, 6):>22}  {mp.nstr(pred_b, 6):>22}  {mp.nstr(pred_c, 6):>20}")
    print("-" * 88)
    print("Note: These scaling diagnostics are exploratory. With N in {16, 20, 24}, Diagnostic A")
    print("is a 3-point interpolation, while Diagnostics B and C have a single degree of information.")
    print("They serve as diagnostics to guide larger-N calculations rather than inferential proof.")
    print()
    print("Strategic observation: If Delta_N = -log|D_0| + log|D_1| behaves smoothly across N")
    print("while -log|D_0| and -log|D_1| grow rapidly, this confirms that the ratio D_0 / D_1")
    print("(the first-jet scale u_1) is the analytically tractable object to target for higher-N")
    print("calculations (e.g. N = 40, 60, 100).")
    print()
    print("=" * 80)
    print("END OF CELL 54")
    print("=" * 80)


if __name__ == "__main__":
    main()
