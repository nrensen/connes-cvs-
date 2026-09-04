"""
CELL 53 — DUAL-SCALE BOUNDARY LAYER DECOUPLING: ENDPOINT JET CANCELLATION
SCALES AND LARGE-DEVIATION RATE FUNCTION

Following the decisive results of Cell 52:
1. Decoupling of Spectral-Edge vs. Cancellation Scales:
       u_edge = 1 / (kappa^2 * N^2) ~ N^-2
   acts as a sharp crossover scale for the resolvent integral (D_BL/D_total
   transitions from 1.0 to 1e-16 across sigma = x/u_edge in [0.01, 100]), but
   normalized heat profiles Theta_N(s) = H_N(s u_edge) / T_N(0) diverge rapidly
   (5900 -> 68700 at s = 1) because T_N(0) vanishes much faster than T_N''(0).
2. Emerging Large-Deviation Structure:
   On the negative axis, -(1 / (kappa*N)) log|D_N| displays stability across
   N in {8, 12, 16, 20, 24} (0.854 -> 0.719 at xi = 1.07), suggesting an
   N-dependent WKB / large-deviation rate function:
       |D_N(-1/r^2)| ~ exp[ -kappa * N * I(r / (kappa * N)) ].

Cell 53 executes the targeted two-prong investigation proposed by the reviewer:

PART 1: THE ENDPOINT JET HIERARCHY & CANCELLATION TIME SCALES
  Computes the even endpoint derivatives D_k(N) = T_{v_N}^{(2k)}(0) for k = 0, ..., 5.
  Defines dimensionalized cancellation time scales:
      u_{k, N} = (|D_0(N)| / |D_k(N)|)^{1/k}
  and compares them against u_edge(N) = 1 / (kappa^2 * N^2).
  Tests the dimensionless cancellation parameter:
      s_cancel(N) = kappa^2 * N^2 * (|D_0| / |D_1|)
  and the dimensionless jet shape invariants:
      beta_N = D_0 * D_2 / D_1^2,   gamma_N = D_0^2 * D_3 / D_1^3.

PART 2: SCALED HEAT PROFILES AT THE CANCELLATION SCALE
  Evaluates H_N(u) under the true cancellation scale:
      u = theta * u_{1, N} = theta * (|D_0| / |D_1|)
  for theta in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0].
  Tests whether Theta_N^{cancel}(theta) = H_N(theta * u_cancel) / T_N(0) collapses
  universally across N.

PART 3: THE LARGE-DEVIATION RATE FUNCTION I_N(xi)
  Samples xi = r / (kappa * N) across a dense grid xi in [0.20, 1.40]
  with automated pole guards ensuring dist(N*xi, Z) >= 0.05 for all N.
  Computes:
      I_N(xi) = - (1 / N) * log|D_N(-1 / (kappa^2 * N^2 * xi^2))|
      J_N(xi) = I_N(xi) / kappa = - (1 / (kappa * N)) * log|D_N|
      S_N(xi) = I_N(xi) / xi   (effective large-xi slope)
  Tests convergence I_N(xi) -> I(xi) and asymptotic linearity I(xi) ~ C * xi.
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
# Evaluators, Derivatives, and Pole Guards
# ---------------------------------------------------------------------------

def D_eval(v, z, kappa):
    """
    D_N(z) = v_0 + sqrt(2) sum_{m=1}^N v_m / (1 + kappa^2 * m^2 * z).
    """
    res = v[0]
    for m in range(1, len(v)):
        denom = 1 + (kappa * m) ** 2 * z
        res += mp.sqrt(2) * v[m] / denom
    return res


def H_heat(v, u, kappa):
    """
    H_N(u) = [e^{-u L} T_N](0) = v_0 + sqrt(2) sum_{m=1}^N v_m * exp(-kappa^2 * m^2 * u).
    """
    res = v[0]
    for m in range(1, len(v)):
        arg = - (kappa * m) ** 2 * u
        res += mp.sqrt(2) * v[m] * mp.exp(arg)
    return res


def compute_endpoint_jets(v, kappa, max_k=5):
    """
    Computes endpoint derivatives D_k = T_{v_N}^{(2k)}(0) for k = 0, ..., max_k.
    D_0 = T_N(0) = v_0 + sqrt(2) * sum_{m=1}^N v_m
    D_k = (-1)^k * sqrt(2) * sum_{m=1}^N (kappa * m)^{2k} * v_m   (k >= 1)
    """
    jets = {}
    # k = 0
    jets[0] = v[0] + mp.sqrt(2) * sum(v[m] for m in range(1, len(v)))

    # k >= 1
    for k in range(1, max_k + 1):
        sign = -1 if (k % 2 == 1) else 1
        s = mp.mpf(0)
        for m in range(1, len(v)):
            am2k = (kappa * m) ** (2 * k)
            s += am2k * v[m]
        jets[k] = sign * mp.sqrt(2) * s

    return jets


def pole_distance(coord):
    """Returns (distance_to_nearest_integer, nearest_integer_pole)."""
    nearest = int(mp.floor(coord + mp.mpf("0.5")))
    dist = abs(coord - nearest)
    return dist, nearest


def guard_pole(coord, min_dist=mp.mpf("0.05")):
    """
    Guarantees that coord is at least min_dist away from any integer pole.
    If dist(coord, Z) < min_dist, shifts outward away from the pole.
    Returns (safe_coord, was_nudged).
    """
    dist, nearest = pole_distance(coord)
    if dist < min_dist:
        shift = min_dist if coord >= nearest else -min_dist
        return nearest + shift, True
    return coord, False


# ---------------------------------------------------------------------------
# MAIN INVESTIGATION
# ---------------------------------------------------------------------------

def main():
    print("=" * 80)
    print("CELL 53 — DUAL-SCALE BOUNDARY LAYER DECOUPLING & LARGE-DEVIATION RATE")
    print("=" * 80)
    print(f"Parameters: c = {C_PARAM}, L = ln(13) = {mp.nstr(L_PARAM, 10)}, kappa = {mp.nstr(KAPPA, 10)}")
    print(f"Working precision: {DPS_RUN} decimal digits")
    print()

    # Step 1: Load ground states via persistent cache (.cell_cache)
    ground_states = {}
    jets_all = {}
    print("Loading ground states via persistent cache (.cell_cache)...")
    for N in N_LIST:
        lam, v, meta = get_ground_state(c=C_PARAM, N=N, T=T_QUAD, dps=DPS_RUN, verbose=False)
        jets = compute_endpoint_jets(v, KAPPA, max_k=5)
        ground_states[N] = (lam, v, jets[0])
        jets_all[N] = jets
        hit_status = "cache hit" if meta.get("cache_hit") else "computed"
        sec_str = f" ({meta.get('total_seconds', 0):.2f}s)" if 'total_seconds' in meta else ""
        print(f"    N = {N:2d} [{hit_status}{sec_str}]: lambda_min = {mp.nstr(lam, 8)}, D_0 = {mp.nstr(jets[0], 6)}, D_1 = {mp.nstr(jets[1], 6)}")

    print()

    # -----------------------------------------------------------------------
    # Part 1: Endpoint Jet Hierarchy and Cancellation Time Scales
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("1. ENDPOINT JET HIERARCHY & CANCELLATION TIME SCALES")
    print("=" * 80)
    print("A. Even endpoint derivatives D_k = T_N^{(2k)}(0) for k = 0, ..., 4:")
    print(f"{'N':>4}  {'D_0':>14}  {'D_1':>14}  {'D_2':>14}  {'D_3':>14}  {'D_4':>14}")
    print("-" * 76)
    for N in N_LIST:
        j = jets_all[N]
        print(f"{N:>4d}  {mp.nstr(j[0], 6):>14}  {mp.nstr(j[1], 6):>14}  {mp.nstr(j[2], 6):>14}  {mp.nstr(j[3], 6):>14}  {mp.nstr(j[4], 6):>14}")
    print("-" * 76)
    print()

    print("B. Dimensionalized cancellation time scales u_{k, N} = (|D_0| / |D_k|)^{1/k}:")
    print(f"{'N':>4}  {'u_edge = 1/(kN)^2':>18}  {'u_{1} = D_0/D_1':>16}  {'u_{2}':>14}  {'u_{3}':>14}  {'u_{4}':>14}")
    print("-" * 82)
    u_cancel_map = {}
    for N in N_LIST:
        j = jets_all[N]
        u_edge = 1 / ((KAPPA * N) ** 2)
        d0_abs = abs(j[0])
        u1 = d0_abs / abs(j[1])
        u2 = mp.sqrt(d0_abs / abs(j[2]))
        u3 = (d0_abs / abs(j[3])) ** (mp.mpf(1) / 3)
        u4 = (d0_abs / abs(j[4])) ** (mp.mpf(1) / 4)
        u_cancel_map[N] = u1
        print(f"{N:>4d}  {mp.nstr(u_edge, 6):>18}  {mp.nstr(u1, 6):>16}  {mp.nstr(u2, 6):>14}  {mp.nstr(u3, 6):>14}  {mp.nstr(u4, 6):>14}")
    print("-" * 82)
    print()

    print("C. Scale ratios against spectral edge: R_{k, N} = u_{k, N} / u_edge = (kappa*N)^2 * u_{k, N}:")
    print("   Note: s_cancel = R_{1, N} = (kappa*N)^2 * (D_0 / D_1).")
    print(f"{'N':>4}  {'s_cancel = R_1':>16}  {'R_2':>14}  {'R_3':>14}  {'R_4':>14}  {'u_1 / u_2':>12}")
    print("-" * 76)
    for N in N_LIST:
        j = jets_all[N]
        aN2 = (KAPPA * N) ** 2
        d0_abs = abs(j[0])
        u1 = d0_abs / abs(j[1])
        u2 = mp.sqrt(d0_abs / abs(j[2]))
        u3 = (d0_abs / abs(j[3])) ** (mp.mpf(1) / 3)
        u4 = (d0_abs / abs(j[4])) ** (mp.mpf(1) / 4)
        r1 = aN2 * u1
        r2 = aN2 * u2
        r3 = aN2 * u3
        r4 = aN2 * u4
        ratio_12 = u1 / u2
        print(f"{N:>4d}  {mp.nstr(r1, 6):>16}  {mp.nstr(r2, 6):>14}  {mp.nstr(r3, 6):>14}  {mp.nstr(r4, 6):>14}  {mp.nstr(ratio_12, 6):>12}")
    print("-" * 76)
    print()

    print("D. Dimensionless jet shape invariants:")
    print("   beta_N = D_0 * D_2 / D_1^2   and   gamma_N = D_0^2 * D_3 / D_1^3:")
    print(f"{'N':>4}  {'beta_N':>16}  {'gamma_N':>16}  {'D_1^2 / (D_0*D_2)':>18}")
    print("-" * 58)
    for N in N_LIST:
        j = jets_all[N]
        d0, d1, d2, d3 = j[0], j[1], j[2], j[3]
        beta = (d0 * d2) / (d1 ** 2)
        gamma = (d0 ** 2 * d3) / (d1 ** 3)
        inv_beta = 1 / beta if beta != 0 else mp.mpf(0)
        print(f"{N:>4d}  {mp.nstr(beta, 6):>16}  {mp.nstr(gamma, 6):>16}  {mp.nstr(inv_beta, 6):>18}")
    print("-" * 58)
    print()

    # -----------------------------------------------------------------------
    # Part 2: Scaled Heat Profiles at the Cancellation Scale u_cancel = D_0 / D_1
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("2. SCALED HEAT PROFILES AT THE CANCELLATION SCALE u = theta * (D_0 / D_1)")
    print("=" * 80)
    print("Testing whether Theta_N^{cancel}(theta) = H_N(theta * u_cancel) / T_N(0)")
    print("collapses universally across N (where linear growth 1 + theta is identical for all N):")
    print()

    theta_test = [
        mp.mpf("0.01"),
        mp.mpf("0.05"),
        mp.mpf("0.1"),
        mp.mpf("0.2"),
        mp.mpf("0.5"),
        mp.mpf("1.0"),
        mp.mpf("2.0"),
        mp.mpf("5.0"),
        mp.mpf("10.0"),
    ]

    print(f"{'theta = u / u_cancel':>22}  {'N = 8':>12}  {'N = 12':>12}  {'N = 16':>12}  {'N = 20':>12}  {'N = 24':>12}")
    print("-" * 88)

    for theta in theta_test:
        row = []
        for N in N_LIST:
            _, v, t0 = ground_states[N]
            u_c = u_cancel_map[N]
            u_phys = theta * u_c
            h_val = H_heat(v, u_phys, KAPPA)
            ratio = h_val / t0
            row.append(f"{mp.nstr(ratio, 6):>12}")
        print(f"{mp.nstr(theta, 4):>22}  {'  '.join(row)}")

    print("-" * 88)
    print()

    # -----------------------------------------------------------------------
    # Part 3: The Large-Deviation Rate Function I_N(xi)
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("3. LARGE-DEVIATION RATE FUNCTION: I_N(xi) = -(1/N) * log|D_N(-1/(kappa^2*N^2*xi^2))|")
    print("=" * 80)
    print("Testing for large-N rate function collapse I_N(xi) -> I(xi) across xi in [0.25, 1.35]:")
    print("Automated pole guard guarantees dist(N*xi, Z) >= 0.05 for every evaluation.")
    print()

    # Selected non-pole target positions across bulk, transition, and exterior
    xi_targets = [
        mp.mpf("0.27"),
        mp.mpf("0.37"),
        mp.mpf("0.47"),
        mp.mpf("0.53"),
        mp.mpf("0.63"),
        mp.mpf("0.71"),
        mp.mpf("0.83"),
        mp.mpf("0.87"),
        mp.mpf("0.95"),
        mp.mpf("1.03"),
        mp.mpf("1.07"),
        mp.mpf("1.15"),
        mp.mpf("1.23"),
        mp.mpf("1.31"),
    ]

    print(f"{'Target xi':>10}  {'I_8(xi)':>12}  {'I_12(xi)':>12}  {'I_16(xi)':>12}  {'I_20(xi)':>12}  {'I_24(xi)':>12}  {'I_24 / xi':>10}")
    print("-" * 86)

    for xi_t in xi_targets:
        row = []
        i24_val = None
        for N in N_LIST:
            _, v_n, _ = ground_states[N]
            raw_m = N * xi_t
            safe_m, _ = guard_pole(raw_m, min_dist=mp.mpf("0.05"))
            r_val = KAPPA * safe_m
            z_val = -1 / (r_val ** 2)
            d_val = abs(D_eval(v_n, z_val, KAPPA))
            i_n = -mp.log(d_val) / N if d_val > 0 else mp.mpf(0)
            row.append(f"{mp.nstr(i_n, 5):>12}")
            if N == 24:
                i24_val = i_n

        slope_str = f"{mp.nstr(i24_val / xi_t, 4):>10}" if i24_val is not None else ""
        print(f"{mp.nstr(xi_t, 3):>10}  {'  '.join(row)}  {slope_str}")

    print("-" * 86)
    print("Note: If I(xi) ~ C * xi for large xi, then -log|D_N| ~ C * N * xi = (C / kappa) * r,")
    print("providing a rigorous large-N derivation of the exponential negative-axis decay envelope.")
    print()

    # -----------------------------------------------------------------------
    # Part 4: Decoupled Two-Scale Asymptotic Comparison
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("4. DECOUPLED TWO-SCALE ASYMPTOTIC SUMMARY")
    print("=" * 80)
    print(f"{'N':>4}  {'u_edge ~ N^-2':>16}  {'u_cancel ~ D_0/D_1':>20}  {'u_cancel / u_edge':>18}  {'Slope I_N(1.07)/1.07':>22}")
    print("-" * 86)
    for N in N_LIST:
        j = jets_all[N]
        u_edge = 1 / ((KAPPA * N) ** 2)
        u_c = abs(j[0]) / abs(j[1])
        ratio_scale = u_c / u_edge
        # Rate function at xi = 1.07
        _, v_n, _ = ground_states[N]
        safe_m, _ = guard_pole(N * mp.mpf("1.07"), min_dist=mp.mpf("0.05"))
        r_v = KAPPA * safe_m
        d_v = abs(D_eval(v_n, -1 / (r_v ** 2), KAPPA))
        i_val = -mp.log(d_v) / N if d_v > 0 else mp.mpf(0)
        slope_val = i_val / mp.mpf("1.07")
        print(f"{N:>4d}  {mp.nstr(u_edge, 6):>16}  {mp.nstr(u_c, 6):>20}  {mp.nstr(ratio_scale, 6):>18}  {mp.nstr(slope_val, 5):>22}")
    print("-" * 86)
    print()
    print("=" * 80)
    print("END OF CELL 53")
    print("=" * 80)


if __name__ == "__main__":
    main()
