"""
CELL 52 — DOUBLE-SCALING BOUNDARY LAYER, HEAT-RESOLVENT SCALING COLLAPSE,
AND POLE-PROTECTED NEGATIVE-AXIS CAUCHY LATTICE PROBE

Following the findings of Cell 51:
1. Rejection of raw e^{-Cr} fits: |D_N(-1/r^2)| displays persistent lattice-scale
   oscillations across pole cells with irregular mode signs b_m = (-1)^m v_m.
2. Heat relaxation: H_N(u) = [e^{-u L} T_N](0) plunges from O(1) down to
   H_N(10^-6) = 1.77e-20 ~= T_N(0), suggesting a shrinking boundary layer
   associated with the spectral edge scale:
       u_N ~ 1 / a_N^2 = 1 / (kappa^2 * N^2).
   For N = 24, u_24 ~ 2.9e-4.

Cell 52 executes the targeted double-scaling boundary-layer experiment proposed
by the reviewer, with four concrete safeguards:
1. Double-Scaled Heat Semigroup:
       H_N(s / (kappa^2 * N^2))
   testing whether curves approach an N-independent limiting profile H_infty(s).
2. Normalized Boundary Profile:
       Theta_N(s) = H_N(s / (kappa^2 * N^2)) / T_N(0)
   testing whether H_N(u) ~ T_N(0) * Theta(kappa^2 * N^2 * u) exhibits universal
   boundary-layer profile collapse (without assuming monotonicity).
3. Scaled Resolvent (Double-Scaling Limit):
       D_N(sigma / (kappa^2 * N^2))
   for sigma = kappa^2 * N^2 * z fixed as N -> infty, z -> 0^+.
4. Pole-Protected Negative-Axis Probe:
   Non-pole spectral offsets (delta = 0.37) and automated pole guards guaranteeing
   dist(r/kappa, Z) >= 0.05, eliminating all accidental pole encounters while
   explicitly isolating lattice oscillations from envelope decay.
5. Two-Scale Resolvent Decomposition (Exact Closed Form):
   Decomposes D_N(x) = D_N^{BL}(x) + D_N^{bulk}(x) at u_N = 1 / (kappa*N)^2,
   testing how the ratio D_BL / D_total transitions across sigma = kappa^2 * N^2 * x.
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
# Evaluators and Pole Guards
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


def boundary_T0(v):
    """Boundary contact T_N(0) = v_0 + sqrt(2) * sum_{m=1}^N v_m."""
    return v[0] + mp.sqrt(2) * sum(v[m] for m in range(1, len(v)))


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
    print("CELL 52 — DOUBLE-SCALING BOUNDARY LAYER & CAUCHY LATTICE ARCHITECTURE")
    print("=" * 80)
    print(f"Parameters: c = {C_PARAM}, L = ln(13) = {mp.nstr(L_PARAM, 10)}, kappa = {mp.nstr(KAPPA, 10)}")
    print(f"Working precision: {DPS_RUN} decimal digits")
    print()

    # Step 1: Load ground states via persistent cache (.cell_cache)
    ground_states = {}
    print("Loading ground states via persistent cache (.cell_cache)...")
    for N in N_LIST:
        lam, v, meta = get_ground_state(c=C_PARAM, N=N, T=T_QUAD, dps=DPS_RUN, verbose=False)
        t0_val = boundary_T0(v)
        ground_states[N] = (lam, v, t0_val)
        hit_status = "cache hit" if meta.get("cache_hit") else "computed"
        sec_str = f" ({meta.get('total_seconds', 0):.2f}s)" if 'total_seconds' in meta else ""
        print(f"    N = {N:2d} [{hit_status}{sec_str}]: lambda_min = {mp.nstr(lam, 8)}, v_0 = {mp.nstr(v[0], 8)}, T_N(0) = {mp.nstr(t0_val, 6)}")

    print()

    # -----------------------------------------------------------------------
    # Part 1: Double-Scaled Heat Semigroup H_N(s / (kappa^2 * N^2))
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("1. DOUBLE-SCALED HEAT SEMIGROUP: H_N(s / (kappa^2 * N^2))")
    print("=" * 80)
    print("Testing whether curves approach an N-independent limiting profile H_infty(s)")
    print("under the boundary-layer variable s = kappa^2 * N^2 * u as N -> infty...")
    print()

    s_test = [
        mp.mpf("1e-4"),
        mp.mpf("1e-3"),
        mp.mpf("1e-2"),
        mp.mpf("0.1"),
        mp.mpf("0.5"),
        mp.mpf("1.0"),
        mp.mpf("5.0"),
        mp.mpf("10.0"),
        mp.mpf("50.0"),
        mp.mpf("100.0"),
    ]

    print(f"{'Scaled time s':>16}  {'N = 8':>12}  {'N = 12':>12}  {'N = 16':>12}  {'N = 20':>12}  {'N = 24':>12}")
    print("-" * 82)

    for s in s_test:
        row = []
        for N in N_LIST:
            _, v, _ = ground_states[N]
            u_phys = s / ((KAPPA * N) ** 2)
            h_val = H_heat(v, u_phys, KAPPA)
            row.append(f"{mp.nstr(h_val, 5):>12}")
        print(f"{mp.nstr(s, 4):>16}  {'  '.join(row)}")

    print("-" * 82)
    print()

    # -----------------------------------------------------------------------
    # Part 2: Normalized Boundary-Layer Profile Theta_N(s) = H_N(s/(kappa^2*N^2)) / T_N(0)
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("2. NORMALIZED BOUNDARY PROFILE: Theta_N(s) = H_N(s / (kappa^2 * N^2)) / T_N(0)")
    print("=" * 80)
    print("Testing whether H_N(u) ~ T_N(0) * Theta(kappa^2 * N^2 * u) exhibits universal")
    print("scaling profile collapse (without assuming monotonic growth)...")
    print("Note: Theta_N(0) = 1.0 identically. Inspecting profile evolution for small s...")
    print()

    s_small = [
        mp.mpf("1e-5"),
        mp.mpf("1e-4"),
        mp.mpf("1e-3"),
        mp.mpf("0.01"),
        mp.mpf("0.05"),
        mp.mpf("0.1"),
        mp.mpf("0.2"),
        mp.mpf("0.5"),
        mp.mpf("1.0"),
    ]

    print(f"{'Scaled time s':>16}  {'N = 8':>12}  {'N = 12':>12}  {'N = 16':>12}  {'N = 20':>12}  {'N = 24':>12}")
    print("-" * 82)

    for s in s_small:
        row = []
        for N in N_LIST:
            _, v, t0 = ground_states[N]
            u_phys = s / ((KAPPA * N) ** 2)
            h_val = H_heat(v, u_phys, KAPPA)
            ratio = h_val / t0
            row.append(f"{mp.nstr(ratio, 5):>12}")
        print(f"{mp.nstr(s, 4):>16}  {'  '.join(row)}")

    print("-" * 82)
    print()

    # -----------------------------------------------------------------------
    # Part 3: Scaled Positive Resolvent D_N(sigma / (kappa^2 * N^2))
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("3. SCALED POSITIVE RESOLVENT: D_N(sigma / (kappa^2 * N^2))")
    print("=" * 80)
    print("Evaluating the double-scaling limit N -> infty, z -> 0^+ with sigma = kappa^2 * N^2 * z fixed...")
    print("Probing the resolvent at the boundary-layer scale x = sigma / (kappa^2 * N^2)...")
    print()

    sigma_test = [
        mp.mpf("1e-3"),
        mp.mpf("1e-2"),
        mp.mpf("0.1"),
        mp.mpf("1.0"),
        mp.mpf("10.0"),
        mp.mpf("100.0"),
    ]

    print(f"{'Scaled res sigma':>18}  {'N = 8':>12}  {'N = 12':>12}  {'N = 16':>12}  {'N = 20':>12}  {'N = 24':>12}")
    print("-" * 84)

    for sigma in sigma_test:
        row = []
        for N in N_LIST:
            _, v, _ = ground_states[N]
            x_phys = sigma / ((KAPPA * N) ** 2)
            d_val = D_eval(v, x_phys, KAPPA)
            row.append(f"{mp.nstr(d_val, 5):>12}")
        print(f"{mp.nstr(sigma, 4):>18}  {'  '.join(row)}")

    print("-" * 84)
    print()

    # -----------------------------------------------------------------------
    # Part 4: Pole-Protected Negative-Axis Probe & Lattice Oscillation Decoupling
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("4. POLE-PROTECTED NEGATIVE-AXIS PROBE & LATTICE DECOUPLING (N = 24)")
    print("=" * 80)
    print("A. Dense spectral sweep: m_coord = 24*xi + 0.37 across xi in [0.10, 1.30] (N = 24).")
    print("   Automated pole guard guarantees dist(m_coord, Z) >= 0.05 (no pole can be hit).")
    print()

    print(f"{'xi':>8}  {'r':>10}  {'m_coord':>10}  {'Pole m*':>8}  {'Dist delta':>12}  {'|D_24|':>14}  {'-log|D|':>12}  {'-log|D| / r':>14}")
    print("-" * 94)

    _, v24, _ = ground_states[24]
    a24 = KAPPA * 24

    # Dense test points: xi from 0.10 to 1.30 in steps of 0.05 with non-pole offset delta = 0.37
    xi_dense = [mp.mpf(step) * mp.mpf("0.05") for step in range(2, 27)]

    for xi in xi_dense:
        raw_m = 24 * xi + mp.mpf("0.37")
        safe_m, was_nudged = guard_pole(raw_m, min_dist=mp.mpf("0.05"))
        dist_delta, m_nearest = pole_distance(safe_m)

        r_val = KAPPA * safe_m
        z_val = -1 / (r_val ** 2)
        b_val = abs(D_eval(v24, z_val, KAPPA))
        log_b = -mp.log(b_val) if b_val > 0 else mp.mpf(0)
        rate = log_b / r_val if r_val > 0 else mp.mpf(0)

        nudge_flag = "*" if was_nudged else " "
        print(f"{mp.nstr(xi, 3):>8}  {mp.nstr(r_val, 5):>10}  {mp.nstr(safe_m, 5):>9}{nudge_flag} {m_nearest:>8d}  {mp.nstr(dist_delta, 4):>12}  {mp.nstr(b_val, 5):>14}  {mp.nstr(log_b, 5):>12}  {mp.nstr(rate, 5):>14}")

    print("-" * 94)
    print("(* indicates automated pole guard adjusted coordinate to maintain delta >= 0.05)")
    print()

    # B. Multi-depth pole-cell delta traversals (extending Cell 51 delta-sampling)
    print("B. Multi-depth pole-cell delta traversals: r = kappa * (m + delta) for N = 24:")
    print("   Testing across low bulk (m=6), mid bulk (m=14), near edge (m=22), and exterior (m=26):")
    print(f"{'Cell m':>8}  {'delta = 0.10':>16}  {'delta = 0.25':>16}  {'delta = 0.50':>16}  {'delta = 0.75':>16}  {'delta = 0.90':>16}")
    print("-" * 94)

    for m_cell in [6, 14, 22, 26]:
        row = []
        for delta_val in [mp.mpf("0.10"), mp.mpf("0.25"), mp.mpf("0.50"), mp.mpf("0.75"), mp.mpf("0.90")]:
            m_probe = mp.mpf(m_cell) + delta_val
            r_probe = KAPPA * m_probe
            z_probe = -1 / (r_probe ** 2)
            d_cell = abs(D_eval(v24, z_probe, KAPPA))
            row.append(f"{mp.nstr(d_cell, 5):>16}")
        print(f"{m_cell:>8d}  {'  '.join(row)}")

    print("-" * 94)
    print()

    # C. Cross-dimension comparison with automated pole protection
    print("C. Cross-dimension comparison of -log|D_N| / (kappa*N) at target scaled locations:")
    print("   Using target locations with automated pole guard ensuring dist(N*xi, Z) >= 0.05:")
    print(f"{'Target xi':>10}  {'N = 8':>14}  {'N = 12':>14}  {'N = 16':>14}  {'N = 20':>14}  {'N = 24':>14}")
    print("-" * 86)

    for xi_target in [mp.mpf("0.53"), mp.mpf("0.83"), mp.mpf("1.07")]:
        row = []
        for N in N_LIST:
            _, v_n, _ = ground_states[N]
            raw_m_n = N * xi_target
            safe_m_n, _ = guard_pole(raw_m_n, min_dist=mp.mpf("0.05"))
            r_target = KAPPA * safe_m_n
            z_target = -1 / (r_target ** 2)
            d_val = abs(D_eval(v_n, z_target, KAPPA))
            a_n = KAPPA * N
            val_scaled = -mp.log(d_val) / a_n if d_val > 0 else mp.mpf(0)
            dist_n, _ = pole_distance(safe_m_n)
            row.append(f"{mp.nstr(val_scaled, 5):>8} (d={mp.nstr(dist_n, 2):>4})")
        print(f"{mp.nstr(xi_target, 3):>10}  {'  '.join(row)}")

    print("-" * 86)
    print()

    # -----------------------------------------------------------------------
    # Part 5: Two-Scale Resolvent Decomposition (Exact Closed Form)
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("5. TWO-SCALE RESOLVENT INTEGRAL SPLITTING (N = 24)")
    print("=" * 80)
    print("Splitting the Laplace integral at the boundary-layer scale u_N = 1 / (kappa*N)^2:")
    print("    D_N(x) = (1/x) int_0^{u_N} e^{-u/x} H_N(u) du  +  (1/x) int_{u_N}^infty e^{-u/x} H_N(u) du")
    print("           = D_N^{BL}(x)  +  D_N^{bulk}(x)")
    print()
    print("Physical interpretation:")
    print("    The relevant dimensionless ratio is sigma = x / u_N = kappa^2 * N^2 * x:")
    print("    - sigma << 1: resolvent kernel samples predominantly the boundary layer")
    print("    - sigma ~  1: transition regime between boundary layer and bulk")
    print("    - sigma >> 1: resolvent kernel samples predominantly the continuum bulk")
    print()
    print("Using exact algebraic closed-form evaluation:")
    print("    Term m:  (1/x) int_0^{u_N} e^{-u/x} e^{-kappa^2 m^2 u} du = [1 - exp(-(1/x + kappa^2 m^2) u_N)] / [1 + kappa^2 m^2 x]")
    print()

    u_scale_24 = 1 / ((KAPPA * 24) ** 2)
    print(f"Boundary layer threshold for N = 24: u_24 = {mp.nstr(u_scale_24, 6)}")
    print()

    def exact_splitting(v, x_val, kappa, u_cut):
        """
        Computes (D_BL, D_bulk, D_total) in exact closed form for all modes.
        """
        # Mode 0 (constant term v_0)
        exp_0 = mp.exp(-u_cut / x_val)
        d_bl = v[0] * (1 - exp_0)
        d_bulk = v[0] * exp_0

        # Modes m = 1, ..., N
        for m in range(1, len(v)):
            km2 = (kappa * m) ** 2
            denom = 1 + km2 * x_val
            exp_factor = mp.exp(-(1 / x_val + km2) * u_cut)
            d_bl += mp.sqrt(2) * v[m] * (1 - exp_factor) / denom
            d_bulk += mp.sqrt(2) * v[m] * exp_factor / denom

        d_total = d_bl + d_bulk
        return d_bl, d_bulk, d_total

    print("A. Fixed physical resolvent points x in {0.1, 1.0, 10.0}:")
    print(f"{'x':>8}  {'D_24(x) exact':>18}  {'D_BL':>22}  {'D_bulk':>22}  {'D_BL / D_total':>16}  {'|D_BL| / |D_bulk|':>18}")
    print("-" * 110)

    for x_val in [mp.mpf("0.1"), mp.mpf("1.0"), mp.mpf("10.0")]:
        d_bl, d_bulk, d_sum = exact_splitting(v24, x_val, KAPPA, u_scale_24)
        d_direct = D_eval(v24, x_val, KAPPA)
        ratio_total = d_bl / d_sum
        ratio_bl_bulk = abs(d_bl) / abs(d_bulk) if abs(d_bulk) > 0 else mp.mpf(0)
        print(f"{mp.nstr(x_val, 3):>8}  {mp.nstr(d_direct, 10):>18}  {mp.nstr(d_bl, 8):>22}  {mp.nstr(d_bulk, 8):>22}  {mp.nstr(ratio_total, 6):>16}  {mp.nstr(ratio_bl_bulk, 6):>18}")

    print("-" * 110)
    print()

    print("B. Scaled boundary-layer resolvent points x = sigma / (kappa*N)^2:")
    print(f"{'sigma':>10}  {'x = sigma*u_24':>16}  {'D_24(x)':>18}  {'D_BL':>20}  {'D_bulk':>20}  {'D_BL / D_total':>16}")
    print("-" * 106)

    for sigma in [mp.mpf("0.01"), mp.mpf("0.1"), mp.mpf("1.0"), mp.mpf("10.0"), mp.mpf("100.0")]:
        x_scaled = sigma * u_scale_24
        d_bl, d_bulk, d_sum = exact_splitting(v24, x_scaled, KAPPA, u_scale_24)
        ratio_total = d_bl / d_sum
        print(f"{mp.nstr(sigma, 4):>10}  {mp.nstr(x_scaled, 6):>16}  {mp.nstr(d_sum, 8):>18}  {mp.nstr(d_bl, 8):>20}  {mp.nstr(d_bulk, 8):>20}  {mp.nstr(ratio_total, 6):>16}")

    print("-" * 106)
    print()
    print("=" * 80)
    print("END OF CELL 52")
    print("=" * 80)


if __name__ == "__main__":
    main()
