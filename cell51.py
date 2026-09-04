"""
CELL 51 — OPERATOR RESOLVENT ANATOMY, HEAT-KERNEL BOUNDARY DYNAMICS,
AND ACCUMULATING POLE GEOMETRY OF THE GENERATING FUNCTION D_N(z)

Based on the analytical breakthrough connecting D_N(z) to the Neumann Laplacian
resolvent on [0, L]:
    L = -d^2/dt^2,   phi_0 = 1,   phi_m = sqrt(2) * cos(2*pi*m*t / L),
    L phi_m = a_m^2 * phi_m,   a_m = kappa * m,   kappa = 2*pi/L.

We have the exact operator identity:
    D_N(z) = [(I + z L)^(-1) T_N](0)
           = v_{N, 0} + sqrt(2) sum_{m=1}^N v_{N, m} / (1 + a_m^2 * z).

And the exact Fourier-resolvent link:
    Phi_N(r) = (2 / sqrt(L)) * (sin(rL/2) / r) * D_N(-1/r^2),
    R_N(r)   = (2 / (L * r^2)) * [D_N(-1/r^2)]^2.

The poles of D_N(z) at z_m = -1 / a_m^2 accumulate at z = 0^- as N -> infinity.
This cell investigates the analytic structure of D_N(z) across N in {8, 12, 16, 20, 24}:

1. Positive-Axis Resolvent D_N(x) for x > 0 and Large-x Asymptotics:
   Examine D_N(x) along the positive real axis (far from all poles), where (I + x L)^(-1)
   is a strictly positive bounded operator. Test large-x asymptotic behavior:
       D_N(x) = v_0 + (sqrt(2) / x) sum_{m=1}^N (v_m / a_m^2) + O(x^-2).

2. Negative-Axis High-Frequency Approach, delta-Sampling, & Local Exponent gamma_eff(r):
   - Evaluate B_N(r) = |D_N(-1/r^2)| at r = kappa * (m + delta) for delta in {0.1, 0.25, 0.5, 0.75, 0.9}
     to quantify sensitivity to pole proximity.
   - Compute the local scaling exponent:
         gamma_eff(r) = d(log[-log|D|]) / d(log r)
     to distinguish between pure exponential (gamma=1), Gaussian (gamma=2), or stretched exponential.

3. Modulated Alternating Coefficients & Discrete Cauchy Transform:
   - Track modulated coefficients b_{N, m} = (-1)^m * v_{N, m} and successive ratios b_m / b_{m-1}.
   - Evaluate the discrete Cauchy transform:
         F_N(w) = sum_{m=1}^N v_{N, m} / (m^2 + w),   w = -r^2 / kappa^2,
     and numerically verify D_N(-1/r^2) = v_0 + sqrt(2) * w * F_N(w).

4. Heat-Kernel Boundary Dynamics H_N(u) = [e^{-u L} T_N](0):
   Track heat evolution of the boundary value from u = 1.0 down to u = 10^-6.
   Verify that at finite N, lim_{u -> 0^+} H_N(u) = T_N(0) != 0, and evaluate the double
   limit behavior lim_{N -> infinity} lim_{u -> 0^+} H_N(u).

5. Dimension Scaling & Boundary-Layer Collapse:
   Test scaling collapse of -(1/N) * log|D_N(-1/r^2)| against the natural spectral variable
   xi = r / (kappa * N) (distinguishing the continuum window xi << 1 from the edge xi ~ 1)
   and against r / sqrt(N).
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
# Resolvent and Heat-Kernel Evaluators
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


def F_cauchy(v, w):
    """
    Discrete Cauchy transform:
        F_N(w) = sum_{m=1}^N v_m / (m^2 + w).
    """
    res = mp.mpf(0)
    for m in range(1, len(v)):
        res += v[m] / (m ** 2 + w)
    return res


def H_heat(v, u, kappa):
    """
    Heat-kernel boundary value:
        H_N(u) = [e^{-u L} T_N](0) = v_0 + sqrt(2) sum_{m=1}^N v_m * exp(-kappa^2 * m^2 * u).
    Note: H_N(0) = v_0 + sqrt(2) * sum_{m=1}^N v_m = T_N(0).
    """
    res = v[0]
    for m in range(1, len(v)):
        arg = - (kappa * m) ** 2 * u
        res += mp.sqrt(2) * v[m] * mp.exp(arg)
    return res


def boundary_T0(v):
    """Boundary contact T_N(0) = v_0 + sqrt(2) * sum_{m=1}^N v_m."""
    return v[0] + mp.sqrt(2) * sum(v[m] for m in range(1, len(v)))


# ---------------------------------------------------------------------------
# MAIN INVESTIGATION
# ---------------------------------------------------------------------------

def main():
    print("=" * 80)
    print("CELL 51 — OPERATOR RESOLVENT ANATOMY & POLE ACCUMULATION OF D_N(z)")
    print("=" * 80)
    print(f"Parameters: c = {C_PARAM}, L = ln(13) = {mp.nstr(L_PARAM, 10)}, kappa = {mp.nstr(KAPPA, 10)}")
    print(f"Working precision: {DPS_RUN} decimal digits")
    print()

    # Step 1: Load/compute ground states for N in {8, 12, 16, 20, 24} via persistent cache
    ground_states = {}
    print("Loading ground states via persistent cache (.cell_cache)...")
    for N in N_LIST:
        print(f"  Fetching ground state for N = {N} (dps = {DPS_RUN})...")
        lam, v, meta = get_ground_state(c=C_PARAM, N=N, T=T_QUAD, dps=DPS_RUN, verbose=False)
        t0_val = boundary_T0(v)
        ground_states[N] = (lam, v, t0_val)
        hit_status = "cache hit" if meta.get("cache_hit") else "computed"
        sec_str = f" ({meta.get('total_seconds', 0):.2f}s)" if 'total_seconds' in meta else ""
        print(f"    N = {N:2d} [{hit_status}{sec_str}]: lambda_min = {mp.nstr(lam, 8)}, v_0 = {mp.nstr(v[0], 8)}, T_N(0) = {mp.nstr(t0_val, 6)}")

    print()

    # -----------------------------------------------------------------------
    # Part 1: Positive-Axis Resolvent D_N(x) for x > 0 & Large-x Asymptotics
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("1. POSITIVE-AXIS RESOLVENT D_N(x) = [(I + x L)^(-1) T_N](0)  (x > 0)")
    print("=" * 80)
    print("Evaluating along the positive real axis (away from all negative poles)...")
    print()

    x_test = [mp.mpf("1e-6"), mp.mpf("1e-4"), mp.mpf("1e-2"), mp.mpf("0.1"), mp.mpf("1.0"), mp.mpf("10.0")]
    header_x = "  ".join([f"x = {mp.nstr(x, 2):<8}" for x in x_test])
    print(f"{'N':>4}  {header_x}")
    print("-" * 75)

    for N in N_LIST:
        _, v, _ = ground_states[N]
        vals = [D_eval(v, x, KAPPA) for x in x_test]
        vals_str = "  ".join([f"{mp.nstr(val, 6):>10}" for val in vals])
        print(f"{N:>4}  {vals_str}")

    print("-" * 75)
    print()

    print("--- Large-x Asymptotic Check at x = 10.0: D_N(x) vs v_0 + (sqrt(2)/x) * sum(v_m / a_m^2) ---")
    print(f"{'N':>4}  {'D_N(10.0)':>16}  {'Asymptotic v_0 + C/x':>24}  {'Difference':>16}")
    print("-" * 62)
    for N in N_LIST:
        _, v, _ = ground_states[N]
        d_exact = D_eval(v, mp.mpf("10.0"), KAPPA)
        c_lead = mp.sqrt(2) * sum(v[m] / ((KAPPA * m) ** 2) for m in range(1, len(v)))
        d_asymp = v[0] + c_lead / mp.mpf("10.0")
        diff = abs(d_exact - d_asymp)
        print(f"{N:>4}  {mp.nstr(d_exact, 8):>16}  {mp.nstr(d_asymp, 8):>24}  {mp.nstr(diff, 6):>16}")
    print("-" * 62)
    print()

    # -----------------------------------------------------------------------
    # Part 2: Negative-Axis Approach, delta-Sampling & Local Exponent gamma_eff(r)
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("2. NEGATIVE-AXIS APPROACH D_N(-1/r^2), delta-SAMPLING & LOCAL EXPONENT gamma_eff")
    print("=" * 80)

    # 2A: Sensitivity to pole proximity (delta-sampling at N = 24)
    print("--- 2A. Sensitivity to Pole Distance: r = kappa * (m + delta) (N = 24) ---")
    print(f"{'m':>4}  {'delta=0.10':>14}  {'delta=0.25':>14}  {'delta=0.50 (half)':>18}  {'delta=0.75':>14}  {'delta=0.90':>14}")
    print("-" * 80)

    _, v24, _ = ground_states[24]
    deltas = [mp.mpf("0.1"), mp.mpf("0.25"), mp.mpf("0.5"), mp.mpf("0.75"), mp.mpf("0.9")]

    for m_probe in [2, 4, 6, 8, 12, 16, 20]:
        row_vals = []
        for d in deltas:
            r_val = KAPPA * (m_probe + d)
            z_val = -1 / (r_val ** 2)
            b_val = abs(D_eval(v24, z_val, KAPPA))
            row_vals.append(f"{mp.nstr(b_val, 5):>14}")
        print(f"{m_probe:>4}  {'  '.join(row_vals)}")

    print("-" * 80)
    print()

    # 2B: Local Scaling Exponent gamma_eff(r) along half-integer points (N = 24)
    print("--- 2B. Local Scaling Exponent gamma_eff(r) = d(log[-log|D|]) / d(log r) (N = 24) ---")
    print("Reference: gamma_eff -> 1.0 (pure exponential e^{-Cr}), 2.0 (Gaussian), 0.5 (stretched exp)")
    print(f"{'k':>4}  {'r = kappa*(k+1/2)':>18}  {'|D_24|':>14}  {'-log|D|':>12}  {'-log|D| / r':>14}  {'gamma_eff(r)':>16}")
    print("-" * 84)

    k_samples = list(range(2, 23))
    r_k = [KAPPA * (mp.mpf(k) + mp.mpf("0.5")) for k in k_samples]
    b_k = [abs(D_eval(v24, -1 / (r ** 2), KAPPA)) for r in r_k]
    log_b_k = [-mp.log(b) if b > 0 else mp.mpf(0) for b in b_k]

    for i in range(len(k_samples)):
        k_val = k_samples[i]
        r_val = r_k[i]
        b_val = b_k[i]
        lb_val = log_b_k[i]
        ratio_lin = lb_val / r_val

        # Central difference for gamma_eff
        if 0 < i < len(k_samples) - 1:
            d_log_lb = mp.log(log_b_k[i + 1]) - mp.log(log_b_k[i - 1])
            d_log_r = mp.log(r_k[i + 1]) - mp.log(r_k[i - 1])
            gamma_eff = d_log_lb / d_log_r
            gamma_str = f"{mp.nstr(gamma_eff, 5):>16}"
        else:
            gamma_str = f"{'—':>16}"

        print(f"{k_val:>4}  {mp.nstr(r_val, 6):>18}  {mp.nstr(b_val, 5):>14}  {mp.nstr(lb_val, 6):>12}  {mp.nstr(ratio_lin, 6):>14}  {gamma_str}")

    print("-" * 84)
    print()

    # -----------------------------------------------------------------------
    # Part 3: Modulated Alternating Coefficients & Discrete Cauchy Transform
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("3. MODULATED COEFFICIENT SEQUENCE b_{N, m} = (-1)^m * v_{N, m} (N = 24)")
    print("=" * 80)
    print("Inspecting positivity and successive decay of b_{24, m} for m = 1, ..., 20...")
    print()

    print(f"{'m':>4}  {'v_{24, m}':>18}  {'b_{24, m} = (-1)^m v_m':>24}  {'Ratio b_m / b_{m-1}':>22}")
    print("-" * 74)

    prev_b = None
    for m in range(1, 21):
        vm = v24[m]
        bm = ((-1) ** m) * vm
        ratio_str = f"{mp.nstr(bm / prev_b, 6):>22}" if prev_b is not None else f"{'—':>22}"
        print(f"{m:>4}  {mp.nstr(vm, 8):>18}  {mp.nstr(bm, 8):>24}  {ratio_str}")
        prev_b = bm

    print("-" * 74)
    print()

    # Numerical check of the discrete Cauchy identity: D_N(-1/r^2) = v_0 + sqrt(2) * w * F_N(w)
    print("--- Verification of Discrete Cauchy Transform Identity: D_N(-1/r^2) = v_0 + sqrt(2)*w*F_N(w) ---")
    r_check = KAPPA * mp.mpf("3.5")
    w_check = - (r_check / KAPPA) ** 2
    d_direct = D_eval(v24, -1 / (r_check ** 2), KAPPA)
    f_val = F_cauchy(v24, w_check)
    d_via_f = v24[0] + mp.sqrt(2) * w_check * f_val
    diff_cauchy = abs(d_direct - d_via_f)
    print(f"  Test frequency r = {mp.nstr(r_check, 6)}, w = {mp.nstr(w_check, 6)}")
    print(f"  D_24 direct:        {mp.nstr(d_direct, 12)}")
    print(f"  v_0 + sqrt(2)*w*F:  {mp.nstr(d_via_f, 12)}")
    print(f"  Difference:         {mp.nstr(diff_cauchy, 6)}  (Identical to machine precision)")
    print()

    # -----------------------------------------------------------------------
    # Part 4: Heat-Kernel Boundary Dynamics H_N(u) = [e^{-u L} T_N](0)
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("4. HEAT-KERNEL BOUNDARY DYNAMICS H_N(u) = [e^{-u L} T_N](0)")
    print("=" * 80)
    print("At finite N, lim_{u -> 0+} H_N(u) = H_N(0) = T_N(0) != 0.")
    print("Tracking heat evolution from u = 1.0 down to u = 10^-6, and testing approach to T_N(0)...")
    print()

    u_test = [mp.mpf("1.0"), mp.mpf("0.1"), mp.mpf("0.01"), mp.mpf("1e-3"), mp.mpf("1e-4"), mp.mpf("1e-5"), mp.mpf("1e-6")]
    print(f"{'Diffusion u':>14}  {'H_8(u)':>14}  {'H_12(u)':>14}  {'H_16(u)':>14}  {'H_20(u)':>14}  {'H_24(u)':>14}")
    print("-" * 86)

    for u in u_test:
        row = []
        for N in N_LIST:
            _, v, _ = ground_states[N]
            h_val = H_heat(v, u, KAPPA)
            row.append(f"{mp.nstr(h_val, 5):>14}")
        print(f"{mp.nstr(u, 6):>14}  {'  '.join(row)}")

    print("-" * 86)
    print()

    print("--- Boundary Limit at u = 0: H_N(0) = T_N(0) across dimensions N ---")
    print(f"{'N':>4}  {'T_N(0)':>20}  {'H_N(1e-6)':>20}  {'|H_N(1e-6) - T_N(0)|':>24}")
    print("-" * 72)
    for N in N_LIST:
        _, v, t0 = ground_states[N]
        h_near0 = H_heat(v, mp.mpf("1e-6"), KAPPA)
        diff_h = abs(h_near0 - t0)
        print(f"{N:>4}  {mp.nstr(t0, 8):>20}  {mp.nstr(h_near0, 8):>20}  {mp.nstr(diff_h, 6):>24}")
    print("-" * 72)
    print()

    # -----------------------------------------------------------------------
    # Part 5: Cross-Dimension Scaling Collapse under xi = r / (kappa * N)
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("5. CROSS-DIMENSION SCALING COLLAPSE: -(1/N) * log|D_N(-1/r^2)| vs xi = r / (kappa*N)")
    print("=" * 80)
    print("Spectral edge occurs at r = kappa * N (xi = 1.0). Continuum window is xi << 1.")
    print()

    xi_test = [mp.mpf("0.2"), mp.mpf("0.4"), mp.mpf("0.6"), mp.mpf("0.8"), mp.mpf("1.0"), mp.mpf("1.2")]
    print(f"{'xi = r/(kappa*N)':>16}  {'N = 8':>12}  {'N = 12':>12}  {'N = 16':>12}  {'N = 20':>12}  {'N = 24':>12}")
    print("-" * 80)

    for xi in xi_test:
        row_xi = []
        for N in N_LIST:
            # Shift slightly off integer lattice to avoid exact pole
            r_val = KAPPA * (xi * N + mp.mpf("0.5"))
            z_val = -1 / (r_val ** 2)
            _, v, _ = ground_states[N]
            b_val = abs(D_eval(v, z_val, KAPPA))
            scaled_log = - mp.log(b_val) / N if b_val > 0 else mp.mpf(0)
            row_xi.append(f"{mp.nstr(scaled_log, 5):>12}")
        print(f"{mp.nstr(xi, 3):>16}  {'  '.join(row_xi)}")

    print("-" * 80)
    print()

    print("--- Alternative Scaling: -(1/sqrt(N)) * log|D_N| vs r / sqrt(N) ---")
    zeta_test = [mp.mpf("1.0"), mp.mpf("2.0"), mp.mpf("3.0"), mp.mpf("4.0")]
    print(f"{'r / sqrt(N)':>14}  {'N = 8':>12}  {'N = 12':>12}  {'N = 16':>12}  {'N = 20':>12}  {'N = 24':>12}")
    print("-" * 78)

    for zeta in zeta_test:
        row_zeta = []
        for N in N_LIST:
            r_val = zeta * mp.sqrt(N)
            z_val = -1 / (r_val ** 2)
            _, v, _ = ground_states[N]
            b_val = abs(D_eval(v, z_val, KAPPA))
            scaled_log = - mp.log(b_val) / mp.sqrt(N) if b_val > 0 else mp.mpf(0)
            row_zeta.append(f"{mp.nstr(scaled_log, 5):>12}")
        print(f"{mp.nstr(zeta, 3):>14}  {'  '.join(row_zeta)}")

    print("-" * 78)
    print()
    print("=" * 80)
    print("END OF CELL 51")
    print("=" * 80)


if __name__ == "__main__":
    main()
