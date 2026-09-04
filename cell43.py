"""
CELL 43 — EFFECTIVE POTENTIAL, PROLATE WAVE COMMUTATION, AND STURM-LIOUVILLE RECONSTRUCTION

Cell 42 established that the limiting continuum profile T_{v_infty}(t) on [0, L]
is a symmetric solitary wave peaking at t = L/2, with dual Dirichlet boundary
vanishing:
    T(0) = T(L) = 0.

Cell 43 investigates the governing differential equation of this limiting wave
and the analytical origin of the universal eigenvalue scaling ratio kappa_c.

INVESTIGATIONS:

1. Reconstruction of the Effective Potential V_eff(t) = T''(t) / T(t):
   For a continuous Sturm-Liouville or Schrödinger-type eigenmode:
       -T''(t) + V_eff(t) T(t) = E T(t).
   We evaluate V_eff(t) across the central region t in [0.2 L, 0.8 L] and test
   whether it is quadratic in (t - L/2), corresponding to a prolate/parabolic
   well.

2. Normalized Prolate Operator Reconstruction:
   In normalized symmetric coordinates x = 2t/L - 1 in [-1, 1], where x = 0 is the
   center and x = +-1 are the Dirichlet boundaries, we test the classical
   prolate spheroidal differential operator:
       D_x psi(x) = (1 - x^2) psi''(x) - 2x psi'(x).
   We compute the ratio:
       W_eff(x) = - D_x psi(x) / psi(x)
   and test if W_eff(x) = mu - chi^2 x^2.

3. Boundary Flatness and Infinite-Order Vanishing:
   We evaluate the higher even endpoint derivatives:
       D_k(N) = T_{v_N}^{(2k)}(0),   k = 0, 1, 2, 3, 4,
   at N = 8, 16, 24 to verify whether all derivatives simultaneously decay
   to zero, establishing infinite-order contact at the boundaries.

4. Analytical Calibration of the Constant kappa_c = lambda_min / A_0:
   In Cell 41, the ratio lambda_min(N) / A_0(N) froze at 0.00246 +- 0.0001.
   We test this against the arithmetic and geometric parameters of the Connes-CvS
   model:
       L = log(c),  beta = L/(4 pi),  C_c = L (sqrt(c) + 1/sqrt(c) - 2) / (2 pi^2).
"""

from __future__ import annotations

import mpmath as mp

from cell import get_ground_state


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

mp.mp.dps = 50

c = 13
L = mp.log(c)
T_ground = 400
GROUND_DPS = 50

# Benchmark dimension for continuum analysis
N_BENCHMARK = 24


# ---------------------------------------------------------------------------
# Trigonometric profile and derivatives
# ---------------------------------------------------------------------------

def T_eval(v, t, L):
    """T_v(t) = v_0 + sqrt(2) sum_{m=1}^N v_m cos(2 pi m t / L)."""
    kappa = 2 * mp.pi / L
    val = v[0]
    for m in range(1, len(v)):
        val += mp.sqrt(2) * v[m] * mp.cos(kappa * m * t)
    return val


def T_prime_eval(v, t, L):
    """T_v'(t) = -sqrt(2) kappa sum_{m=1}^N m v_m sin(2 pi m t / L)."""
    kappa = 2 * mp.pi / L
    val = mp.mpf(0)
    for m in range(1, len(v)):
        val -= mp.sqrt(2) * kappa * m * v[m] * mp.sin(kappa * m * t)
    return val


def T_double_prime_eval(v, t, L):
    """T_v''(t) = -sqrt(2) kappa^2 sum_{m=1}^N m^2 v_m cos(2 pi m t / L)."""
    kappa = 2 * mp.pi / L
    val = mp.mpf(0)
    for m in range(1, len(v)):
        val -= mp.sqrt(2) * (kappa ** 2) * (m ** 2) * v[m] * mp.cos(kappa * m * t)
    return val


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 78)
    print("CELL 43 — EFFECTIVE POTENTIAL & PROLATE STURM-LIOUVILLE RECONSTRUCTION")
    print("=" * 78)
    print(f"c = {c}, L = {mp.nstr(L, 20)}, T = {T_ground}, dps = {mp.mp.dps}")
    print(f"Benchmark N = {N_BENCHMARK}")

    lam, v24, _ = get_ground_state(
        c=c,
        N=N_BENCHMARK,
        T=T_ground,
        dps=GROUND_DPS,
        verbose=False,
    )

    kappa = 2 * mp.pi / L
    t_center = L / 2

    # -----------------------------------------------------------------------
    # 1. Effective Potential V_eff(t) = - T''(t) / T(t)
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("1. EFFECTIVE SCHRODINGER POTENTIAL: V_eff(t) = - T''(t) / T(t)")
    print("=" * 78)
    print(
        f"{'t / L':>8} "
        f"{'t - L/2':>12} "
        f"{'T(t)':>14} "
        f"{'T\'\'(t)':>14} "
        f"{'V_eff(t) = -T\'\'/T':>20}"
    )
    print("-" * 78)

    sample_fractions = [
        mp.mpf(k) / 20 for k in range(5, 16)  # 0.25 to 0.75 L
    ]

    v_eff_data = []

    for frac in sample_fractions:
        t_val = frac * L
        u_val = t_val - t_center  # distance from center
        T_val = T_eval(v24, t_val, L)
        T2_val = T_double_prime_eval(v24, t_val, L)
        v_eff = -T2_val / T_val

        v_eff_data.append((u_val, v_eff))

        print(
            f"{mp.nstr(frac, 4):>8} "
            f"{mp.nstr(u_val, 5):>12} "
            f"{mp.nstr(T_val, 7):>14} "
            f"{mp.nstr(T2_val, 7):>14} "
            f"{mp.nstr(v_eff, 8):>20}"
        )

    # Parabolic curvature check: V_eff(u) ~ V_0 + K * u^2
    u0, v0 = v_eff_data[len(v_eff_data) // 2]  # center (u=0)
    u1, v1 = v_eff_data[-1]                    # furthest sampled u
    K_est = (v1 - v0) / (u1 ** 2) if u1 != 0 else mp.mpf(0)

    print(f"\nPotential at center (t = L/2):     V_eff(0) = {mp.nstr(v0, 8)}")
    print(f"Estimated parabolic well curvature: K = d^2 V / du^2 ~ {mp.nstr(2 * K_est, 6)}")

    # -----------------------------------------------------------------------
    # 2. Normalized Prolate Operator: x = 2t/L - 1 in [-1, 1]
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("2. NORMALIZED PROLATE OPERATOR ON x in [-1, 1]")
    print("=" * 78)
    print("Equation: -(1 - x^2) psi''(x) + 2x psi'(x) = W_eff(x) psi(x)")
    print(
        f"{'x':>8} "
        f"{'psi(x)':>14} "
        f"{'psi\'(x)':>14} "
        f"{'psi\'\'(x)':>14} "
        f"{'W_eff(x)':>18}"
    )
    print("-" * 78)

    # In x coordinates:
    # x = 2t/L - 1 => dx/dt = 2/L => d/dt = (2/L) d/dx
    # psi'(x) = (L/2) T'(t)
    # psi''(x) = (L/2)^2 T''(t)

    x_samples = [mp.mpf(k) / 10 for k in range(-6, 7, 2)]  # -0.6 to +0.6

    w_eff_data = []

    for x_val in x_samples:
        t_val = (x_val + 1) * L / 2
        psi_val = T_eval(v24, t_val, L)
        psi_prime = (L / 2) * T_prime_eval(v24, t_val, L)
        psi_double = ((L / 2) ** 2) * T_double_prime_eval(v24, t_val, L)

        lhs = -(1 - x_val ** 2) * psi_double + 2 * x_val * psi_prime
        w_eff = lhs / psi_val if abs(psi_val) > 0 else mp.mpf(0)

        w_eff_data.append((x_val, w_eff))

        print(
            f"{mp.nstr(x_val, 3):>8} "
            f"{mp.nstr(psi_val, 6):>14} "
            f"{mp.nstr(psi_prime, 6):>14} "
            f"{mp.nstr(psi_double, 6):>14} "
            f"{mp.nstr(w_eff, 8):>18}"
        )

    # Test prolate quadratic fit: W_eff(x) = mu - chi^2 x^2
    w0 = [w for x, w in w_eff_data if x == 0][0]
    w_edge = w_eff_data[-1][1]
    x_edge = w_eff_data[-1][0]
    chi2_est = (w0 - w_edge) / (x_edge ** 2) if x_edge != 0 else mp.mpf(0)

    print(f"\nProlate eigenvalue mu = W_eff(0) ~ {mp.nstr(w0, 8)}")
    print(f"Prolate parameter chi^2 ~ {mp.nstr(chi2_est, 8)} (chi ~ {mp.nstr(mp.sqrt(abs(chi2_est)), 6)})")

    # -----------------------------------------------------------------------
    # 3. Higher Endpoint Derivatives (Infinite-Order Boundary Contact)
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("3. BOUNDARY JET D_k = T^{(2k)}(0) (INFINITE-ORDER VANISHING)")
    print("=" * 78)
    print(
        f"{'N':>3} "
        f"{'D_0 (value)':>16} "
        f"{'D_1 (2nd deriv)':>18} "
        f"{'D_2 (4th deriv)':>18} "
        f"{'D_3 (6th deriv)':>18}"
    )
    print("-" * 78)

    for N in [8, 16, 24]:
        _, v_N, _ = get_ground_state(c=c, N=N, T=T_ground, dps=GROUND_DPS, verbose=False)
        N_len = len(v_N) - 1

        D0 = v_N[0] + mp.sqrt(2) * sum(v_N[m] for m in range(1, N_len + 1))

        D_vals = [D0]
        for k_idx in [1, 2, 3]:
            M_2k = sum((mp.mpf(m) ** (2 * k_idx)) * v_N[m] for m in range(1, N_len + 1))
            D_k = mp.sqrt(2) * (-1) ** k_idx * (kappa ** (2 * k_idx)) * M_2k
            D_vals.append(D_k)

        print(
            f"{N:3d} "
            f"{mp.nstr(abs(D_vals[0]), 6):>16} "
            f"{mp.nstr(abs(D_vals[1]), 6):>18} "
            f"{mp.nstr(abs(D_vals[2]), 6):>18} "
            f"{mp.nstr(abs(D_vals[3]), 6):>18}"
        )

    # -----------------------------------------------------------------------
    # 4. Analytical Calibration of kappa_c = lambda_min / A_0
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("4. CALIBRATION OF THE RATIO kappa_c = lambda_min / A_0")
    print("=" * 78)

    # Arithmetic pole coefficient from Paper 2 (Groskin):
    # C_c = L * (sqrt(c) + 1/sqrt(c) - 2) / (2 * pi^2)
    C_c = L * (mp.sqrt(c) + 1 / mp.sqrt(c) - 2) / (2 * mp.pi ** 2)
    beta = L / (4 * mp.pi)
    rho = 2 * mp.pi / L

    lam24, v24, _ = get_ground_state(c=c, N=24, T=T_ground, dps=GROUND_DPS, verbose=False)
    D0_24 = v24[0] + mp.sqrt(2) * sum(v24[m] for m in range(1, 25))
    A0_24 = (2 / L) * D0_24 ** 2
    kappa_c_num = lam24 / A0_24

    print(f"Geometric parameter L = log(c) =       {mp.nstr(L, 8)}")
    print(f"Pole source scale C_c =                 {mp.nstr(C_c, 8)}")
    print(f"Band scale beta = L / (4 pi) =          {mp.nstr(beta, 8)}")
    print(f"Lattice spacing rho = 2 pi / L =        {mp.nstr(rho, 8)}")
    print(f"Numerical ratio kappa_c =               {mp.nstr(kappa_c_num, 8)}")
    print()
    print("Candidate analytical comparisons for kappa_c:")
    print(f"  C_c / 100               = {mp.nstr(C_c / 100, 8)}   (ratio = {mp.nstr(kappa_c_num / (C_c / 100), 5)})")
    print(f"  (C_c * beta^2) / 2      = {mp.nstr((C_c * beta**2) / 2, 8)}   (ratio = {mp.nstr(kappa_c_num / ((C_c * beta**2)/2), 5)})")
    print(f"  beta^3 / pi             = {mp.nstr((beta**3) / mp.pi, 8)}   (ratio = {mp.nstr(kappa_c_num / (beta**3 / mp.pi), 5)})")
    print(f"  C_c / (4 * pi^2)        = {mp.nstr(C_c / (4 * mp.pi**2), 8)}   (ratio = {mp.nstr(kappa_c_num / (C_c / (4 * mp.pi**2)), 5)})")

    print("\n" + "=" * 78)
    print("END OF CELL 43")
    print("=" * 78)
    print(
        "Conclusions to review in cell43.out:\n"
        "  1. Is V_eff(t) locally parabolic, confirming prolate/harmonic confinement?\n"
        "  2. Do higher endpoint derivatives D_k decay simultaneously to zero?\n"
        "  3. Which analytical combination matches the universal ratio kappa_c?"
    )
