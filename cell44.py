"""
CELL 44 — WKB TUNNELING BARRIER, EXACT LEGENDRE SPECTRUM, AND SLEPIAN PROLATE COUPLING

Cell 43 established three fundamental properties of the continuum ground state:
1. An effective confining Schrödinger well V_conf(t) = E - V_eff(t) with minimum at t = L/2.
2. Infinite-order flat boundary contact: T^{(k)}(0) = T^{(k)}(L) = 0 for all k >= 0.
3. Universal eigenvalue proportionality lambda_min(N) ~ kappa_c * c^{-N}.

Cell 44 tests the physical and mathematical mechanism of this confinement:

1. Log-Barrier Potential & Boundary Divergence:
   We evaluate S(t) = -log(T(t)) and compute the effective divergence index
   p_eff(t) = -t S'(t) / S(t) = d log S / d log(1/t) as t -> 0.

2. WKB Quantum Tunneling Barrier Penetration:
   In quantum tunneling, boundary decay is governed by the classically forbidden
   barrier action:
       S_WKB = int_0^{t_turn} sqrt(T''(t) / T(t)) dt.
   We test whether this WKB tunneling integral quantitatively accounts for the
   observed 20 orders of boundary suppression (log(T(L/2)/T(0)) ~ 46.55).

3. Exact Legendre Multipole Spectrum:
   On x = 2t/L - 1 in [-1, 1], the normalized wave psi(x) = T((x+1)L/2) is even.
   Its Legendre series psi(x) = sum_{k=0}^K c_{2k} P_{2k}(x) has exact closed-form
   coefficients via Bauer's spherical Bessel identity:
       c_0 = v_0 + sqrt(2) sum_{m=1}^N (-1)^m v_m j_0(pi m),
       c_{2k} = (4k + 1) sqrt(2) (-1)^k sum_{m=1}^N (-1)^m v_m j_{2k}(pi m).
   We evaluate the exact Legendre multipole spectrum c_{2k} to 50 decimal digits.

4. Slepian Prolate Recurrence Residual:
   For a true prolate spheroidal wave function S_{00}(c_0, x), the Legendre coefficients
   satisfy the classical Slepian-Bouwkamp three-term recurrence:
       A_k c_{2k+2} + (B_k - mu) c_{2k} + C_k c_{2k-2} = 0.
   We optimize the prolate parameter c_0 to minimize the spectral recurrence residual,
   testing the exact prolate spheroidal identity of the Connes-CvS ground state.
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
N_BENCHMARK = 24


# ---------------------------------------------------------------------------
# Spherical Bessel function j_n(z)
# ---------------------------------------------------------------------------

def spherical_jn(n, z):
    """j_n(z) = sqrt(pi / (2z)) J_{n + 1/2}(z)."""
    if z == 0:
        return mp.mpf(1) if n == 0 else mp.mpf(0)
    order = mp.mpf(n) + mp.mpf("0.5")
    return mp.sqrt(mp.pi / (2 * z)) * mp.besselj(order, z)


# ---------------------------------------------------------------------------
# Trigonometric profile and derivatives
# ---------------------------------------------------------------------------

def T_eval(v, t, L):
    kappa = 2 * mp.pi / L
    val = v[0]
    for m in range(1, len(v)):
        val += mp.sqrt(2) * v[m] * mp.cos(kappa * m * t)
    return val


def T_prime_eval(v, t, L):
    kappa = 2 * mp.pi / L
    val = mp.mpf(0)
    for m in range(1, len(v)):
        val -= mp.sqrt(2) * kappa * m * v[m] * mp.sin(kappa * m * t)
    return val


def T_double_prime_eval(v, t, L):
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
    print("CELL 44 — WKB TUNNELING BARRIER & EXACT SLEPIAN PROLATE SPECTRUM")
    print("=" * 78)
    print(f"c = {c}, L = {mp.nstr(L, 20)}, T = {T_ground}, dps = {mp.mp.dps}")
    print(f"Benchmark N = {N_BENCHMARK}")

    lam24, v24, _ = get_ground_state(
        c=c,
        N=N_BENCHMARK,
        T=T_ground,
        dps=GROUND_DPS,
        verbose=False,
    )

    t_center = L / 2

    # -----------------------------------------------------------------------
    # 1. Log-Barrier Potential S(t) = -log(T(t)) and Singularity Index
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("1. LOG-BARRIER POTENTIAL S(t) = -log(T(t)) AND DIVERGENCE INDEX")
    print("=" * 78)
    print(
        f"{'t / L':>8} "
        f"{'t':>10} "
        f"{'T(t)':>14} "
        f"{'S(t) = -log T':>16} "
        f"{'S\'(t) = -T\'/T':>16} "
        f"{'p_eff = -t S\'/S':>16}"
    )
    print("-" * 78)

    sample_fractions = [
        mp.mpf("0.005"),
        mp.mpf("0.01"),
        mp.mpf("0.02"),
        mp.mpf("0.05"),
        mp.mpf("0.10"),
        mp.mpf("0.15"),
        mp.mpf("0.20"),
        mp.mpf("0.30"),
        mp.mpf("0.40"),
        mp.mpf("0.50"),
    ]

    for frac in sample_fractions:
        t_val = frac * L
        T_val = T_eval(v24, t_val, L)
        Tp_val = T_prime_eval(v24, t_val, L)

        if T_val > 0:
            S_val = -mp.log(T_val)
            Sp_val = -Tp_val / T_val
            p_eff = -t_val * Sp_val / S_val if S_val != 0 else mp.mpf(0)
        else:
            S_val = mp.inf
            Sp_val = mp.inf
            p_eff = mp.nan

        print(
            f"{mp.nstr(frac, 4):>8} "
            f"{mp.nstr(t_val, 5):>10} "
            f"{mp.nstr(T_val, 7):>14} "
            f"{mp.nstr(S_val, 7):>16} "
            f"{mp.nstr(Sp_val, 7):>16} "
            f"{mp.nstr(p_eff, 6):>16}"
        )

    # -----------------------------------------------------------------------
    # 2. WKB Quantum Tunneling Barrier Action
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("2. WKB QUANTUM TUNNELING INTEGRAL")
    print("=" * 78)

    # Robust high-precision bisection for turning point T''(t_turn) = 0
    # From Cell 43: T''(0.4 L) = +5.38 > 0, T''(0.45 L) = -30.37 < 0
    a_turn = mp.mpf("0.40") * L
    b_turn = mp.mpf("0.45") * L
    fa_turn = T_double_prime_eval(v24, a_turn, L)
    for _ in range(250):
        mid_turn = (a_turn + b_turn) / 2
        fmid_turn = T_double_prime_eval(v24, mid_turn, L)
        if (fa_turn > 0 and fmid_turn > 0) or (fa_turn < 0 and fmid_turn < 0):
            a_turn = mid_turn
            fa_turn = fmid_turn
        else:
            b_turn = mid_turn
    t_turn = (a_turn + b_turn) / 2
    print(f"Classical turning point (inflection point): t_turn = {mp.nstr(t_turn, 8)} ({mp.nstr(t_turn / L, 5)} L)")

    # WKB integrand: k(s) = sqrt(T''(s) / T(s)) for s in [0, t_turn]
    def wkb_integrand(s):
        T_val = T_eval(v24, s, L)
        T2_val = T_double_prime_eval(v24, s, L)
        ratio = T2_val / T_val if T_val > 0 else mp.mpf(0)
        return mp.sqrt(max(mp.mpf(0), ratio))

    # Integrate from 0 to t_turn
    S_wkb = mp.quad(wkb_integrand, [0, t_turn])

    T_peak = T_eval(v24, t_center, L)
    T_edge = T_eval(v24, mp.mpf(0), L)
    actual_suppression = mp.log(T_peak / T_edge)

    print(f"Peak value T(L/2):                       {mp.nstr(T_peak, 8)}")
    print(f"Boundary value T(0):                     {mp.nstr(T_edge, 8)}")
    print(f"Actual suppression log(T_peak / T_edge): {mp.nstr(actual_suppression, 8)}")
    print(f"WKB tunneling action S_WKB:             {mp.nstr(S_wkb, 8)}")
    print(f"WKB predicted ratio exp(-S_WKB):        {mp.nstr(mp.exp(-S_wkb), 8)}")
    print(f"Actual boundary ratio T(0)/T_peak:      {mp.nstr(T_edge / T_peak, 8)}")
    print(f"WKB action ratio (actual / WKB):        {mp.nstr(actual_suppression / S_wkb, 6)}")

    # -----------------------------------------------------------------------
    # 3. Exact Legendre Multipole Spectrum of the Limiting Wave
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("3. EXACT LEGENDRE MULTIPOLE SPECTRUM (BAUER-BESSEL EXPANSION)")
    print("=" * 78)
    print("Decomposition: psi(x) = sum_{k=0}^K c_{2k} P_{2k}(x) on x in [-1, 1]")
    print(
        f"{'2k':>4} "
        f"{'Legendre coef c_{2k}':>24} "
        f"{'Fractional Energy c_{2k}^2/(4k+1)':>32}"
    )
    print("-" * 78)

    N_len = len(v24) - 1
    K_max = 10
    c_legendre = []

    # c_0 = v_0 + sqrt(2) sum_{m=1}^N (-1)^m v_m j_0(pi m)
    c0_val = v24[0]
    for m in range(1, N_len + 1):
        j0_val = spherical_jn(0, mp.pi * m)
        c0_val += mp.sqrt(2) * ((-1) ** m) * v24[m] * j0_val
    c_legendre.append(c0_val)

    # c_{2k} = (4k + 1) sqrt(2) (-1)^k sum_{m=1}^N (-1)^m v_m j_{2k}(pi m)
    for k in range(1, K_max + 1):
        deg = 2 * k
        sum_m = mp.mpf(0)
        for m in range(1, N_len + 1):
            j_val = spherical_jn(deg, mp.pi * m)
            sum_m += ((-1) ** m) * v24[m] * j_val
        ck_val = (2 * deg + 1) * mp.sqrt(2) * ((-1) ** k) * sum_m
        c_legendre.append(ck_val)

    # Print spectrum and energy fractions
    # Total energy: int_{-1}^1 psi^2 dx = sum_{k=0}^K 2/(4k+1) c_{2k}^2
    total_energy = sum(2 * (c_legendre[k] ** 2) / (4 * k + 1) for k in range(K_max + 1))

    for k in range(K_max + 1):
        deg = 2 * k
        val = c_legendre[k]
        energy_k = 2 * (val ** 2) / (2 * deg + 1)
        fraction = energy_k / total_energy
        print(
            f"{deg:4d} "
            f"{mp.nstr(val, 12):>24} "
            f"{mp.nstr(fraction, 8):>32}"
        )

    print(f"\nTotal reconstructed energy on [-1, 1]:  {mp.nstr(total_energy, 12)}")
    print(f"Expected L^2 norm (2/L) int_0^L T^2 dt: {mp.nstr(2, 12)}")

    # -----------------------------------------------------------------------
    # 4. Slepian Prolate Recurrence Residual Optimization
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("4. SLEPIAN-BOUWKAMP PROLATE RECURRENCE RESIDUAL OPTIMIZATION")
    print("=" * 78)
    print("Recurrence: A_k c_{2k+2} + (B_k(c_0) - mu) c_{2k} + C_k c_{2k-2} = 0")

    # For each candidate prolate bandwidth c_0, we can eliminate mu using k = 0:
    # A_0 c_2 + (B_0 - mu) c_0 = 0 => mu(c_0) = B_0(c_0) + A_0(c_0) * (c_2 / c_0)
    # Then evaluate the residual of the recurrence at k = 1, 2, 3...

    def recurrence_residual(c0_val):
        c0_sq = c0_val ** 2

        # k = 0:
        # A_0 = (1 * 2) / (3 * 5) * c0_sq = (2/15) c0_sq
        # B_0 = 0 + (-1) / ((-1) * 3) * c0_sq = (1/3) c0_sq
        A0 = (mp.mpf(2) / 15) * c0_sq
        B0 = (mp.mpf(1) / 3) * c0_sq
        mu = B0 + A0 * (c_legendre[1] / c_legendre[0])

        residuals = []
        for k in range(1, 6):
            deg = 2 * k
            Ak = mp.mpf((deg + 1) * (deg + 2)) / ((2 * deg + 3) * (2 * deg + 5)) * c0_sq
            Bk = mp.mpf(deg * (deg + 1)) + mp.mpf(2 * deg ** 2 + 2 * deg - 1) / ((2 * deg - 1) * (2 * deg + 3)) * c0_sq
            Ck = mp.mpf(deg * (deg - 1)) / ((2 * deg - 3) * (2 * deg - 1)) * c0_sq

            res = Ak * c_legendre[k + 1] + (Bk - mu) * c_legendre[k] + Ck * c_legendre[k - 1]
            # Normalize by c_{2k}
            rel_res = res / c_legendre[k] if c_legendre[k] != 0 else mp.mpf(0)
            residuals.append(rel_res)

        return mu, residuals

    print(
        f"{'c_0':>8} "
        f"{'c_0^2':>10} "
        f"{'mu(c_0)':>14} "
        f"{'Residual k=1':>16} "
        f"{'Residual k=2':>16} "
        f"{'Residual k=3':>16}"
    )
    print("-" * 78)

    candidate_c0 = [
        mp.mpf("2.0"),
        mp.mpf("3.0"),
        mp.mpf("4.0"),
        mp.mpf("5.0"),
        mp.mpf("6.0"),
        mp.mpf("7.0"),
        mp.mpf("8.0"),
        mp.mpf("9.0"),
    ]

    best_c0 = None
    min_norm = mp.inf

    for c0_test in candidate_c0:
        mu_test, res_list = recurrence_residual(c0_test)
        res_norm = mp.sqrt(sum(r ** 2 for r in res_list[:3]))
        if res_norm < min_norm:
            min_norm = res_norm
            best_c0 = c0_test

        print(
            f"{mp.nstr(c0_test, 3):>8} "
            f"{mp.nstr(c0_test**2, 4):>10} "
            f"{mp.nstr(mu_test, 6):>14} "
            f"{mp.nstr(res_list[0], 6):>16} "
            f"{mp.nstr(res_list[1], 6):>16} "
            f"{mp.nstr(res_list[2], 6):>16}"
        )

    print(f"\nBest discrete candidate: c_0 ~ {best_c0} (residual norm ~ {mp.nstr(min_norm, 5)})")

    print("\n" + "=" * 78)
    print("END OF CELL 44")
    print("=" * 78)
    print(
        "Conclusions to review in cell44.out:\n"
        "  1. How does the divergence index p_eff behave near the boundary t -> 0?\n"
        "  2. Does the WKB quantum tunneling action match the boundary suppression?\n"
        "  3. How rapidly does the Legendre multipole spectrum c_{2k} decay?\n"
        "  4. What is the optimal prolate bandwidth c_0 matching the Slepian recurrence?"
    )
