"""
CELL 47 — MULTI-c SCALING OF THE WEIL GROUND STATE, WKB TUNNELING, AND ARITHMETIC ENERGY DISTRIBUTION

In Cells 41–46, for the canonical scaling cutoff c = 13, we established:
1. Geometric eigenvalue scaling: lambda_min(N) ~ kappa_c * c^{-N}.
2. WKB quantum tunneling law: S_WKB quantitatively predicts boundary suppression within 5.6%.
3. Super-polynomial resolvent decay and tail hierarchy extinction: A_k -> 0.
4. Exact tri-partite zero-energy balance: Q_pole(infty) + Q_prime(infty) + Q_arch(infty) = 0.

Cell 47 tests the universality and scaling of these fundamental laws across multiple cutoffs:
    c in {5, 7, 11, 13, 17}

Objectives:
1. Multi-c Ground-State Eigenvalue Scaling Law:
   For each cutoff c, evaluate lambda_min(N; c) across N in {4, 8, 12, 16, 20} to test:
       lambda_min(N; c) ~ kappa_c(c) * c^{-N}.
   Verify that the effective decay base b_eff = (lambda_min(N-4) / lambda_min(N))^{1/4} converges to c.

2. Scaling of the Universal Ratio kappa_c(c) = lambda_min / A_0:
   Track kappa_c(N) = lambda_min(N) / A_0(N) as N increases, and calibrate the stabilized limit
   kappa_c(c) against candidate analytical models:
       - C_c / 100, where C_c = L (sqrt(c) + 1/sqrt(c) - 2) / (2 pi^2)
       - beta^3 / pi, where beta = L / (4 pi)
       - C_c * beta^2 / 2
       - C_c / (4 pi^2)

3. Multi-c WKB Quantum Barrier Penetration:
   For each cutoff c at N = 20:
   Locate the classical turning point t_turn(c) where T''(t_turn) = 0 and evaluate the WKB action:
       S_WKB(c) = int_0^{t_turn} sqrt(T''(t) / T(t)) dt.
   Compare S_WKB(c) against the actual boundary suppression log(T_max / T(0)) across all c.
   Test the scaling hypothesis S_WKB(c) ~ (N/2) log(c) = (N/2) L.

4. Multi-c Tri-Partite Weil Form Equilibrium & Prime Energy Share:
   For each cutoff c at N = 20, compute the three arithmetic pieces:
       Q_pole(c), Q_prime(c), Q_arch(c).
   Verify that Q_total(c) = Q_pole + Q_prime + Q_arch matches lambda_min(N=20; c).
   Determine the energy partition fractions:
       f_prime(c) = |Q_prime(c)| / Q_pole(c),
       f_arch(c)  = |Q_arch(c)|  / Q_pole(c).
   Trace how the negative energy burden shifts from the continuous Archimedean place
   to the discrete prime powers as c grows.
"""

from __future__ import annotations

import mpmath as mp

from cell import (
    canonical_to_full,
    get_ground_state,
    h_plus,
    prime_power_terms,
)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

mp.mp.dps = 50

CUTOFFS = [5, 7, 11, 13, 17]
N_LIST = [4, 8, 12, 16, 20]
T_GROUND = 400
GROUND_DPS = 50


# ---------------------------------------------------------------------------
# Trigonometric profile and derivatives
# ---------------------------------------------------------------------------

def T_eval(v, t, L):
    kappa = 2 * mp.pi / L
    val = v[0]
    for m in range(1, len(v)):
        val += mp.sqrt(2) * v[m] * mp.cos(kappa * m * t)
    return val


def T_double_prime_eval(v, t, L):
    kappa = 2 * mp.pi / L
    val = mp.mpf(0)
    for m in range(1, len(v)):
        val -= mp.sqrt(2) * (kappa * m) ** 2 * v[m] * mp.cos(kappa * m * t)
    return val


# ---------------------------------------------------------------------------
# Entire amplitude Phi_v(r)
# ---------------------------------------------------------------------------

def Phi_eval(v, r, L):
    if r == 0:
        return v[0] * mp.sqrt(L)
    sin_term = mp.sin(r * L / 2)
    kappa = 2 * mp.pi / L
    sum_m = mp.mpf(0)
    for m in range(1, len(v)):
        am = kappa * m
        denom = r ** 2 - am ** 2
        sum_m += mp.sqrt(2) * v[m] * r * sin_term / denom
    return (2 / mp.sqrt(L)) * (v[0] * sin_term / r + sum_m)


# ---------------------------------------------------------------------------
# Prime and Pole Matrix Builders
# ---------------------------------------------------------------------------

def psi_prime_val(x, q, Lambda_q, L):
    a = 1 - mp.log(q) / L
    prefactor = -1 / mp.pi * Lambda_q / mp.sqrt(q)
    return prefactor * mp.sin(2 * mp.pi * x * a)


def psi_prime_deriv_val(x, q, Lambda_q, L):
    a = 1 - mp.log(q) / L
    prefactor = -1 / mp.pi * Lambda_q / mp.sqrt(q)
    return prefactor * 2 * mp.pi * a * mp.cos(2 * mp.pi * x * a)


def build_prime_matrix(N, c, L):
    size = 2 * N + 1
    Q_prime = mp.matrix(size, size)
    terms = prime_power_terms(c)

    for q, Lambda_q in terms:
        vals = [psi_prime_val(x, q, Lambda_q, L) for x in range(-N, N + 1)]
        ders = [psi_prime_deriv_val(x, q, Lambda_q, L) for x in range(-N, N + 1)]

        for i, m in enumerate(range(-N, N + 1)):
            for j, n in enumerate(range(-N, N + 1)):
                if m != n:
                    Q_prime[i, j] += (vals[i] - vals[j]) / mp.mpf(m - n)
                else:
                    Q_prime[i, j] += ders[i]
    return Q_prime


def psi_pole_val(x, L):
    if x == 0:
        return mp.mpf(0)
    integrand = lambda y: 2 * mp.cosh(y / 2) * mp.sin(2 * mp.pi * x * (1 - y / L))
    return (1 / mp.pi) * mp.quad(integrand, [0, L])


def psi_pole_deriv_val(x, L):
    integrand = lambda y: 2 * mp.cosh(y / 2) * (2 * mp.pi * (1 - y / L)) * mp.cos(2 * mp.pi * x * (1 - y / L))
    return (1 / mp.pi) * mp.quad(integrand, [0, L])


def build_pole_matrix(N, L):
    size = 2 * N + 1
    Q_pole = mp.matrix(size, size)
    vals = [psi_pole_val(x, L) for x in range(-N, N + 1)]
    ders = [psi_pole_deriv_val(x, L) for x in range(-N, N + 1)]

    for i, m in enumerate(range(-N, N + 1)):
        for j, n in enumerate(range(-N, N + 1)):
            if m != n:
                Q_pole[i, j] = (vals[i] - vals[j]) / mp.mpf(m - n)
            else:
                Q_pole[i, j] = ders[i]
    return Q_pole


# ---------------------------------------------------------------------------
# Main Analysis
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 78)
    print("CELL 47 — MULTI-c SCALING OF THE WEIL GROUND STATE, WKB TUNNELING,")
    print("         AND ARITHMETIC ENERGY DISTRIBUTION")
    print("=" * 78)
    print(f"Cutoffs: {CUTOFFS}")
    print(f"Dimensions: {N_LIST}")
    print(f"T = {T_GROUND}, dps = {mp.mp.dps}")

    # Data store: ground_states[c][N] = (lambda_min, vec)
    ground_states = {}

    # -----------------------------------------------------------------------
    # 1. Multi-c Ground-State Eigenvalue Scaling Law
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("1. MULTI-c GROUND-STATE EIGENVALUE SCALING LAW: lambda_min(N; c) ~ c^{-N}")
    print("=" * 78)

    stabilized_kappa = {}

    for c_val in CUTOFFS:
        L_val = mp.log(c_val)
        ground_states[c_val] = {}
        print(f"\n--- Cutoff c = {c_val} (L = {mp.nstr(L_val, 8)}) ---")
        print(
            f"{'N':>3} "
            f"{'lambda_min(N)':>20} "
            f"{'b_eff':>12} "
            f"{'|T(0)|':>16} "
            f"{'A_0(N)':>16} "
            f"{'kappa_c(N)':>14}"
        )
        print("-" * 85)

        prev_lam = None
        for N in N_LIST:
            lam, vec, _ = get_ground_state(
                c=c_val,
                N=N,
                T=T_GROUND,
                dps=GROUND_DPS,
                verbose=False,
            )
            ground_states[c_val][N] = (lam, vec)

            # Effective base b_eff = (lam_{N-4} / lam_N)^(1/4)
            b_eff_str = "---"
            if prev_lam is not None and lam > 0:
                b_eff = (prev_lam / lam) ** mp.mpf("0.25")
                b_eff_str = mp.nstr(b_eff, 6)
            prev_lam = lam

            # Endpoint values
            T_0 = T_eval(vec, mp.mpf(0), L_val)
            A_0 = (2 / L_val) * T_0 ** 2
            kappa_N = lam / A_0 if A_0 > 0 else mp.mpf(0)

            print(
                f"{N:3d} "
                f"{mp.nstr(lam, 12):>20} "
                f"{b_eff_str:>12} "
                f"{mp.nstr(abs(T_0), 8):>16} "
                f"{mp.nstr(A_0, 8):>16} "
                f"{mp.nstr(kappa_N, 8):>14}"
            )

        stabilized_kappa[c_val] = kappa_N

    # -----------------------------------------------------------------------
    # 2. Scaling of the Universal Ratio kappa_c(c)
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("2. SCALING OF THE UNIVERSAL RATIO kappa_c(c) = lambda_min / A_0")
    print("=" * 78)
    print(
        f"{'c':>3} "
        f"{'L':>8} "
        f"{'kappa_c':>12} "
        f"{'C_c/100':>12} "
        f"{'ratio(C_c/100)':>16} "
        f"{'beta^3/pi':>12} "
        f"{'ratio(beta^3)':>14}"
    )
    print("-" * 82)

    for c_val in CUTOFFS:
        L_val = mp.log(c_val)
        C_c = L_val * (mp.sqrt(c_val) + 1 / mp.sqrt(c_val) - 2) / (2 * mp.pi ** 2)
        beta = L_val / (4 * mp.pi)
        model_C = C_c / 100
        model_beta = (beta ** 3) / mp.pi

        k_c = stabilized_kappa[c_val]
        r_C = k_c / model_C if model_C > 0 else mp.mpf(0)
        r_beta = k_c / model_beta if model_beta > 0 else mp.mpf(0)

        print(
            f"{c_val:3d} "
            f"{mp.nstr(L_val, 6):>8} "
            f"{mp.nstr(k_c, 6):>12} "
            f"{mp.nstr(model_C, 6):>12} "
            f"{mp.nstr(r_C, 6):>16} "
            f"{mp.nstr(model_beta, 6):>12} "
            f"{mp.nstr(r_beta, 6):>14}"
        )

    # -----------------------------------------------------------------------
    # 3. Multi-c WKB Quantum Tunneling Barrier Action
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("3. MULTI-c WKB QUANTUM TUNNELING BARRIER PENETRATION (N = 20)")
    print("=" * 78)
    print(
        f"{'c':>3} "
        f"{'t_turn/L':>10} "
        f"{'S_WKB':>14} "
        f"{'Actual Suppr':>16} "
        f"{'Suppr/S_WKB':>14} "
        f"{'S_WKB/L':>12} "
        f"{'S_WKB/(NL/2)':>14}"
    )
    print("-" * 88)

    wkb_actions = {}

    for c_val in CUTOFFS:
        L_val = mp.log(c_val)
        _, v20 = ground_states[c_val][20]

        # Scan s in [0.10, 0.49] with step 0.01 to find sign change of T''
        a_bracket, b_bracket = None, None
        s_prev = mp.mpf("0.10")
        f_prev = T_double_prime_eval(v20, s_prev * L_val, L_val)

        for step in range(11, 50):
            s_curr = mp.mpf(step) / 100
            f_curr = T_double_prime_eval(v20, s_curr * L_val, L_val)
            if (f_prev > 0 and f_curr < 0) or (f_prev < 0 and f_curr > 0):
                a_bracket = s_prev * L_val
                b_bracket = s_curr * L_val
                break
            s_prev = s_curr
            f_prev = f_curr

        if a_bracket is None:
            # Fallback default bracket
            a_bracket = mp.mpf("0.35") * L_val
            b_bracket = mp.mpf("0.45") * L_val

        # Robust 250-iteration bisection
        fa = T_double_prime_eval(v20, a_bracket, L_val)
        a_cur, b_cur = a_bracket, b_bracket
        for _ in range(250):
            mid = (a_cur + b_cur) / 2
            fmid = T_double_prime_eval(v20, mid, L_val)
            if (fa > 0 and fmid > 0) or (fa < 0 and fmid < 0):
                a_cur = mid
                fa = fmid
            else:
                b_cur = mid
        t_turn = (a_cur + b_cur) / 2

        # WKB integrand: sqrt(T''(t) / T(t))
        def wkb_integrand(t):
            T_val = T_eval(v20, t, L_val)
            T2_val = T_double_prime_eval(v20, t, L_val)
            ratio = T2_val / T_val if T_val > 0 else mp.mpf(0)
            return mp.sqrt(max(mp.mpf(0), ratio))

        S_wkb = mp.quad(wkb_integrand, [0, t_turn])
        wkb_actions[c_val] = S_wkb

        T_peak = T_eval(v20, L_val / 2, L_val)
        T_edge = T_eval(v20, mp.mpf(0), L_val)
        actual_suppr = mp.log(T_peak / abs(T_edge))

        acc_ratio = actual_suppr / S_wkb if S_wkb > 0 else mp.mpf(0)
        s_per_L = S_wkb / L_val
        s_per_half_NL = S_wkb / (10 * L_val)  # N/2 * L = 10 L

        print(
            f"{c_val:3d} "
            f"{mp.nstr(t_turn / L_val, 6):>10} "
            f"{mp.nstr(S_wkb, 8):>14} "
            f"{mp.nstr(actual_suppr, 8):>16} "
            f"{mp.nstr(acc_ratio, 6):>14} "
            f"{mp.nstr(s_per_L, 6):>12} "
            f"{mp.nstr(s_per_half_NL, 6):>14}"
        )

    # -----------------------------------------------------------------------
    # 4. Multi-c Tri-Partite Weil Form Equilibrium & Prime Energy Share
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("4. MULTI-c TRI-PARTITE WEIL ENERGY DISTRIBUTION (N = 20)")
    print("=" * 78)
    print(
        f"{'c':>3} "
        f"{'Q_pole':>14} "
        f"{'Q_prime':>14} "
        f"{'Q_arch':>14} "
        f"{'Q_total':>16} "
        f"{'f_prime (%)':>14} "
        f"{'f_arch (%)':>14}"
    )
    print("-" * 92)

    for c_val in CUTOFFS:
        L_val = mp.log(c_val)
        lam20, v20 = ground_states[c_val][20]
        u20 = canonical_to_full(v20)

        Q_pr = build_prime_matrix(20, c_val, L_val)
        Q_po = build_pole_matrix(20, L_val)

        pole_val = mp.fdot(u20, Q_po * u20)
        prime_val = mp.fdot(u20, Q_pr * u20)

        arch_val = (1 / mp.pi) * mp.quad(
            lambda r: h_plus(r) * Phi_eval(v20, r, L_val) ** 2,
            [0, 80],
        )

        total_val = pole_val + prime_val + arch_val
        f_prime = 100 * abs(prime_val) / pole_val
        f_arch = 100 * abs(arch_val) / pole_val

        print(
            f"{c_val:3d} "
            f"{mp.nstr(pole_val, 8):>14} "
            f"{mp.nstr(prime_val, 8):>14} "
            f"{mp.nstr(arch_val, 8):>14} "
            f"{mp.nstr(total_val, 6):>16} "
            f"{mp.nstr(f_prime, 5):>14} "
            f"{mp.nstr(f_arch, 5):>14}"
        )

    print("\n" + "=" * 78)
    print("END OF CELL 47")
    print("=" * 78)
    print(
        "Conclusions to review in cell47.out:\n"
        "  1. Does b_eff converge to c for every cutoff c in {5, 7, 11, 13, 17}?\n"
        "  2. How does kappa_c scale with c and L = log(c)?\n"
        "  3. Does WKB barrier action S_WKB predict actual boundary decay within 5-10% for all c?\n"
        "  4. How does the prime energy share f_prime evolve as c increases from 5 to 17?"
    )
