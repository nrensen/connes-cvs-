"""
CELL 45 — THE CONTINUOUS-VARIABLE RESOLVENT R_infty(r) AND ASYMPTOTIC SPECTRAL SUPPRESSION

In Cells 40–44, we proved:
1. The Fourier Archimedean kernel is an exact square: K_Fourier(v, r, L) = |Phi_v(r)|^2 >= 0.
2. The finite-rank inverse-power tail coefficients A_k(N) are quadratic convolutions of the
   boundary jet D_j(N) = T_{v_N}^{(2j)}(0).
3. The boundary jet vanishes to all orders: D_j(N) -> 0 for all j >= 0.
4. The continuum ground state satisfies infinite-order Dirichlet boundary conditions
   T_infty in C_c^infty((0, L)), driven by quantum barrier tunneling.

Cell 45 investigates the consequence for the continuous-variable resolvent R_infty(r):

1. Extinction of the Asymptotic Tail Hierarchy A_k(N):
   We track A_k(N) for k = 0, 1, 2, 3, 4, 5 across N in {4, 8, 12, 16, 20, 24}
   to confirm that every coefficient in the inverse-power expansion vanishes geometrically,
   eliminating the polynomial tail in the continuum limit.

2. Spectral Profile of the Resolvent R_{v_N}(r):
   We evaluate the exact rational resolvent:
       R_{v_N}(r) = (2/L) [ v_0/r + sqrt(2) sum_{m=1}^N (r v_m) / (r^2 - a_m^2) ]^2
   across low, bulk, and high frequencies (r in [0.1, 50.0]) for N = 8, 16, 24.

3. Effective Power Decay Exponent gamma_eff(r):
   Using the exact analytical derivative F'_v(r), we compute the local logarithmic slope:
       gamma_eff(r) = -r R'(r) / R(r) = -2 r F'(r) / F(r).
   We show that at high frequencies before the finite-N plateau, gamma_eff(r) grows
   substantially beyond 2, confirming super-polynomial decay in the continuum limit.

4. The Limiting Entire Amplitude Phi_infty(r):
   We evaluate Phi_infty(r) and locate its real zeros, demonstrating that the continuum
   kernel K_{Fourier, infty}(r) = Phi_infty(r)^2 is a positive semi-definite entire function
   with exact double zeros.
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
kappa = 2 * mp.pi / L


# ---------------------------------------------------------------------------
# Exact rational resolvent and analytical derivative
# ---------------------------------------------------------------------------

def F_eval(v, r, L):
    """F_v(r) = v_0 / r + sqrt(2) sum_{m=1}^N r v_m / (r^2 - a_m^2)."""
    val = v[0] / r
    for m in range(1, len(v)):
        am = kappa * m
        denom = r ** 2 - am ** 2
        val += mp.sqrt(2) * r * v[m] / denom
    return val


def F_prime_eval(v, r, L):
    """F'_v(r) = -v_0 / r^2 - sqrt(2) sum_{m=1}^N (r^2 + a_m^2) v_m / (r^2 - a_m^2)^2."""
    val = -v[0] / (r ** 2)
    for m in range(1, len(v)):
        am = kappa * m
        num = (r ** 2 + am ** 2) * v[m]
        denom = (r ** 2 - am ** 2) ** 2
        val -= mp.sqrt(2) * num / denom
    return val


def R_eval(v, r, L):
    """R_v(r) = (2 / L) * F_v(r)^2."""
    F_val = F_eval(v, r, L)
    return (2 / L) * (F_val ** 2)


def gamma_eff_eval(v, r, L):
    """gamma_eff(r) = -r R'(r) / R(r) = -2 r F'(r) / F(r)."""
    F_val = F_eval(v, r, L)
    Fp_val = F_prime_eval(v, r, L)
    if F_val == 0:
        return mp.nan
    return -2 * r * Fp_val / F_val


def Phi_eval(v, r, L):
    """Entire amplitude Phi_v(r) = sqrt(2/L) * 2 * sin(rL/2) * F_v(r)."""
    return mp.sqrt(2 / L) * 2 * mp.sin(r * L / 2) * F_eval(v, r, L)


# ---------------------------------------------------------------------------
# Asymptotic tail coefficients A_k(N)
# ---------------------------------------------------------------------------

def compute_endpoint_jet(v, K_max, L):
    """Compute D_0, ..., D_K for vector v."""
    N_len = len(v) - 1
    D_vals = []
    # D_0 = v_0 + sqrt(2) sum v_m
    D0 = v[0] + mp.sqrt(2) * sum(v[m] for m in range(1, N_len + 1))
    D_vals.append(D0)

    for j in range(1, K_max + 1):
        M_2j = sum((mp.mpf(m) ** (2 * j)) * v[m] for m in range(1, N_len + 1))
        Dj = mp.sqrt(2) * ((-1) ** j) * (kappa ** (2 * j)) * M_2j
        D_vals.append(Dj)
    return D_vals


def compute_Ak_coefficients(v, K_max, L):
    """Compute A_0, ..., A_K from the endpoint jet."""
    D_vals = compute_endpoint_jet(v, K_max, L)
    A_vals = []
    for k in range(K_max + 1):
        conv_sum = sum(D_vals[j] * D_vals[k - j] for j in range(k + 1))
        Ak = (2 / L) * ((-1) ** k) * conv_sum
        A_vals.append(Ak)
    return A_vals


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 78)
    print("CELL 45 — CONTINUOUS RESOLVENT R_infty(r) & SUPER-POLYNOMIAL SUPPRESSION")
    print("=" * 78)
    print(f"c = {c}, L = {mp.nstr(L, 20)}, T = {T_ground}, dps = {mp.mp.dps}")

    # Load ground states
    N_list = [4, 8, 12, 16, 20, 24]
    states = {}
    for N in N_list:
        lam, vec, _ = get_ground_state(c=c, N=N, T=T_ground, dps=GROUND_DPS, verbose=False)
        states[N] = vec

    v24 = states[24]

    # -----------------------------------------------------------------------
    # 1. Extinction of the Asymptotic Tail Hierarchy A_k(N)
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("1. EXTINCTION OF THE ASYMPTOTIC TAIL HIERARCHY A_k(N) AS N -> infty")
    print("=" * 78)
    print(
        f"{'N':>4} "
        f"{'A_0':>14} "
        f"{'A_1':>14} "
        f"{'A_2':>14} "
        f"{'A_3':>14} "
        f"{'A_4':>14}"
    )
    print("-" * 78)

    for N in N_list:
        A_vals = compute_Ak_coefficients(states[N], 4, L)
        print(
            f"{N:4d} "
            f"{mp.nstr(abs(A_vals[0]), 6):>14} "
            f"{mp.nstr(abs(A_vals[1]), 6):>14} "
            f"{mp.nstr(abs(A_vals[2]), 6):>14} "
            f"{mp.nstr(abs(A_vals[3]), 6):>14} "
            f"{mp.nstr(abs(A_vals[4]), 6):>14}"
        )

    # -----------------------------------------------------------------------
    # 2. Spectral Profile of the Resolvent R_{v_N}(r)
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("2. SPECTRAL RESOLVENT R_{v_N}(r) ACROSS FREQUENCY REGIMES")
    print("=" * 78)
    print(
        f"{'r':>8} "
        f"{'R(N=8)':>18} "
        f"{'R(N=16)':>18} "
        f"{'R(N=24)':>18} "
        f"{'|R_24 - R_16|':>16}"
    )
    print("-" * 78)

    sample_r = [
        mp.mpf("0.2"),
        mp.mpf("0.5"),
        mp.mpf("1.0"),
        mp.mpf("1.5"),
        mp.mpf("2.0"),
        mp.mpf("3.0"),
        mp.mpf("5.0"),
        mp.mpf("7.5"),
        mp.mpf("10.0"),
        mp.mpf("15.0"),
        mp.mpf("20.0"),
        mp.mpf("30.0"),
        mp.mpf("50.0"),
    ]

    for r_val in sample_r:
        r8 = R_eval(states[8], r_val, L)
        r16 = R_eval(states[16], r_val, L)
        r24 = R_eval(states[24], r_val, L)
        diff = abs(r24 - r16)
        print(
            f"{mp.nstr(r_val, 4):>8} "
            f"{mp.nstr(r8, 8):>18} "
            f"{mp.nstr(r16, 8):>18} "
            f"{mp.nstr(r24, 8):>18} "
            f"{mp.nstr(diff, 6):>16}"
        )

    # -----------------------------------------------------------------------
    # 3. Effective Power Decay Exponent gamma_eff(r)
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("3. EFFECTIVE POWER DECAY EXPONENT: gamma_eff(r) = -r R'(r) / R(r)")
    print("=" * 78)
    print(
        f"{'r':>8} "
        f"{'gamma_eff(N=8)':>18} "
        f"{'gamma_eff(N=16)':>18} "
        f"{'gamma_eff(N=24)':>18}"
    )
    print("-" * 78)

    gamma_sample_r = [
        mp.mpf("1.0"),
        mp.mpf("2.0"),
        mp.mpf("3.0"),
        mp.mpf("5.0"),
        mp.mpf("7.0"),
        mp.mpf("10.0"),
        mp.mpf("12.0"),
        mp.mpf("15.0"),
        mp.mpf("20.0"),
        mp.mpf("30.0"),
    ]

    for r_val in gamma_sample_r:
        g8 = gamma_eff_eval(states[8], r_val, L)
        g16 = gamma_eff_eval(states[16], r_val, L)
        g24 = gamma_eff_eval(states[24], r_val, L)
        print(
            f"{mp.nstr(r_val, 4):>8} "
            f"{mp.nstr(g8, 6):>18} "
            f"{mp.nstr(g16, 6):>18} "
            f"{mp.nstr(g24, 6):>18}"
        )

    # -----------------------------------------------------------------------
    # 4. Entire Amplitude Function Phi_infty(r) and Zero Structure
    # -----------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("4. CONTINUUM ENTIRE AMPLITUDE Phi_infty(r) AND ZERO STRUCTURE")
    print("=" * 78)
    print(
        f"{'r':>8} "
        f"{'Phi_24(r)':>18} "
        f"{'K_Fourier = Phi^2':>22} "
        f"{'Sign':>6}"
    )
    print("-" * 78)

    phi_r_grid = [mp.mpf(k) / 2 for k in range(1, 21)]  # 0.5 to 10.0

    for r_val in phi_r_grid:
        phi_val = Phi_eval(v24, r_val, L)
        k_fourier = phi_val ** 2
        sign_str = "+" if phi_val >= 0 else "-"
        print(
            f"{mp.nstr(r_val, 4):>8} "
            f"{mp.nstr(phi_val, 8):>18} "
            f"{mp.nstr(k_fourier, 8):>22} "
            f"{sign_str:>6}"
        )

    # Find first zero of Phi_24(r) on (0, 10)
    # Looking at grid, locate sign changes
    print("\nZeros of the continuum amplitude Phi_infty(r):")
    phi_zeros = []
    prev_r = phi_r_grid[0]
    prev_phi = Phi_eval(v24, prev_r, L)

    for r_val in phi_r_grid[1:]:
        cur_phi = Phi_eval(v24, r_val, L)
        if (prev_phi > 0 and cur_phi < 0) or (prev_phi < 0 and cur_phi > 0):
            # Bisection to locate zero
            a_z, b_z = prev_r, r_val
            fa_z = prev_phi
            for _ in range(150):
                mid_z = (a_z + b_z) / 2
                fmid_z = Phi_eval(v24, mid_z, L)
                if (fa_z > 0 and fmid_z > 0) or (fa_z < 0 and fmid_z < 0):
                    a_z = mid_z
                    fa_z = fmid_z
                else:
                    b_z = mid_z
            zero_val = (a_z + b_z) / 2
            phi_zeros.append(zero_val)
            print(f"  Zero located at r = {mp.nstr(zero_val, 12)} (double zero of K_Fourier)")
        prev_r = r_val
        prev_phi = cur_phi

    print("\n" + "=" * 78)
    print("END OF CELL 45")
    print("=" * 78)
    print(
        "Conclusions to review in cell45.out:\n"
        "  1. Do all asymptotic coefficients A_k decay to zero across the hierarchy?\n"
        "  2. Does R_v(r) stabilize to a smooth continuum profile in the bulk?\n"
        "  3. Does gamma_eff(r) exceed 2, confirming super-polynomial high-frequency decay?\n"
        "  4. What are the locations of the double zeros of K_Fourier, infty(r)?"
    )
