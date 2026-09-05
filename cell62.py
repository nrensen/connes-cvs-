#!/usr/bin/env python3
"""
================================================================================
CELL 62 — OPERATOR DECOMPOSITION OF THE FIRST-JET SCALE & EVEN RESOLVENT PROFILE
================================================================================

PURPOSE:
--------
Investigate the analytical operator mechanism governing the first-jet ratio
D_1 / D_0 using the exact representation of Theorem 7.3:
    D_1 / D_0 = kappa^2 [ <d, R_even s_2> - D_0^2 M_2 ]
where s_2 = K psi + M_1 d.

We decompose <d, R_even s_2> into two exact operator matrix elements:
    T_diag  = M_1 * <d, R_even d>
    T_cross = <d, R_even K psi>

and analyze:
1. The balance between T_diag and T_cross across N in {8, 12, 16, 20, 24}.
2. The mode-by-mode uniformity of the even overlap ratio d_k^2 / Delta_{even, k}.
3. The mode-by-mode uniformity of the odd overlap ratio a_j^2 / Delta_{odd, j}.
4. The spatial profile of the resolvent solution vector w_d = R_even d in c^perp.
5. The large-N scaling exponents of M_1, <d, R_even d>, and D_1 / D_0.

TERMINATION:
------------
Terminates strictly with: CELL 62 EXECUTION COMPLETE
================================================================================
"""

from __future__ import annotations

import time
import mpmath as mp

from connes_cvs import build_galerkin_matrix
from cell import (
    get_ground_state,
    canonical_to_full,
)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

mp.mp.dps = 50

C_PARAM = 13
L_PARAM = mp.log(C_PARAM)
T_PARAM = 400
GROUND_DPS = 50
KAPPA = 2 * mp.pi / L_PARAM

N_LIST = [8, 12, 16, 20, 24]


# -----------------------------------------------------------------------------
# Basis and Matrix Operations
# -----------------------------------------------------------------------------

def full_parity_basis(N: int) -> tuple[mp.matrix, mp.matrix]:
    """Construct orthonormal basis matrices E (even) and O (odd) for R^{2N+1}."""
    dim = 2 * N + 1
    E = mp.matrix(dim, N + 1)
    O = mp.matrix(dim, N)

    E[N, 0] = mp.mpf(1)
    for m in range(1, N + 1):
        inv_sqrt2 = 1 / mp.sqrt(2)
        E[N + m, m] = inv_sqrt2
        E[N - m, m] = inv_sqrt2
        O[N + m, m - 1] = inv_sqrt2
        O[N - m, m - 1] = -inv_sqrt2

    return E, O


def psi_vector(Q: mp.matrix, N: int) -> mp.matrix:
    """Extract symbol generator vector psi(m) = m * Q_{0, m} from row N of Q."""
    dim = 2 * N + 1
    psi = mp.matrix(dim, 1)
    for j in range(dim):
        m = j - N
        if m != 0:
            psi[j, 0] = mp.mpf(m) * Q[N, j]
        else:
            psi[j, 0] = mp.mpf(0)
    return psi


def constant_vector(N: int) -> mp.matrix:
    """Vector d = (1, 1, ..., 1)^T in R^{2N+1}."""
    return mp.matrix([[mp.mpf(1)] for _ in range(2 * N + 1)])


def K_apply(x: mp.matrix, N: int) -> mp.matrix:
    """Apply Fourier index operator K = diag(-N, ..., N) to column vector x."""
    dim = 2 * N + 1
    out = mp.matrix(dim, 1)
    for i in range(dim):
        out[i, 0] = mp.mpf(i - N) * x[i, 0]
    return out


def dot(x: mp.matrix, y: mp.matrix) -> mp.mpf:
    """Inner product <x, y> for column vectors."""
    return (x.T * y)[0, 0]


def norm2(x: mp.matrix) -> mp.mpf:
    """Squared Euclidean norm ||x||^2."""
    return dot(x, x)


def eigensystem_symmetric(A: mp.matrix) -> tuple[list[mp.mpf], mp.matrix]:
    """Symmetric eigendecomposition returning (sorted_eigenvalues, eigenvector_matrix)."""
    vals, vecs = mp.eigsy(A)
    return list(vals), vecs


# -----------------------------------------------------------------------------
# Main Execution Suite
# -----------------------------------------------------------------------------

def main():
    print("=" * 106)
    print("CELL 62 — OPERATOR DECOMPOSITION OF THE FIRST-JET SCALE & EVEN RESOLVENT PROFILE")
    print("=" * 106)
    print(f"Parameters: c = {C_PARAM}, L = {mp.nstr(L_PARAM, 15)}, T = {T_PARAM}, dps = {mp.mp.dps}")
    print("Analyzing the exact decomposition <d, R_even s_2> = T_cross + M_1 * <d, R_even d>.\n")

    results_across_N = {}

    for N in N_LIST:
        t0 = time.perf_counter()
        lam0, v_can, _ = get_ground_state(
            c=C_PARAM,
            N=N,
            T=T_PARAM,
            dps=GROUND_DPS,
            verbose=False,
        )
        c = canonical_to_full(v_can)
        Q = build_galerkin_matrix(
            c=C_PARAM,
            N=N,
            T=T_PARAM,
            dps=GROUND_DPS,
        )

        E, O = full_parity_basis(N)
        Q_even = E.T * Q * E
        Q_odd = O.T * Q * O

        even_vals, even_vecs = eigensystem_symmetric(Q_even)
        odd_vals, odd_vecs = eigensystem_symmetric(Q_odd)

        psi = psi_vector(Q, N)
        psi_odd = O.T * psi

        d = constant_vector(N)
        d_even = E.T * d
        D0 = dot(d, c)

        # First endpoint derivative D_1
        # In canonical coordinates, D_1 = - kappa^2 * sqrt(2) * sum_{m=1}^N m^2 v_m
        sum_m2_v = mp.mpf(0)
        for m in range(1, len(v_can)):
            sum_m2_v += (m ** 2) * v_can[m]
        D1 = - (KAPPA ** 2) * mp.sqrt(2) * sum_m2_v
        D1_over_D0 = D1 / D0

        # Odd sector resolvent sums: M1 and M2
        M1 = mp.mpf(0)
        M2 = mp.mpf(0)
        odd_mode_data = []

        for j in range(len(odd_vals)):
            e_j = odd_vecs[:, j]
            aj = dot(e_j, psi_odd)
            denom_j = odd_vals[j] - lam0
            term_M1 = (aj ** 2) / denom_j
            term_M2 = (aj ** 2) / (denom_j ** 2)
            M1 += term_M1
            M2 += term_M2
            odd_mode_data.append({
                "j": j + 1,
                "mu_j": odd_vals[j],
                "Delta_j": denom_j,
                "a_j": aj,
                "ratio_a2_Delta": term_M1,
                "cum_M1": M1,
            })

        # Second-jet source vector s_2 = K psi + M1 * d
        Kpsi = K_apply(psi, N)
        Kpsi_even = E.T * Kpsi
        s2 = Kpsi + M1 * d
        s2_even = E.T * s2

        # Resolvent decomposition on even subspace:
        # We project d, Kpsi, and s2 onto the excited even modes u_k (k >= 1)
        even_mode_data = []
        d_R_d = mp.mpf(0)
        d_R_Kpsi = mp.mpf(0)
        d_R_s2 = mp.mpf(0)

        for k in range(1, len(even_vals)):
            u_k = even_vecs[:, k]
            denom_k = even_vals[k] - lam0
            dk = dot(u_k, d_even)
            kpsi_k = dot(u_k, Kpsi_even)
            s2_k = dot(u_k, s2_even)

            term_dRd = (dk ** 2) / denom_k
            term_dRKpsi = (dk * kpsi_k) / denom_k
            term_dRs2 = (dk * s2_k) / denom_k

            d_R_d += term_dRd
            d_R_Kpsi += term_dRKpsi
            d_R_s2 += term_dRs2

            even_mode_data.append({
                "k": k,
                "E_k": even_vals[k],
                "Delta_k": denom_k,
                "d_k": dk,
                "kpsi_k": kpsi_k,
                "s2_k": s2_k,
                "ratio_d2_Delta": term_dRd,
                "cum_dRd": d_R_d,
                "term_dRs2": term_dRs2,
                "cum_dRs2": d_R_s2,
            })

        T_diag = M1 * d_R_d
        T_cross = d_R_Kpsi
        T_total = T_diag + T_cross

        subtraction_norm = (D0 ** 2) * M2
        resolvent_D1_over_D0 = (KAPPA ** 2) * (T_total - subtraction_norm)
        abs_error = abs(resolvent_D1_over_D0 - D1_over_D0)

        # Spatial resolvent vector w_d in c^perp:
        # w_d = sum_{k=1}^N (d_k / Delta_k) * u_k
        w_d = mp.matrix(len(even_vals), 1)
        for k in range(1, len(even_vals)):
            u_k = even_vecs[:, k]
            denom_k = even_vals[k] - lam0
            dk = dot(u_k, d_even)
            w_d += (dk / denom_k) * u_k

        elapsed = time.perf_counter() - t0

        results_across_N[N] = {
            "N": N,
            "lam0": lam0,
            "D0": D0,
            "D1": D1,
            "D1_over_D0": D1_over_D0,
            "M1": M1,
            "M2": M2,
            "subtraction_norm": subtraction_norm,
            "d_R_d": d_R_d,
            "T_diag": T_diag,
            "T_cross": T_cross,
            "T_total": T_total,
            "resolvent_D1_over_D0": resolvent_D1_over_D0,
            "abs_error": abs_error,
            "odd_mode_data": odd_mode_data,
            "even_mode_data": even_mode_data,
            "w_d": w_d,
            "elapsed": elapsed,
        }

    # =========================================================================
    # PART 1: THE OPERATOR DECOMPOSITION ACROSS DIMENSIONS
    # =========================================================================
    print("=" * 106)
    print("PART 1: EXACT OPERATOR RESOLVENT DECOMPOSITION: D_1 / D_0 vs T_diag & T_cross")
    print("=" * 106)
    print("Theorem 7.3: D_1 / D_0 = kappa^2 [ T_diag + T_cross - D_0^2 M_2 ]")
    print("where T_diag = M_1 * <d, R_even d>  and  T_cross = <d, R_even K psi>\n")

    print(
        f"{'N':>3} "
        f"{'D_1 / D_0':>18} "
        f"{'kappa^2 T_diag':>18} "
        f"{'kappa^2 T_cross':>18} "
        f"{'kappa^2 ||Kc||^2':>16} "
        f"{'T_diag / Total':>14} "
        f"{'Abs Error':>12}"
    )
    print("-" * 106)

    for N in N_LIST:
        res = results_across_N[N]
        d1_d0 = res["D1_over_D0"]
        t_diag_scaled = (KAPPA ** 2) * res["T_diag"]
        t_cross_scaled = (KAPPA ** 2) * res["T_cross"]
        sub_scaled = (KAPPA ** 2) * res["subtraction_norm"]
        ratio_diag = res["T_diag"] / res["T_total"]
        err = res["abs_error"]

        print(
            f"{N:3d} "
            f"{mp.nstr(d1_d0, 12):>18} "
            f"{mp.nstr(t_diag_scaled, 12):>18} "
            f"{mp.nstr(t_cross_scaled, 12):>18} "
            f"{mp.nstr(sub_scaled, 8):>16} "
            f"{mp.nstr(ratio_diag, 8):>14} "
            f"{mp.nstr(err, 6):>12}"
        )

    # =========================================================================
    # PART 2: EVEN-SECTOR MODE-BY-MODE UNIFORMITY (N = 24)
    # =========================================================================
    print("\n" + "=" * 106)
    print("PART 2: EVEN-SECTOR MODE-BY-MODE RESOLVENT SPECTRUM (N = 24)")
    print("=" * 106)
    print("Auditing d_k = <u_k, d> and the rescaled overlap d_k^2 / Delta_{even, k}:")
    print(
        f"\n{'Mode k':>6} "
        f"{'E_k':>18} "
        f"{'Delta_k = E_k - lam0':>22} "
        f"{'d_k':>18} "
        f"{'d_k^2 / Delta_k':>18} "
        f"{'Cum <d, R d>':>18}"
    )
    print("-" * 106)

    even_data_24 = results_across_N[24]["even_mode_data"]
    for item in even_data_24[:12]:
        k = item["k"]
        Ek = item["E_k"]
        Dk = item["Delta_k"]
        dk = item["d_k"]
        r_d = item["ratio_d2_Delta"]
        cum_d = item["cum_dRd"]
        print(
            f"{k:6d} "
            f"{mp.nstr(Ek, 10):>18} "
            f"{mp.nstr(Dk, 12):>22} "
            f"{mp.nstr(dk, 10):>18} "
            f"{mp.nstr(r_d, 10):>18} "
            f"{mp.nstr(cum_d, 12):>18}"
        )

    # =========================================================================
    # PART 3: ODD-SECTOR MODE-BY-MODE UNIFORMITY (N = 24)
    # =========================================================================
    print("\n" + "=" * 106)
    print("PART 3: ODD-SECTOR MODE-BY-MODE RESOLVENT SPECTRUM (N = 24)")
    print("=" * 106)
    print("Auditing a_j = <e_j, psi> and the rescaled overlap a_j^2 / Delta_{odd, j}:")
    print(
        f"\n{'Mode j':>6} "
        f"{'mu_j':>18} "
        f"{'Delta_j = mu_j - lam0':>22} "
        f"{'a_j':>18} "
        f"{'a_j^2 / Delta_j':>18} "
        f"{'Cum M_1':>18}"
    )
    print("-" * 106)

    odd_data_24 = results_across_N[24]["odd_mode_data"]
    for item in odd_data_24[:12]:
        j = item["j"]
        muj = item["mu_j"]
        Dj = item["Delta_j"]
        aj = item["a_j"]
        r_a = item["ratio_a2_Delta"]
        cum_M1 = item["cum_M1"]
        print(
            f"{j:6d} "
            f"{mp.nstr(muj, 10):>18} "
            f"{mp.nstr(Dj, 12):>22} "
            f"{mp.nstr(aj, 10):>18} "
            f"{mp.nstr(r_a, 10):>18} "
            f"{mp.nstr(cum_M1, 12):>18}"
        )

    # =========================================================================
    # PART 4: SPATIAL PROFILE OF THE RESOLVENT VECTOR w_d = R_even d (N = 24)
    # =========================================================================
    print("\n" + "=" * 106)
    print("PART 4: SPATIAL COORDINATE PROFILE OF w_d = R_even d (N = 24)")
    print("=" * 106)
    print("Inspecting the components w_d(m) in canonical basis for m = 0, ..., N:")
    print(
        f"\n{'m':>4} "
        f"{'Normalized s = m/N':>20} "
        f"{'w_d(m)':>24} "
        f"{'w_d(m) / ||w_d||':>22}"
    )
    print("-" * 74)

    w_d_24 = results_across_N[24]["w_d"]
    norm_wd = mp.sqrt(dot(w_d_24, w_d_24))

    for m in range(25):
        s_norm = mp.mpf(m) / 24
        wm = w_d_24[m, 0]
        wm_norm = wm / norm_wd
        print(
            f"{m:4d} "
            f"{mp.nstr(s_norm, 6):>20} "
            f"{mp.nstr(wm, 14):>24} "
            f"{mp.nstr(wm_norm, 10):>22}"
        )

    # =========================================================================
    # PART 5: LARGE-N ASYMPTOTIC SCALING EXCLUDED FROM EXPONENTIAL SMALLNESS
    # =========================================================================
    print("\n" + "=" * 106)
    print("PART 5: ASYMPTOTIC SCALING OF RESOLVENT ELEMENTS ACROSS DIMENSIONS")
    print("=" * 106)
    print("Tracking M_1(N), <d, R_even d>(N), and effective power-law growth:")
    print(
        f"\n{'N':>3} "
        f"{'M_1(N)':>16} "
        f"{'<d, R_even d>':>18} "
        f"{'Product M_1 <d,Rd>':>20} "
        f"{'D_1 / D_0':>18} "
        f"{'Ratio (D_1/D_0) / (kappa^2 M_1 <d,Rd>)':>36}"
    )
    print("-" * 116)

    for N in N_LIST:
        res = results_across_N[N]
        m1 = res["M1"]
        drd = res["d_R_d"]
        prod = m1 * drd
        d1_d0 = res["D1_over_D0"]
        ratio_eff = d1_d0 / ((KAPPA ** 2) * prod)

        print(
            f"{N:3d} "
            f"{mp.nstr(m1, 8):>16} "
            f"{mp.nstr(drd, 10):>18} "
            f"{mp.nstr(prod, 12):>20} "
            f"{mp.nstr(d1_d0, 10):>18} "
            f"{mp.nstr(ratio_eff, 8):>36}"
        )

    print("\n" + "=" * 106)
    print("CELL 62 EXECUTION COMPLETE")
    print("=" * 106)


if __name__ == "__main__":
    main()
