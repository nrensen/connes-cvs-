#!/usr/bin/env python3
"""
================================================================================
CELL 72 — TARGETED AUDIT OF SEMICLASSICAL ACTION COMPETITION AND EVEN-RESOLVENT
         POLE CANCELLATION IN THE FIRST EXCITED TAIL MODE (j = 2)
================================================================================

PURPOSE:
--------
Perform a dedicated numerical and asymptotic reconnaissance of the primary excited
wavepacket tail mode j = 2 (Paper NR2 §8.23, Milestone M19).

From Cell 70, mode j = 2 constitutes 99.997% of the entire wavepacket tail
sum_{j >= 2} T_j at N = 24. From Cell 71, the overlap ratio a_2^2 / a_1^2 is not
polynomial, but grows from 1.66e4 at N = 8 to 1.62e5 at N = 24 due to the steep
near-pole growth of G_d'(mu_1).

This script audits the precise mechanism that defeats this overlap growth:

1. Mode j = 2 Three-Factor Semiclassical Decomposition:
       T_2 = K_2 * G_2 * S_2 = K_2 * Q_2,
   where:
       K_2 = ||K u_2||^2 / ||K u_1||^2             (kinetic energy ratio, O(1)),
       G_2 = G_d'(mu_1) / G_d'(mu_2)               (Stieltjes derivative ratio),
       S_2 = ((mu_1 - lam_0) / (mu_2 - lam_0))^2   (squared gap suppression),
       Q_2 = G_2 * S_2                             (combined filtering product).

2. Exact Even-Resolvent Pole-Cancellation Identity:
   Verify the exact finite-N algebraic identity (Eq. 8.23.4):
       Q_2 = [ D_0^2 + (mu_1 - lam_0)^2 * E_even(mu_1) ] /
             [ D_0^2 + (mu_2 - lam_0)^2 * E_even(mu_2) ],
   where E_even(mu) = sum_{k=1}^N d_k^2 / (mu - E_k)^2 is the excited even-resolvent sum.
   Confirm that the ground-state pole D_0^2 cancels symmetrically, and determine the
   exact magnitude of the numerator and denominator components.

3. Semiclassical Action Slopes and Rate Inequality:
   Fit consecutive logarithmic rates across N in {8, 12, 16, 20, 24}:
       gamma_2  =  Delta log(G_2) / Delta N     (overlap growth rate),
       tau_2    = -Delta log(S_2) / Delta N     (gap suppression rate),
       sigma_Q2 = -Delta log(Q_2) / Delta N     (net filtering rate),
       sigma_T2 = -Delta log(T_2) / Delta N     (net tail extinction rate).
   Audit the fundamental action inequality:
       tau_2 > gamma_2   <===>   Delta S_2 = tau_2 - gamma_2 > 0.

4. Mode j = 3 Comparison:
   Audit the same quantities for j = 3 to verify whether the action gap
   Delta S_j = tau_j - gamma_j widens further up the bound-state ladder.

OUTPUT CONSTRAINTS:
-------------------
- Dispassionate, dry, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 72 EXECUTION COMPLETE.
================================================================================
"""

from __future__ import annotations

import time
import mpmath as mp

from connes_cvs import build_galerkin_matrix
from cell import get_galerkin_matrix

# -----------------------------------------------------------------------------
# Precision and Parameter Configuration
# -----------------------------------------------------------------------------

mp.mp.dps = 50

C_PARAM = 13
L_PARAM = mp.log(C_PARAM)
T_PARAM = 400
GROUND_DPS = 50

N_LIST = [8, 12, 16, 20, 24]
BARRIER_V_STAR = mp.mpf('1.0')


# -----------------------------------------------------------------------------
# Parity Basis & Operator Construction
# -----------------------------------------------------------------------------

def full_parity_basis(N: int) -> tuple[mp.matrix, mp.matrix]:
    """Construct orthonormal basis matrices E (even) and O (odd) for R^{2N+1}."""
    dim = 2 * N + 1
    E = mp.matrix(dim, N + 1)
    O = mp.matrix(dim, N)

    E[N, 0] = mp.mpf(1)

    inv_sqrt2 = 1 / mp.sqrt(2)
    for m in range(1, N + 1):
        E[N + m, m] = inv_sqrt2
        E[N - m, m] = inv_sqrt2
        O[N + m, m - 1] = inv_sqrt2
        O[N - m, m - 1] = -inv_sqrt2

    return E, O


def solve_parity_eigensystems(
    Q_full: mp.matrix, N: int
) -> tuple[
    mp.mpf, mp.matrix, mp.matrix,
    list[mp.mpf], mp.matrix,
    list[mp.mpf], mp.matrix
]:
    """Compute even and odd spectra and eigenvectors."""
    E, O = full_parity_basis(N)

    Q_even = E.T * Q_full * E
    Q_even = mp.mpf('0.5') * (Q_even + Q_even.T)

    Q_odd = O.T * Q_full * O
    Q_odd = mp.mpf('0.5') * (Q_odd + Q_odd.T)

    evals_e, V_e = mp.eigsy(Q_even)
    evals_o, V_o = mp.eigsy(Q_odd)

    # Sort even
    idx_e = sorted(range(N + 1), key=lambda i: evals_e[i])
    sorted_evals_e = [evals_e[i] for i in idx_e]
    sorted_V_e = mp.matrix(N + 1, N + 1)
    for col_idx, orig_col in enumerate(idx_e):
        for row_idx in range(N + 1):
            sorted_V_e[row_idx, col_idx] = V_e[row_idx, orig_col]

    # Sort odd
    idx_o = sorted(range(N), key=lambda i: evals_o[i])
    sorted_evals_o = [evals_o[i] for i in idx_o]
    sorted_V_o = mp.matrix(N, N)
    for col_idx, orig_col in enumerate(idx_o):
        for row_idx in range(N):
            sorted_V_o[row_idx, col_idx] = V_o[row_idx, orig_col]

    lam_0 = sorted_evals_e[0]
    return lam_0, E, O, sorted_evals_e, sorted_V_e, sorted_evals_o, sorted_V_o


# -----------------------------------------------------------------------------
# Main Execution Function
# -----------------------------------------------------------------------------

def run_cell72() -> None:
    print("=" * 105)
    print("CELL 72 — TARGETED AUDIT OF SEMICLASSICAL ACTION COMPETITION AND EVEN-RESOLVENT")
    print("         POLE CANCELLATION IN THE FIRST EXCITED TAIL MODE (j = 2)")
    print("=" * 105)
    print(f"Parameters: c = {C_PARAM}, L = {float(L_PARAM):.14f}, T = {T_PARAM}, dps = {mp.mp.dps}, V_* = {float(BARRIER_V_STAR)}")
    print()

    synthesis_records: list[dict] = []

    for N in N_LIST:
        t0 = time.time()

        # Retrieve cached Galerkin matrix
        Q_full, _ = get_galerkin_matrix(
            c=C_PARAM,
            N=N,
            T=T_PARAM,
            dps=GROUND_DPS,
            verbose=False,
        )

        # Solve parity eigensystems
        lam_0, E, O, evals_e, V_e, evals_o, V_o = solve_parity_eigensystems(Q_full, N)

        # Ground state c in R^{N+1}
        c_even = V_e[:, 0]

        # Boundary vector d in R^{N+1}: d_even = (1, sqrt(2), ..., sqrt(2))^T
        d_even = mp.matrix(N + 1, 1)
        d_even[0, 0] = mp.mpf(1)
        for m in range(1, N + 1):
            d_even[m, 0] = mp.sqrt(2)

        D_0 = sum(d_even[m, 0] * c_even[m, 0] for m in range(N + 1))
        D_0_sq = D_0 ** 2

        # Boundary overlaps in even sector: d_k = <u_k^{even}, d>
        d_k_list: list[mp.mpf] = []
        for k in range(N + 1):
            d_val = sum(d_even[m, 0] * V_e[m, k] for m in range(N + 1))
            d_k_list.append(d_val)

        # Source vector psi = K Q e_0 on R^{2N+1}, mapped to H_odd via O
        # psi_full[N + m] = m * Q_full[N + m, N], psi_full[N - m] = -m * Q_full[N - m, N]
        # In odd basis: psi_odd[m-1] = sqrt(2) * m * Q_full[N + m, N]
        psi_odd = mp.matrix(N, 1)
        for m in range(1, N + 1):
            psi_odd[m - 1, 0] = mp.sqrt(2) * m * Q_full[N + m, N]

        # Mode overlaps in odd sector: a_j = <psi, u_j>
        a_j_list: list[mp.mpf] = []
        for j in range(N):
            a_val = sum(psi_odd[m - 1, 0] * V_o[m - 1, j] for m in range(1, N + 1))
            a_j_list.append(a_val)

        # Coordinate transition dipoles b_{0j} = <c, K u_j>
        b_0j_list: list[mp.mpf] = []
        for j in range(N):
            b_j = sum(c_even[m, 0] * m * V_o[m - 1, j] for m in range(1, N + 1))
            b_0j_list.append(b_j)

        # Coordinate kinetic energies: ||K u_j||^2 = sum_{m=1}^N m^2 (u_{j, m-1})^2
        ke_list: list[mp.mpf] = []
        for j in range(N):
            ke_val = sum((m ** 2) * (V_o[m - 1, j] ** 2) for m in range(1, N + 1))
            ke_list.append(ke_val)

        # Stieltjes derivatives and excited even-resolvent sums
        g_prime_list: list[mp.mpf] = []
        e_even_list: list[mp.mpf] = []
        for j in range(N):
            mu_j = evals_o[j]
            gp = sum((d_k_list[k] ** 2) / ((mu_j - evals_e[k]) ** 2) for k in range(N + 1))
            g_prime_list.append(gp)
            e_ev = sum((d_k_list[k] ** 2) / ((mu_j - evals_e[k]) ** 2) for k in range(1, N + 1))
            e_even_list.append(e_ev)

        # Mode 1 reference values
        gap_1 = evals_o[1] - lam_0
        ke_1 = ke_list[1]
        gp_1 = g_prime_list[1]
        e_even_1 = e_even_list[1]
        b_01 = b_0j_list[1]
        b_01_sq = b_01 ** 2

        # Numerator of pole cancellation identity
        num_pole_1 = D_0_sq + (gap_1 ** 2) * e_even_1

        # Evaluate modes j = 2 and j = 3
        mode_records: dict[int, dict] = {}
        for j in [2, 3]:
            if j >= N:
                continue

            gap_j = evals_o[j] - lam_0
            ke_j = ke_list[j]
            gp_j = g_prime_list[j]
            e_even_j = e_even_list[j]
            b_0j = b_0j_list[j]
            b_0j_sq = b_0j ** 2

            # Direct tail ratio
            T_j_direct = b_0j_sq / b_01_sq

            # Three-factor components
            K_j = ke_j / ke_1
            G_j = gp_1 / gp_j
            S_j = (gap_1 / gap_j) ** 2
            Q_j = G_j * S_j
            T_j_fact = K_j * Q_j

            # Even-resolvent pole cancellation formula
            den_pole_j = D_0_sq + (gap_j ** 2) * e_even_j
            Q_j_pole = num_pole_1 / den_pole_j
            T_j_pole = K_j * Q_j_pole

            # Residuals
            pole_id_res = abs(Q_j - Q_j_pole)
            tail_id_res = abs(T_j_direct - T_j_fact)

            # Component shares in denominator
            den_term_D0 = D_0_sq
            den_term_exc = (gap_j ** 2) * e_even_j

            mode_records[j] = {
                "j": j,
                "gap_j": gap_j,
                "ke_j": ke_j,
                "gp_j": gp_j,
                "e_even_j": e_even_j,
                "b_0j_sq": b_0j_sq,
                "T_j_direct": T_j_direct,
                "K_j": K_j,
                "G_j": G_j,
                "S_j": S_j,
                "Q_j": Q_j,
                "Q_j_pole": Q_j_pole,
                "T_j_fact": T_j_fact,
                f"K_{j}": K_j,
                f"G_{j}": G_j,
                f"S_{j}": S_j,
                f"Q_{j}": Q_j,
                f"Q_{j}_pole": Q_j_pole,
                "pole_id_res": pole_id_res,
                "tail_id_res": tail_id_res,
                "den_term_D0": den_term_D0,
                "den_term_exc": den_term_exc,
            }

        elapsed = time.time() - t0

        record = {
            "N": N,
            "elapsed": elapsed,
            "D_0_sq": D_0_sq,
            "gap_1": gap_1,
            "ke_1": ke_1,
            "gp_1": gp_1,
            "e_even_1": e_even_1,
            "num_pole_1": num_pole_1,
            "num_term_exc": (gap_1 ** 2) * e_even_1,
            "b_01_sq": b_01_sq,
            "modes": mode_records,
        }
        synthesis_records.append(record)

        print(f"--- DIMENSION N = {N:2d} (Elapsed: {elapsed:.2f}s) ---")
        print(f"  Tunneling gap:       mu_1 - lam_0                   = {mp.nstr(gap_1, 9)}")
        print(f"  Boundary squared:    D_0^2                          = {mp.nstr(D_0_sq, 9)}")
        print(f"  Num pole num:        D_0^2 + (mu_1-lam)^2*E_even(1) = {mp.nstr(num_pole_1, 9)}")
        print(f"  Num exc term:        (mu_1-lam)^2*E_even(1)         = {mp.nstr(record['num_term_exc'], 4)}")

        if 2 in mode_records:
            m2 = mode_records[2]
            print(f"  Mode 2 kinetic ratio:   K_2 = ||Ku_2||^2 / ||Ku_1||^2 = {mp.nstr(m2['K_j'], 6)}")
            print(f"  Mode 2 Stieltjes ratio: G_2 = G'(1) / G'(2)           = {mp.nstr(m2['G_j'], 6)}")
            print(f"  Mode 2 gap ratio:       S_2 = ((mu_1-lam)/(mu_2-lam))^2= {mp.nstr(m2['S_j'], 6)}")
            print(f"  Mode 2 product:         Q_2 = G_2 * S_2               = {mp.nstr(m2['Q_j'], 6)}")
            print(f"  Mode 2 pole formula:    Q_2^{{pole}}                   = {mp.nstr(m2['Q_j_pole'], 6)}")
            print(f"  Pole id residual:       |Q_2 - Q_2^{{pole}}|           = {mp.nstr(m2['pole_id_res'], 4)}")
            print(f"  Mode 2 tail ratio:      T_2 = b_{{02}}^2 / b_{{01}}^2   = {mp.nstr(m2['T_j_direct'], 6)}")
            print(f"  Denom D_0^2 term:       D_0^2                         = {mp.nstr(m2['den_term_D0'], 4)}")
            print(f"  Denom exc term:         (mu_2-lam)^2 * E_even(2)      = {mp.nstr(m2['den_term_exc'], 4)}")
        print()

    # -------------------------------------------------------------------------
    # Multi-Dimension Synthesis Tables
    # -------------------------------------------------------------------------

    # Table 1: Mode j = 2 Semiclassical Action Decomposition
    print("=" * 115)
    print("SYNTHESIS TABLE 1: MODE j = 2 SEMICLASSICAL ACTION DECOMPOSITION (T_2 = K_2 * G_2 * S_2 = K_2 * Q_2)")
    print("=" * 115)
    header1 = (
        f"{'N':>3} | {'K_2':>8} | {'G_2 = G\'(1)/G\'(2)':>20} | {'S_2 (gap ratio)':>16} | "
        f"{'Q_2 = G_2 * S_2':>16} | {'T_2 (direct)':>14} | {'T_2 (fact)':>14}"
    )
    print(header1)
    print("-" * len(header1))

    for rec in synthesis_records:
        if 2 in rec["modes"]:
            m = rec["modes"][2]
            row = (
                f"{rec['N']:3d} | "
                f"{mp.nstr(m['K_j'], 5):>8} | "
                f"{mp.nstr(m['G_j'], 6):>20} | "
                f"{mp.nstr(m['S_j'], 6):>16} | "
                f"{mp.nstr(m['Q_j'], 6):>16} | "
                f"{mp.nstr(m['T_j_direct'], 6):>14} | "
                f"{mp.nstr(m['T_j_fact'], 6):>14}"
            )
            print(row)

    print()

    # Table 2: Exact Even-Resolvent Pole Cancellation Architecture
    print("=" * 125)
    print("SYNTHESIS TABLE 2: EXACT EVEN-RESOLVENT POLE CANCELLATION: Q_2 = [D_0^2 + (mu_1-lam)^2*E_1] / [D_0^2 + (mu_2-lam)^2*E_2]")
    print("=" * 125)
    header2 = (
        f"{'N':>3} | {'D_0^2':>14} | {'(mu_1-lam)^2*E_1':>18} | {'Numerator':>14} | "
        f"{'(mu_2-lam)^2*E_2':>18} | {'Denominator':>14} | {'Q_2^{pole}':>14} | {'Pole Id Res':>11}"
    )
    print(header2)
    print("-" * len(header2))

    for rec in synthesis_records:
        if 2 in rec["modes"]:
            m = rec["modes"][2]
            row = (
                f"{rec['N']:3d} | "
                f"{mp.nstr(rec['D_0_sq'], 5):>14} | "
                f"{mp.nstr(rec['num_term_exc'], 5):>18} | "
                f"{mp.nstr(rec['num_pole_1'], 5):>14} | "
                f"{mp.nstr(m['den_term_exc'], 5):>18} | "
                f"{mp.nstr(m['den_term_D0'] + m['den_term_exc'], 5):>14} | "
                f"{mp.nstr(m['Q_j_pole'], 5):>14} | "
                f"{mp.nstr(m['pole_id_res'], 4):>11}"
            )
            print(row)

    print()

    # Table 3: Consecutive Logarithmic Scaling Slopes and Action Competition
    print("=" * 115)
    print("SYNTHESIS TABLE 3: CONSECUTIVE LOGARITHMIC SLOPES AND ACTION COMPETITION (tau_2 > gamma_2)")
    print("=" * 115)
    header3 = (
        f"{'Interval':>10} | {'gamma_2 (overlap)':>18} | {'tau_2 (gap)':>16} | "
        f"{'Delta S_2 = tau-gamma':>22} | {'sigma_Q2':>14} | {'sigma_T2':>14}"
    )
    print(header3)
    print("-" * len(header3))

    for idx in range(len(synthesis_records) - 1):
        rA = synthesis_records[idx]
        rB = synthesis_records[idx + 1]
        NA, NB = rA["N"], rB["N"]
        dN = NB - NA

        if 2 in rA["modes"] and 2 in rB["modes"]:
            mA = rA["modes"][2]
            mB = rB["modes"][2]

            gamma_2 = (mp.log(mB["G_j"]) - mp.log(mA["G_j"])) / dN
            tau_2 = -(mp.log(mB["S_j"]) - mp.log(mA["S_j"])) / dN
            delta_S2 = tau_2 - gamma_2
            sigma_Q2 = -(mp.log(mB["Q_j"]) - mp.log(mA["Q_j"])) / dN
            sigma_T2 = -(mp.log(mB["T_j_direct"]) - mp.log(mA["T_j_direct"])) / dN

            interval_str = f"N={NA}->{NB}"
            row = (
                f"{interval_str:>10} | "
                f"{mp.nstr(gamma_2, 5):>18} | "
                f"{mp.nstr(tau_2, 5):>16} | "
                f"{mp.nstr(delta_S2, 5):>22} | "
                f"{mp.nstr(sigma_Q2, 5):>14} | "
                f"{mp.nstr(sigma_T2, 5):>14}"
            )
            print(row)

    print()

    # Table 4: Mode j = 3 Comparison
    print("=" * 115)
    print("SYNTHESIS TABLE 4: HIGHER-MODE COMPARISON: MODE j = 3 (T_3 = K_3 * G_3 * S_3 = K_3 * Q_3)")
    print("=" * 115)
    header4 = (
        f"{'N':>3} | {'K_3':>8} | {'G_3':>14} | {'S_3':>16} | "
        f"{'Q_3':>16} | {'T_3 (direct)':>14} | {'Pole Id Res':>11}"
    )
    print(header4)
    print("-" * len(header4))

    for rec in synthesis_records:
        if 3 in rec["modes"]:
            m = rec["modes"][3]
            row = (
                f"{rec['N']:3d} | "
                f"{mp.nstr(m['K_j'], 5):>8} | "
                f"{mp.nstr(m['G_j'], 5):>14} | "
                f"{mp.nstr(m['S_j'], 5):>16} | "
                f"{mp.nstr(m['Q_j'], 5):>16} | "
                f"{mp.nstr(m['T_j_direct'], 5):>14} | "
                f"{mp.nstr(m['pole_id_res'], 4):>11}"
            )
            print(row)

    print()
    print("=" * 80)
    print("CELL 72 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell72()
