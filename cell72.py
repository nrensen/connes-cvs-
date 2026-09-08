#!/usr/bin/env python3
"""
================================================================================
CELL 72 — TARGETED AUDIT OF SEMICLASSICAL ACTION COMPETITION AND EVEN-RESOLVENT
         POLE CANCELLATION IN THE FIRST EXCITED TAIL MODE (j = 2)
================================================================================

PURPOSE:
--------
Perform a dedicated numerical and asymptotic reconnaissance of the primary excited
wavepacket tail mode j = 2 (Paper 4B §8.23, Milestone M19).

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
    E_basis, O_basis = full_parity_basis(N)

    # Project full Galerkin matrix onto parity subspaces
    Q_even = E_basis.T * Q_full * E_basis
    Q_odd = O_basis.T * Q_full * O_basis

    # Symmetrize to eliminate any numerical asymmetry
    Q_even = mp.mpf('0.5') * (Q_even + Q_even.T)
    Q_odd = mp.mpf('0.5') * (Q_odd + Q_odd.T)

    # Even eigensystem
    E_vals_raw, E_vecs_raw = mp.eigsy(Q_even)
    even_pairs = sorted(zip([mp.mpf(x) for x in E_vals_raw], range(N + 1)), key=lambda p: p[0])
    even_evals = [p[0] for p in even_pairs]
    even_evecs = mp.matrix(N + 1, N + 1)
    for new_col, (_, old_col) in enumerate(even_pairs):
        for row in range(N + 1):
            even_evecs[row, new_col] = E_vecs_raw[row, old_col]

    # Odd eigensystem
    O_vals_raw, O_vecs_raw = mp.eigsy(Q_odd)
    odd_pairs = sorted(zip([mp.mpf(x) for x in O_vals_raw], range(N)), key=lambda p: p[0])
    odd_evals = [p[0] for p in odd_pairs]
    odd_evecs = mp.matrix(N, N)
    for new_col, (_, old_col) in enumerate(odd_pairs):
        for row in range(N):
            odd_evecs[row, new_col] = O_vecs_raw[row, old_col]

    # Ground state and eigenvalue
    lam_0 = even_evals[0]
    c_even = mp.matrix(N + 1, 1)
    for row in range(N + 1):
        c_even[row, 0] = even_evecs[row, 0]

    # Full ground state vector in R^{2N+1}
    c_full = E_basis * c_even
    if c_full[N, 0] < 0:
        c_full = -c_full
        c_even = -c_even
        for row in range(N + 1):
            even_evecs[row, 0] = -even_evecs[row, 0]

    return lam_0, c_full, Q_even, even_evals, even_evecs, odd_evals, odd_evecs


def compute_coordinate_derivative_norms(
    odd_evecs: mp.matrix, N: int
) -> list[mp.mpf]:
    """
    Compute discrete coordinate kinetic energy ||K u_j||^2 for each odd mode.
    In the odd basis e_{m-1}^O = (e_m - e_{-m})/sqrt(2), K is diagonal:
    K e_{m-1}^O = m e_{m-1}^O.
    Thus ||K u_j||^2 = sum_{m=1}^N m^2 (u_{j, m-1})^2.
    """
    ke_norms = []
    for j in range(N):
        norm_sq = mp.mpf(0)
        for m in range(1, N + 1):
            coord_val = odd_evecs[m - 1, j]
            norm_sq += (mp.mpf(m) ** 2) * (coord_val ** 2)
        ke_norms.append(norm_sq)
    return ke_norms


def compute_stieltjes_derivatives(
    even_evals: list[mp.mpf],
    d_coeffs: list[mp.mpf],
    mu_val: mp.mpf,
    N: int
) -> tuple[mp.mpf, mp.mpf]:
    """
    Compute full Stieltjes derivative:
        G_d'(mu) = sum_{k=0}^N d_k^2 / (mu - E_k)^2
    and excited even-resolvent sum (k >= 1):
        E_even(mu) = sum_{k=1}^N d_k^2 / (mu - E_k)^2.
    """
    gp_full = mp.mpf(0)
    e_even = mp.mpf(0)
    for k in range(N + 1):
        term = (d_coeffs[k] ** 2) / ((mu_val - even_evals[k]) ** 2)
        gp_full += term
        if k >= 1:
            e_even += term
    return gp_full, e_even


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
        Q_full = get_galerkin_matrix(C_PARAM, L_PARAM, N, T_PARAM, GROUND_DPS)

        # Solve parity eigensystems
        lam_0, c_full, Q_even, even_evals, even_evecs, odd_evals, odd_evecs = solve_parity_eigensystems(Q_full, N)

        # Dirichlet constant vector d = (1, ..., 1)^T in R^{2N+1}
        dim = 2 * N + 1
        d_full = mp.matrix([mp.mpf(1)] * dim)
        E_basis, O_basis = full_parity_basis(N)
        d_even = E_basis.T * d_full

        # Even boundary overlaps d_k = <d, u_k^{even}>
        d_coeffs = []
        for k in range(N + 1):
            u_k = even_evecs[:, k]
            val = mp.mpf(0)
            for row in range(N + 1):
                val += d_even[row, 0] * u_k[row, 0]
            d_coeffs.append(val)

        D_0 = d_coeffs[0]
        D_0_sq = D_0 ** 2

        # Odd source overlaps a_j = <psi, u_j>
        e_0 = mp.matrix(dim, 1)
        e_0[N, 0] = mp.mpf(1)
        Q_e0 = Q_full * e_0
        psi_full = mp.matrix(dim, 1)
        for idx in range(dim):
            m = idx - N
            psi_full[idx, 0] = mp.mpf(m) * Q_e0[idx, 0]
        psi_odd = O_basis.T * psi_full

        a_coeffs = []
        for j in range(N):
            u_j = odd_evecs[:, j]
            val = mp.mpf(0)
            for row in range(N):
                val += psi_odd[row, 0] * u_j[row, 0]
            a_coeffs.append(val)

        # Coordinate kinetic energy norms
        ke_norms = compute_coordinate_derivative_norms(odd_evecs, N)

        # Mode 1 reference values
        gap_1 = odd_evals[1] - lam_0
        ke_1 = ke_norms[1]
        gp_1, e_even_1 = compute_stieltjes_derivatives(even_evals, d_coeffs, odd_evals[1], N)
        b_01 = a_coeffs[1] * D_0 / gap_1
        b_01_sq = b_01 ** 2

        # Numerator of pole cancellation identity
        num_pole_1 = D_0_sq + (gap_1 ** 2) * e_even_1

        # Evaluate modes j = 2 and j = 3
        mode_records: dict[int, dict] = {}
        for j in [2, 3]:
            if j >= N:
                continue

            gap_j = odd_evals[j] - lam_0
            ke_j = ke_norms[j]
            gp_j, e_even_j = compute_stieltjes_derivatives(even_evals, d_coeffs, odd_evals[j], N)
            b_0j = a_coeffs[j] * D_0 / gap_j
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
            print(f"  Mode 2 kinetic ratio:   K_2 = ||Ku_2||^2 / ||Ku_1||^2 = {mp.nstr(m2['K_2'], 6)}")
            print(f"  Mode 2 Stieltjes ratio: G_2 = G'(1) / G'(2)           = {mp.nstr(m2['G_2'], 6)}")
            print(f"  Mode 2 gap ratio:       S_2 = ((mu_1-lam)/(mu_2-lam))^2= {mp.nstr(m2['S_2'], 6)}")
            print(f"  Mode 2 product:         Q_2 = G_2 * S_2               = {mp.nstr(m2['Q_2'], 6)}")
            print(f"  Mode 2 pole formula:    Q_2^{{pole}}                   = {mp.nstr(m2['Q_2_pole'], 6)}")
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
                f"{mp.nstr(m['K_2'], 5):>8} | "
                f"{mp.nstr(m['G_2'], 6):>20} | "
                f"{mp.nstr(m['S_2'], 6):>16} | "
                f"{mp.nstr(m['Q_2'], 6):>16} | "
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

            gamma_2 = (mp.log(mB["G_2"]) - mp.log(mA["G_2"])) / dN
            tau_2 = -(mp.log(mB["S_2"]) - mp.log(mA["S_2"])) / dN
            delta_S2 = tau_2 - gamma_2
            sigma_Q2 = -(mp.log(mB["Q_2"]) - mp.log(mA["Q_2"])) / dN
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
                f"{mp.nstr(m['K_3'], 5):>8} | "
                f"{mp.nstr(m['G_3'], 5):>14} | "
                f"{mp.nstr(m['S_3'], 5):>16} | "
                f"{mp.nstr(m['Q_3'], 5):>16} | "
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
