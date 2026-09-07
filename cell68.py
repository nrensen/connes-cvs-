#!/usr/bin/env python3
"""
================================================================================
CELL 68 — AUDIT OF THE GLOBAL WEIGHTED RESOLVENT TRACE T_N, DISCRETE
COORDINATE KINETIC ENERGIES, AND EXPONENTIAL QUENCHING D_0^2 T_N
================================================================================

PURPOSE:
--------
Investigate the global weighted resolvent trace:
    T_N = Tr_{H_odd}[ K^2 (Q_odd - lambda I)^{-1} P_{perp u_0} ]
        = sum_{m=1}^N m^2 [R_{odd, perp}(lambda)]_{mm}
        = sum_{j=1}^{N-1} ||K u_j||^2 / (mu_j - lambda)
across discrete Galerkin dimensions N in {8, 12, 16, 20, 24} (Paper 4B §8.17 & §8.18):

1. Exact Identity Verification:
   Verify the exact equality between the Fourier-basis diagonal trace:
       T_N^{trace} = sum_{m=1}^N m^2 [R_{odd, perp}]_{mm}
   and the odd-eigenbasis coordinate-energy sum:
       T_N^{spec} = sum_{j=1}^{N-1} ||K u_j||^2 / (mu_j - lambda).

2. Bound-State Coordinate Kinetic Energy Localization (Bridge B2):
   Compute the discrete coordinate kinetic energy:
       ||K u_j||^2 = sum_{m=1}^N m^2 u_{j, m}^2
   for bound and low-lying excited odd modes j in {0, 1, 2, 3, 4, 5}.
   Verify that ||K u_j||^2 remains O(1) across dimensions N, while the total trace
   satisfies the exact identity:
       sum_{j=0}^{N-1} ||K u_j||^2 = N(N+1)(2N+1) / 6.

3. Semiclassical Trace Partitioning (Bound vs. High Energy):
   Decompose T_N into bound-state and high-energy (above-barrier) sectors:
       T_N = T_{bound}(N) + T_{high}(N),
   where modes with mu_j < V_* (V_* = 1.0) constitute the bound ladder.
   Demonstrate that T_{bound} grows exponentially like 1 / (mu_1 - lambda) ~ e^{+sigma_1 N},
   proving that the bare trace T_N is exponentially divergent (not polynomially bounded).

4. Exponential Quenching of the Physical Scaled Trace D_0^2 T_N:
   Multiply T_N by the ground-state tunneling amplitude D_0^2 ~ e^{-2 sigma_0 N}.
   Verify that because of barrier thinning (sigma_0 > sigma_1), the scaled trace:
       D_0^2 T_N = D_0^2 T_{bound} + D_0^2 T_{high}
   is exponentially quenched to zero:
       D_0^2 T_N ~ exp(- (sigma_0 + Delta sigma) N) -> 0.
   Extract the effective decay action sigma_T(N) = - (1 / N) log(D_0^2 T_N).

5. Coordinate Residual Sum S_{coord}(N):
   Compute the exact excited coordinate-energy sum:
       S_{coord}(N) = sum_{j=1}^{N-1} [D_0^2 / (mu_j - lambda)] ||P_{perp u_0^{even}} K u_j||^2
   and confirm S_{coord}(N) <= D_0^2 T_N, verifying that Route B continuum decoupling
   holds unconditionally with respect to the polynomial trace hypothesis.

OUTPUT CONSTRAINTS:
-------------------
- Strictly dispassionate, objective numerical tables.
- Zero narrative editorializing.
- Terminating sentinel: CELL 68 EXECUTION COMPLETE.
================================================================================
"""

from __future__ import annotations

import time
import mpmath as mp

from connes_cvs import build_galerkin_matrix

# -----------------------------------------------------------------------------
# Precision and Parameter Configuration
# -----------------------------------------------------------------------------

mp.mp.dps = 50

C_PARAM = 13
L_PARAM = mp.log(C_PARAM)
T_PARAM = 400
GROUND_DPS = 50
KAPPA = 2 * mp.pi / L_PARAM

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
# Main Diagnostic Execution
# -----------------------------------------------------------------------------

def run_cell68() -> None:
    print("=" * 96)
    print("CELL 68 — AUDIT OF THE GLOBAL WEIGHTED RESOLVENT TRACE T_N, DISCRETE")
    print("         COORDINATE KINETIC ENERGIES, AND EXPONENTIAL QUENCHING D_0^2 T_N")
    print("=" * 96)
    print(
        f"Parameters: c = {C_PARAM}, L = {mp.nstr(L_PARAM, 15)}, "
        f"T = {T_PARAM}, dps = {mp.mp.dps}, V_* = {mp.nstr(BARRIER_V_STAR, 4)}"
    )
    print()

    synthesis_records: list[dict] = []

    for N in N_LIST:
        t0 = time.time()

        # Build Galerkin operator
        Q_full = build_galerkin_matrix(N, C_PARAM, L_PARAM, T_PARAM, dps=GROUND_DPS)
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

        # Odd eigenvalues mu_j (j = 0, ..., N-1)
        mu_0 = evals_o[0]
        mu_1 = evals_o[1]
        gap_0 = mu_0 - lam_0
        gap_1 = mu_1 - lam_0

        # Exact total trace of K^2 on H_odd
        exact_trace_K2 = mp.mpf(N * (N + 1) * (2 * N + 1)) / 6

        # Coordinate kinetic energies: ||K u_j||^2 = sum_{m=1}^N m^2 (V_o[m-1, j])^2
        kinetic_energies: list[mp.mpf] = []
        for j in range(N):
            ke_j = sum((m ** 2) * (V_o[m - 1, j] ** 2) for m in range(1, N + 1))
            kinetic_energies.append(ke_j)

        sum_kinetic = sum(kinetic_energies)
        trace_residual = abs(sum_kinetic - exact_trace_K2)

        # Even-ground state overlap b_{0j} = <u_0^{even}, K u_j>
        # In coordinates: (K_{eo} u_j)[m] = m * u_{j, m} for m >= 1; 0 at m = 0.
        # So b_{0j} = sum_{m=1}^N c_even[m] * m * V_o[m-1, j]
        b_0j_sq_list: list[mp.mpf] = []
        for j in range(N):
            b_val = sum(c_even[m, 0] * m * V_o[m - 1, j] for m in range(1, N + 1))
            b_0j_sq_list.append(b_val ** 2)

        # Spectral expansion of T_N
        T_spec = mp.mpf(0)
        T_bound = mp.mpf(0)
        T_high = mp.mpf(0)
        S_coord = mp.mpf(0)
        n_bound = 0

        for j in range(1, N):
            mu_j = evals_o[j]
            gap_j = mu_j - lam_0
            term = kinetic_energies[j] / gap_j
            T_spec += term

            # Projected coordinate energy: ||P_{perp u_0^{even}} K u_j||^2 = ||K u_j||^2 - b_{0j}^2
            proj_ke = kinetic_energies[j] - b_0j_sq_list[j]
            S_coord += (D_0_sq / gap_j) * proj_ke

            if mu_j < BARRIER_V_STAR:
                T_bound += term
                n_bound += 1
            else:
                T_high += term

        # Fourier trace computation of T_N:
        # R_{odd, perp} = sum_{j=1}^{N-1} (1 / (mu_j - lam_0)) u_j u_j^T
        # [R_{odd, perp}]_{mm} = sum_{j=1}^{N-1} (V_o[m-1, j])^2 / (mu_j - lam_0)
        # T_trace = sum_{m=1}^N m^2 [R_{odd, perp}]_{mm}
        T_trace = mp.mpf(0)
        for m in range(1, N + 1):
            diag_R_m = sum((V_o[m - 1, j] ** 2) / (evals_o[j] - lam_0) for j in range(1, N))
            T_trace += (m ** 2) * diag_R_m

        identity_error = abs(T_spec - T_trace)

        # Scaled trace D_0^2 T_N
        D0_sq_T_N = D_0_sq * T_spec
        D0_sq_T_bound = D_0_sq * T_bound
        D0_sq_T_high = D_0_sq * T_high

        # Semiclassical actions
        sigma_0 = - (mp.mpf(1) / N) * mp.log(D_0_sq) / 2
        sigma_1 = - (mp.mpf(1) / N) * mp.log(gap_1)
        sigma_T = - (mp.mpf(1) / N) * mp.log(D0_sq_T_N)

        elapsed = time.time() - t0

        record = {
            "N": N,
            "elapsed": elapsed,
            "D_0_sq": D_0_sq,
            "gap_1": gap_1,
            "T_spec": T_spec,
            "T_trace": T_trace,
            "identity_error": identity_error,
            "T_bound": T_bound,
            "T_high": T_high,
            "n_bound": n_bound,
            "D0_sq_T_N": D0_sq_T_N,
            "D0_sq_T_bound": D0_sq_T_bound,
            "D0_sq_T_high": D0_sq_T_high,
            "S_coord": S_coord,
            "ke_0": kinetic_energies[0],
            "ke_1": kinetic_energies[1],
            "ke_2": kinetic_energies[2],
            "ke_3": kinetic_energies[3] if N > 3 else mp.mpf(0),
            "trace_residual": trace_residual,
            "sigma_0": sigma_0,
            "sigma_1": sigma_1,
            "sigma_T": sigma_T,
        }
        synthesis_records.append(record)

        print(f"--- DIMENSION N = {N:2d} (Elapsed: {elapsed:.2f}s) ---")
        print(f"  Ground scale:      D_0^2                = {mp.nstr(D_0_sq, 9)}")
        print(f"  Excited gap:       mu_1 - lam_0         = {mp.nstr(gap_1, 9)}")
        print(f"  Exact K^2 Trace:   sum_j ||Ku_j||^2     = {mp.nstr(sum_kinetic, 10)} (Residual: {mp.nstr(trace_residual, 4)})")
        print(f"  Ground KE:         ||Ku_0||^2           = {mp.nstr(kinetic_energies[0], 9)}")
        print(f"  Excited Mode 1 KE: ||Ku_1||^2           = {mp.nstr(kinetic_energies[1], 9)}")
        print(f"  Excited Mode 2 KE: ||Ku_2||^2           = {mp.nstr(kinetic_energies[2], 9)}")
        print(f"  Bare Trace T_N:    T_spec               = {mp.nstr(T_spec, 9)}")
        print(f"  Fourier Trace T_N: T_trace              = {mp.nstr(T_trace, 9)}")
        print(f"  Trace Identity:    |T_spec - T_trace|   = {mp.nstr(identity_error, 4)}")
        print(f"  Bound modes (<1):  N_bound              = {n_bound}")
        print(f"  Bound trace:       T_bound              = {mp.nstr(T_bound, 9)}")
        print(f"  High-energy trace: T_high               = {mp.nstr(T_high, 9)}")
        print(f"  Scaled Trace:      D_0^2 * T_N          = {mp.nstr(D0_sq_T_N, 9)}")
        print(f"  Scaled Bound:      D_0^2 * T_bound      = {mp.nstr(D0_sq_T_bound, 9)}")
        print(f"  Scaled High:       D_0^2 * T_high       = {mp.nstr(D0_sq_T_high, 9)}")
        print(f"  Coord Residual:    S_coord(N)           = {mp.nstr(S_coord, 9)}")
        print(f"  Actions:           sigma_0 = {mp.nstr(sigma_0, 6)}, sigma_1 = {mp.nstr(sigma_1, 6)}, sigma_T = {mp.nstr(sigma_T, 6)}")
        print()

    # -------------------------------------------------------------------------
    # Synthesis Tables
    # -------------------------------------------------------------------------
    print("=" * 116)
    print("SYNTHESIS TABLE 1: EXACT TRACE IDENTITY, COORDINATE KINETIC ENERGIES & BARE TRACE T_N")
    print("=" * 116)
    print(f"{'N':>4} | {'||Ku_0||^2':>12} | {'||Ku_1||^2':>12} | {'||Ku_2||^2':>12} | {'Tr(K^2) Error':>14} | {'Bare T_N':>14} | {'Identity Error':>14}")
    print("-" * 116)
    for r in synthesis_records:
        print(
            f"{r['N']:4d} | "
            f"{mp.nstr(r['ke_0'], 6):>12} | "
            f"{mp.nstr(r['ke_1'], 6):>12} | "
            f"{mp.nstr(r['ke_2'], 6):>12} | "
            f"{mp.nstr(r['trace_residual'], 4):>14} | "
            f"{mp.nstr(r['T_spec'], 6):>14} | "
            f"{mp.nstr(r['identity_error'], 4):>14}"
        )
    print("=" * 116)
    print()

    print("=" * 116)
    print("SYNTHESIS TABLE 2: SEMICLASSICAL SECTOR PARTITIONING, SCALED TRACE & EXPONENTIAL QUENCHING")
    print("=" * 116)
    print(f"{'N':>4} | {'T_bound':>13} | {'T_high':>13} | {'D_0^2 * T_N':>14} | {'S_coord(N)':>14} | {'sigma_0':>8} | {'sigma_1':>8} | {'sigma_T':>8}")
    print("-" * 116)
    for r in synthesis_records:
        print(
            f"{r['N']:4d} | "
            f"{mp.nstr(r['T_bound'], 6):>13} | "
            f"{mp.nstr(r['T_high'], 6):>13} | "
            f"{mp.nstr(r['D0_sq_T_N'], 6):>14} | "
            f"{mp.nstr(r['S_coord'], 6):>14} | "
            f"{mp.nstr(r['sigma_0'], 5):>8} | "
            f"{mp.nstr(r['sigma_1'], 5):>8} | "
            f"{mp.nstr(r['sigma_T'], 5):>8}"
        )
    print("=" * 116)
    print()

    print("=" * 80)
    print("CELL 68 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell68()
