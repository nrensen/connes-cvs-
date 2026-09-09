#!/usr/bin/env python3
# =============================================================================
# cell95.py — Archimedean Cutoff Calibration & High-T Spectral Recovery Sweep
#             (T in {200, 400, 600, 800, 1000, 1200, 1600} at N in {160, 176, 192})
#
# Milestone M41 — Connes–van Suijlekom Galerkin Positivity Programme
# =============================================================================
#
# RATIONALE & CONTEXT:
# In Cell 94 (Milestone M40, T = 400), the 11|12 boundary gap g_11 = E_12 - E_11
# exhibited remarkable stability across N in {64, ..., 160}, declining by only
# 1.56% (from 0.42046 to 0.41390) with local logarithmic slopes s_11 <= 0.030
# and Cauchy projector convergence ||Delta P_11||_op -> 0.00328.
#
# However, at N = 176, E_11 and E_12 abruptly collapsed by 24 orders of magnitude,
# and at N = 192 a spurious negative block appeared (E_11 = -0.185), accompanied
# by complete decoupling of the projector geometry (||Delta P_11||_op = 1.0,
# sin theta_max = 1.0).
#
# MATHEMATICAL DIAGNOSIS (THE RESONANCE FRONTIER):
# The Archimedean Mellin integral in psi_arch(x) splits at kinks tau = +/- alpha_x,
# where alpha_x = 2*pi*x / L. For cutoff c = 13 and L = log(13) ~ 2.56495,
# the mode frequency is alpha_x ~ 2.44963 * x.
#
# For fixed T = 400, the maximum resolvable dimension is:
#     N_Nyquist(T=400) = 400 * log(13) / (2*pi) ~ 163.29.
#
# For N <= 160, all alpha_x <= 391.94 < 400 lie strictly INSIDE [-T, T].
# For N = 176, alpha_176 ~ 431.14 > 400 lies OUTSIDE [-T, T].
# The integration interval [-400, 400] truncates before reaching the kernel's
# resonance, generating severe cutoff leakage that pollutes the deep spectrum.
#
# EXPERIMENTAL OBJECTIVES OF CELL 95:
# 1. TEST A (Control Stability at N = 160):
#    Audit whether E_11 ~ 0.003975 and E_12 ~ 0.417871 (g_11 ~ 0.4139) remain
#    stable as T increases from 400 to 1600. Test whether T = 200 (< alpha_160)
#    exhibits cutoff artifact.
# 2. TEST B (High-T Recovery at N = 176):
#    Audit whether increasing T >= 600 past alpha_176 ~ 431.14 eliminates the
#    collapse and recovers the 0.41-scale boundary separation g_11 ~ 0.41.
# 3. TEST C (High-T Recovery at N = 192):
#    Audit whether increasing T >= 600 past alpha_192 ~ 470.33 eliminates the
#    spurious negative block (E_11 = -0.185) and restores positive gap g_11 ~ 0.41.
# 4. TEST D (Resolution Frontier Scaling T_req(N)):
#    Map the empirical resolution threshold T_min(N) against the theoretical
#    resonance threshold alpha_N = 2*pi*N / log(13).
#
# =============================================================================

from __future__ import annotations

import time
import mpmath as mp

from cell import get_galerkin_matrix

# -----------------------------------------------------------------------------
# Precision and Parameter Configuration
# -----------------------------------------------------------------------------

# Working precision of 70 decimal digits
mp.mp.dps = 70

C_PARAM = 13
L_PARAM = mp.log(C_PARAM)
GROUND_DPS = 50

# Sweep configurations
T_LIST = [200, 400, 600, 800, 1000, 1200, 1600]
N_LIST = [160, 176, 192]
MODES_EXTRACT = [9, 10, 11, 12, 13, 14]


# -----------------------------------------------------------------------------
# Resonance Frequency Helper
# -----------------------------------------------------------------------------

def resonance_frequency(N: int) -> mp.mpf:
    """Return alpha_N = 2*pi*N / log(c), the critical Archimedean frequency."""
    return 2 * mp.pi * N / L_PARAM


# -----------------------------------------------------------------------------
# Parity Basis & Eigensystem Solver
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


def solve_parity_eigensystem(
    Q_full: mp.matrix, N: int
) -> tuple[mp.matrix, list[mp.mpf], mp.matrix]:
    """Compute even spectrum, eigenvectors, and projected matrix."""
    E, _ = full_parity_basis(N)

    Q_even = E.T * Q_full * E
    Q_even = mp.mpf('0.5') * (Q_even + Q_even.T)

    evals_e, V_e = mp.eigsy(Q_even)

    idx_e = sorted(range(N + 1), key=lambda i: evals_e[i])
    sorted_evals_e = [evals_e[i] for i in idx_e]
    sorted_V_e = mp.matrix(N + 1, N + 1)
    for col_idx, orig_col in enumerate(idx_e):
        for row_idx in range(N + 1):
            sorted_V_e[row_idx, col_idx] = V_e[row_idx, orig_col]

    return Q_even, sorted_evals_e, sorted_V_e


# -----------------------------------------------------------------------------
# Main Execution Function
# -----------------------------------------------------------------------------

def run_cell95() -> None:
    print("=" * 140)
    print("CELL 95 — ARCHIMEDEAN CUTOFF CALIBRATION & HIGH-T SPECTRAL RECOVERY SWEEP")
    print("          (T in {200, 400, 600, 800, 1000, 1200, 1600} at N in {160, 176, 192})")
    print("=" * 140)
    print(f"Configuration: c = {C_PARAM}, L = log(c) = {mp.nstr(L_PARAM, 12)}")
    print(f"Working Precision: {mp.mp.dps} decimal digits (generation dps = {GROUND_DPS})")
    print(f"Archimedean Cutoff Sweep: T in {T_LIST}")
    print(f"Target Dimensions: N in {N_LIST}")
    print(f"Resonance Frequencies alpha_N = 2*pi*N/L:")
    for N in N_LIST:
        alpha = resonance_frequency(N)
        print(f"   N = {N:3d} : alpha_{N} = {mp.nstr(alpha, 7)}  (T_min ~ {int(mp.ceil(alpha))})")
    print("=" * 140)

    # Dictionary: (T, N) -> result record
    sweep_data: dict[tuple[int, int], dict] = {}
    total_start_time = time.time()

    # Outer loop over T ensures that for each T, the basis evaluations
    # psi(n) are cached incrementally from N = 160 -> 176 -> 192.
    for T in T_LIST:
        print(f"\n####################################################################################################")
        print(f">>> STARTING ARCHIMEDEAN CUTOFF T = {T} ...")
        print(f"####################################################################################################")
        t_block_start = time.time()

        for N in N_LIST:
            step_start = time.time()
            alpha_N = resonance_frequency(N)
            is_above_resonance = (T > alpha_N)

            print(f"\n   >>> Processing (T = {T}, N = {N}) [dim = {N+1}, alpha_{N} = {mp.nstr(alpha_N, 6)}, T > alpha: {is_above_resonance}] ...")

            Q_full, _ = get_galerkin_matrix(
                c=C_PARAM,
                N=N,
                T=T,
                dps=GROUND_DPS,
                verbose=False,
            )

            Q_even, evals_e, V_e = solve_parity_eigensystem(Q_full, N)
            step_elapsed = time.time() - step_start

            # Extract low-energy spectrum
            extracted_evals = {k: evals_e[k] for k in MODES_EXTRACT if k < len(evals_e)}

            E_10 = extracted_evals.get(10, mp.mpf('nan'))
            E_11 = extracted_evals.get(11, mp.mpf('nan'))
            E_12 = extracted_evals.get(12, mp.mpf('nan'))
            E_13 = extracted_evals.get(13, mp.mpf('nan'))
            E_14 = extracted_evals.get(14, mp.mpf('nan'))

            g_10 = E_11 - E_10
            g_11 = E_12 - E_11
            g_12 = E_13 - E_12
            g_13 = E_14 - E_13

            # Classification of spectral status
            has_negative = any(evals_e[k] < -1e-12 for k in range(min(15, len(evals_e))))
            is_deep_collapsed = (abs(E_12) < 1e-10)
            is_recovered = (g_11 > mp.mpf('0.35')) and not has_negative

            if has_negative:
                status = "NEGATIVE_BLOCK_POLLUTED"
            elif is_deep_collapsed:
                status = "DEEP_SPECTRUM_COLLAPSED"
            elif is_recovered:
                if N == 160:
                    status = "CONTROL_STABLE"
                else:
                    status = "RECOVERED (g11 ~ 0.41)"
            else:
                status = "UNRESOLVED_TRANSIENT"

            # Compute prominence ratio if denominators permit
            denom = max(abs(g_10), abs(g_12))
            gamma_11 = (g_11 / denom) if denom > 1e-60 else mp.mpf('nan')

            rec = {
                "T": T,
                "N": N,
                "alpha_N": alpha_N,
                "is_above_resonance": is_above_resonance,
                "evals": extracted_evals,
                "min_eval": evals_e[0],
                "E_10": E_10,
                "E_11": E_11,
                "E_12": E_12,
                "E_13": E_13,
                "g_10": g_10,
                "g_11": g_11,
                "g_12": g_12,
                "g_13": g_13,
                "gamma_11": gamma_11,
                "status": status,
                "time": step_elapsed,
            }
            sweep_data[(T, N)] = rec

            print(f"      Completed (T={T}, N={N}) in {step_elapsed:.2f}s | E_11 = {mp.nstr(E_11, 8)} | E_12 = {mp.nstr(E_12, 8)} | g_11 = {mp.nstr(g_11, 8)} | Status: {status}")

        t_block_elapsed = time.time() - t_block_start
        print(f"   Finished T = {T} in {t_block_elapsed:.2f}s.")

    total_elapsed = time.time() - total_start_time
    print(f"\nAll sweeps completed in {total_elapsed:.2f}s ({total_elapsed/60:.2f}m).")

    # =========================================================================
    # PRIMARY EXECUTIVE SUMMARY TABLE: (T, N) GRID OF BOUNDARY RESOLUTION
    # =========================================================================
    print("\n" + "#" * 140)
    print("PRIMARY EXECUTIVE SUMMARY TABLE: ARCHIMEDEAN CUTOFF T vs BOUNDARY RESOLUTION")
    print("Columns: [T, N, alpha_N, T > alpha_N, E_11, E_12, g_10, g_11, g_12, Status]")
    print("#" * 140)
    header = f"{'T':>5} | {'N':>4} | {'alpha_N':>7} | {'> alpha?':>8} | {'E_{11}':>14} | {'E_{12}':>14} | {'g_{10}':>14} | {'g_{11}':>14} | {'g_{12}':>14} | {'Status'}"
    print(header)
    print("-" * 140)

    for T in T_LIST:
        for N in N_LIST:
            rec = sweep_data[(T, N)]
            e11_s = mp.nstr(rec['E_11'], 7)
            e12_s = mp.nstr(rec['E_12'], 7)
            g10_s = mp.nstr(rec['g_10'], 7)
            g11_s = mp.nstr(rec['g_11'], 7)
            g12_s = mp.nstr(rec['g_12'], 7)
            alpha_s = mp.nstr(rec['alpha_N'], 5)
            above_s = "YES" if rec['is_above_resonance'] else "NO"

            print(f"{T:5d} | {N:4d} | {alpha_s:>7} | {above_s:>8} | {e11_s:>14} | {e12_s:>14} | {g10_s:>14} | {g11_s:>14} | {g12_s:>14} | {rec['status']}")
        print("-" * 140)
    print("#" * 140)

    # =========================================================================
    # TEST A: Control Stability at N = 160 across Archimedean Cutoffs T
    # =========================================================================
    print("\n" + "=" * 140)
    print("TEST A: CONTROL STABILITY AT N = 160 ACROSS ARCHIMEDEAN CUTOFFS T in {200, ..., 1600}")
    print("Resonance Frequency alpha_160 ~ 391.94 (T = 200 is UNDER-RESOLVED; T >= 400 is RESOLVED)")
    print("=" * 140)
    print(f"{'T':>5} | {'T > alpha?':>10} | {'E_{10}':>14} | {'E_{11}':>14} | {'E_{12}':>14} | {'g_{10}':>14} | {'g_{11}':>14} | {'g_{12}':>14} | {'Status'}")
    print("-" * 140)

    ref_g11_160 = sweep_data.get((400, 160), {}).get('g_11', mp.mpf('nan'))

    for T in T_LIST:
        rec = sweep_data[(T, 160)]
        e10_s = mp.nstr(rec['E_10'], 7)
        e11_s = mp.nstr(rec['E_11'], 7)
        e12_s = mp.nstr(rec['E_12'], 7)
        g10_s = mp.nstr(rec['g_10'], 7)
        g11_s = mp.nstr(rec['g_11'], 7)
        g12_s = mp.nstr(rec['g_12'], 7)
        above_s = "YES" if rec['is_above_resonance'] else "NO"

        print(f"{T:5d} | {above_s:>10} | {e10_s:>14} | {e11_s:>14} | {e12_s:>14} | {g10_s:>14} | {g11_s:>14} | {g12_s:>14} | {rec['status']}")

    # =========================================================================
    # TEST B: High-T Recovery at N = 176 across Archimedean Cutoffs T
    # =========================================================================
    print("\n" + "=" * 140)
    print("TEST B: HIGH-T SPECTRAL RECOVERY AT N = 176 ACROSS ARCHIMEDEAN CUTOFFS T in {200, ..., 1600}")
    print("Resonance Frequency alpha_176 ~ 431.14 (T in {200, 400} are UNDER-RESOLVED; T >= 600 is RESOLVED)")
    print("=" * 140)
    print(f"{'T':>5} | {'T > alpha?':>10} | {'E_{10}':>14} | {'E_{11}':>14} | {'E_{12}':>14} | {'g_{10}':>14} | {'g_{11}':>14} | {'g_{12}':>14} | {'Status'}")
    print("-" * 140)

    for T in T_LIST:
        rec = sweep_data[(T, 176)]
        e10_s = mp.nstr(rec['E_10'], 7)
        e11_s = mp.nstr(rec['E_11'], 7)
        e12_s = mp.nstr(rec['E_12'], 7)
        g10_s = mp.nstr(rec['g_10'], 7)
        g11_s = mp.nstr(rec['g_11'], 7)
        g12_s = mp.nstr(rec['g_12'], 7)
        above_s = "YES" if rec['is_above_resonance'] else "NO"

        print(f"{T:5d} | {above_s:>10} | {e10_s:>14} | {e11_s:>14} | {e12_s:>14} | {g10_s:>14} | {g11_s:>14} | {g12_s:>14} | {rec['status']}")

    # =========================================================================
    # TEST C: High-T Recovery at N = 192 across Archimedean Cutoffs T
    # =========================================================================
    print("\n" + "=" * 140)
    print("TEST C: HIGH-T SPECTRAL RECOVERY AT N = 192 ACROSS ARCHIMEDEAN CUTOFFS T in {200, ..., 1600}")
    print("Resonance Frequency alpha_192 ~ 470.33 (T in {200, 400} are UNDER-RESOLVED; T >= 600 is RESOLVED)")
    print("=" * 140)
    print(f"{'T':>5} | {'T > alpha?':>10} | {'E_{10}':>14} | {'E_{11}':>14} | {'E_{12}':>14} | {'g_{10}':>14} | {'g_{11}':>14} | {'g_{12}':>14} | {'Status'}")
    print("-" * 140)

    for T in T_LIST:
        rec = sweep_data[(T, 192)]
        e10_s = mp.nstr(rec['E_10'], 7)
        e11_s = mp.nstr(rec['E_11'], 7)
        e12_s = mp.nstr(rec['E_12'], 7)
        g10_s = mp.nstr(rec['g_10'], 7)
        g11_s = mp.nstr(rec['g_11'], 7)
        g12_s = mp.nstr(rec['g_12'], 7)
        above_s = "YES" if rec['is_above_resonance'] else "NO"

        print(f"{T:5d} | {above_s:>10} | {e10_s:>14} | {e11_s:>14} | {e12_s:>14} | {g10_s:>14} | {g11_s:>14} | {g12_s:>14} | {rec['status']}")

    # =========================================================================
    # TEST D: Resolution Scaling Frontier T_req(N) vs Resonance Threshold alpha_N
    # =========================================================================
    print("\n" + "=" * 140)
    print("TEST D: RESOLUTION FRONTIER T_req(N) vs RESONANCE THRESHOLD alpha_N = 2*pi*N / log(13)")
    print("Auditing Minimum Archimedean Cutoff T_min Required for g_11 >= 0.35 and Positive Spectrum")
    print("=" * 140)
    print(f"{'N':>4} | {'alpha_N':>10} | {'Min T Tested':>12} | {'T_resolved':>12} | {'Safety Margin (T - alpha)/alpha':>32} | {'g_11 at T_resolved':>20}")
    print("-" * 140)

    for N in N_LIST:
        alpha_N = resonance_frequency(N)
        # Find smallest T where status is resolved
        resolved_T = None
        resolved_g11 = None
        for T in T_LIST:
            rec = sweep_data[(T, N)]
            if "RECOVERED" in rec['status'] or "CONTROL_STABLE" in rec['status']:
                resolved_T = T
                resolved_g11 = rec['g_11']
                break

        alpha_s = mp.nstr(alpha_N, 6)
        min_t_s = str(T_LIST[0])
        if resolved_T is not None:
            res_t_s = str(resolved_T)
            margin = (resolved_T - alpha_N) / alpha_N
            margin_s = f"{mp.nstr(margin * 100, 4)}%"
            g11_s = mp.nstr(resolved_g11, 7)
        else:
            res_t_s = "NOT_RECOVERED"
            margin_s = "—"
            g11_s = "—"

        print(f"{N:4d} | {alpha_s:>10} | {min_t_s:>12} | {res_t_s:>12} | {margin_s:>32} | {g11_s:>20}")

    print("\n" + "=" * 80)
    print("CELL 95 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell95()
