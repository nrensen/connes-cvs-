#!/usr/bin/env python3
"""
CELL 147 — Continuum Threshold Stability Audit, Asymptotic Box Offset Robustness,
and Variational Lower Bounds.

Target Gate: Gate 1 (Finite-N Spectral Mechanism & Asymptotic Tail Extinction)
Milestone: Milestone M-G1.6 / Boundary Quantization & Resolvent Coercivity
Investigating:
  1. High-N Extension: Dimension sweep extended to N in [48, 56, 64, 72, 80, 88, 96].
  2. Archimedean Cutoff Robustness: Dual audit across T in {600, 800}.
  3. Threshold Stability: Testing whether candidate threshold E_infty ~ 0.002461
     survives T-variation and higher N, or if it drifts toward zero.
  4. Local Exponent Descent: Tracking whether a_loc(E_11) continues to descend
     below 0.6194 toward zero, confirming crossover to a positive constant offset.
  5. Resolvent Denominator Floor: Tracking D_10(N) = (E_11 - E_2)(E_11 - E_3)
     against the candidate asymptotic floor E_infty^2.

Pre-Flight Invariants (N = 64, T = 600, 50 dps):
  omega_0 = 2.9315259531 (residual < 1e-8)
  nu_0    = 4.2604954421 (residual < 1e-8)
  mu_0    = -0.4869792197 (residual < 1e-8)
"""

import time
import mpmath as mp

from cell import (
    h_plus,
    get_galerkin_matrix,
    canonical_even_projector,
    load_prime_powers_table,
    build_W_tilde,
    build_D_tilde_per,
    build_Delta_D_tilde_closed,
    symmetric_eigendecomposition,
    linear_fit,
    build_continuum_operators,
)


# =============================================================================
# Precision & Physical Parameters
# =============================================================================
mp.mp.dps = 50

C_PARAM = 13
L_PARAM = mp.log(mp.mpf(C_PARAM))
N_BOUND = 11  # 12th state is index 11

SWEEP_N = [48, 56, 64, 72, 80, 88, 96]
CUTOFF_T_LIST = [600, 800]


# =============================================================================
# Pre-Flight Continuum Hard Invariant Check (N = 64, T = 600)
# =============================================================================
def preflight_continuum_audit_64(sys_64: dict, prime_data: list) -> tuple:
    """
    Audit continuum competition invariants at N = 64 against certified baselines.
    """
    ops = build_continuum_operators(
        sys_64["V_Qeven"],
        prime_data,
        L=L_PARAM,
        n_bound=N_BOUND,
    )
    return ops["omega_0"], ops["nu_0"], ops["mu_0"]


# =============================================================================
# Eigensystem Builder for Given (N, T)
# =============================================================================
def solve_system_at_N_T(N: int, T_val: int) -> dict:
    t0 = time.time()
    Q_full, _ = get_galerkin_matrix(
        c=C_PARAM,
        N=N,
        T=T_val,
        dps=50,
        verbose=False,
    )
    V_even = canonical_even_projector(N)
    Q_even = V_even.T * Q_full * V_even
    Q_even = mp.mpf("0.5") * (Q_even + Q_even.T)
    evals_e, V_e = symmetric_eigendecomposition(Q_even)

    op_norm = max(abs(e) for e in evals_e)
    t_solve = time.time() - t0

    return {
        "N": N,
        "T": T_val,
        "dim_even": N + 1,
        "evals_Qeven": evals_e,
        "V_Qeven": V_e,
        "op_norm": op_norm,
        "solve_time": t_solve,
    }


# =============================================================================
# Main Execution Suite
# =============================================================================
def run_cell147_suite():
    t_suite_start = time.time()

    print("=" * 80)
    print("CELL 147 — CONTINUUM THRESHOLD STABILITY & ASYMPTOTIC BOX OFFSET AUDIT")
    print(f"Parameters: c = {C_PARAM}, L = {float(L_PARAM):.11f}, dps = 50")
    print(f"Testing Dimensions: N in {SWEEP_N}")
    print(f"Testing Archimedean Cutoffs: T in {CUTOFF_T_LIST}")
    print("Auditing E_infty Stability, Local Exponent Descent, and Variational Floor")
    print("=" * 80)

    prime_data = load_prime_powers_table(C_PARAM)
    print(f"Loaded {len(prime_data)} prime powers up to c = {C_PARAM}.\n")

    # =========================================================================
    # PRE-FLIGHT HARD REGRESSION AUDIT (N = 64, T = 600, 50 dps)
    # =========================================================================
    print("=" * 80)
    print("PRE-FLIGHT HARD REGRESSION AUDIT AGAINST CERTIFIED INVARIANTS (N = 64, T = 600)")
    print("=" * 80)
    t_pre_start = time.time()

    sys_64_600 = solve_system_at_N_T(64, 600)
    omega_0_64, nu_0_64, mu_0_64 = preflight_continuum_audit_64(sys_64_600, prime_data)

    expected_omega_0 = mp.mpf("2.9315259531")
    expected_nu_0 = mp.mpf("4.2604954421")
    expected_mu_0 = mp.mpf("-0.4869792197")

    err_omega = abs(omega_0_64 - expected_omega_0)
    err_nu = abs(nu_0_64 - expected_nu_0)
    err_mu = abs(mu_0_64 - expected_mu_0)

    print(f"  omega_0 = {float(omega_0_64):.10f} (Expected: {float(expected_omega_0):.10f}, Residual: {float(err_omega):.2e})")
    print(f"  nu_0    = {float(nu_0_64):.10f} (Expected: {float(expected_nu_0):.10f}, Residual: {float(err_nu):.2e})")
    print(f"  mu_0    = {float(mu_0_64):.10f} (Expected: {float(expected_mu_0):.10f}, Residual: {float(err_mu):.2e})")

    if err_omega > mp.mpf("1e-8") or err_nu > mp.mpf("1e-8") or err_mu > mp.mpf("1e-8"):
        print("FATAL: REGRESSION AUDIT FAILED AGAINST CERTIFIED INVARIANTS! ABORTING.")
        raise RuntimeError("Operator regression failure between certified baselines and Cell 147.")

    t_pre = time.time() - t_pre_start
    print(f"  REGRESSION AUDIT PASSED in {t_pre:.2f}s: Operators match certified invariants.\n")

    # =========================================================================
    # MULTI-T DIMENSION SWEEPS
    # =========================================================================
    results = {}
    for T_val in CUTOFF_T_LIST:
        print("=" * 80)
        print(f"EXECUTING DIMENSION SWEEP FOR ARCHIMEDEAN CUTOFF T = {T_val}...")
        print("=" * 80)
        results[T_val] = {}
        for N in SWEEP_N:
            t_solve_start = time.time()
            if T_val == 600 and N == 64:
                sys_N = sys_64_600
                label = "reused from pre-flight"
            else:
                sys_N = solve_system_at_N_T(N, T_val)
                label = f"solved in {time.time() - t_solve_start:.2f}s"
            results[T_val][N] = sys_N
            E11 = sys_N["evals_Qeven"][11]
            E2 = sys_N["evals_Qeven"][2]
            E3 = sys_N["evals_Qeven"][3]
            D10 = (E11 - E2) * (E11 - E3)
            print(f"  T = {T_val}, N = {N:2d} | dim = {sys_N['dim_even']:2d} | E_11 = {mp.nstr(E11, 6)} | D_10 = {mp.nstr(D10, 5)} ({label})")
        print()

    # =========================================================================
    # TABLE 1: Dual-T High-N Continuum Base Trajectory
    # =========================================================================
    print("=" * 80)
    print("TABLE 1: DUAL-T HIGH-N CONTINUUM BASE TRAJECTORY (N in [48..96])")
    print("Tracking E_11(N; T=600) vs E_11(N; T=800) and Relative Cutoff Discrepancy delta_T")
    print("=" * 80)
    print(
        f"{'N':>3} | "
        f"{'E_11(T=600)':>14} | "
        f"{'E_11(T=800)':>14} | "
        f"{'delta_T (Rel Diff)':>18} | "
        f"{'D_10(T=600)':>14} | "
        f"{'D_10(T=800)':>14}"
    )
    print("-" * 85)

    table1_summary = []
    for N in SWEEP_N:
        e11_600 = results[600][N]["evals_Qeven"][11]
        e11_800 = results[800][N]["evals_Qeven"][11]
        rel_diff = abs(e11_800 - e11_600) / e11_600

        e2_600 = results[600][N]["evals_Qeven"][2]
        e3_600 = results[600][N]["evals_Qeven"][3]
        d10_600 = (e11_600 - e2_600) * (e11_600 - e3_600)

        e2_800 = results[800][N]["evals_Qeven"][2]
        e3_800 = results[800][N]["evals_Qeven"][3]
        d10_800 = (e11_800 - e2_800) * (e11_800 - e3_800)

        print(
            f"{N:3d} | "
            f"{mp.nstr(e11_600, 6):>14} | "
            f"{mp.nstr(e11_800, 6):>14} | "
            f"{mp.nstr(rel_diff, 4):>18} | "
            f"{mp.nstr(d10_600, 5):>14} | "
            f"{mp.nstr(d10_800, 5):>14}"
        )
        table1_summary.append({
            "N": N,
            "e11_600": e11_600,
            "e11_800": e11_800,
            "rel_diff": rel_diff,
            "d10_600": d10_600,
            "d10_800": d10_800,
        })
    print("-" * 85)
    print("Observation: Check whether relative discrepancy delta_T remains small and stable as N -> oo.\n")

    # =========================================================================
    # TABLE 2: Extended Local Exponent Deceleration Table
    # =========================================================================
    print("=" * 80)
    print("TABLE 2: EXTENDED LOCAL EXPONENT DECELERATION TABLE (N in [48..96])")
    print("Tracking a_loc(E_11) and p_D_loc across Consecutive Dimensions for T = 600 and T = 800")
    print("=" * 80)
    print(
        f"{'N_1 -> N_2':>11} | "
        f"{'a_loc (T=600)':>14} | "
        f"{'a_loc (T=800)':>14} | "
        f"{'p_D_loc (T=600)':>16} | "
        f"{'p_D_loc (T=800)':>16}"
    )
    print("-" * 85)

    last_a_loc_600 = None
    last_a_loc_800 = None

    for i in range(len(SWEEP_N) - 1):
        N1 = SWEEP_N[i]
        N2 = SWEEP_N[i + 1]
        log_ratio = mp.log(mp.mpf(N2) / mp.mpf(N1))

        e11_1_600 = results[600][N1]["evals_Qeven"][11]
        e11_2_600 = results[600][N2]["evals_Qeven"][11]
        a_loc_600 = -mp.log(e11_2_600 / e11_1_600) / log_ratio

        e11_1_800 = results[800][N1]["evals_Qeven"][11]
        e11_2_800 = results[800][N2]["evals_Qeven"][11]
        a_loc_800 = -mp.log(e11_2_800 / e11_1_800) / log_ratio

        d10_1_600 = table1_summary[i]["d10_600"]
        d10_2_600 = table1_summary[i + 1]["d10_600"]
        p_d_600 = -mp.log(d10_2_600 / d10_1_600) / log_ratio

        d10_1_800 = table1_summary[i]["d10_800"]
        d10_2_800 = table1_summary[i + 1]["d10_800"]
        p_d_800 = -mp.log(d10_2_800 / d10_1_800) / log_ratio

        pair_label = f"{N1:2d} -> {N2:2d}"
        print(
            f"{pair_label:>11} | "
            f"{mp.nstr(a_loc_600, 4):>14} | "
            f"{mp.nstr(a_loc_800, 4):>14} | "
            f"{mp.nstr(p_d_600, 4):>16} | "
            f"{mp.nstr(p_d_800, 4):>16}"
        )
        last_a_loc_600 = a_loc_600
        last_a_loc_800 = a_loc_800

    print("-" * 85)
    print("Observation: Confirm whether a_loc drops below 0.6194 toward zero,")
    print("as predicted when the eigenvalue sequence crosses over to a positive offset E_infty > 0.\n")

    # =========================================================================
    # TABLE 3: Candidate Asymptotic Models & Threshold Stability
    # =========================================================================
    print("=" * 80)
    print("TABLE 3: ASYMPTOTIC MODEL SELECTION & THRESHOLD STABILITY (N in [48..96])")
    print("Comparing Free Power Law vs Offset Box Model across T in {600, 800}")
    print("=" * 80)

    log_N = [mp.log(N) for N in SWEEP_N]
    inv_N2 = [1 / (mp.mpf(N) ** 2) for N in SWEEP_N]

    models_extracted = {}
    for T_val in CUTOFF_T_LIST:
        e11_vals = [results[T_val][N]["evals_Qeven"][11] for N in SWEEP_N]
        log_e11 = [mp.log(e) for e in e11_vals]

        # Model 1: Free Power Law E_11 ~ C * N^(-a)
        slope_pow, inter_pow, r2_pow = linear_fit(log_N, log_e11)
        a_pow = -slope_pow
        C_pow = mp.exp(inter_pow)

        # Model 2: Offset Box Model E_11 ~ E_infty + C_box / N^2
        slope_box, inter_box, r2_box = linear_fit(inv_N2, e11_vals)
        C_box = slope_box
        E_infty = inter_box

        models_extracted[T_val] = {
            "a_pow": a_pow,
            "C_pow": C_pow,
            "r2_pow": r2_pow,
            "C_box": C_box,
            "E_infty": E_infty,
            "r2_box": r2_box,
        }

    print(f"{'Cutoff T':<10} | {'Model':<25} | {'Fitted Parameters':<32} | {'Fit Quality R^2':>15}")
    print("-" * 90)
    for T_val in CUTOFF_T_LIST:
        m = models_extracted[T_val]
        pow_str = f"E_11 ~ {mp.nstr(m['C_pow'], 4)} * N^{mp.nstr(-m['a_pow'], 4)}"
        box_str = f"E_11 ~ {mp.nstr(m['E_infty'], 5)} + {mp.nstr(m['C_box'], 4)}/N^2"

        print(
            f"T = {T_val:<6} | "
            f"{'1. Free Power Law':<25} | "
            f"{pow_str:<32} | "
            f"{mp.nstr(m['r2_pow'], 6):>15}"
        )
        print(
            f"T = {T_val:<6} | "
            f"{'2. Offset Box Model':<25} | "
            f"{box_str:<32} | "
            f"{mp.nstr(m['r2_box'], 6):>15}"
        )
        print("-" * 90)

    print()

    # =========================================================================
    # TABLE 4: Resolvent Denominator & Gate 1 Product Floor
    # =========================================================================
    print("=" * 80)
    print("TABLE 4: RESOLVENT DENOMINATOR AND GATE 1 PRODUCT FLOOR")
    print("Tracking D_10(N), ||Q_even||_op, and R_spec(N) against Asymptotic Candidate Floor E_infty^2")
    print("=" * 80)
    print(
        f"{'T':>4} | {'N':>3} | "
        f"{'D_10(N)':>12} | "
        f"{'E_infty^2 Floor':>16} | "
        f"{'D_10 / E_infty^2':>16} | "
        f"{'||Q_even||_op':>14} | "
        f"{'R_spec, 10':>12}"
    )
    print("-" * 85)

    for T_val in CUTOFF_T_LIST:
        e_infty = models_extracted[T_val]["E_infty"]
        floor_val = e_infty ** 2 if e_infty > 0 else mp.mpf("0")

        for N in SWEEP_N:
            sys_N = results[T_val][N]
            evals_e = sys_N["evals_Qeven"]
            e11 = evals_e[11]
            e2 = evals_e[2]
            e3 = evals_e[3]
            d10 = (e11 - e2) * (e11 - e3)
            ratio_floor = d10 / floor_val if floor_val > 0 else mp.mpf("0")
            op_norm = sys_N["op_norm"]
            r_spec = op_norm / d10 if d10 > 0 else mp.mpf("0")

            print(
                f"{T_val:4d} | {N:3d} | "
                f"{mp.nstr(d10, 5):>12} | "
                f"{mp.nstr(floor_val, 5):>16} | "
                f"{mp.nstr(ratio_floor, 4):>16} | "
                f"{mp.nstr(op_norm, 5):>14} | "
                f"{mp.nstr(r_spec, 4):>12}"
            )
        print("-" * 85)

    print()

    # =========================================================================
    # SYNTHESIS & KEY OBSERVATIONS
    # =========================================================================
    print("=" * 80)
    print("SYNTHESIS & KEY OBSERVATIONS:")
    print(f"  1. T-Cutoff Robustness: Compare E_infty(600) = {mp.nstr(models_extracted[600]['E_infty'], 5)} vs E_infty(800) = {mp.nstr(models_extracted[800]['E_infty'], 5)}.")
    print(f"  2. Box Coefficient Stability: C_box(600) = {mp.nstr(models_extracted[600]['C_box'], 4)} vs C_box(800) = {mp.nstr(models_extracted[800]['C_box'], 4)}.")
    print(f"  3. Local Exponent Descent: a_loc reaches {mp.nstr(last_a_loc_600, 4)} (T=600) and {mp.nstr(last_a_loc_800, 4)} (T=800) by N=96.")
    print("  4. Epistemic Verdict: If E_infty > 0 is verified to survive cutoff variation, Rayleigh-Ritz monotonicity")
    print("     E_11(N) >= E_11(oo) = E_infty unconditionally guarantees a non-vanishing floor D(N) >= E_infty^2 > 0.")
    print("=" * 80)

    t_suite = time.time() - t_suite_start
    print(f"Total suite execution time: {t_suite:.2f}s")
    print("=" * 80)
    print("CELL 147 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell147_suite()
