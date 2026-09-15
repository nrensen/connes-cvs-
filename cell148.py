#!/usr/bin/env python3
"""
CELL 148 — Variational Coercivity, Component Energy Dissection of the Threshold State,
and Analytical Lower Bounds on B_11^perp.

Target Gate: Gate 1 (Finite-N Spectral Mechanism & Asymptotic Tail Extinction)
Milestone: Milestone M-G1.6 / Variational Lower Bound & Resolvent Coercivity
Investigating:
  1. Exact Arithmetic Component Dissection:
     E_11(N) = E_arch[v_11] + E_prime[v_11] + E_pole[v_11], with E_prime = -W + D_trans.
  2. Competition Energy and Pareto Deficit Geometry:
     K_rest[v_11], W_perp[v_11], H_1[v_11], and Pareto deficit pair (Delta K, Delta W).
  3. Modal Mass Distribution on W_perp Eigenmodes:
     Cluster mass P_<=3 vs bulk mass P_>=4, testing where the threshold state lives.
  4. Coercivity Components & Scaling:
     Tracking individual component scaling exponents to isolate the analytical mechanism
     preserving E_11 >= E_infty > 0.

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
    build_Q_prime_even,
    build_Q_pole_even,
    symmetric_eigendecomposition,
    extract_continuum_projector,
    build_continuum_operators,
)


# =============================================================================
# Precision & Physical Parameters
# =============================================================================
mp.mp.dps = 50

C_PARAM = 13
L_PARAM = mp.log(mp.mpf(C_PARAM))
N_BOUND = 11  # 12th state is index 11
T_PARAM = 600

SWEEP_N = [48, 56, 64, 72, 80, 88, 96]


# =============================================================================
# Pre-Flight Hard Regression Audit (N = 64, T = 600)
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
# Eigensystem & Component Builder
# =============================================================================
def build_full_dissection_at_N(N: int, T_val: int, prime_data: list) -> dict:
    t0 = time.time()
    dim_even = N + 1
    q_cont = N - N_BOUND + 1
    PI = mp.pi

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

    # Component matrices
    W_tilde = build_W_tilde(N, prime_data, L_PARAM)
    D_per = build_D_tilde_per(N, prime_data, L_PARAM)
    Delta_D = build_Delta_D_tilde_closed(N, prime_data, L_PARAM)
    D_trans = D_per + Delta_D
    Q_prime_even = build_Q_prime_even(N, prime_data, L_PARAM)

    Q_pole_even = build_Q_pole_even(N, L_PARAM)
    Q_arch_even = Q_even - Q_prime_even - Q_pole_even
    Q_arch_even = mp.mpf("0.5") * (Q_arch_even + Q_arch_even.T)

    # Continuum subspace projector
    U_cont = extract_continuum_projector(V_e, n_bound=N_BOUND)

    W_hat_perp = U_cont.T * W_tilde * U_cont
    W_hat_perp = mp.mpf("0.5") * (W_hat_perp + W_hat_perp.T)

    Omega_diag = mp.matrix(dim_even, dim_even)
    for m in range(dim_even):
        a_m = mp.mpf("2") * PI * mp.mpf(m) / L_PARAM if m > 0 else mp.mpf("0")
        h_val = h_plus(a_m)
        Omega_diag[m, m] = h_val + D_per[m, m]

    K_rest = U_cont.T * (Omega_diag + Delta_D) * U_cont
    K_rest = mp.mpf("0.5") * (K_rest + K_rest.T)

    Q_hat_comp = K_rest - W_hat_perp
    Q_hat_comp = mp.mpf("0.5") * (Q_hat_comp + Q_hat_comp.T)

    # W_hat_perp eigensystem (sorted descending for near-degenerate cluster analysis)
    evals_W_desc, V_W_desc = symmetric_eigendecomposition(W_hat_perp, sort_descending=True)

    # Threshold state v_11
    v_11 = [V_e[r, N_BOUND] for r in range(dim_even)]
    v_11_mat = mp.matrix(dim_even, 1)
    for r in range(dim_even):
        v_11_mat[r, 0] = v_11[r]

    # Quadratic forms on v_11
    E_11 = evals_e[N_BOUND]
    E_arch = float((v_11_mat.T * Q_arch_even * v_11_mat)[0, 0])
    E_prime = float((v_11_mat.T * Q_prime_even * v_11_mat)[0, 0])
    E_pole = float((v_11_mat.T * Q_pole_even * v_11_mat)[0, 0])
    W_val = float((v_11_mat.T * W_tilde * v_11_mat)[0, 0])
    D_trans_val = float((v_11_mat.T * D_trans * v_11_mat)[0, 0])

    # Competition and Pareto coordinates
    # In U_cont coordinates, v_11 is precisely the first basis vector e_0 = (1, 0, ..., 0)^T
    K_rest_val = float(K_rest[0, 0])
    W_perp_val = float(W_hat_perp[0, 0])
    H1_val = float(Q_hat_comp[0, 0])

    evals_K, _ = symmetric_eigendecomposition(K_rest)
    omega_0 = float(evals_K[0])
    nu_0 = float(evals_W_desc[0])
    evals_Qcomp, _ = symmetric_eigendecomposition(Q_hat_comp)
    mu_0 = float(evals_Qcomp[0])

    delta_K = K_rest_val - omega_0
    delta_W = nu_0 - W_perp_val
    delta_tot = delta_K + delta_W
    c_gain = mu_0 - (omega_0 - nu_0)
    pareto_margin = delta_tot - c_gain

    # Projections of v_11 onto W_perp eigenbasis
    # In U_cont, v_11 is e_0. Its inner product with column k of V_W_desc is V_W_desc[0, k]
    P_W = [float(V_W_desc[0, k] ** 2) for k in range(q_cont)]
    P_cluster = sum(P_W[:4])  # top 4 modes (0, 1, 2, 3)
    P_bulk = sum(P_W[4:])     # deeper modes

    t_solve = time.time() - t0

    return {
        "N": N,
        "dim_even": dim_even,
        "q_cont": q_cont,
        "E_11": float(E_11),
        "E_arch": E_arch,
        "E_prime": E_prime,
        "E_pole": E_pole,
        "W_val": W_val,
        "D_trans_val": D_trans_val,
        "closure_res": abs(float(E_11) - (E_arch + E_prime + E_pole)),
        "K_rest_val": K_rest_val,
        "W_perp_val": W_perp_val,
        "H1_val": H1_val,
        "omega_0": omega_0,
        "nu_0": nu_0,
        "mu_0": mu_0,
        "delta_K": delta_K,
        "delta_W": delta_W,
        "delta_tot": delta_tot,
        "c_gain": c_gain,
        "pareto_margin": pareto_margin,
        "P_cluster": P_cluster,
        "P_bulk": P_bulk,
        "P_top": P_W[0],
        "t_solve": t_solve,
    }


# =============================================================================
# Execution Suite
# =============================================================================
def run_cell148_suite():
    t_suite_start = time.time()

    print("=" * 80)
    print("CELL 148 — VARIATIONAL COERCIVITY & THRESHOLD STATE COMPONENT DISSECTION")
    print(f"Parameters: c = {C_PARAM}, L = {float(L_PARAM):.11f}, T = {T_PARAM}, dps = 50")
    print("Auditing Arithmetic Partition, Competition Balance, and Pareto Deficit on B_11^perp")
    print("=" * 80)

    prime_data = load_prime_powers_table(C_PARAM)
    print(f"Loaded {len(prime_data)} prime powers up to c = {C_PARAM}.\n")

    # =========================================================================
    # PRE-FLIGHT HARD REGRESSION AUDIT (N = 64, T = 600)
    # =========================================================================
    print("=" * 80)
    print("PRE-FLIGHT HARD REGRESSION AUDIT AGAINST CERTIFIED INVARIANTS (N = 64, T = 600)")
    print("=" * 80)
    t_pre_start = time.time()

    # Fast build at N=64 for pre-flight
    Q_full_64, _ = get_galerkin_matrix(c=C_PARAM, N=64, T=600, dps=50, verbose=False)
    V_even_64 = canonical_even_projector(64)
    Q_even_64 = V_even_64.T * Q_full_64 * V_even_64
    Q_even_64 = mp.mpf("0.5") * (Q_even_64 + Q_even_64.T)
    evals_e_64, V_e_64 = symmetric_eigendecomposition(Q_even_64)

    sys_64 = {"V_Qeven": V_e_64}
    omega_0_64, nu_0_64, mu_0_64 = preflight_continuum_audit_64(sys_64, prime_data)

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
        raise RuntimeError("Operator regression failure between certified baselines and Cell 148.")

    t_pre = time.time() - t_pre_start
    print(f"  REGRESSION AUDIT PASSED in {t_pre:.2f}s: Operators match certified invariants.\n")

    # =========================================================================
    # EXECUTE FULL COMPONENT DISSECTION SWEEP ACROSS N
    # =========================================================================
    print(f"Executing operator dissection sweeps across N in {SWEEP_N} at T = {T_PARAM}...")
    results = []
    for N in SWEEP_N:
        res = build_full_dissection_at_N(N, T_PARAM, prime_data)
        results.append(res)
        print(f"  Completed N = {N:2d} (dim = {res['dim_even']:2d}, q = {res['q_cont']:2d}) in {res['t_solve']:5.2f}s | E_11 = {res['E_11']:.6f}")

    print()

    # =========================================================================
    # TABLE 1: EXACT OPERATOR COMPONENT DISSECTION OF E_11(N)
    # =========================================================================
    print("=" * 80)
    print("TABLE 1: EXACT OPERATOR COMPONENT DISSECTION OF THE THRESHOLD STATE E_11(N)")
    print("Tracking E_11 = E_arch + E_prime + E_pole, with E_prime = -W + D_trans")
    print("=" * 80)
    print(
        f"{'N':>3} | "
        f"{'E_11':>10} | "
        f"{'E_arch':>10} | "
        f"{'E_prime':>10} | "
        f"{'E_pole':>10} | "
        f"{'-W[v_11]':>10} | "
        f"{'D_trans':>10} | "
        f"{'Closure Res':>11}"
    )
    print("-" * 80)

    for r in results:
        print(
            f"{r['N']:3d} | "
            f"{r['E_11']:10.6f} | "
            f"{r['E_arch']:10.6f} | "
            f"{r['E_prime']:10.6f} | "
            f"{r['E_pole']:10.6f} | "
            f"{-r['W_val']:10.6f} | "
            f"{r['D_trans_val']:10.6f} | "
            f"{r['closure_res']:11.2e}"
        )

    print("-" * 80)
    print("Observation: Confirm whether E_arch + E_prime < 0 with E_pole providing positive rescue,\n"
          "or if E_arch + E_prime itself provides a strictly positive floor.\n")

    # =========================================================================
    # TABLE 2: COMPETITION OPERATOR & PARETO DEFICIT PAIR ON v_11^(N)
    # =========================================================================
    print("=" * 80)
    print("TABLE 2: COMPETITION OPERATOR & PARETO DEFICIT PAIR ON v_11^(N)")
    print("Tracking K_rest[v_11], W_perp[v_11], H_1[v_11], Deficits Delta K, Delta W, and Pareto Margin")
    print("=" * 80)
    print(
        f"{'N':>3} | "
        f"{'K_rest':>9} | "
        f"{'W_perp':>9} | "
        f"{'H_1[v_11]':>10} | "
        f"{'Delta K':>9} | "
        f"{'Delta W':>9} | "
        f"{'Delta Tot':>10} | "
        f"{'C_gain':>9} | "
        f"{'Margin':>9}"
    )
    print("-" * 80)

    for r in results:
        print(
            f"{r['N']:3d} | "
            f"{r['K_rest_val']:9.5f} | "
            f"{r['W_perp_val']:9.5f} | "
            f"{r['H1_val']:10.5f} | "
            f"{r['delta_K']:9.5f} | "
            f"{r['delta_W']:9.5f} | "
            f"{r['delta_tot']:10.5f} | "
            f"{r['c_gain']:9.5f} | "
            f"{r['pareto_margin']:9.5f}"
        )

    print("-" * 80)
    print("Observation: Check whether v_11 pays a substantial kinetic excess Delta K or well sacrifice Delta W,\n"
          "confirming that the threshold state is strongly penalized away from the individual optima.\n")

    # =========================================================================
    # TABLE 3: SPECTRAL MASS DISTRIBUTION OF v_11 ACROSS W_perp EIGENMODES
    # =========================================================================
    print("=" * 80)
    print("TABLE 3: SPECTRAL MASS DISTRIBUTION OF v_11 ACROSS W_perp EIGENMODES")
    print("Tracking Top Mode Mass P_top, Cluster Mass P_<=3 (Modes 0-3), and Deep Bulk Mass P_>=4")
    print("=" * 80)
    print(
        f"{'N':>3} | "
        f"{'q_cont':>6} | "
        f"{'P_top (Mode 0)':>16} | "
        f"{'P_<=3 (Cluster)':>16} | "
        f"{'P_>=4 (Bulk)':>14} | "
        f"{'W_perp / nu_0':>14}"
    )
    print("-" * 80)

    for r in results:
        ratio_W = r['W_perp_val'] / r['nu_0'] if r['nu_0'] > 0 else 0
        print(
            f"{r['N']:3d} | "
            f"{r['q_cont']:6d} | "
            f"{r['P_top']:16.4%} | "
            f"{r['P_cluster']:16.4%} | "
            f"{r['P_bulk']:14.4%} | "
            f"{ratio_W:14.4f}"
        )

    print("-" * 80)
    print("Observation: Contrast with v_phys (which had 49.3% in cluster, 50.7% in bulk).\n"
          "Examine whether v_11 is concentrated in the bulk or localized in the top modes.\n")

    # =========================================================================
    # TABLE 4: COERCIVITY SURPLUS AND SCALING ANALYSIS
    # =========================================================================
    print("=" * 80)
    print("TABLE 4: COERCIVITY SURPLUS & SCALING COMPONENTS")
    print("Tracking E_11 vs Potential Residual R_pot = D_trans - W and Pole Surplus E_pole")
    print("=" * 80)
    print(
        f"{'N':>3} | "
        f"{'E_11':>10} | "
        f"{'E_arch':>10} | "
        f"{'D_trans - W':>12} | "
        f"{'E_arch + E_prime':>16} | "
        f"{'E_pole':>10} | "
        f"{'H_1 + 0.50':>10}"
    )
    print("-" * 80)

    for r in results:
        r_pot = r['D_trans_val'] - r['W_val']
        e_arch_prime = r['E_arch'] + r['E_prime']
        h1_offset = r['H1_val'] + 0.50
        print(
            f"{r['N']:3d} | "
            f"{r['E_11']:10.6f} | "
            f"{r['E_arch']:10.6f} | "
            f"{r_pot:12.6f} | "
            f"{e_arch_prime:16.6f} | "
            f"{r['E_pole']:10.6f} | "
            f"{h1_offset:10.6f}"
        )

    print("-" * 80)
    print("Observation: Identify which component acts as the definitive positive barrier preventing E_11 -> 0.\n")

    # =========================================================================
    # SYNTHESIS & PRELIMINARY VERDICT
    # =========================================================================
    print("=" * 80)
    print("SYNTHESIS & KEY OBSERVATIONS:")
    print(f"  1. Algebraic Closure: All component decompositions verify with residual < {max(r['closure_res'] for r in results):.2e}.")
    print(f"  2. Threshold Energy Range: E_11 traverses {results[0]['E_11']:.6f} (N=48) down to {results[-1]['E_11']:.6f} (N=96).")
    print(f"  3. Pole vs Archimedean Role: E_pole ranges from {results[0]['E_pole']:.6f} to {results[-1]['E_pole']:.6f}.")
    print(f"  4. Competition Operator: H_1[v_11] evaluated against ground state mu_0 = {results[0]['mu_0']:.6f}.")
    print("=" * 80)

    t_suite = time.time() - t_suite_start
    print(f"Total suite execution time: {t_suite:.2f}s")
    print("=" * 80)
    print("CELL 148 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell148_suite()
