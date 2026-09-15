#!/usr/bin/env python3
"""
CELL 149 — Spectral Deficit Inequality, Master Operator Decomposition,
and Variational Coercivity Certificate on B_11^perp.

Target Gate: Gate 1 (Finite-N Spectral Mechanism & Asymptotic Tail Extinction)
Milestone: Milestone M-G1.6 / Variational Lower Bound & Resolvent Coercivity
Investigating:
  1. Complete Deficit Spectrum of W_perp: Computing descending eigenvalues
     nu_0 >= nu_1 >= ... and modal deficits delta_nu_k = nu_0 - nu_k.
  2. Cluster Gap Floor & Step Gap: Evaluating delta_*^(4) = nu_0 - nu_4 and
     Delta_step = nu_3 - nu_4 across N in [48, 56, 64, 72, 80, 88, 96].
  3. Exact Modal Deficit Conservation: Verifying Delta W[v_11] === Delta W_spec[v_11]
     and establishing the rigorous bulk inequality Delta W[v_11] >= delta_*^(4) P_>=4.
  4. Master Operator Tri-Partition: Verifying the exact operator identity
     Q_even = Q_comp + Delta_arch + Q_pole, where Delta_arch = Q_arch - diag(h_+).
  5. Variational Coercivity Certificate: Evaluating the analytical lower bound
     E_11 >= mu_0 + delta_*^(4) P_>=4 + Delta_arch[v_11] + E_pole[v_11] > 0.

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
# System Solver & Spectral Deficit Dissection at Dimension N
# =============================================================================
def solve_deficit_system_at_N(N: int, T_val: int, prime_data: list) -> dict:
    t0 = time.time()
    dim_even = N + 1
    q_cont = N - N_BOUND + 1
    PI = mp.pi

    # Full Galerkin matrix and canonical even projection
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

    # Diagonal Fourier Archimedean multiplier and competition operator
    Omega_arch = mp.matrix(dim_even, dim_even)
    for m in range(dim_even):
        a_m = mp.mpf("2") * PI * mp.mpf(m) / L_PARAM if m > 0 else mp.mpf("0")
        Omega_arch[m, m] = h_plus(a_m)

    Q_comp_even = Omega_arch + Q_prime_even
    Q_comp_even = mp.mpf("0.5") * (Q_comp_even + Q_comp_even.T)

    # Non-local Archimedean kernel remainder Delta_arch = Q_arch - Omega_arch
    Delta_arch_even = Q_arch_even - Omega_arch
    Delta_arch_even = mp.mpf("0.5") * (Delta_arch_even + Delta_arch_even.T)

    # Threshold state v_11
    v_11_vec = [V_e[r, N_BOUND] for r in range(dim_even)]
    v_11_mat = mp.matrix(dim_even, 1)
    for r in range(dim_even):
        v_11_mat[r, 0] = v_11_vec[r]

    E_11 = evals_e[N_BOUND]
    E_arch = float((v_11_mat.T * Q_arch_even * v_11_mat)[0, 0])
    E_prime = float((v_11_mat.T * Q_prime_even * v_11_mat)[0, 0])
    E_pole = float((v_11_mat.T * Q_pole_even * v_11_mat)[0, 0])

    # Master tri-partition forms on v_11
    H1_val = float((v_11_mat.T * Q_comp_even * v_11_mat)[0, 0])
    Delta_arch_val = float((v_11_mat.T * Delta_arch_even * v_11_mat)[0, 0])
    tri_closure_res = abs(float(E_11) - (H1_val + Delta_arch_val + E_pole))

    # Continuum subspace restriction
    U_cont = extract_continuum_projector(V_e, n_bound=N_BOUND)
    W_hat_perp = U_cont.T * W_tilde * U_cont
    W_hat_perp = mp.mpf("0.5") * (W_hat_perp + W_hat_perp.T)

    H1_perp = U_cont.T * Q_comp_even * U_cont
    H1_perp = mp.mpf("0.5") * (H1_perp + H1_perp.T)
    evals_H1, _ = symmetric_eigendecomposition(H1_perp)
    mu_0 = evals_H1[0]

    # Eigensystem of W_hat_perp (sorted descending: nu_0 >= nu_1 >= ...)
    evals_W_desc, V_W_desc = symmetric_eigendecomposition(W_hat_perp, sort_descending=True)
    nu_0 = evals_W_desc[0]

    # Modal deficits delta_nu_k = nu_0 - nu_k
    delta_nu = [float(nu_0 - evals_W_desc[k]) for k in range(q_cont)]

    # Spectral gap definitions
    nu_top_list = [float(evals_W_desc[k]) for k in range(min(8, q_cont))]
    delta_nu_top = delta_nu[:min(8, q_cont)]
    delta_star_4 = delta_nu[4] if q_cont > 4 else 0.0
    delta_step = float(evals_W_desc[3] - evals_W_desc[4]) if q_cont > 4 else 0.0

    # In continuum coordinates U_cont, v_11 is basis vector e_0 = (1, 0, ..., 0)^T
    # The modal mass along eigenvector y_k of W_perp is P_W(k) = (V_W_desc[0, k])^2
    P_W = [float(V_W_desc[0, k] ** 2) for k in range(q_cont)]
    P_cluster = sum(P_W[:4])
    P_bulk = sum(P_W[4:])

    # Modal deficit conservation
    delta_W_cluster = sum(delta_nu[k] * P_W[k] for k in range(1, min(4, q_cont)))
    delta_W_bulk = sum(delta_nu[k] * P_W[k] for k in range(4, q_cont))
    delta_W_spec = delta_W_cluster + delta_W_bulk

    W_perp_val = float(W_hat_perp[0, 0])
    delta_W_exact = float(nu_0) - W_perp_val
    deficit_res = abs(delta_W_exact - delta_W_spec)

    # Bulk lower bound: Delta W >= delta_*^(4) * P_bulk
    bulk_bound_floor = delta_star_4 * P_bulk
    bulk_bound_satisfied = (delta_W_exact >= bulk_bound_floor - 1e-12)

    # Coercivity certificate floor: E_11 >= mu_0 + delta_*^(4) P_bulk + Delta_arch + E_pole
    # Note: on B_11^perp, H1[v] = H1[v_phys] + excess >= mu_0 + delta_W_excess
    E_11_cert = float(mu_0) + bulk_bound_floor + Delta_arch_val + E_pole

    t_solve = time.time() - t0

    return {
        "N": N,
        "dim_even": dim_even,
        "q_cont": q_cont,
        "E_11": float(E_11),
        "E_arch": E_arch,
        "E_prime": E_prime,
        "E_pole": E_pole,
        "H1_val": H1_val,
        "Delta_arch_val": Delta_arch_val,
        "tri_closure_res": tri_closure_res,
        "mu_0": float(mu_0),
        "nu_0": float(nu_0),
        "nu_top_list": nu_top_list,
        "delta_nu_top": delta_nu_top,
        "delta_star_4": delta_star_4,
        "delta_step": delta_step,
        "P_top": P_W[0],
        "P_cluster": P_cluster,
        "P_bulk": P_bulk,
        "delta_W_exact": delta_W_exact,
        "delta_W_cluster": delta_W_cluster,
        "delta_W_bulk": delta_W_bulk,
        "delta_W_spec": delta_W_spec,
        "deficit_res": deficit_res,
        "bulk_bound_floor": bulk_bound_floor,
        "bulk_bound_satisfied": bulk_bound_satisfied,
        "E_11_cert": E_11_cert,
        "t_solve": t_solve,
    }


# =============================================================================
# Main Execution Suite
# =============================================================================
def run_cell149_suite():
    t_suite_start = time.time()

    print("=" * 80)
    print("CELL 149 — SPECTRAL DEFICIT INEQUALITY & MASTER OPERATOR COERCIVITY")
    print(f"Parameters: c = {C_PARAM}, L = {float(L_PARAM):.11f}, T = {T_PARAM}, dps = 50")
    print("Target: Potential Well Deficit Spectrum, Modal Conservation, and Coercivity Bounds")
    print("=" * 80)

    prime_data = load_prime_powers_table(C_PARAM)
    print(f"Loaded {len(prime_data)} prime powers up to c = {C_PARAM}.\n")

    # -------------------------------------------------------------------------
    # Pre-Flight Hard Regression Audit (N = 64, T = 600)
    # -------------------------------------------------------------------------
    print("=" * 80)
    print("PRE-FLIGHT HARD REGRESSION AUDIT AGAINST CERTIFIED INVARIANTS (N = 64, T = 600)")
    print("=" * 80)

    t_pf0 = time.time()
    Q_full_64, _ = get_galerkin_matrix(c=C_PARAM, N=64, T=600, dps=50, verbose=False)
    V_even_64 = canonical_even_projector(64)
    Q_even_64 = V_even_64.T * Q_full_64 * V_even_64
    Q_even_64 = mp.mpf("0.5") * (Q_even_64 + Q_even_64.T)
    evals_e_64, V_e_64 = symmetric_eigendecomposition(Q_even_64)

    sys_64 = {"V_Qeven": V_e_64}
    omega_0_64, nu_0_64, mu_0_64 = preflight_continuum_audit_64(sys_64, prime_data)
    t_pf = time.time() - t_pf0

    target_omega = mp.mpf("2.9315259531")
    target_nu = mp.mpf("4.2604954421")
    target_mu = mp.mpf("-0.4869792197")

    res_omega = abs(omega_0_64 - target_omega)
    res_nu = abs(nu_0_64 - target_nu)
    res_mu = abs(mu_0_64 - target_mu)

    print(f"  omega_0 = {float(omega_0_64):.10f} (Expected: 2.9315259531, Residual: {float(res_omega):.2e})")
    print(f"  nu_0    = {float(nu_0_64):.10f} (Expected: 4.2604954421, Residual: {float(res_nu):.2e})")
    print(f"  mu_0    = {float(mu_0_64):.10f} (Expected: -0.4869792197, Residual: {float(res_mu):.2e})")

    assert res_omega < mp.mpf("1e-8"), f"omega_0 regression failed: residual {res_omega}"
    assert res_nu < mp.mpf("1e-8"), f"nu_0 regression failed: residual {res_nu}"
    assert res_mu < mp.mpf("1e-8"), f"mu_0 regression failed: residual {res_mu}"
    print(f"  REGRESSION AUDIT PASSED in {t_pf:.2f}s: Operators match certified invariants.\n")

    # -------------------------------------------------------------------------
    # Dimension Sweeps Across Sweep Range
    # -------------------------------------------------------------------------
    print(f"Executing spectral deficit sweeps across N in {SWEEP_N} at T = {T_PARAM}...")
    results = []
    for N_val in SWEEP_N:
        res = solve_deficit_system_at_N(N_val, T_PARAM, prime_data)
        results.append(res)
        print(
            f"  Completed N = {res['N']:2d} (dim = {res['dim_even']:2d}, q = {res['q_cont']:2d}) "
            f"in {res['t_solve']:.2f}s | E_11 = {res['E_11']:.6f} | Delta W = {res['delta_W_exact']:.6f}"
        )
    print()

    # -------------------------------------------------------------------------
    # TABLE 1: WELL DEFICIT SPECTRUM & CLUSTER GAP FLOOR
    # -------------------------------------------------------------------------
    print("=" * 80)
    print("TABLE 1: WELL DEFICIT SPECTRUM delta_nu_k AND CLUSTER GAP FLOOR (N in [48..96])")
    print("Tracking nu_0, Cluster Eigenvalues nu_1..nu_3, Bulk Boundary nu_4, and Floor delta_*^(4)")
    print("=" * 80)
    print("  N |       nu_0 |       nu_1 |       nu_2 |       nu_3 |       nu_4 | delta_*^(4) |  Delta_step")
    print("-" * 80)
    for r in results:
        top = r["nu_top_list"]
        print(
            f"{r['N']:3d} | "
            f"{top[0]:10.6f} | "
            f"{top[1]:10.6f} | "
            f"{top[2]:10.6f} | "
            f"{top[3]:10.6f} | "
            f"{top[4]:10.6f} | "
            f"{r['delta_star_4']:11.6f} | "
            f"{r['delta_step']:10.6f}"
        )
    print("-" * 80)
    print("Observation: Verify whether the cluster gap floor delta_*^(4) = nu_0 - nu_4 remains")
    print("macroscopically bounded away from zero (delta_* > 0) as N -> oo.\n")

    # -------------------------------------------------------------------------
    # TABLE 2: EXACT MODAL DEFICIT CONSERVATION & BULK INEQUALITY
    # -------------------------------------------------------------------------
    print("=" * 80)
    print("TABLE 2: EXACT MODAL DEFICIT CONSERVATION & BULK PENALTY (N in [48..96])")
    print("Tracking Exact Delta W, Cluster Deficit, Bulk Deficit, Bulk Bound, and Deficit Residual")
    print("=" * 80)
    print("  N |    Delta W |  DeltaW_cls |  DeltaW_blk | delta_*4*P4 | Blk/Tot (%) | Deficit Res")
    print("-" * 80)
    for r in results:
        pct_bulk = (r["delta_W_bulk"] / r["delta_W_exact"]) * 100.0 if r["delta_W_exact"] > 0 else 0.0
        print(
            f"{r['N']:3d} | "
            f"{r['delta_W_exact']:10.6f} | "
            f"{r['delta_W_cluster']:11.6f} | "
            f"{r['delta_W_bulk']:11.6f} | "
            f"{r['bulk_bound_floor']:11.6f} | "
            f"{pct_bulk:10.2f}% | "
            f"{r['deficit_res']:11.2e}"
        )
    print("-" * 80)
    print("Observation: Confirm that bulk modes k >= 4 generate virtually 100% of the well sacrifice,")
    print("and verify the analytical bound Delta W >= delta_*^(4) * P_>=4.\n")

    # -------------------------------------------------------------------------
    # TABLE 3: MASTER OPERATOR TRI-PARTITION OF E_11(N)
    # -------------------------------------------------------------------------
    print("=" * 80)
    print("TABLE 3: MASTER OPERATOR TRI-PARTITION: E_11 = H_1 + Delta_arch + E_pole")
    print("Tracking Competition Energy H_1, Archimedean Non-Local Boost Delta_arch, and Pole Rescue")
    print("=" * 80)
    print("  N |       E_11 |   H_1[v_11] |  Delta_arch |      E_pole |  H_1 + D_ar | Tri-Res")
    print("-" * 80)
    for r in results:
        h1_plus_darch = r["H1_val"] + r["Delta_arch_val"]
        print(
            f"{r['N']:3d} | "
            f"{r['E_11']:10.6f} | "
            f"{r['H1_val']:11.6f} | "
            f"{r['Delta_arch_val']:11.6f} | "
            f"{r['E_pole']:11.6f} | "
            f"{h1_plus_darch:11.6f} | "
            f"{r['tri_closure_res']:9.2e}"
        )
    print("-" * 80)
    print("Observation: Verify exact algebraic closure of the Master Tri-Partition, and track")
    print("the non-local Archimedean kernel boost Delta_arch[v_11] > 0.\n")

    # -------------------------------------------------------------------------
    # TABLE 4: VARIATIONAL COERCIVITY CERTIFICATE & LOWER BOUND FLOOR
    # -------------------------------------------------------------------------
    print("=" * 80)
    print("TABLE 4: VARIATIONAL COERCIVITY CERTIFICATE & LOWER BOUND COMPARISON")
    print("Tracking Actual E_11 vs Analytical Lower Bound E_11^cert = mu_0 + delta_*4*P4 + D_ar + E_pol")
    print("=" * 80)
    print("  N |       E_11 |  mu_0(H_1) | delta_*4*P4 |  Delta_arch |      E_pole |   E_11^cert | Cert Margin")
    print("-" * 80)
    for r in results:
        cert_margin = r["E_11"] - r["E_11_cert"]
        print(
            f"{r['N']:3d} | "
            f"{r['E_11']:10.6f} | "
            f"{r['mu_0']:10.6f} | "
            f"{r['bulk_bound_floor']:11.6f} | "
            f"{r['Delta_arch_val']:11.6f} | "
            f"{r['E_pole']:11.6f} | "
            f"{r['E_11_cert']:11.6f} | "
            f"{cert_margin:11.6f}"
        )
    print("-" * 80)
    print("Observation: Confirm whether the certified lower bound E_11^cert strictly preserves")
    print("positivity across all tested dimensions.\n")

    # -------------------------------------------------------------------------
    # SYNTHESIS & KEY OBSERVATIONS
    # -------------------------------------------------------------------------
    print("=" * 80)
    print("SYNTHESIS & KEY OBSERVATIONS:")
    max_tri_res = max(r["tri_closure_res"] for r in results)
    max_def_res = max(r["deficit_res"] for r in results)
    print(f"  1. Master Closure: Identity Q_even = Q_comp + Delta_arch + Q_pole holds with residual < {max_tri_res:.2e}.")
    print(f"  2. Deficit Conservation: Modal deficit identity Delta W === Delta W_spec holds with residual < {max_def_res:.2e}.")
    print(f"  3. Cluster Gap Floor: delta_*^(4) ranges from {results[0]['delta_star_4']:.6f} to {results[-1]['delta_star_4']:.6f}.")
    print(f"  4. Bulk Deficit Dominance: Bulk modes k >= 4 account for {results[-1]['delta_W_bulk'] / results[-1]['delta_W_exact'] * 100:.2f}% of total well sacrifice at N=96.")
    print(f"  5. Coercivity Lower Bound: E_11^cert ranges from {results[0]['E_11_cert']:.6f} to {results[-1]['E_11_cert']:.6f}.")
    print("=" * 80)

    t_suite = time.time() - t_suite_start
    print(f"Total suite execution time: {t_suite:.2f}s")
    print("=" * 80)
    print("CELL 149 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell149_suite()
