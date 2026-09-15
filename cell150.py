#!/usr/bin/env python3
"""
CELL 150 — Spectrum, Definiteness, and Operator Dominance of the Combined
Non-Local Positive Form A_pos = Delta_arch + Q_pole on B_11^perp.

Target Gate: Gate 1 (Finite-N Spectral Mechanism & Asymptotic Tail Extinction)
Milestone: Milestone M-G1.6 / Operator Dominance & Non-Local Positivity on B_11^perp
Investigating:
  1. Complete Spectrum of Constituent Positive Operators:
     Evaluating eigenvalues of Q_pole,perp, Delta_arch,perp, and A_pos,perp
     on the continuum subspace B_11^perp across N in [48, 56, 64, 72, 80, 88, 96].
  2. Definiteness & Coercivity Lower Bound:
     Testing whether alpha_0(N) = lambda_min(A_pos,perp) >= alpha_* > 0.
  3. Overlap Anatomy of the Threshold State v_11:
     Projecting v_11 onto the eigenmodes z_j of A_pos,perp and verifying
     spectral expectation reconstruction.
  4. Operator Dominance Quotient & Surplus Margin:
     Evaluating R_dom[v_11] = A_pos[v_11] / (-H_1[v_11]) and surplus E_11.
  5. Non-Local Archimedean Kernel Structure (Delta_arch):
     Decomposing Delta_arch into diagonal shifts vs off-diagonal Frobenius coupling.

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
# System Solver & Non-Local Positivity Dissection at Dimension N
# =============================================================================
def solve_positivity_system_at_N(N: int, T_val: int, prime_data: list) -> dict:
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

    # Prime, pole, and Archimedean components
    Q_prime_even = build_Q_prime_even(N, prime_data, L_PARAM)
    Q_pole_even = build_Q_pole_even(N, L_PARAM)
    Q_arch_even = Q_even - Q_prime_even - Q_pole_even
    Q_arch_even = mp.mpf("0.5") * (Q_arch_even + Q_arch_even.T)

    # Diagonal Archimedean multiplier Omega_arch and competition operator Q_comp
    Omega_arch = mp.matrix(dim_even, dim_even)
    for m in range(dim_even):
        a_m = mp.mpf("2") * PI * mp.mpf(m) / L_PARAM if m > 0 else mp.mpf("0")
        Omega_arch[m, m] = h_plus(a_m)

    Q_comp_even = Omega_arch + Q_prime_even
    Q_comp_even = mp.mpf("0.5") * (Q_comp_even + Q_comp_even.T)

    # Non-local Archimedean remainder Delta_arch = Q_arch - Omega_arch
    Delta_arch_even = Q_arch_even - Omega_arch
    Delta_arch_even = mp.mpf("0.5") * (Delta_arch_even + Delta_arch_even.T)

    # Combined non-local positive operator A_pos = Delta_arch + Q_pole
    A_pos_even = Delta_arch_even + Q_pole_even
    A_pos_even = mp.mpf("0.5") * (A_pos_even + A_pos_even.T)

    # Threshold state v_11
    v_11_vec = [V_e[r, N_BOUND] for r in range(dim_even)]
    v_11_mat = mp.matrix(dim_even, 1)
    for r in range(dim_even):
        v_11_mat[r, 0] = v_11_vec[r]

    E_11 = evals_e[N_BOUND]
    E_arch = float((v_11_mat.T * Q_arch_even * v_11_mat)[0, 0])
    E_prime = float((v_11_mat.T * Q_prime_even * v_11_mat)[0, 0])
    E_pole = float((v_11_mat.T * Q_pole_even * v_11_mat)[0, 0])
    H1_val = float((v_11_mat.T * Q_comp_even * v_11_mat)[0, 0])
    Delta_arch_val = float((v_11_mat.T * Delta_arch_even * v_11_mat)[0, 0])
    A_pos_val = float((v_11_mat.T * A_pos_even * v_11_mat)[0, 0])
    tri_closure_res = abs(float(E_11) - (H1_val + Delta_arch_val + E_pole))
    A_pos_sum_res = abs(A_pos_val - (Delta_arch_val + E_pole))

    # Continuum subspace restriction onto B_11^perp
    U_cont = extract_continuum_projector(V_e, n_bound=N_BOUND)

    # Restricted constituent operators
    Q_pole_perp = U_cont.T * Q_pole_even * U_cont
    Q_pole_perp = mp.mpf("0.5") * (Q_pole_perp + Q_pole_perp.T)

    Delta_arch_perp = U_cont.T * Delta_arch_even * U_cont
    Delta_arch_perp = mp.mpf("0.5") * (Delta_arch_perp + Delta_arch_perp.T)

    A_pos_perp = U_cont.T * A_pos_even * U_cont
    A_pos_perp = mp.mpf("0.5") * (A_pos_perp + A_pos_perp.T)

    # Eigensystem of Q_pole,perp
    evals_pole, _ = symmetric_eigendecomposition(Q_pole_perp, sort_descending=False)
    lambda_min_pole = float(evals_pole[0])
    lambda_max_pole = float(evals_pole[-1])

    # Eigensystem of Delta_arch,perp
    evals_darch, _ = symmetric_eigendecomposition(Delta_arch_perp, sort_descending=False)
    lambda_min_darch = float(evals_darch[0])
    lambda_max_darch = float(evals_darch[-1])

    # Eigensystem of A_pos,perp (sorted ascending: alpha_0 <= alpha_1 <= ...)
    evals_Apos, V_Apos = symmetric_eigendecomposition(A_pos_perp, sort_descending=False)
    alpha_0 = float(evals_Apos[0])
    alpha_1 = float(evals_Apos[1]) if q_cont > 1 else alpha_0
    alpha_2 = float(evals_Apos[2]) if q_cont > 2 else alpha_1
    alpha_3 = float(evals_Apos[3]) if q_cont > 3 else alpha_2
    alpha_4 = float(evals_Apos[4]) if q_cont > 4 else alpha_3
    alpha_max = float(evals_Apos[-1])
    gap_alpha_10 = alpha_1 - alpha_0
    cond_Apos = (alpha_max / alpha_0) if alpha_0 > 0 else float("inf")

    # Module 2: Overlap Anatomy of Threshold State v_11 in A_pos,perp basis
    # In continuum coordinates U_cont, v_11 is basis vector e_0 = (1, 0, ..., 0)^T
    # Overlap with eigenvector z_j of A_pos,perp is V_Apos[0, j]
    P_A = [float(V_Apos[0, j] ** 2) for j in range(q_cont)]
    P_A_0 = P_A[0]
    P_A_le4 = sum(P_A[:min(5, q_cont)])
    P_A_ge5 = sum(P_A[5:]) if q_cont > 5 else 0.0

    # Spectral expectation reconstruction
    A_pos_spec = sum(float(evals_Apos[j]) * P_A[j] for j in range(q_cont))
    spec_recon_res = abs(A_pos_val - A_pos_spec)

    # Module 3: Operator Dominance Quotient & Surplus Margin
    # Competition energy on v_11 is H1_val < 0; -H_1[v_11] > 0
    neg_H1 = -H1_val
    R_dom = (A_pos_val / neg_H1) if neg_H1 > 0 else float("nan")
    surplus_margin = float(E_11)  # by identity A_pos - (-H_1) = A_pos + H_1 = E_11
    rel_margin = (surplus_margin / A_pos_val) if A_pos_val > 0 else 0.0

    # Module 4: Non-Local Archimedean Kernel Structure (Delta_arch)
    # Decompose Delta_arch into diagonal shift and off-diagonal coupling
    Delta_diag_mat = mp.matrix(dim_even, dim_even)
    Delta_off_mat = mp.matrix(dim_even, dim_even)
    diag_frob_sq = mp.mpf("0")
    off_frob_sq = mp.mpf("0")
    for r in range(dim_even):
        for c_idx in range(dim_even):
            val = Delta_arch_even[r, c_idx]
            if r == c_idx:
                Delta_diag_mat[r, c_idx] = val
                diag_frob_sq += val ** 2
            else:
                Delta_off_mat[r, c_idx] = val
                off_frob_sq += val ** 2

    frob_diag = float(mp.sqrt(diag_frob_sq))
    frob_off = float(mp.sqrt(off_frob_sq))
    frob_tot = float(mp.sqrt(diag_frob_sq + off_frob_sq))
    rho_off = (frob_off / frob_tot * 100.0) if frob_tot > 0 else 0.0

    Delta_diag_val = float((v_11_mat.T * Delta_diag_mat * v_11_mat)[0, 0])
    Delta_off_val = float((v_11_mat.T * Delta_off_mat * v_11_mat)[0, 0])
    darch_decomp_res = abs(Delta_arch_val - (Delta_diag_val + Delta_off_val))

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
        "neg_H1": neg_H1,
        "Delta_arch_val": Delta_arch_val,
        "A_pos_val": A_pos_val,
        "tri_closure_res": tri_closure_res,
        "A_pos_sum_res": A_pos_sum_res,
        # Module 1
        "lambda_min_pole": lambda_min_pole,
        "lambda_max_pole": lambda_max_pole,
        "lambda_min_darch": lambda_min_darch,
        "lambda_max_darch": lambda_max_darch,
        "alpha_0": alpha_0,
        "alpha_1": alpha_1,
        "alpha_2": alpha_2,
        "alpha_3": alpha_3,
        "alpha_4": alpha_4,
        "alpha_max": alpha_max,
        "gap_alpha_10": gap_alpha_10,
        "cond_Apos": cond_Apos,
        # Module 2
        "P_A_0": P_A_0,
        "P_A_le4": P_A_le4,
        "P_A_ge5": P_A_ge5,
        "A_pos_spec": A_pos_spec,
        "spec_recon_res": spec_recon_res,
        # Module 3
        "R_dom": R_dom,
        "surplus_margin": surplus_margin,
        "rel_margin": rel_margin,
        # Module 4
        "Delta_diag_val": Delta_diag_val,
        "Delta_off_val": Delta_off_val,
        "frob_diag": frob_diag,
        "frob_off": frob_off,
        "frob_tot": frob_tot,
        "rho_off": rho_off,
        "darch_decomp_res": darch_decomp_res,
        "t_solve": t_solve,
    }


# =============================================================================
# Main Execution Suite
# =============================================================================
def run_cell150_suite():
    t_suite_start = time.time()

    print("=" * 80)
    print("CELL 150 — NON-LOCAL POSITIVITY & OPERATOR DOMINANCE ON B_11^perp")
    print(f"Parameters: c = {C_PARAM}, L = {float(L_PARAM):.11f}, T = {T_PARAM}, dps = 50")
    print("Target: Combined Positive Operator A_pos = Delta_arch + Q_pole on B_11^perp")
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
    print(f"Executing non-local positivity sweeps across N in {SWEEP_N} at T = {T_PARAM}...")
    results = []
    for N_val in SWEEP_N:
        res = solve_positivity_system_at_N(N_val, T_PARAM, prime_data)
        results.append(res)
        print(
            f"  Completed N = {res['N']:2d} (dim = {res['dim_even']:2d}, q = {res['q_cont']:2d}) "
            f"in {res['t_solve']:.2f}s | E_11 = {res['E_11']:.6f} | "
            f"alpha_0 = {res['alpha_0']:.6f} | R_dom = {res['R_dom']:.4f}"
        )
    print()

    # -------------------------------------------------------------------------
    # TABLE 1: SPECTRUM & DEFINITENESS OF CONSTITUENT POSITIVE OPERATORS
    # -------------------------------------------------------------------------
    print("=" * 80)
    print("TABLE 1: CONSTITUENT POSITIVE OPERATOR SPECTRUM ON B_11^perp (N in [48..96])")
    print("Tracking lambda_min(Q_pole), lambda_min(Delta_arch), alpha_0(A_pos), and Cond(A_pos)")
    print("=" * 80)
    print("  N | dim | q_cont | min(Q_pole) | min(D_arch) | alpha_0(A_pos) | alpha_max(A_pos) | Cond(A_pos)")
    print("-" * 80)
    for r in results:
        print(
            f"{r['N']:3d} | "
            f"{r['dim_even']:3d} | "
            f"{r['q_cont']:6d} | "
            f"{r['lambda_min_pole']:11.6f} | "
            f"{r['lambda_min_darch']:11.6f} | "
            f"{r['alpha_0']:14.6f} | "
            f"{r['alpha_max']:16.6f} | "
            f"{r['cond_Apos']:11.2f}"
        )
    print("-" * 80)
    print("Observation: Check whether Q_pole >= 0, Delta_arch > 0, and alpha_0(A_pos) >= alpha_* > 0")
    print("strictly across all tested discrete Galerkin dimensions.\n")

    # -------------------------------------------------------------------------
    # TABLE 2: LOWEST EIGENVALUES OF A_pos,perp (alpha_0 <= ... <= alpha_4)
    # -------------------------------------------------------------------------
    print("=" * 80)
    print("TABLE 2: LOWEST EIGENVALUES OF A_pos,perp ON B_11^perp (alpha_0 <= ... <= alpha_4)")
    print("Tracking the lowest 5 eigenvalues and fundamental spectral gap Delta alpha_10")
    print("=" * 80)
    print("  N |    alpha_0 |    alpha_1 |    alpha_2 |    alpha_3 |    alpha_4 | Gap (a_1-a_0)")
    print("-" * 80)
    for r in results:
        print(
            f"{r['N']:3d} | "
            f"{r['alpha_0']:10.6f} | "
            f"{r['alpha_1']:10.6f} | "
            f"{r['alpha_2']:10.6f} | "
            f"{r['alpha_3']:10.6f} | "
            f"{r['alpha_4']:10.6f} | "
            f"{r['gap_alpha_10']:13.6f}"
        )
    print("-" * 80)
    print("Observation: Verify whether the bottom eigenvalue alpha_0 converges to a strictly")
    print("positive continuum limit, or exhibits clustering/gap closure.\n")

    # -------------------------------------------------------------------------
    # TABLE 3: THRESHOLD STATE EXPECTATION & OPERATOR DOMINANCE RATIO
    # -------------------------------------------------------------------------
    print("=" * 80)
    print("TABLE 3: THRESHOLD STATE EXPECTATION & OPERATOR DOMINANCE RATIO R_dom[v_11]")
    print("Tracking E_11, A_pos[v_11], -H_1[v_11], R_dom[v_11], P_A(0), and P_A(<=4)")
    print("=" * 80)
    print("  N |       E_11 | A_pos[v_11] |  -H_1[v_11] | R_dom[v_11] | Surplus (%) |  P_A(0) (%) | P_A(<=4) (%) | Spec Res")
    print("-" * 80)
    for r in results:
        pct_surplus = r["rel_margin"] * 100.0
        pct_p0 = r["P_A_0"] * 100.0
        pct_ple4 = r["P_A_le4"] * 100.0
        print(
            f"{r['N']:3d} | "
            f"{r['E_11']:10.6f} | "
            f"{r['A_pos_val']:11.6f} | "
            f"{r['neg_H1']:11.6f} | "
            f"{r['R_dom']:11.6f} | "
            f"{pct_surplus:10.2f}% | "
            f"{pct_p0:10.2f}% | "
            f"{pct_ple4:11.2f}% | "
            f"{r['spec_recon_res']:8.2e}"
        )
    print("-" * 80)
    print("Observation: Confirm whether R_dom[v_11] > 1 strictly, ensuring A_pos dominates -H_1,")
    print("and evaluate whether v_11 is concentrated in the lowest modes of A_pos,perp.\n")

    # -------------------------------------------------------------------------
    # TABLE 4: NON-LOCAL ARCHIMEDEAN KERNEL STRUCTURE (Delta_arch)
    # -------------------------------------------------------------------------
    print("=" * 80)
    print("TABLE 4: NON-LOCAL ARCHIMEDEAN KERNEL ANATOMY: DIAGONAL SHIFT VS OFF-DIAGONAL")
    print("Tracking Delta_arch[v_11], Diagonal Component, Off-Diagonal Coupling, and Frob Ratios")
    print("=" * 80)
    print("  N | Delta_arch |  Delta_diag |   Delta_off | ||D_diag||_F |  ||D_off||_F | Off/Tot (%) | Decomp Res")
    print("-" * 80)
    for r in results:
        print(
            f"{r['N']:3d} | "
            f"{r['Delta_arch_val']:10.6f} | "
            f"{r['Delta_diag_val']:11.6f} | "
            f"{r['Delta_off_val']:11.6f} | "
            f"{r['frob_diag']:12.6f} | "
            f"{r['frob_off']:12.6f} | "
            f"{r['rho_off']:10.2f}% | "
            f"{r['darch_decomp_res']:10.2e}"
        )
    print("-" * 80)
    print("Observation: Determine whether the positive Archimedean remainder Delta_arch[v_11] ~ +0.196")
    print("originates primarily from diagonal multiplier correction or off-diagonal coherence.\n")

    # -------------------------------------------------------------------------
    # SYNTHESIS & KEY OBSERVATIONS
    # -------------------------------------------------------------------------
    print("=" * 80)
    print("SYNTHESIS & KEY OBSERVATIONS:")
    max_tri_res = max(r["tri_closure_res"] for r in results)
    max_pos_res = max(r["A_pos_sum_res"] for r in results)
    max_spec_res = max(r["spec_recon_res"] for r in results)
    max_decomp_res = max(r["darch_decomp_res"] for r in results)
    print(f"  1. Algebraic Closure: Tri-partition residual < {max_tri_res:.2e}; A_pos decomposition residual < {max_pos_res:.2e}.")
    print(f"  2. Spectral Reconstruction: Eigenmode projection reconstructs A_pos[v_11] to < {max_spec_res:.2e}.")
    print(f"  3. Archimedean Kernel Decomposition: Frobenius split verified to residual < {max_decomp_res:.2e}.")
    print(f"  4. Definiteness Audits: min(Q_pole) = {results[-1]['lambda_min_pole']:.6f}, min(Delta_arch) = {results[-1]['lambda_min_darch']:.6f}, alpha_0 = {results[-1]['alpha_0']:.6f} at N=96.")
    print(f"  5. Operator Dominance: R_dom[v_11] ranges from {results[0]['R_dom']:.6f} to {results[-1]['R_dom']:.6f} (margin E_11 > 0).")
    print("=" * 80)

    t_suite = time.time() - t_suite_start
    print(f"Total suite execution time: {t_suite:.2f}s")
    print("=" * 80)
    print("CELL 150 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell150_suite()
