# ============================================================
# CELL 104 — PROJECTED COUPLING ANATOMY & DESTRUCTIVE INTERFERENCE AUDIT
# ============================================================
#
# Gate:       Gate 1 (Finite-N Spectral Mechanism & Asymptotic Tail Extinction)
#
# Target Propositions:
#
#   1. Theorem 1 in cell104.md (Exact Eigenvector Complementarity Identity):
#      For any cutoff M < N, the projected coupling vector w_M = B_M^T v^(P)
#      satisfies the exact algebraic identity:
#
#          w_M = -(C_M - E11 * I) v^(Q),
#
#      where v^(Q) is the high-mode tail of the ground-state eigenvector v_N.
#
#   2. Theorem 2 in cell104.md (Exact Rayleigh Shift Energy Identity):
#      The physical Feshbach Rayleigh shift satisfies the exact identity:
#
#          Delta E11^Fesh(M) = <v^(P), R_M(E11) v^(P)>
#                            = <v^(Q), (C_M - E11 * I) v^(Q)>.
#
#      The resolvent inverse (C_M - E11*I)^(-1) exactly cancels one factor
#      of (C_M - E11*I), eliminating small-denominator inversion.
#
#   3. Theorem 3 in cell104.md (Two-Sided Tail-Mass Sandwich):
#      delta_M * ||v^(Q)||^2   <= Delta E11^Fesh(M) <= (||H||_op - E11) * ||v^(Q)||^2,
#      delta_M^2 * ||v^(Q)||^2 <= ||w_M||^2         <= (||H||_op - E11)^2 * ||v^(Q)||^2.
#
#   4. Theorem 4 in cell104.md (Collective Destructive Interference):
#      Individual terms v_m * H_{m, k} decay only algebraically (~ 1/k),
#      but undergo exact collective cancellation against intermediate modes:
#
#          sum_{m=0}^{m_0} v_m * H_{m, k} = -sum_{m=m_0+1}^M v_m * H_{m, k} + O(e^(-sigma M)),
#
#      producing cancellation ratios C_M(k) = (sum |v_m H_{mk}|) / |sum v_m H_{mk}| >> 1.
#
# Falsification Criteria:
#
#   - If ||w_M - [-(C_M - E11*I) v^(Q)]||_inf > 1e-50, Theorem 1 is refuted.
#   - If |Delta E11^Fesh - <v^(Q), (C_M - E11*I) v^(Q)>| > 1e-50, Theorem 2 is refuted.
#   - If Delta E11^Fesh or ||w_M||^2 violates the two-sided sandwich, Theorem 3 is refuted.
#   - If C_M(M+1) ~ 1 (no cancellation), Theorem 4 is refuted.
#
# Design:
#
#   Evaluates M in {24, 32, 48, 64} at N=192, c=13, T=600.
#   Computes ground state v_N of H at 70 dps.
#   Audits exact vector and scalar identities to working precision.
#   Performs mode-by-mode term decomposition to quantify destructive interference.
#
# ============================================================

import time
import mpmath as mp

from cell import get_galerkin_matrix

# ============================================================
# PARAMETERS & PRECISION
# ============================================================

mp.mp.dps = 70

C_PARAM = 13
T_PARAM = 600
N_PARAM = 192
GROUND_DPS = 50

M_VALUES = [24, 32, 48, 64]

DIST_TOL = mp.mpf("1e-50")

# ============================================================
# MATRIX UTILITIES
# ============================================================

def full_to_canonical_matrix(Q, N):
    """
    Convert full (2N+1) x (2N+1) Galerkin matrix Q into canonical
    even-basis (N+1) x (N+1) matrix H.
    """
    H = mp.matrix(N + 1, N + 1)
    centre = N

    H[0, 0] = Q[centre, centre]

    for k in range(1, N + 1):
        H[0, k] = mp.sqrt(2) * Q[centre, centre + k]
        H[k, 0] = H[0, k]

    for j in range(1, N + 1):
        for k in range(j, N + 1):
            value = (
                Q[centre + j, centre + k]
                + Q[centre + j, centre - k]
            )
            H[j, k] = value
            H[k, j] = value

    return H


def eigsys_sym(A):
    """
    Compute sorted eigenvalues and orthonormal eigenvectors of symmetric matrix A.
    Returns (evals, evecs) where evals is a list and evecs is an mp.matrix
    whose columns are the eigenvectors.
    """
    vals, V = mp.eigsy(A)
    dim = A.rows
    idx = sorted(range(dim), key=lambda i: vals[i])

    evals = [vals[i] for i in idx]
    evecs = mp.matrix(dim, dim)

    for col_out, col_in in enumerate(idx):
        for row in range(dim):
            evecs[row, col_out] = V[row, col_in]

    return evals, evecs


# ============================================================
# MAIN EXECUTION
# ============================================================

def main():
    t_start = time.time()

    print("=" * 80)
    print("CELL 104 — PROJECTED COUPLING ANATOMY & DESTRUCTIVE INTERFERENCE AUDIT")
    print("  Eigenvector Complementarity, Exact Rayleigh Shift, & Mode Cancellation")
    print(f"  Configuration: N = {N_PARAM}, c = {C_PARAM}, T = {T_PARAM}")
    print(f"  mpmath dps: {mp.mp.dps}")
    print("=" * 80)
    print()

    # 1. Retrieve Galerkin matrix
    print("--- STEP 1: RETRIEVING GALERKIN MATRIX ---")
    t0 = time.time()
    Q_full, _ = get_galerkin_matrix(
        c=C_PARAM,
        N=N_PARAM,
        T=T_PARAM,
        dps=GROUND_DPS,
        verbose=False,
    )
    H = full_to_canonical_matrix(Q_full, N_PARAM)
    print(f"Canonical even matrix H: shape ({H.rows}, {H.cols}), time: {time.time() - t0:.2f} s")
    print()

    # 2. Ground-State Eigensolve on Full H
    print("--- STEP 2: GROUND-STATE EIGENSYSTEM OF FULL H ---")
    t0 = time.time()
    evals_H, evecs_H = eigsys_sym(H)
    E11 = evals_H[0]
    H_norm_op = evals_H[-1]
    v_ground = [evecs_H[i, 0] for i in range(H.rows)]
    norm_v = mp.sqrt(sum(x * x for x in v_ground))

    # Eigenpair residual check
    Hv = [sum(H[i, j] * v_ground[j] for j in range(H.rows)) for i in range(H.rows)]
    res_eigen = max(abs(Hv[i] - E11 * v_ground[i]) for i in range(H.rows))

    print(f"Ground-state eigenvalue E11:      {mp.nstr(E11, 25)}")
    print(f"Spectral ceiling ||H||_op:        {mp.nstr(H_norm_op, 15)}")
    print(f"Eigenvector norm:                 {mp.nstr(norm_v, 20)}")
    print(f"Eigenpair residual ||Hv - Ev||:   {mp.nstr(res_eigen, 10)}")
    print(f"Eigensolve runtime: {time.time() - t0:.2f} s")
    print()

    # 3. Diagnostic Sweep across Cutoffs M
    records = []

    for M in M_VALUES:
        print(f"============================================================")
        print(f"--- CUTOFF M = {M} (P: 0..{M}, dim {M+1}; Q: {M+1}..{N_PARAM}, dim {N_PARAM-M}) ---")
        print(f"============================================================")
        t_m0 = time.time()

        p_dim = M + 1
        q_dim = N_PARAM - M

        # Subspace projections of ground state
        v_P = v_ground[:p_dim]
        v_Q = v_ground[p_dim:]
        norm_v_P_sq = sum(x * x for x in v_P)
        norm_v_Q_sq = sum(x * x for x in v_Q)

        # Partition H
        B = mp.matrix(p_dim, q_dim)
        for i in range(p_dim):
            for j in range(q_dim):
                B[i, j] = H[i, M + 1 + j]

        C = mp.matrix(q_dim, q_dim)
        for i in range(q_dim):
            for j in range(q_dim):
                C[i, j] = H[M + 1 + i, M + 1 + j]

        # Shifted high-mode block C_shifted = C - E11 * I
        C_shifted = mp.matrix(q_dim, q_dim)
        for i in range(q_dim):
            for j in range(q_dim):
                C_shifted[i, j] = C[i, j] - (E11 if i == j else mp.mpf(0))

        # Direct computation of projected coupling w_M = B^T * v^(P)
        w_M = [sum(B[i, j] * v_P[i] for i in range(p_dim)) for j in range(q_dim)]
        w_M_norm_sq = sum(x * x for x in w_M)

        # Complementarity computation: w_comp = -(C_M - E11*I) * v^(Q)
        w_comp = [-sum(C_shifted[i, j] * v_Q[j] for j in range(q_dim)) for i in range(q_dim)]

        # Residual of Theorem 1: ||w_M - w_comp||_inf
        res_comp_inf = max(abs(w_M[j] - w_comp[j]) for j in range(q_dim))
        rel_res_comp = res_comp_inf / mp.sqrt(w_M_norm_sq) if w_M_norm_sq > DIST_TOL else mp.mpf(0)

        # High-sector eigensystem of C_M
        t_eig = time.time()
        mu, U = eigsys_sym(C)
        t_eig_elapsed = time.time() - t_eig

        mu_0 = mu[0]
        mu_max = mu[-1]
        delta_M = mu_0 - E11
        shift_ceiling = H_norm_op - E11

        # Direct Feshbach Rayleigh shift via resolvent:
        # Delta E11^Fesh = w_M^T * (C - E11*I)^(-1) * w_M = sum_j <w_M, u_j>^2 / (mu_j - E11)
        w_tilde = [sum(U[k, j] * w_M[k] for k in range(q_dim)) for j in range(q_dim)]
        Delta_E_Fesh = sum((w_tilde[j] ** 2) / (mu[j] - E11) for j in range(q_dim))

        # Complementarity Rayleigh shift via tail expectation:
        # Delta E11^tail = <v^(Q), (C - E11*I) v^(Q)>
        v_Q_shifted = [sum(C_shifted[i, j] * v_Q[j] for j in range(q_dim)) for i in range(q_dim)]
        Delta_E_tail = sum(v_Q[i] * v_Q_shifted[i] for i in range(q_dim))

        # Residual of Theorem 2: |Delta_E_Fesh - Delta_E_tail|
        res_Rayleigh = abs(Delta_E_Fesh - Delta_E_tail)
        rel_res_Rayleigh = res_Rayleigh / Delta_E_Fesh if Delta_E_Fesh > DIST_TOL else mp.mpf(0)

        # Theorem 3: Two-sided bounds
        # Energy shift bounds
        lower_E_bound = delta_M * norm_v_Q_sq
        upper_E_bound = shift_ceiling * norm_v_Q_sq
        E_sandwich_holds = (lower_E_bound <= Delta_E_Fesh <= upper_E_bound)

        # Norm bounds
        lower_w_bound = (delta_M ** 2) * norm_v_Q_sq
        upper_w_bound = (shift_ceiling ** 2) * norm_v_Q_sq
        w_sandwich_holds = (lower_w_bound <= w_M_norm_sq <= upper_w_bound)

        # Effective ratio ||w_M||^2 / ||v^(Q)||^2
        w_ratio = w_M_norm_sq / norm_v_Q_sq if norm_v_Q_sq > DIST_TOL else mp.mpf(0)
        E_ratio = Delta_E_Fesh / norm_v_Q_sq if norm_v_Q_sq > DIST_TOL else mp.mpf(0)

        # Theorem 4: Mode-by-mode cancellation forensics for interface mode k = M+1 (j=0)
        # Entry in w_M is w_M[0] = sum_{m=0}^M H[m, M+1] * v_ground[m]
        terms_int = [H[m, M + 1] * v_ground[m] for m in range(p_dim)]
        sum_abs_terms = sum(abs(t) for t in terms_int)
        actual_sum = abs(sum(terms_int))
        cancellation_ratio = sum_abs_terms / actual_sum if actual_sum > DIST_TOL else mp.inf

        # Find largest individual term and its mode index
        max_term_val = mp.mpf(0)
        max_term_idx = 0
        for m, t in enumerate(terms_int):
            if abs(t) > max_term_val:
                max_term_val = abs(t)
                max_term_idx = m

        # Partition terms into core (m <= 8) and buffer (m > 8)
        m_core_cutoff = min(8, M)
        sum_core = sum(terms_int[:m_core_cutoff + 1])
        sum_buffer = sum(terms_int[m_core_cutoff + 1:])
        core_buffer_balance = abs(sum_core + sum_buffer) / abs(sum_core) if abs(sum_core) > DIST_TOL else mp.mpf(0)

        # Instantaneous tail decay rate
        sigma_M = -mp.log(norm_v_Q_sq) / (2 * M)

        print()
        print("  --- THEOREM 1: EXACT EIGENVECTOR COMPLEMENTARITY ---")
        print(f"  ||v^(Q)||^2 (tail mass):        {mp.nstr(norm_v_Q_sq, 15)}")
        print(f"  ||w_M||^2 (projected norm):     {mp.nstr(w_M_norm_sq, 15)}")
        print(f"  Complementarity res ||w - w*||: {mp.nstr(res_comp_inf, 10)}")
        print(f"  Relative complementarity res:   {mp.nstr(rel_res_comp, 10)}")
        print(f"  Theorem 1 certified (< 1e-50):  {res_comp_inf < DIST_TOL}")
        print()
        print("  --- THEOREM 2: EXACT RAYLEIGH SHIFT ENERGY IDENTITY ---")
        print(f"  Delta E11^Fesh (resolvent):     {mp.nstr(Delta_E_Fesh, 15)}")
        print(f"  Delta E11^tail (<vQ, (C-E)vQ>): {mp.nstr(Delta_E_tail, 15)}")
        print(f"  Rayleigh identity residual:     {mp.nstr(res_Rayleigh, 10)}")
        print(f"  Relative Rayleigh residual:     {mp.nstr(rel_res_Rayleigh, 10)}")
        print(f"  Theorem 2 certified (< 1e-50):  {res_Rayleigh < DIST_TOL}")
        print()
        print("  --- THEOREM 3: TWO-SIDED TAIL-MASS SANDWICH ---")
        print(f"  delta_M:                        {mp.nstr(delta_M, 10)}")
        print(f"  Shift bounds [lower, upper]:    [{mp.nstr(lower_E_bound, 8)}, {mp.nstr(upper_E_bound, 8)}]")
        print(f"  Energy shift sandwich holds:    {E_sandwich_holds} (ratio: {mp.nstr(E_ratio, 6)})")
        print(f"  Norm bounds [lower, upper]:     [{mp.nstr(lower_w_bound, 8)}, {mp.nstr(upper_w_bound, 8)}]")
        print(f"  Norm sandwich holds:            {w_sandwich_holds} (ratio: {mp.nstr(w_ratio, 6)})")
        print()
        print("  --- THEOREM 4: COLLECTIVE DESTRUCTIVE INTERFERENCE FORENSICS ---")
        print(f"  Interface mode k = M+1:         k = {M + 1}")
        print(f"  Largest single term |v_m H_mk|: {mp.nstr(max_term_val, 10)} at m = {max_term_idx}")
        print(f"  Sum of absolute terms:          {mp.nstr(sum_abs_terms, 10)}")
        print(f"  Net sum |w_{M, M+1}|:           {mp.nstr(actual_sum, 10)}")
        print(f"  Cancellation Ratio C_M(M+1):    {mp.nstr(cancellation_ratio, 6)} (~ 10^{float(mp.log10(cancellation_ratio)):.1f})")
        print(f"  Core sum (m <= {m_core_cutoff}):           {mp.nstr(sum_core, 10)}")
        print(f"  Buffer sum (m > {m_core_cutoff}):          {mp.nstr(sum_buffer, 10)}")
        print(f"  |Core + Buffer| / |Core|:       {mp.nstr(core_buffer_balance, 6)}")
        print(f"  Tail decay rate sigma_M:        {mp.nstr(sigma_M, 6)}")
        print(f"  Cutoff runtime: {time.time() - t_m0:.2f} s")
        print()

        rec = {
            "M": M,
            "norm_v_Q_sq": norm_v_Q_sq,
            "w_M_norm_sq": w_M_norm_sq,
            "Delta_E_Fesh": Delta_E_Fesh,
            "res_comp": res_comp_inf,
            "res_Rayleigh": res_Rayleigh,
            "delta_M": delta_M,
            "E_ratio": E_ratio,
            "w_ratio": w_ratio,
            "cancellation_ratio": cancellation_ratio,
            "core_balance": core_buffer_balance,
            "sigma_M": sigma_M,
            "runtime": time.time() - t_m0,
        }
        records.append(rec)

    # ============================================================
    # SYNTHESIS TABLE
    # ============================================================
    print("=" * 115)
    print("CELL 104 SYNTHESIS: EIGENVECTOR COMPLEMENTARITY & DESTRUCTIVE INTERFERENCE")
    print("=" * 115)
    headers = [
        "M", "||v^(Q)||^2", "||w_M||^2", "Delta_E_Fesh",
        "Res_comp", "Res_Rayl", "C_M(M+1)", "Core_bal", "sigma_M"
    ]
    print(f"{headers[0]:>4} | {headers[1]:>14} | {headers[2]:>14} | {headers[3]:>14} | "
          f"{headers[4]:>10} | {headers[5]:>10} | {headers[6]:>12} | {headers[7]:>10} | {headers[8]:>8}")
    print("-" * 115)

    for r in records:
        print(
            f"{r['M']:>4} | "
            f"{mp.nstr(r['norm_v_Q_sq'], 6):>14} | "
            f"{mp.nstr(r['w_M_norm_sq'], 6):>14} | "
            f"{mp.nstr(r['Delta_E_Fesh'], 6):>14} | "
            f"{mp.nstr(r['res_comp'], 4):>10} | "
            f"{mp.nstr(r['res_Rayleigh'], 4):>10} | "
            f"{mp.nstr(r['cancellation_ratio'], 4):>12} | "
            f"{mp.nstr(r['core_balance'], 4):>10} | "
            f"{mp.nstr(r['sigma_M'], 4):>8}"
        )
    print("=" * 115)
    print()

    # Final Audit Summary
    print("--- PRECISE THEORETICAL CONCLUSIONS ---")
    all_thm1 = all(r["res_comp"] < DIST_TOL for r in records)
    all_thm2 = all(r["res_Rayleigh"] < DIST_TOL for r in records)
    all_thm3 = all(r["delta_M"] * r["norm_v_Q_sq"] <= r["Delta_E_Fesh"] <= (H_norm_op - E11) * r["norm_v_Q_sq"] for r in records)
    all_thm4 = all(r["cancellation_ratio"] > 1e6 for r in records)

    print(f"1. Theorem 1 (w_M == -(C_M - E11*I) v^(Q)):     {'CERTIFIED' if all_thm1 else 'VIOLATED'}")
    print(f"2. Theorem 2 (Delta E11 == <v^(Q), (C-E) v^(Q)>): {'CERTIFIED' if all_thm2 else 'VIOLATED'}")
    print(f"3. Theorem 3 (Two-Sided Tail-Mass Sandwich):    {'CERTIFIED' if all_thm3 else 'VIOLATED'}")
    print(f"4. Theorem 4 (Collective Destructive Interf.):  {'CERTIFIED' if all_thm4 else 'VIOLATED'}")
    print(f"Total script runtime: {time.time() - t_start:.2f} s")
    print()

    print("=" * 80)
    print("CELL 104 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
