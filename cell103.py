# ============================================================
# CELL 103 — ASYMPTOTIC SCALING AUDIT: GEOMETRIC PREFACTOR & ISOTROPIC BASELINE
# ============================================================
#
# Gate:       Gate 1 (Finite-N Spectral Mechanism & Asymptotic Tail Extinction)
#
# Target Propositions:
#
#   1. Theorem 1 in cell103.md (Universal Spectral Bandwidth Bound):
#      For any cutoff M < N, the geometric prefactor C_geom(M) satisfies
#      the exact unconditional upper bound:
#
#          C_geom(M) <= (mu_{q-1} - mu_0) / delta_M = diam(sigma(C_M)) / delta_M.
#
#      Because ||H||_op <= M(c, T) and delta_M >= delta_infty > 0,
#      C_geom(M) is uniformly bounded for all M, rigorously establishing
#      that C_geom(M) * D_KS(M) -> 0 is equivalent to D_KS(M) -> 0.
#
#   2. Theorem 2 in cell103.md (Isotropic Baseline Factorization):
#      The isotropic baseline factors as:
#
#          S_M^iso(E11) = ||B_M||_F^2 * G_bar_M(E11),
#
#      where G_bar_M(E11) = (1/q_M) tr((C_M - E11*I)^(-1)) in [1/mu_max, 1/delta_M]
#      is strictly Theta(1), so S_M^iso -> 0 iff ||B_M||_F^2 -> 0.
#
#   3. Theorem 3 in cell103.md (Master Operator Norm Enclosure):
#      The Feshbach trace and relative excess satisfy:
#
#          S_M(E11) <= ||B_M||_F^2 * [ G_bar_M(E11) + D_KS(M) / delta_M ],
#          |epsilon_M| <= C_geom(M) * D_KS(M) <= [diam(sigma(C_M)) / delta_M] * D_KS(M).
#
#   4. Ground-State Low-Rank Projection (Mechanism 2):
#      The physical Feshbach Rayleigh shift of the ground state satisfies:
#
#          Delta E11^Fesh(M) = <v^(P), R_M(E11) v^(P)> <= E_proj(M) / delta_M,
#
#      where E_proj(M) = ||B_M^T v^(P)||^2 plummets due to exponential
#      spatial localization of the ground state v_N.
#
# Falsification Criteria:
#
#   - If C_geom(M) > diam(sigma(C_M)) / delta_M for any M, Theorem 1 is refuted.
#   - If G_bar_M(E11) falls outside [1/(mu_max - E11), 1/delta_M], Proposition 2.1 is refuted.
#   - If S_M(E11) exceeds the Master Bound ||B||_F^2 * [G_bar + D_KS/delta_M], Theorem 3 is refuted.
#   - If Delta E11^Fesh(M) exceeds E_proj(M) / delta_M, projection bounding is refuted.
#
# Design:
#
#   Evaluates M in {24, 32, 48, 64} at N=192, c=13, T=600.
#   Vectorized projection V = B * U for rapid evaluation of coupling norms.
#   Ground-state eigenvector v_N computed directly via eigsys_sym(H).
#   Audits both full P-space trace S_M and projected ground-state shift Delta E11^Fesh.
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
    print("CELL 103 — ASYMPTOTIC SCALING AUDIT")
    print("  Geometric Prefactor C_geom, Isotropic Baseline S_M^iso, & Master Enclosure")
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
    v_ground = [evecs_H[i, 0] for i in range(H.rows)]
    norm_v = mp.sqrt(sum(x * x for x in v_ground))

    print(f"Ground-state eigenvalue E11: {mp.nstr(E11, 20)}")
    print(f"First excited eigenvalue E1: {mp.nstr(evals_H[1], 20)}")
    print(f"Spectral gap E1 - E11:       {mp.nstr(evals_H[1] - E11, 20)}")
    print(f"Ground-state eigenvector norm: {mp.nstr(norm_v, 20)}")
    print(f"Eigensolve runtime: {time.time() - t0:.2f} s")
    print()

    # 3. Diagnostic Sweep across M
    records = []

    for M in M_VALUES:
        print(f"============================================================")
        print(f"--- CUTOFF M = {M} (P-space: 0..{M}, dim {M+1}; Q-space: {M+1}..{N_PARAM}, dim {N_PARAM-M}) ---")
        print(f"============================================================")
        t_m0 = time.time()

        p_dim = M + 1
        q_dim = N_PARAM - M

        # Partition H
        B = mp.matrix(p_dim, q_dim)
        for i in range(p_dim):
            for j in range(q_dim):
                B[i, j] = H[i, M + 1 + j]

        C = mp.matrix(q_dim, q_dim)
        for i in range(q_dim):
            for j in range(q_dim):
                C[i, j] = H[M + 1 + i, M + 1 + j]

        # Frobenius norm of B
        frob_B_sq = sum(B[i, j] ** 2 for i in range(p_dim) for j in range(q_dim))
        frob_B = mp.sqrt(frob_B_sq)

        # High-sector eigensystem C_M
        t_eig = time.time()
        mu, U = eigsys_sym(C)
        t_eig_elapsed = time.time() - t_eig

        mu_0 = mu[0]
        mu_max = mu[-1]
        delta_M = mu_0 - E11
        spec_diam = mu_max - mu_0

        print(f"  ||B_M||_F^2:                 {mp.nstr(frob_B_sq, 12)}")
        print(f"  mu_0 (lowest high-mode):     {mp.nstr(mu_0, 12)}")
        print(f"  mu_max (highest high-mode):  {mp.nstr(mu_max, 12)}")
        print(f"  delta_M = mu_0 - E11:        {mp.nstr(delta_M, 12)}")
        print(f"  diam(sigma(C_M)):            {mp.nstr(spec_diam, 12)}")
        print(f"  Eigensolve Q-sector ({q_dim}x{q_dim}): {t_eig_elapsed:.2f} s")

        # Resolvent summands g_j = 1 / (mu_j - E11)
        g = [1 / (mu[j] - E11) for j in range(q_dim)]
        g_0 = g[0]
        g_max = g[-1]
        g_mean = sum(g) / q_dim  # G_bar_M(E11)

        # Vectorized projection V = B * U
        t_proj = time.time()
        V = B * U
        a = [sum(V[i, j] ** 2 for i in range(p_dim)) for j in range(q_dim)]
        t_proj_elapsed = time.time() - t_proj

        # Parseval check
        sum_a = sum(a)
        parseval_err = abs(sum_a - frob_B_sq)
        a_bar = frob_B_sq / q_dim
        r = [a[j] / a_bar for j in range(q_dim)]

        # Actual trace S_M and Isotropic trace S_M^iso
        S_M = sum(a[j] * g[j] for j in range(q_dim))
        S_M_iso = a_bar * sum(g)  # = frob_B_sq * g_mean
        Delta_S = abs(S_M - S_M_iso)
        epsilon_M = (S_M / S_M_iso) - mp.mpf(1)

        # Geometric prefactor C_geom
        C_geom_exact = (g_0 - g_max) / g_mean
        C_geom_bound = spec_diam / delta_M
        geom_slack = C_geom_bound / C_geom_exact

        # Kolmogorov-Smirnov distance
        F_coup = mp.mpf(0)
        D_KS = mp.mpf(0)
        for j in range(q_dim):
            F_coup += a[j] / frob_B_sq
            F_unif = mp.mpf(j + 1) / mp.mpf(q_dim)
            diff = abs(F_coup - F_unif)
            if diff > D_KS:
                D_KS = diff

        # Bound on S_M - S_M^iso via KS
        B_KS = frob_B_sq * D_KS * (g_0 - g_max)
        ks_slack = B_KS / Delta_S if Delta_S > DIST_TOL else mp.inf

        # Bound on relative excess epsilon_M
        eps_bound_exact = C_geom_exact * D_KS
        eps_bound_univ = C_geom_bound * D_KS

        # Master Operator Norm Bound (Theorem 3)
        B_master = frob_B_sq * (g_mean + D_KS / delta_M)
        master_slack = B_master / S_M

        # Normalized resolvent trace bounds (Proposition 2.1)
        G_bar_lower = 1 / (mu_max - E11)
        G_bar_upper = 1 / delta_M
        g_mean_valid = (G_bar_lower <= g_mean <= G_bar_upper)

        # Mechanism 2: Ground-State Projection
        v_P = v_ground[:p_dim]
        v_Q = v_ground[p_dim:]
        norm_v_P_sq = sum(x * x for x in v_P)
        norm_v_Q_sq = sum(x * x for x in v_Q)

        # w^(Q) = B^T * v^(P) in R^{q_dim}
        w_Q = [sum(B[i, j] * v_P[i] for i in range(p_dim)) for j in range(q_dim)]
        E_proj = sum(x * x for x in w_Q)

        # Rayleigh shift: Delta E11^Fesh = <v^(P), R_M(E11) v^(P)>
        # In terms of U-coordinates: w_tilde = U^T * w_Q
        w_tilde = [sum(U[k, j] * w_Q[k] for k in range(q_dim)) for j in range(q_dim)]
        Delta_E_Fesh = sum((w_tilde[j] ** 2) * g[j] for j in range(q_dim))
        B_proj = E_proj / delta_M
        proj_slack = B_proj / Delta_E_Fesh if Delta_E_Fesh > DIST_TOL else mp.inf

        # Per-mode average trace in P-space
        S_bar_P = S_M / p_dim

        print()
        print("  --- THEOREM 1: UNIVERSAL BANDWIDTH BOUND ON C_geom ---")
        print(f"  C_geom (exact):              {mp.nstr(C_geom_exact, 10)}")
        print(f"  C_geom bound (diam/delta):   {mp.nstr(C_geom_bound, 10)}")
        print(f"  Theorem 1 verified:          {C_geom_exact <= C_geom_bound} (slack: {mp.nstr(geom_slack, 4)}x)")
        print()
        print("  --- THEOREM 2: ISOTROPIC BASELINE STRUCTURE ---")
        print(f"  G_bar_M(E11) = <g>_unif:     {mp.nstr(g_mean, 10)}")
        print(f"  G_bar bounds:                [{mp.nstr(G_bar_lower, 8)}, {mp.nstr(G_bar_upper, 8)}]")
        print(f"  G_bar bounds satisfied:      {g_mean_valid}")
        print(f"  S_M^iso = ||B||_F^2 * G_bar: {mp.nstr(S_M_iso, 10)}")
        print(f"  Actual Feshbach trace S_M:   {mp.nstr(S_M, 10)}")
        print(f"  Trace discrepancy Delta S:   {mp.nstr(Delta_S, 10)}")
        print(f"  Relative excess epsilon_M:   {mp.nstr(epsilon_M, 10)}")
        print(f"  D_KS distance:               {mp.nstr(D_KS, 10)}")
        print(f"  KS discrepancy bound B_KS:   {mp.nstr(B_KS, 10)} (slack: {mp.nstr(ks_slack, 4)}x)")
        print()
        print("  --- THEOREM 3: MASTER ENCLOSURE & RELATIVE BOUNDS ---")
        print(f"  Master Bound B_master:       {mp.nstr(B_master, 10)} (slack: {mp.nstr(master_slack, 4)}x)")
        print(f"  Master Bound verified:       {S_M <= B_master}")
        print(f"  |eps_M| <= C_geom * D_KS:    {abs(epsilon_M) <= eps_bound_exact} (bound: {mp.nstr(eps_bound_exact, 8)})")
        print(f"  |eps_M| <= Univ * D_KS:      {abs(epsilon_M) <= eps_bound_univ} (bound: {mp.nstr(eps_bound_univ, 8)})")
        print()
        print("  --- MECHANISM 2: GROUND-STATE PROJECTION AUDIT ---")
        print(f"  ||v^(P)||^2 (ground in P):   {mp.nstr(norm_v_P_sq, 10)}")
        print(f"  ||v^(Q)||^2 (ground in Q):   {mp.nstr(norm_v_Q_sq, 10)}")
        print(f"  Projected energy E_proj:     {mp.nstr(E_proj, 10)}")
        print(f"  Rayleigh shift Delta E_Fesh: {mp.nstr(Delta_E_Fesh, 10)}")
        print(f"  Projection bound E_proj/del: {mp.nstr(B_proj, 10)} (slack: {mp.nstr(proj_slack, 4)}x)")
        print(f"  Projection bound verified:   {Delta_E_Fesh <= B_proj}")
        print(f"  Per-mode average trace S/p:  {mp.nstr(S_bar_P, 10)}")
        print(f"  Cutoff runtime: {time.time() - t_m0:.2f} s")
        print()

        rec = {
            "M": M,
            "p": p_dim,
            "q": q_dim,
            "frob_B_sq": frob_B_sq,
            "delta_M": delta_M,
            "diam": spec_diam,
            "G_bar": g_mean,
            "C_geom_exact": C_geom_exact,
            "C_geom_bound": C_geom_bound,
            "D_KS": D_KS,
            "S_M": S_M,
            "S_M_iso": S_M_iso,
            "epsilon_M": epsilon_M,
            "B_master": B_master,
            "norm_v_P_sq": norm_v_P_sq,
            "E_proj": E_proj,
            "Delta_E_Fesh": Delta_E_Fesh,
            "B_proj": B_proj,
            "S_bar_P": S_bar_P,
            "runtime": time.time() - t_m0,
        }
        records.append(rec)

    # ============================================================
    # SYNTHESIS TABLE
    # ============================================================
    print("=" * 110)
    print("CELL 103 SYNTHESIS: UNIVERSAL SCALING, ISOTROPY, & GROUND-STATE PROJECTION")
    print("=" * 110)
    headers = [
        "M", "delta_M", "||B||_F^2", "G_bar", "C_exact", "C_bound",
        "D_KS", "S_M", "S_iso", "eps_M", "E_proj", "Delta_E_Fesh"
    ]
    print(f"{headers[0]:>4} | {headers[1]:>8} | {headers[2]:>9} | {headers[3]:>8} | {headers[4]:>8} | {headers[5]:>8} | "
          f"{headers[6]:>7} | {headers[7]:>8} | {headers[8]:>8} | {headers[9]:>8} | {headers[10]:>10} | {headers[11]:>12}")
    print("-" * 110)

    for r in records:
        print(
            f"{r['M']:>4} | "
            f"{mp.nstr(r['delta_M'], 4):>8} | "
            f"{mp.nstr(r['frob_B_sq'], 5):>9} | "
            f"{mp.nstr(r['G_bar'], 4):>8} | "
            f"{mp.nstr(r['C_geom_exact'], 4):>8} | "
            f"{mp.nstr(r['C_geom_bound'], 4):>8} | "
            f"{mp.nstr(r['D_KS'], 4):>7} | "
            f"{mp.nstr(r['S_M'], 4):>8} | "
            f"{mp.nstr(r['S_M_iso'], 4):>8} | "
            f"{mp.nstr(r['epsilon_M'], 4):>8} | "
            f"{mp.nstr(r['E_proj'], 4):>10} | "
            f"{mp.nstr(r['Delta_E_Fesh'], 5):>12}"
        )
    print("=" * 110)
    print()

    # Final Audit Summary
    print("--- PRECISE THEORETICAL CONCLUSIONS ---")
    all_thm1 = all(r["C_geom_exact"] <= r["C_geom_bound"] for r in records)
    all_master = all(r["S_M"] <= r["B_master"] for r in records)
    all_proj = all(r["Delta_E_Fesh"] <= r["B_proj"] for r in records)

    print(f"1. Theorem 1 (C_geom <= diam / delta_M):        {'CERTIFIED' if all_thm1 else 'VIOLATED'}")
    print(f"2. Theorem 3 (S_M <= Master Enclosure):        {'CERTIFIED' if all_master else 'VIOLATED'}")
    print(f"3. Mechanism 2 (Delta E_Fesh <= E_proj/delta):  {'CERTIFIED' if all_proj else 'VIOLATED'}")
    print(f"Total script runtime: {time.time() - t_start:.2f} s")
    print()

    print("=" * 80)
    print("CELL 103 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
