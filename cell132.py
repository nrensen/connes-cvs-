"""
CELL 132 — GEOMETRIC DISSECTION OF THE COMPETITION MINIMUM & COUPLED CANCELLATION ON PHI^PERP
=============================================================================================

Target Gate: Gate 1 (Finite-N Spectral Mechanism & Asymptotic Tail Extinction,
             Milestone M-G1.5 / Continuum Spectral Threshold & Mode Quenching)

Target Propositions & Tested Hypotheses:
  1. Full Matrix Operator Norm Residual Audit:
       R_{Phi^perp} === Q_hat_{even} - (Omega_hat + Delta_Q_hat_{arch} - K_hat_{neg} + Q_hat_{pole})
     Certify that ||R_{Phi^perp}||_{op} < 10^{-45} dynamically to 50 decimal digits, confirming
     exact operator conservation on Phi^perp without hidden cancellations.
  2. Algebraic Rayleigh Identity Audit:
       Delta_{Rayleigh} = |R_{net}(w) - (R_{comp}(w) + R_{arch}(w) + R_{pole}(w))| < 10^{-45}
     for all negative eigenvectors w of Q_hat_{comp}.
  3. Negative Eigenspace Dimension (k_{neg}) & Spectrum of Q_hat_{comp}:
       Track the number of negative eigenvalues k_{neg} and the lowest eigenvalues
       mu_0 <= mu_1 <= mu_2 across N in [16, 64].
  4. The Coupled Positivity Mechanism on Negative Directions:
       Evaluate the Rayleigh quotient budget along the principal negative directions
       w_0, w_1, ... of Q_hat_{comp} = Omega_hat - K_hat_{neg}, testing how Delta_Q_hat_{arch}
       and Q_hat_{pole} restore positivity (rho_restore = (R_arch + R_pole) / |R_comp| > 1).
  5. Physical Modal Anatomy of v_{bad} = U_{cont} * w_0:
       Profile the Fourier energy distribution |(v_{bad})_m|^2 across m in {0, ..., N}.

Falsification Criteria:
  - If ||R_{Phi^perp}||_{op} > 10^{-40}, the operator decomposition audit FAILS.
  - If Delta_{Rayleigh} > 10^{-40}, the Rayleigh quotient budget is algebraically inconsistent.
"""

import time
import mpmath as mp
from connes_cvs.operator import (
    psi_prime,
    psi_prime_deriv,
    psi_pole,
    psi_pole_deriv,
    prime_powers_up_to,
    h_plus,
)
from cell import get_galerkin_matrix

# ============================================================
# CONFIGURATION & PARAMETERS
# ============================================================

mp.mp.dps = 50

C_PARAM = 13
L_PARAM = mp.log(C_PARAM)
T_PARAM = 600
GROUND_DPS = 50

# Bound-state core dimension at c = 13
N_BOUND = 11

# Discrete dimension grid for continuum subspace sweep
N_GRID = [16, 20, 24, 28, 32, 40, 48, 64]


# ============================================================
# MATHEMATICAL UTILITIES & OPERATOR BUILDERS
# ============================================================

def canonical_even_projector(N: int) -> mp.matrix:
    """
    Construct the (2N+1) x (N+1) orthonormal isometry V_even mapping the canonical
    even basis v in R^{N+1} to the full exponential basis c in R^{2N+1}.
    """
    dim_full = 2 * N + 1
    dim_even = N + 1
    V = mp.matrix(dim_full, dim_even)
    V[N, 0] = mp.mpf("1")
    inv_sqrt2 = mp.mpf("1") / mp.sqrt(mp.mpf("2"))
    for m in range(1, dim_even):
        V[N + m, m] = inv_sqrt2
        V[N - m, m] = inv_sqrt2
    return V


def assemble_divided_difference_matrix(psi_vals: list, psi_deriv_vals: list, N: int) -> mp.matrix:
    """
    Assemble the (2N+1) x (2N+1) symmetric Galerkin matrix from basis functional
    values psi(n) and psi'(n) for n in [0, N], using exact parity identities.
    """
    dim = 2 * N + 1
    full_psi = [mp.mpf("0")] * dim
    full_psi_d = [mp.mpf("0")] * dim

    for n in range(N + 1):
        p = psi_vals[n]
        pd = psi_deriv_vals[n]
        full_psi[N + n] = p
        full_psi[N - n] = -p
        full_psi_d[N + n] = pd
        full_psi_d[N - n] = pd

    Q = mp.matrix(dim, dim)
    for i in range(dim):
        m = i - N
        p_m = full_psi[i]
        for j in range(i, dim):
            n = j - N
            if m == n:
                val = full_psi_d[j]
            else:
                val = (p_m - full_psi[j]) / mp.mpf(m - n)
            Q[i, j] = val
            Q[j, i] = val
    return Q


def build_W_tilde(N: int, prime_data: list) -> mp.matrix:
    """
    Construct the (N+1) x (N+1) step-potential matrix W_tilde using exact Fourier moments:
      W_tilde[0, 0] = -psi_prime'(0)
      W_tilde[0, n] = -(sqrt(2)/n) * psi_prime(n)
      W_tilde[m, n] = -psi_prime(m-n)/(m-n) - psi_prime(m+n)/(m+n)   (m != n, m, n >= 1)
      W_tilde[m, m] = -psi_prime'(0) - psi_prime(2m)/(2m)            (m >= 1)
    """
    dim = N + 1
    W = mp.matrix(dim, dim)
    psi_0_d = psi_prime_deriv(0, L_PARAM, prime_data)
    sqrt2 = mp.sqrt(mp.mpf("2"))

    psi_vals = [psi_prime(k, L_PARAM, prime_data) for k in range(2 * N + 1)]

    W[0, 0] = -psi_0_d

    for n in range(1, dim):
        val = -sqrt2 * psi_vals[n] / mp.mpf(n)
        W[0, n] = val
        W[n, 0] = val

    for m in range(1, dim):
        for n in range(m, dim):
            if m == n:
                val = -psi_0_d - psi_vals[2 * m] / mp.mpf(2 * m)
            else:
                diff_idx = abs(m - n)
                sum_idx = m + n
                val = -psi_vals[diff_idx] / mp.mpf(diff_idx) - psi_vals[sum_idx] / mp.mpf(sum_idx)
            W[m, n] = val
            W[n, m] = val

    return mp.mpf("0.5") * (W + W.T)


def build_D_tilde_per(N: int, prime_data: list) -> mp.matrix:
    """
    Construct the (N+1) x (N+1) diagonal periodic translation-defect matrix:
      D_tilde^{per} = diag(0, 4*M(1), 4*M(2), ..., 4*M(N))
      where M(m) = sum_{q <= c} w_q * sin^2(pi * m * log(q) / L).
    """
    dim = N + 1
    D_per = mp.matrix(dim, dim)
    PI = mp.pi

    for m in range(1, dim):
        M_m = mp.mpf("0")
        for (_, logq, w_q) in prime_data:
            theta = PI * mp.mpf(m) * logq / L_PARAM
            M_m += w_q * (mp.sin(theta) ** 2)
        D_per[m, m] = mp.mpf("4") * M_m

    return D_per


def build_Delta_D_tilde_closed(N: int, prime_data: list) -> mp.matrix:
    """
    Construct the (N+1) x (N+1) boundary truncation matrix Delta_D_tilde using the
    exact closed-form trigonometric formulas from Theorem 3.3.
    """
    dim = N + 1
    Delta_D = mp.matrix(dim, dim)
    PI = mp.pi

    for m in range(1, dim):
        for n in range(m, dim):
            entry = mp.mpf("0")
            for (_, logq, w_q) in prime_data:
                theta_m = PI * mp.mpf(m) * logq / L_PARAM
                theta_n = PI * mp.mpf(n) * logq / L_PARAM
                sin_prod = mp.sin(theta_m) * mp.sin(theta_n)

                if m == n:
                    J_val = mp.mpf("0.5") * logq - (L_PARAM / (mp.mpf("4") * PI * mp.mpf(m))) * mp.sin(mp.mpf("2") * theta_m)
                else:
                    diff_m_n = mp.mpf(m - n)
                    sum_m_n = mp.mpf(m + n)
                    sin_diff = mp.sin(PI * diff_m_n * logq / L_PARAM)
                    sin_sum = mp.sin(PI * sum_m_n * logq / L_PARAM)
                    J_val = (L_PARAM / (mp.mpf("2") * PI)) * (sin_diff / diff_m_n - sin_sum / sum_m_n)

                term = -(mp.mpf("8") / L_PARAM) * w_q * sin_prod * J_val
                entry += term

            Delta_D[m, n] = entry
            Delta_D[n, m] = entry

    return mp.mpf("0.5") * (Delta_D + Delta_D.T)


def symmetric_eigenvalues(A: mp.matrix) -> list[mp.mpf]:
    """Compute and return sorted eigenvalues of symmetric matrix A in ascending order."""
    vals, _ = mp.eigsy(A)
    return sorted(vals)


def symmetric_eigendecomposition(A: mp.matrix) -> tuple[list[mp.mpf], mp.matrix]:
    """
    Compute sorted eigenvalues and corresponding orthonormal eigenvector matrix V
    such that A = V * diag(evals) * V^T, with evals[0] <= evals[1] <= ...
    """
    dim = A.rows
    vals, V_raw = mp.eigsy(A)
    idx = sorted(range(dim), key=lambda i: vals[i])
    evals = [vals[i] for i in idx]
    V = mp.matrix(dim, dim)
    for col_out, col_in in enumerate(idx):
        for row in range(dim):
            V[row, col_out] = V_raw[row, col_in]
    return evals, V


def matrix_max_norm(A: mp.matrix) -> mp.mpf:
    """Return max_{i, j} |A_{ij}|."""
    m = mp.mpf("0")
    for i in range(A.rows):
        for j in range(A.cols):
            v = abs(A[i, j])
            if v > m:
                m = v
    return m


def matrix_frobenius_norm(A: mp.matrix) -> mp.mpf:
    """Return sqrt(sum_{i, j} A_{ij}^2)."""
    s = mp.mpf("0")
    for i in range(A.rows):
        for j in range(A.cols):
            s += A[i, j] ** 2
    return mp.sqrt(s)


def operator_norm_sym(A: mp.matrix) -> mp.mpf:
    """Return ||A||_{op} = max_i |lambda_i(A)| for symmetric matrix A."""
    evals = symmetric_eigenvalues(A)
    return max(abs(evals[0]), abs(evals[-1]))


# ============================================================
# MAIN AUDIT SUITE
# ============================================================

def run_cell132_audit():
    t_start = time.time()

    print("=" * 80)
    print("CELL 132 — GEOMETRIC DISSECTION OF THE COMPETITION MINIMUM & COUPLED CANCELLATION")
    print(f"Parameters: c = {C_PARAM}, L = {mp.nstr(L_PARAM, 12)}, T = {T_PARAM}, dps = {mp.mp.dps}")
    print(f"Bound-State Cutoff: N_bound = {N_BOUND} (continuum subspace dimension q = N - 10)")
    print("=" * 80)

    prime_data, primes_list = prime_powers_up_to(C_PARAM)
    print(f"Loaded {len(prime_data)} prime powers up to c = {C_PARAM}.")
    print("-" * 80)

    table1_rows = []  # Residual audit
    table2_rows = []  # Competition spectrum
    table3_rows = []  # Rayleigh budget of w_bad
    table4_rows = []  # Physical modal anatomy

    max_global_r_op = mp.mpf("0")
    max_global_rayleigh_err = mp.mpf("0")

    for N in N_GRID:
        t_n_start = time.time()
        dim_even = N + 1
        q_cont = N - N_BOUND + 1  # dimension of Phi^perp: (N+1) - 11 = N - 10

        # ----------------------------------------------------
        # 1. Assembling Full Galerkin Operator & Constituents
        # ----------------------------------------------------
        # Full Galerkin matrix from persistent cache
        Q_full, _ = get_galerkin_matrix(
            c=C_PARAM,
            N=N,
            T=T_PARAM,
            dps=GROUND_DPS,
            verbose=False,
        )
        V_even = canonical_even_projector(N)
        Q_even = V_even.T * Q_full * V_even
        Q_even = mp.mpf("0.5") * (Q_even + Q_even.T)

        # Prime constituent pieces
        psi_pr_vals = [psi_prime(n, L_PARAM, prime_data) for n in range(dim_even)]
        psi_pr_derivs = [psi_prime_deriv(n, L_PARAM, prime_data) for n in range(dim_even)]
        Q_prime_full = assemble_divided_difference_matrix(psi_pr_vals, psi_pr_derivs, N)
        Q_even_prime = V_even.T * Q_prime_full * V_even
        Q_even_prime = mp.mpf("0.5") * (Q_even_prime + Q_even_prime.T)

        # Pole piece
        psi_po_vals = [psi_pole(n, L_PARAM) for n in range(dim_even)]
        psi_po_derivs = [psi_pole_deriv(n, L_PARAM) for n in range(dim_even)]
        Q_pole_full = assemble_divided_difference_matrix(psi_po_vals, psi_po_derivs, N)
        Q_even_pole = V_even.T * Q_pole_full * V_even
        Q_even_pole = mp.mpf("0.5") * (Q_even_pole + Q_even_pole.T)

        # Archimedean piece by exact difference from certified Galerkin operator
        Q_even_arch = Q_even - Q_even_prime - Q_even_pole
        Q_even_arch = mp.mpf("0.5") * (Q_even_arch + Q_even_arch.T)

        # ----------------------------------------------------
        # 2. Decomposed Operators: W_tilde, D_per, Delta_D, K_neg
        # ----------------------------------------------------
        W_tilde = build_W_tilde(N, prime_data)
        D_per = build_D_tilde_per(N, prime_data)
        Delta_D = build_Delta_D_tilde_closed(N, prime_data)

        # Total negative potential K_neg = W_tilde - Delta_D >= 0
        K_neg = W_tilde - Delta_D
        K_neg = mp.mpf("0.5") * (K_neg + K_neg.T)

        # Diagonal multiplier D_mult and backbone Omega_diag
        PI = mp.pi
        D_mult = mp.matrix(dim_even, dim_even)
        Omega_diag = mp.matrix(dim_even, dim_even)
        for m in range(dim_even):
            a_m = mp.mpf("2") * PI * mp.mpf(m) / L_PARAM if m > 0 else mp.mpf("0")
            h_val = h_plus(a_m, GROUND_DPS)
            D_mult[m, m] = h_val
            Omega_diag[m, m] = h_val + D_per[m, m]

        # Archimedean off-diagonal defect
        Delta_Q_arch = Q_even_arch - D_mult
        Delta_Q_arch = mp.mpf("0.5") * (Delta_Q_arch + Delta_Q_arch.T)

        # Full-space exact identity check:
        # Q_even === Omega_diag + Delta_Q_arch - K_neg + Q_even_pole
        Q_full_space_reconstructed = Omega_diag + Delta_Q_arch - K_neg + Q_even_pole
        r_full_space = operator_norm_sym(Q_even - Q_full_space_reconstructed)

        # ----------------------------------------------------
        # 3. Spectral Decomposition & Continuum Isometry U_cont
        # ----------------------------------------------------
        evals_even, V_even_eigs = symmetric_eigendecomposition(Q_even)

        # Continuum subspace isometry U_cont in R^{(N+1) x q_cont}
        U_cont = mp.matrix(dim_even, q_cont)
        for col in range(q_cont):
            orig_col = N_BOUND + col  # indices 11, 12, ..., N
            for row in range(dim_even):
                U_cont[row, col] = V_even_eigs[row, orig_col]

        # ----------------------------------------------------
        # 4. Compress Constituent Operators to Phi^perp
        # ----------------------------------------------------
        Q_hat_even = U_cont.T * Q_even * U_cont
        Q_hat_even = mp.mpf("0.5") * (Q_hat_even + Q_hat_even.T)

        Omega_hat = U_cont.T * Omega_diag * U_cont
        Omega_hat = mp.mpf("0.5") * (Omega_hat + Omega_hat.T)

        Delta_Q_hat_arch = U_cont.T * Delta_Q_arch * U_cont
        Delta_Q_hat_arch = mp.mpf("0.5") * (Delta_Q_hat_arch + Delta_Q_hat_arch.T)

        K_hat_neg = U_cont.T * K_neg * U_cont
        K_hat_neg = mp.mpf("0.5") * (K_hat_neg + K_hat_neg.T)

        Q_hat_comp = Omega_hat - K_hat_neg
        Q_hat_comp = mp.mpf("0.5") * (Q_hat_comp + Q_hat_comp.T)

        Q_hat_pole = U_cont.T * Q_even_pole * U_cont
        Q_hat_pole = mp.mpf("0.5") * (Q_hat_pole + Q_hat_pole.T)

        # ----------------------------------------------------
        # 5. Table 1: Full Matrix Residual Audit on Phi^perp
        # ----------------------------------------------------
        R_mat = Q_hat_even - (Omega_hat + Delta_Q_hat_arch - K_hat_neg + Q_hat_pole)
        R_mat = mp.mpf("0.5") * (R_mat + R_mat.T)
        r_max = matrix_max_norm(R_mat)
        r_frob = matrix_frobenius_norm(R_mat)
        r_op = operator_norm_sym(R_mat)

        if r_op > max_global_r_op:
            max_global_r_op = r_op

        table1_rows.append({
            "N": N, "q": q_cont,
            "r_max": r_max, "r_frob": r_frob, "r_op": r_op,
            "r_full": r_full_space,
        })

        # ----------------------------------------------------
        # 6. Table 2: Eigendecomposition of Q_hat_comp
        # ----------------------------------------------------
        mu_vals, W_comp = symmetric_eigendecomposition(Q_hat_comp)
        k_neg = sum(1 for m in mu_vals if m < 0)
        mu_0 = mu_vals[0]
        mu_1 = mu_vals[1] if q_cont > 1 else mp.mpf("0")
        mu_2 = mu_vals[2] if q_cont > 2 else mp.mpf("0")

        table2_rows.append({
            "N": N, "q": q_cont, "k_neg": k_neg,
            "mu_0": mu_0, "mu_1": mu_1, "mu_2": mu_2,
        })

        # ----------------------------------------------------
        # 7. Table 3: Rayleigh Quotient Budget along w_bad
        # ----------------------------------------------------
        w_bad = mp.matrix(q_cont, 1)
        for i in range(q_cont):
            w_bad[i, 0] = W_comp[i, 0]

        R_comp = (w_bad.T * Q_hat_comp * w_bad)[0, 0]
        R_omega = (w_bad.T * Omega_hat * w_bad)[0, 0]
        R_neg = (w_bad.T * K_hat_neg * w_bad)[0, 0]
        R_arch = (w_bad.T * Delta_Q_hat_arch * w_bad)[0, 0]
        R_pole = (w_bad.T * Q_hat_pole * w_bad)[0, 0]
        R_net = (w_bad.T * Q_hat_even * w_bad)[0, 0]

        # Algebraic Rayleigh identity residual
        delta_rayleigh = abs(R_net - (R_comp + R_arch + R_pole))
        if delta_rayleigh > max_global_rayleigh_err:
            max_global_rayleigh_err = delta_rayleigh

        if R_comp < 0:
            rho_restore = (R_arch + R_pole) / abs(R_comp)
        else:
            rho_restore = mp.mpf("nan")

        table3_rows.append({
            "N": N,
            "R_comp": R_comp,
            "R_omega": R_omega,
            "R_neg": R_neg,
            "R_arch": R_arch,
            "R_pole": R_pole,
            "R_net": R_net,
            "delta_rayleigh": delta_rayleigh,
            "rho": rho_restore,
        })

        # ----------------------------------------------------
        # 8. Table 4: Physical Modal Anatomy of v_bad
        # ----------------------------------------------------
        v_bad = U_cont * w_bad  # in R^{N+1}
        e_0 = (v_bad[0, 0]) ** 2
        e_low = sum((v_bad[m, 0]) ** 2 for m in range(1, min(4, dim_even)))
        e_mid = sum((v_bad[m, 0]) ** 2 for m in range(4, min(11, dim_even)))
        e_tail = sum((v_bad[m, 0]) ** 2 for m in range(11, dim_even))

        peak_idx = 0
        peak_val = mp.mpf("0")
        for m in range(dim_even):
            v_m_sq = (v_bad[m, 0]) ** 2
            if v_m_sq > peak_val:
                peak_val = v_m_sq
                peak_idx = m

        table4_rows.append({
            "N": N,
            "e_0": e_0,
            "e_low": e_low,
            "e_mid": e_mid,
            "e_tail": e_tail,
            "peak_m": peak_idx,
            "peak_val": peak_val,
        })

        t_elapsed = time.time() - t_n_start
        print(f"Completed N = {N:2d} (dim = {dim_even:2d}, q = {q_cont:2d}) in {t_elapsed:.2f}s | "
              f"||R||_op = {float(r_op):.2e}, mu_0 = {float(mu_0):+.6f}, "
              f"R_arch = {float(R_arch):+.6f}, R_net = {float(R_net):+.6f}, "
              f"|Delta_Rayleigh| = {float(delta_rayleigh):.2e}")

    # ============================================================
    # FORMATTED DIAGNOSTIC REPORT
    # ============================================================

    print("\n" + "-" * 80)
    print("TABLE 1: FULL MATRIX RESIDUAL AUDIT ON PHI^PERP")
    print("Residual: R = Q_hat - (Omega_hat + Delta_Q_hat_arch - K_hat_neg + Q_hat_pole)")
    print("-" * 80)
    print(f"{'N':>4} | {'dim':>4} | {'q':>3} | {'||R||_max':>14} | {'||R||_F':>14} | {'||R||_op':>14}")
    print("-" * 80)
    for r in table1_rows:
        print(f"{r['N']:4d} | {r['N']+1:4d} | {r['q']:3d} | {float(r['r_max']):14.4e} | {float(r['r_frob']):14.4e} | {float(r['r_op']):14.4e}")
    print("-" * 80)

    # Dynamic status verification
    if max_global_r_op < mp.mpf("1e-45"):
        print(f"VERIFICATION STATUS: PASSED (max ||R||_op = {float(max_global_r_op):.4e} < 1e-45).")
    else:
        print(f"VERIFICATION STATUS: FAILED (max ||R||_op = {float(max_global_r_op):.4e} >= 1e-45)!")
    print("-" * 80)

    print("\n" + "-" * 80)
    print("TABLE 2: SPECTRUM OF COMPETITION OPERATOR Q_hat_comp = Omega_hat - K_hat_neg")
    print("-" * 80)
    print(f"{'N':>4} | {'q':>3} | {'k_neg':>5} | {'mu_0 (Floor)':>14} | {'mu_1':>14} | {'mu_2':>14}")
    print("-" * 80)
    for r in table2_rows:
        print(f"{r['N']:4d} | {r['q']:3d} | {r['k_neg']:5d} | {float(r['mu_0']):14.6f} | {float(r['mu_1']):14.6f} | {float(r['mu_2']):14.6f}")
    print("-" * 80)

    print("\n" + "-" * 80)
    print("TABLE 3: RAYLEIGH QUOTIENT BUDGET ALONG PRINCIPAL VULNERABLE STATE w_bad")
    print("Net: R_net = R_comp + R_arch + R_pole | Restoring Ratio: rho = (R_arch + R_pole) / |R_comp|")
    print("-" * 80)
    print(f"{'N':>4} | {'R_comp':>10} | {'R_Omega':>10} | {'R_neg':>10} | {'R_arch':>10} | {'R_pole':>10} | {'R_net':>10} | {'rho_restore':>11}")
    print("-" * 80)
    for r in table3_rows:
        rho_str = f"{float(r['rho']):11.6f}" if not mp.isnan(r['rho']) else f"{'N/A (pos)':>11}"
        print(f"{r['N']:4d} | {float(r['R_comp']):10.6f} | {float(r['R_omega']):10.6f} | {float(r['R_neg']):10.6f} | {float(r['R_arch']):10.6f} | {float(r['R_pole']):10.6f} | {float(r['R_net']):10.6f} | {rho_str}")
    print("-" * 80)

    # Dynamic Rayleigh balance verification
    if max_global_rayleigh_err < mp.mpf("1e-45"):
        print(f"RAYLEIGH IDENTITY AUDIT: PASSED (max |Delta_Rayleigh| = {float(max_global_rayleigh_err):.4e} < 1e-45).")
    else:
        print(f"RAYLEIGH IDENTITY AUDIT: FAILED (max |Delta_Rayleigh| = {float(max_global_rayleigh_err):.4e} >= 1e-45)!")
    print("-" * 80)

    print("\n" + "-" * 80)
    print("TABLE 4: PHYSICAL MODAL ANATOMY OF VULNERABLE STATE v_bad = U_cont * w_bad")
    print("-" * 80)
    print(f"{'N':>4} | {'|v_0|^2':>10} | {'Low (1..3)':>10} | {'Mid (4..10)':>11} | {'Tail (>=11)':>11} | {'Peak m*':>7} | {'|v_{m*}|^2':>10}")
    print("-" * 80)
    for r in table4_rows:
        print(f"{r['N']:4d} | {float(r['e_0']):10.6f} | {float(r['e_low']):10.6f} | {float(r['e_mid']):11.6f} | {float(r['e_tail']):11.6f} | {r['peak_m']:7d} | {float(r['peak_val']):10.6f}")
    print("-" * 80)

    # Synthesis at N = 64
    r1 = table1_rows[-1]
    r2 = table2_rows[-1]
    r3 = table3_rows[-1]
    r4 = table4_rows[-1]
    print(f"\nSYNTHESIS AT MAXIMAL RESOLUTION N = {r1['N']}:")
    print(f"  Operator Norm Residual ||R||_op:     {float(r1['r_op']):.4e} (Max norm: {float(r1['r_max']):.4e})")
    print(f"  Negative Eigenspace Dimension k_neg: {r2['k_neg']} / {r2['q']}")
    print(f"  Lowest Competition Eigenvalue mu_0:  {float(r2['mu_0']):+.8f}")
    print(f"  Second Competition Eigenvalue mu_1:  {float(r2['mu_1']):+.8f}")
    print(f"  Third Competition Eigenvalue mu_2:   {float(r2['mu_2']):+.8f}")
    print(f"  Competition Deficit R_comp:          {float(r3['R_comp']):+.8f}")
    print(f"  Archimedean Form Contribution R_arch:{float(r3['R_arch']):+.8f}")
    print(f"  Pole Regularization R_pole:          {float(r3['R_pole']):+.8f}")
    print(f"  Net Quadratic Form R_net:            {float(r3['R_net']):+.8f}")
    print(f"  Algebraic Residual |Delta_Rayleigh|: {float(r3['delta_rayleigh']):.4e}")
    if not mp.isnan(r3['rho']):
        if r3['rho'] > 1:
            print(f"  Restoring Force Efficiency rho:      {float(r3['rho']):.6f} (> 1 confirms net positivity restoration)")
        else:
            print(f"  Restoring Force Efficiency rho:      {float(r3['rho']):.6f} (<= 1: restoration incomplete along w_bad)")
    print(f"  Physical Modal Distribution:         Zero: {float(r4['e_0']):.4f}, Low: {float(r4['e_low']):.4f}, Mid: {float(r4['e_mid']):.4f}, Tail: {float(r4['e_tail']):.4f}")
    print(f"  Peak Modal Concentration:            m* = {r4['peak_m']} (weight = {float(r4['peak_val']):.6f})")

    t_total = time.time() - t_start
    print(f"\nTotal execution time: {t_total:.2f}s")
    print("=" * 80)
    print("CELL 132 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell132_audit()
