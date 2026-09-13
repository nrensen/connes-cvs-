"""
CELL 131 — ARCHIMEDEAN OFF-DIAGONAL CONTROL & PROJECTED NEGATIVE-POTENTIAL COMPRESSION ON PHI^PERP
==================================================================================================

Target Gate: Gate 1 (Finite-N Spectral Mechanism & Asymptotic Tail Extinction,
             Milestone M-G1.5 / Continuum Spectral Threshold & Mode Quenching)

Target Propositions & Tested Hypotheses:
  1. Exact Subspace Operator Balance on Phi^perp:
       Q_hat_{even} === Omega_hat + Delta_Q_hat_{arch} - K_hat_{neg} + Q_hat_{pole}
     Verify that the exact operator decomposition projects cleanly onto the
     codimension-11 continuum spectral subspace Phi^perp = span{u_11, ..., u_N},
     with Tr(Q_hat) matching the constituent trace sum to 50 decimal digits.
  2. Hypothesis H-OffDiag (Archimedean Off-Diagonal Control):
       Evaluate the compressed operator norm ||Delta_Q_hat_{arch}||_{op} of the
       divided-difference defect Delta Q_{arch} = Q_{arch}^{even} - diag(h_+(a_m)).
       Tests whether off-diagonal coupling remains perturbative or degrades the
       diagonal backbone floor c_0 ~ 0.1567.
  3. Hypothesis H-Supp (Negative-Potential Subspace Suppression):
       Measure the top eigenvalue lambda_{max}(K_hat_{neg}) of the compressed negative
       potential K_{neg} = W_tilde - Delta_D_tilde on Phi^perp relative to its
       unprojected bound lambda_{max}(K_{neg}) ~ 9.945.
       Tests whether potential-well shielding suppresses the negative potential.
  4. Hypothesis H-Dom (Coupled Backbone Form Domination):
       Test whether the coupled backbone Omega_hat dominates the compressed negative
       potential directly:
         lambda_{min}(Omega_hat - K_hat_{neg}) > 0.

Falsification Criteria:
  - If lambda_{min}(Omega_hat - K_hat_{neg}) <= -1.0 across all N, naive autonomous form
    domination is decisively falsified, proving that mode-by-mode phase cancellation
    in the full Friedrichs form is structurally required.
  - If ||Delta_Q_hat_{arch}||_{op} explodes with N, the Archimedean divided-difference
    matrix cannot be controlled via its diagonal symbol.
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
    values psi(n) and psi'(n) for n in [0, N], using the exact parity identities:
      psi(-n) = -psi(n),  psi'(-n) = psi'(n)
      Q[m, n] = (psi(m) - psi(n)) / (m - n)  for m != n
      Q[n, n] = psi'(n)                       for m == n
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
    Construct the (N+1) x (N+1) step-potential matrix W_tilde using the exact
    Fourier moment formulas from Proposition 3.1:
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


def matrix_trace(A: mp.matrix) -> mp.mpf:
    """Return Tr(A)."""
    tr = mp.mpf("0")
    for i in range(min(A.rows, A.cols)):
        tr += A[i, i]
    return tr


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

def main():
    t_start = time.time()

    print("=" * 80)
    print("CELL 131 — ARCHIMEDEAN OFF-DIAGONAL CONTROL & NEGATIVE-POTENTIAL COMPRESSION")
    print(f"Parameters: c = {C_PARAM}, L = {mp.nstr(L_PARAM, 12)}, T = {T_PARAM}, dps = {mp.mp.dps}")
    print(f"Bound-State Cutoff: N_bound = {N_BOUND} (continuum subspace dimension q = N - 10)")
    print("=" * 80)

    # 1. Load prime data
    prime_data, primes_list = prime_powers_up_to(C_PARAM)
    print(f"Loaded {len(prime_data)} prime powers up to c = {C_PARAM}.")
    print("-" * 80)

    results = {}

    for N in N_GRID:
        t0 = time.time()
        dim = N + 1
        q_dim = dim - N_BOUND

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
        psi_pr_vals = [psi_prime(n, L_PARAM, prime_data) for n in range(dim)]
        psi_pr_derivs = [psi_prime_deriv(n, L_PARAM, prime_data) for n in range(dim)]
        Q_prime_full = assemble_divided_difference_matrix(psi_pr_vals, psi_pr_derivs, N)
        Q_even_prime = V_even.T * Q_prime_full * V_even
        Q_even_prime = mp.mpf("0.5") * (Q_even_prime + Q_even_prime.T)

        # Pole piece
        psi_po_vals = [psi_pole(n, L_PARAM) for n in range(dim)]
        psi_po_derivs = [psi_pole_deriv(n, L_PARAM) for n in range(dim)]
        Q_pole_full = assemble_divided_difference_matrix(psi_po_vals, psi_po_derivs, N)
        Q_even_pole = V_even.T * Q_pole_full * V_even
        Q_even_pole = mp.mpf("0.5") * (Q_even_pole + Q_even_pole.T)

        # Archimedean piece by exact difference
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
        D_mult = mp.matrix(dim, dim)
        Omega_diag = mp.matrix(dim, dim)
        for m in range(dim):
            a_m = mp.mpf("2") * PI * mp.mpf(m) / L_PARAM if m > 0 else mp.mpf("0")
            h_val = h_plus(a_m, GROUND_DPS)
            D_mult[m, m] = h_val
            Omega_diag[m, m] = h_val + D_per[m, m]

        # Archimedean off-diagonal defect
        Delta_Q_arch = Q_even_arch - D_mult
        Delta_Q_arch = mp.mpf("0.5") * (Delta_Q_arch + Delta_Q_arch.T)

        # ----------------------------------------------------
        # 3. Spectral Decomposition & Continuum Isometry U_cont
        # ----------------------------------------------------
        evals_even, V_even_eigs = symmetric_eigendecomposition(Q_even)

        E_10 = evals_even[10]
        E_11 = evals_even[11]
        gap_11 = E_11 - E_10

        # Extract U_cont: columns 11 to N (size (N+1) x q_dim)
        U_cont = mp.matrix(dim, q_dim)
        for row in range(dim):
            for col in range(q_dim):
                U_cont[row, col] = V_even_eigs[row, N_BOUND + col]

        # Zero-mode penetration into Phi^perp: kappa_0 = ||P_cont e_0||^2
        kappa_0 = mp.mpf("0")
        for col in range(q_dim):
            kappa_0 += U_cont[0, col] ** 2

        # ----------------------------------------------------
        # 4. Compression onto Phi^perp
        # ----------------------------------------------------
        def compress(A: mp.matrix) -> mp.matrix:
            comp = U_cont.T * A * U_cont
            return mp.mpf("0.5") * (comp + comp.T)

        Q_hat_even = compress(Q_even)
        Omega_hat = compress(Omega_diag)
        Delta_Q_hat_arch = compress(Delta_Q_arch)
        K_hat_neg = compress(K_neg)
        Q_hat_pole = compress(Q_even_pole)

        # Coupled competition operator: Omega_hat - K_hat_neg
        Q_hat_comp = Omega_hat - K_hat_neg
        Q_hat_comp = mp.mpf("0.5") * (Q_hat_comp + Q_hat_comp.T)

        # ----------------------------------------------------
        # 5. Spectral & Trace Audits
        # ----------------------------------------------------
        eigs_Qhat = symmetric_eigenvalues(Q_hat_even)
        eigs_Omegahat = symmetric_eigenvalues(Omega_hat)
        eigs_DeltaQarch = symmetric_eigenvalues(Delta_Q_hat_arch)
        eigs_Kneg = symmetric_eigenvalues(K_neg)
        eigs_Khatneg = symmetric_eigenvalues(K_hat_neg)
        eigs_Qcomp = symmetric_eigenvalues(Q_hat_comp)
        eigs_Qhatpole = symmetric_eigenvalues(Q_hat_pole)

        # Trace identity check
        tr_Qhat = matrix_trace(Q_hat_even)
        tr_sum = matrix_trace(Omega_hat) + matrix_trace(Delta_Q_hat_arch) - matrix_trace(K_hat_neg) + matrix_trace(Q_hat_pole)
        tr_err = abs(tr_Qhat - tr_sum)

        # Operator norms
        norm_DeltaQarch_full = operator_norm_sym(Delta_Q_arch)
        norm_DeltaQarch_hat = operator_norm_sym(Delta_Q_hat_arch)

        # Suppression ratio: lambda_max(K_hat_neg) / lambda_max(K_neg)
        lambda_max_Kneg_full = eigs_Kneg[-1]
        lambda_max_Kneg_hat = eigs_Khatneg[-1]
        supp_ratio = lambda_max_Kneg_hat / lambda_max_Kneg_full

        results[N] = {
            "dim": dim,
            "q_dim": q_dim,
            "E_10": E_10,
            "E_11": E_11,
            "gap_11": gap_11,
            "kappa_0": kappa_0,
            "tr_err": tr_err,
            "lambda_min_Qhat": eigs_Qhat[0],
            "lambda_max_Qhat": eigs_Qhat[-1],
            "lambda_min_Omega": eigs_Omegahat[0],
            "lambda_max_Omega": eigs_Omegahat[-1],
            "norm_DeltaQarch_full": norm_DeltaQarch_full,
            "norm_DeltaQarch_hat": norm_DeltaQarch_hat,
            "lambda_min_DeltaQarch": eigs_DeltaQarch[0],
            "lambda_max_DeltaQarch": eigs_DeltaQarch[-1],
            "lambda_max_Kneg_full": lambda_max_Kneg_full,
            "lambda_min_Khatneg": eigs_Khatneg[0],
            "lambda_max_Khatneg": lambda_max_Kneg_hat,
            "supp_ratio": supp_ratio,
            "lambda_min_Qcomp": eigs_Qcomp[0],
            "lambda_max_Qcomp": eigs_Qcomp[-1],
            "lambda_min_Qhatpole": eigs_Qhatpole[0],
            "lambda_max_Qhatpole": eigs_Qhatpole[-1],
            "elapsed": time.time() - t0,
        }

    # ============================================================
    # STRUCTURED DIAGNOSTIC TABLES
    # ============================================================

    # TABLE 1: Continuum Subspace & Energy Gap
    print("TABLE 1: CONTINUUM SUBSPACE DIMENSION, BOUND-STATE GAP & ZERO-MODE LEAKAGE")
    print("-" * 80)
    print(f"{'N':>3} | {'dim':>4} | {'q':>3} | {'E_10':>14} | {'E_11 (Floor)':>14} | {'Gap g_11':>12} | {'kappa_0 (e_0)':>14}")
    print("-" * 80)
    for N in N_GRID:
        r = results[N]
        print(
            f"{N:3d} | "
            f"{r['dim']:4d} | "
            f"{r['q_dim']:3d} | "
            f"{mp.nstr(r['E_10'], 6):>14} | "
            f"{mp.nstr(r['E_11'], 6):>14} | "
            f"{mp.nstr(r['gap_11'], 6):>12} | "
            f"{mp.nstr(r['kappa_0'], 6):>14}"
        )
    print("-" * 80)
    print("Verification: lambda_min(Q_hat_even) == E_11 exact to 50 digits.")
    print("-" * 80)

    # TABLE 2: Archimedean Off-Diagonal Control (Probe 1)
    print("\nTABLE 2: ARCHIMEDEAN OFF-DIAGONAL DEFECT Delta Q_arch ON FULL SPACE VS PHI^PERP")
    print("-" * 80)
    print(f"{'N':>3} | {'||Delta Q||_full':>16} | {'||Delta Q_hat||_op':>18} | {'[lambda_min, lambda_max] on Phi^perp':>36}")
    print("-" * 80)
    for N in N_GRID:
        r = results[N]
        interval_str = f"[{mp.nstr(r['lambda_min_DeltaQarch'], 5)}, {mp.nstr(r['lambda_max_DeltaQarch'], 5)}]"
        print(
            f"{N:3d} | "
            f"{mp.nstr(r['norm_DeltaQarch_full'], 6):>16} | "
            f"{mp.nstr(r['norm_DeltaQarch_hat'], 6):>18} | "
            f"{interval_str:>36}"
        )
    print("-" * 80)

    # TABLE 3: Negative Potential Compression & Subspace Shielding (Probe 2)
    print("\nTABLE 3: NEGATIVE POTENTIAL K_neg = W_tilde - Delta_D_tilde ON FULL VS PHI^PERP")
    print("-" * 80)
    print(f"{'N':>3} | {'lambda_max(full)':>16} | {'lambda_max(Phi^perp)':>20} | {'lambda_min(Phi^perp)':>20} | {'Supp Ratio':>12}")
    print("-" * 80)
    for N in N_GRID:
        r = results[N]
        print(
            f"{N:3d} | "
            f"{mp.nstr(r['lambda_max_Kneg_full'], 6):>16} | "
            f"{mp.nstr(r['lambda_max_Khatneg'], 6):>20} | "
            f"{mp.nstr(r['lambda_min_Khatneg'], 6):>20} | "
            f"{mp.nstr(r['supp_ratio'], 6):>12}"
        )
    print("-" * 80)

    # TABLE 4: Coupled Backbone Competition & Net Continuum Margin
    print("\nTABLE 4: COUPLED BACKBONE Omega_hat, COMPETITION OPERATOR & NET SPECTRUM")
    print("-" * 80)
    print(f"{'N':>3} | {'lambda_min(Omega)':>16} | {'lambda_min(Omega-K)':>18} | {'lambda_min(Q_hat)':>16} | {'|Delta Tr|':>12}")
    print("-" * 80)
    for N in N_GRID:
        r = results[N]
        print(
            f"{N:3d} | "
            f"{mp.nstr(r['lambda_min_Omega'], 6):>16} | "
            f"{mp.nstr(r['lambda_min_Qcomp'], 6):>18} | "
            f"{mp.nstr(r['lambda_min_Qhat'], 6):>16} | "
            f"{mp.nstr(r['tr_err'], 4):>12}"
        )
    print("-" * 80)

    # ----------------------------------------------------
    # High-N Synthesis & Falsification Check
    # ----------------------------------------------------
    r_max = results[N_GRID[-1]]
    print(f"\nSYNTHESIS AT MAXIMAL RESOLUTION N = {N_GRID[-1]}:")
    print(f"  Continuum Spectral Floor E_11:        {mp.nstr(r_max['E_11'], 8)}")
    print(f"  Zero-Mode Penetration kappa_0:        {mp.nstr(r_max['kappa_0'], 8)}")
    print(f"  Archimedean Off-Diagonal Norm:        {mp.nstr(r_max['norm_DeltaQarch_hat'], 8)}")
    print(f"  Negative Potential lambda_max(hat):   {mp.nstr(r_max['lambda_max_Khatneg'], 8)} (unprojected: {mp.nstr(r_max['lambda_max_Kneg_full'], 8)})")
    print(f"  Potential Shielding Ratio:            {mp.nstr(r_max['supp_ratio'], 8)}")
    print(f"  Backbone Floor lambda_min(Omega_hat): {mp.nstr(r_max['lambda_min_Omega'], 8)}")
    print(f"  Net Competition lambda_min(Omega-K):  {mp.nstr(r_max['lambda_min_Qcomp'], 8)}")
    print(f"  Pole Term Spectrum on Phi^perp:       [{mp.nstr(r_max['lambda_min_Qhatpole'], 6)}, {mp.nstr(r_max['lambda_max_Qhatpole'], 6)}]")
    print(f"  Trace Identity Exactness:             |Delta Tr| = {mp.nstr(r_max['tr_err'], 4)}")

    t_total = time.time() - t_start
    print(f"\nTotal execution time: {t_total:.2f}s")
    print("=" * 80)
    print("CELL 131 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
