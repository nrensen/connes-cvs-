"""
CELL 130 — EXACT COMPONENT DECOMPOSITION AUDIT: MATRIX EQUIVALENCE & SPECTRUM OF THE PRIME FORM
================================================================================================

Target Gate: Gate 1 (Finite-N Spectral Mechanism & Asymptotic Tail Extinction,
             Milestone M-G1.5 / Semiclassical Potential Barrier & Mode Quenching)

Target Propositions & Tested Hypotheses:
  1. Tripartite Component Decomposition Identity:
       Q_{prime}^{even} === -W_tilde + D_tilde^{per} + Delta_D_tilde
     Audit whether the even Galerkin truncation of André Weil's prime operator decomposes
     identically into the negative step potential -W_tilde, the periodic translation-defect
     multiplier D_tilde^{per}, and the finite-interval boundary truncation matrix Delta_D_tilde.
  2. Exact Closed-Form Boundary Integration:
       Delta_D_tilde[m, n] evaluated via exact closed-form trigonometric formulas
       derived from midpoint symmetry u = t - (1/2)log(q), cross-validated against
       50-digit numerical quadrature.
  3. Spectral Signatures of Components:
       - W_tilde >= 0 (positive semi-definite step-potential Gram matrix)
       - D_tilde^{per} >= 0 (diagonal positive semi-definite periodic translation matrix)
       - Delta_D_tilde <= 0 (negative semi-definite finite-interval boundary penalty)
  4. Trace and Norm Conservation:
       Tr(Q_{prime}^{even}) === -Tr(W_tilde) + Tr(D_tilde^{per}) + Tr(Delta_D_tilde)
       ||Q_{prime}^{even} - (-W_tilde + D_tilde^{per} + Delta_D_tilde)||_{max} < 10^-45

Falsification Criteria:
  - If ||Q_{prime}^{even} - (-W_tilde + D_tilde^{per} + Delta_D_tilde)||_{max} >= 10^-40,
    the algebraic decomposition is decisively falsified.
  - If ||Delta_D_tilde^{closed} - Delta_D_tilde^{quad}||_{max} >= 10^-40, the closed-form
    boundary integration formula is mathematically incorrect.
"""

import time
import mpmath as mp
from connes_cvs.operator import (
    psi_prime,
    psi_prime_deriv,
    prime_powers_up_to,
)

# ============================================================
# CONFIGURATION & PARAMETERS
# ============================================================

mp.mp.dps = 50

C_PARAM = 13
L_PARAM = mp.log(C_PARAM)
GROUND_DPS = 50

# Discrete dimensions N for full algebraic audit
N_GRID = [4, 8, 12, 16, 20, 24]

# Dimensions for independent quadrature validation of boundary integrals
N_QUAD_VALIDATION = [4, 8]


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


def build_Q_prime_even(N: int, prime_data: list) -> mp.matrix:
    """
    Construct Q_{prime}^{even} on R^{N+1} via the canonical divided-difference matrix.
    """
    psi_pr_vals = [psi_prime(n, L_PARAM, prime_data) for n in range(N + 1)]
    psi_pr_derivs = [psi_prime_deriv(n, L_PARAM, prime_data) for n in range(N + 1)]
    Q_full = assemble_divided_difference_matrix(psi_pr_vals, psi_pr_derivs, N)
    V_even = canonical_even_projector(N)
    Q_even = V_even.T * Q_full * V_even
    return mp.mpf("0.5") * (Q_even + Q_even.T)


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

    # Precompute psi_prime(k) for k in [0, 2N]
    psi_vals = [psi_prime(k, L_PARAM, prime_data) for k in range(2 * N + 1)]

    # (0, 0) entry
    W[0, 0] = -psi_0_d

    # (0, n) and (n, 0) entries
    for n in range(1, dim):
        val = -sqrt2 * psi_vals[n] / mp.mpf(n)
        W[0, n] = val
        W[n, 0] = val

    # (m, n) entries for m, n >= 1
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
    exact closed-form trigonometric formulas derived in Theorem 3.3.
    """
    dim = N + 1
    Delta_D = mp.matrix(dim, dim)
    PI = mp.pi

    # Entries with m = 0 or n = 0 are identically zero
    for m in range(1, dim):
        for n in range(m, dim):
            entry = mp.mpf("0")
            for (_, logq, w_q) in prime_data:
                theta_m = PI * mp.mpf(m) * logq / L_PARAM
                theta_n = PI * mp.mpf(n) * logq / L_PARAM
                sin_prod = mp.sin(theta_m) * mp.sin(theta_n)

                if m == n:
                    # J_{mm}(q) = (1/2)*log(q) - (L / (4*pi*m)) * sin(2*pi*m*log(q)/L)
                    J_val = mp.mpf("0.5") * logq - (L_PARAM / (mp.mpf("4") * PI * mp.mpf(m))) * mp.sin(mp.mpf("2") * theta_m)
                else:
                    # J_{mn}(q) = (L / (2*pi)) * [ sin(pi*(m-n)*log(q)/L)/(m-n) - sin(pi*(m+n)*log(q)/L)/(m+n) ]
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


def build_Delta_D_tilde_quad(N: int, prime_data: list) -> mp.matrix:
    """
    Construct the (N+1) x (N+1) boundary truncation matrix Delta_D_tilde by direct
    numerical quadrature over the boundary strip [0, log(q)]:
      (Delta_D)_{mn} = -(1/L) * sum_{q <= c} w_q * int_0^{log(q)} [phi_m(t) - phi_m(t - log(q))]
                                                                * [phi_n(t) - phi_n(t - log(q))] dt
    Used as an independent cross-check of the closed-form formulas.
    """
    dim = N + 1
    Delta_D = mp.matrix(dim, dim)
    sqrt2 = mp.sqrt(mp.mpf("2"))
    PI = mp.pi

    for m in range(1, dim):
        a_m = mp.mpf("2") * PI * mp.mpf(m) / L_PARAM
        for n in range(m, dim):
            a_n = mp.mpf("2") * PI * mp.mpf(n) / L_PARAM
            total_integral = mp.mpf("0")

            for (_, logq, w_q) in prime_data:
                def integrand(t, a_m=a_m, a_n=a_n, logq=logq):
                    diff_m = sqrt2 * (mp.cos(a_m * t) - mp.cos(a_m * (t - logq)))
                    diff_n = sqrt2 * (mp.cos(a_n * t) - mp.cos(a_n * (t - logq)))
                    return diff_m * diff_n

                val = mp.quad(integrand, [mp.mpf("0"), logq])
                total_integral += -(w_q / L_PARAM) * val

            Delta_D[m, n] = total_integral
            Delta_D[n, m] = total_integral

    return mp.mpf("0.5") * (Delta_D + Delta_D.T)


def matrix_max_norm(A: mp.matrix) -> mp.mpf:
    """Return max_{i, j} |A_{ij}|."""
    max_val = mp.mpf("0")
    for i in range(A.rows):
        for j in range(A.cols):
            v = abs(A[i, j])
            if v > max_val:
                max_val = v
    return max_val


def matrix_frobenius_norm(A: mp.matrix) -> mp.mpf:
    """Return sqrt(sum_{i, j} A_{ij}^2)."""
    s = mp.mpf("0")
    for i in range(A.rows):
        for j in range(A.cols):
            s += A[i, j] ** 2
    return mp.sqrt(s)


def matrix_trace(A: mp.matrix) -> mp.mpf:
    """Return Tr(A)."""
    tr = mp.mpf("0")
    for i in range(min(A.rows, A.cols)):
        tr += A[i, i]
    return tr


def symmetric_eigenvalues(A: mp.matrix) -> list[mp.mpf]:
    """Compute and return sorted eigenvalues of symmetric matrix A in ascending order."""
    vals, _ = mp.eigsy(A)
    return sorted(vals)


# ============================================================
# MAIN AUDIT SUITE
# ============================================================

def main():
    t_start = time.time()

    print("=" * 80)
    print("CELL 130 — EXACT COMPONENT DECOMPOSITION AUDIT")
    print(f"Parameters: c = {C_PARAM}, L = {mp.nstr(L_PARAM, 12)}, dps = {mp.mp.dps}")
    print("=" * 80)

    # 1. Load prime data
    prime_data, primes_list = prime_powers_up_to(C_PARAM)
    print(f"Loaded {len(prime_data)} prime powers up to c = {C_PARAM}:")
    for (q, logq, w) in prime_data:
        print(f"  q = {q:2d}: log(q) = {mp.nstr(logq, 8)}, weight w_q = {mp.nstr(w, 8)}")
    print("-" * 80)

    # 2. Quadrature vs Closed-Form Boundary Cross-Validation
    print("PHASE 1: QUADRATURE VS CLOSED-FORM BOUNDARY CROSS-VALIDATION")
    print("-" * 80)
    quad_passed = True
    for N in N_QUAD_VALIDATION:
        t0 = time.time()
        D_closed = build_Delta_D_tilde_closed(N, prime_data)
        D_quad = build_Delta_D_tilde_quad(N, prime_data)
        diff_mat = D_closed - D_quad
        max_err = matrix_max_norm(diff_mat)
        elapsed = time.time() - t0
        print(f"  N = {N:2d}: ||Delta_D_closed - Delta_D_quad||_max = {mp.nstr(max_err, 6)}  ({elapsed:.2f}s)")
        if max_err > mp.mpf("1e-45"):
            quad_passed = False
    print(f"Quadrature cross-validation status: {'PASSED (err < 10^-45)' if quad_passed else 'FAILED'}")
    print("-" * 80)

    # 3. Component Decomposition Audit across N_GRID
    print("PHASE 2: EXACT COMPONENT DECOMPOSITION AUDIT ACROSS DIMENSIONS N")
    print("-" * 80)
    print(f"{'N':>3} | {'||Res||_max':>14} | {'||Res||_F':>14} | {'Tr(Q_prime)':>14} | {'Tr(RHS)':>14} | {'|Delta Tr|':>12}")
    print("-" * 80)

    audit_records = {}
    decomposition_passed = True

    for N in N_GRID:
        t0 = time.time()

        # Build Pathway 1: Q_{prime}^{even}
        Q_prime_even = build_Q_prime_even(N, prime_data)

        # Build Pathway 2: -W_tilde + D_tilde^{per} + Delta_D_tilde
        W_tilde = build_W_tilde(N, prime_data)
        D_per = build_D_tilde_per(N, prime_data)
        Delta_D = build_Delta_D_tilde_closed(N, prime_data)

        RHS = -W_tilde + D_per + Delta_D
        Residual = Q_prime_even - RHS

        res_max = matrix_max_norm(Residual)
        res_frob = matrix_frobenius_norm(Residual)

        tr_Q = matrix_trace(Q_prime_even)
        tr_RHS = matrix_trace(RHS)
        tr_diff = abs(tr_Q - tr_RHS)

        if res_max > mp.mpf("1e-45"):
            decomposition_passed = False

        # Spectral analysis of components (sorted in ascending order)
        eigs_Q = symmetric_eigenvalues(Q_prime_even)
        eigs_W = symmetric_eigenvalues(W_tilde)
        eigs_Dper = sorted([D_per[i, i] for i in range(D_per.rows)])
        eigs_DeltaD = symmetric_eigenvalues(Delta_D)

        audit_records[N] = {
            "res_max": res_max,
            "res_frob": res_frob,
            "tr_Q": tr_Q,
            "tr_RHS": tr_RHS,
            "tr_diff": tr_diff,
            "lambda_min_Q": eigs_Q[0],
            "lambda_max_Q": eigs_Q[-1],
            "lambda_min_W": eigs_W[0],
            "lambda_max_W": eigs_W[-1],
            "lambda_min_Dper": eigs_Dper[0],
            "lambda_max_Dper": eigs_Dper[-1],
            "lambda_min_DeltaD": eigs_DeltaD[0],
            "lambda_max_DeltaD": eigs_DeltaD[-1],
            "norm_W": matrix_frobenius_norm(W_tilde),
            "norm_Dper": matrix_frobenius_norm(D_per),
            "norm_DeltaD": matrix_frobenius_norm(Delta_D),
            "elapsed": time.time() - t0,
        }

        print(
            f"{N:3d} | "
            f"{mp.nstr(res_max, 6):>14} | "
            f"{mp.nstr(res_frob, 6):>14} | "
            f"{mp.nstr(tr_Q, 6):>14} | "
            f"{mp.nstr(tr_RHS, 6):>14} | "
            f"{mp.nstr(tr_diff, 4):>12}"
        )

    print("-" * 80)
    print(f"Exact decomposition identity status: {'CERTIFIED (||Res||_max < 10^-45)' if decomposition_passed else 'FALSIFIED'}")
    print("-" * 80)

    # 4. Detailed Component Spectral Summary
    print("PHASE 3: COMPONENT SPECTRAL & NORM SIGNATURES")
    print("-" * 80)
    for N in [4, 8, 16, 24]:
        rec = audit_records[N]
        print(f"Dimension N = {N:2d} (dim = {N+1}):")
        print(f"  Q_prime_even:    spectrum in [{mp.nstr(rec['lambda_min_Q'], 6)}, {mp.nstr(rec['lambda_max_Q'], 6)}]")
        print(f"  Step Potential W: spectrum in [{mp.nstr(rec['lambda_min_W'], 6)}, {mp.nstr(rec['lambda_max_W'], 6)}]  (||W||_F = {mp.nstr(rec['norm_W'], 6)})")
        print(f"  Periodic Defect:  spectrum in [{mp.nstr(rec['lambda_min_Dper'], 6)}, {mp.nstr(rec['lambda_max_Dper'], 6)}]  (||D_per||_F = {mp.nstr(rec['norm_Dper'], 6)})")
        print(f"  Boundary Defect:  spectrum in [{mp.nstr(rec['lambda_min_DeltaD'], 6)}, {mp.nstr(rec['lambda_max_DeltaD'], 6)}]  (||Delta_D||_F = {mp.nstr(rec['norm_DeltaD'], 6)})")
        print(f"  Step Potential W >= 0: {'TRUE' if rec['lambda_min_W'] >= -mp.mpf('1e-45') else 'FALSE'}")
        print(f"  Boundary Defect <= 0:  {'TRUE' if rec['lambda_max_DeltaD'] <= mp.mpf('1e-45') else 'FALSE'}")
        print("-" * 40)

    t_total = time.time() - t_start
    print(f"Total execution time: {t_total:.2f}s")
    print("=" * 80)
    print("CELL 130 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
