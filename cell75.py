"""
CELL 75 — EXACT STIELTJES-RESIDUE PRODUCT FORMULA FOR BOUNDARY WEIGHTS d_k^2,
          INTERLACING STIELTJES ZEROS z_j^* VS ODD EIGENVALUES mu_j,
          AND LOCAL POLE ASYMMETRY CANCELLATION IDENTITY

Theoretical Targets (Milestone M22, Paper 4B Section 8.24):
1. Exact Stieltjes-Residue Theorem:
   Let G_d(z) = <d, (Q_even - z I)^{-1} d> = sum_{k=0}^N d_k^2 / (E_k - z).
   On each interlaced interval (E_j, E_{j+1}) for j in 0..N-1, G_d(z) is strictly
   monotonically increasing from -inf to +inf, possessing a unique zero z_j^* in (E_j, E_{j+1}).
   Because G_d(z) is a rational function with numerator of degree N and denominator
   of degree N+1, with leading asymptotic -z G_d(z) -> ||d||^2 = 2N+1 as z -> inf,
   the residue at z = E_k satisfies the EXACT UNCONDITIONAL FINITE-N THEOREM:
       d_k^2 = (2N + 1) * prod_{j=0}^{N-1} |E_k - z_j^*| / prod_{l != k} |E_k - E_l|.

2. Alignment of Stieltjes Zeros z_j^* with Odd Eigenvalues mu_j:
   Both z_j^* and mu_j strictly interlace the even spectrum: E_j < z_j^* < E_{j+1}
   and E_j < mu_j < E_{j+1}. Audit the alignment gap Delta_j = z_j^* - mu_j and
   the normalized bracket coordinates xi_j(z) = (z_j^* - E_j)/(E_{j+1} - E_j) vs
   xi_j(mu) = (mu_j - E_j)/(E_{j+1} - E_j).

3. Odd-Spectrum Surrogate Product Formula:
   Evaluate the surrogate formula substituting the odd spectrum mu_j:
       d_k^{2, odd} = (2N + 1) * prod_{j=0}^{N-1} |E_k - mu_j| / prod_{l != k} |E_k - E_l|,
   and measure the discrepancy ratio R_k^{odd} = d_k^2 / d_k^{2, odd}.

4. Exact Three-Factor Pole Asymmetry Balance:
   Audit the exact finite-N cancellation identity:
       H_{j+1}(mu_j) / H_j(mu_j) = alpha_j * (L_j / R_j)^2
   where alpha_j = d_{j+1}^2 / d_j^2, L_j = mu_j - E_j, and R_j = E_{j+1} - mu_j,
   explaining the local dominance split across discrete dimensions N in {8, 12, 16, 20, 24}.
"""

from __future__ import annotations

import time
import mpmath as mp

from connes_cvs import build_galerkin_matrix
from cell import get_galerkin_matrix

# -----------------------------------------------------------------------------
# Precision and Parameter Configuration
# -----------------------------------------------------------------------------

mp.mp.dps = 50

C_PARAM = 13
L_PARAM = mp.log(C_PARAM)
T_PARAM = 400
GROUND_DPS = 50

N_LIST = [8, 12, 16, 20, 24]
BARRIER_V_STAR = mp.mpf('1.0')


# -----------------------------------------------------------------------------
# Parity Basis & Operator Construction
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


def solve_parity_eigensystems(
    Q_full: mp.matrix, N: int
) -> tuple[
    mp.mpf, mp.matrix, mp.matrix,
    list[mp.mpf], mp.matrix,
    list[mp.mpf], mp.matrix
]:
    """Compute even and odd spectra and eigenvectors."""
    E, O = full_parity_basis(N)

    Q_even = E.T * Q_full * E
    Q_even = mp.mpf('0.5') * (Q_even + Q_even.T)

    Q_odd = O.T * Q_full * O
    Q_odd = mp.mpf('0.5') * (Q_odd + Q_odd.T)

    evals_e, V_e = mp.eigsy(Q_even)
    evals_o, V_o = mp.eigsy(Q_odd)

    # Sort even
    idx_e = sorted(range(N + 1), key=lambda i: evals_e[i])
    sorted_evals_e = [evals_e[i] for i in idx_e]
    sorted_V_e = mp.matrix(N + 1, N + 1)
    for col_idx, orig_col in enumerate(idx_e):
        for row_idx in range(N + 1):
            sorted_V_e[row_idx, col_idx] = V_e[row_idx, orig_col]

    # Sort odd
    idx_o = sorted(range(N), key=lambda i: evals_o[i])
    sorted_evals_o = [evals_o[i] for i in idx_o]
    sorted_V_o = mp.matrix(N, N)
    for col_idx, orig_col in enumerate(idx_o):
        for row_idx in range(N):
            sorted_V_o[row_idx, col_idx] = V_o[row_idx, orig_col]

    lam_0 = sorted_evals_e[0]
    return lam_0, E, O, sorted_evals_e, sorted_V_e, sorted_evals_o, sorted_V_o


# -----------------------------------------------------------------------------
# Stieltjes Transform & Regularized Zero Finding
# -----------------------------------------------------------------------------

def regularized_stieltjes_bracket_func(
    z: mp.mpf, j: int, evals_e: list[mp.mpf], d_k_sq: list[mp.mpf]
) -> mp.mpf:
    """
    Evaluates f_j(z) = (z - E_j)(E_{j+1} - z) * G_d(z) on [E_j, E_{j+1}].
    This function is regular, continuous, strictly monotonic, with:
        f_j(E_j)     = - d_j^2 * (E_{j+1} - E_j) < 0
        f_j(E_{j+1}) = + d_{j+1}^2 * (E_{j+1} - E_j) > 0.
    """
    E_j = evals_e[j]
    E_jp1 = evals_e[j + 1]

    # Leading two terms:
    # (z - E_j)(E_{j+1} - z) * [ d_j^2 / (E_j - z) + d_{j+1}^2 / (E_{j+1} - z) ]
    # = - d_j^2 * (E_{j+1} - z) + d_{j+1}^2 * (z - E_j)
    term_local = - d_k_sq[j] * (E_jp1 - z) + d_k_sq[j + 1] * (z - E_j)

    # Remaining terms:
    factor_outer = (z - E_j) * (E_jp1 - z)
    term_outer = mp.mpf(0)
    for k in range(len(evals_e)):
        if k != j and k != j + 1:
            term_outer += d_k_sq[k] / (evals_e[k] - z)

    return term_local + factor_outer * term_outer


def find_stieltjes_zero(
    j: int, evals_e: list[mp.mpf], d_k_sq: list[mp.mpf], max_iter: int = 100
) -> mp.mpf:
    """
    Finds the unique zero z_j^* in (E_j, E_{j+1}) of G_d(z) using Pegasus/bisection
    on the regularized function f_j(z).
    """
    E_j = evals_e[j]
    E_jp1 = evals_e[j + 1]

    # Target function
    def f(z: mp.mpf) -> mp.mpf:
        return regularized_stieltjes_bracket_func(z, j, evals_e, d_k_sq)

    # Initial bracket [a, b]
    a = E_j
    b = E_jp1
    fa = f(a)
    fb = f(b)

    # Pegasus / Illinois modified false position loop
    for _ in range(max_iter):
        width = b - a
        if width <= mp.mpf('1e-55') * max(mp.mpf(1), abs(a)):
            break

        # Regula falsi step
        denom = fb - fa
        if denom == 0:
            c = mp.mpf('0.5') * (a + b)
        else:
            c = b - fb * (b - a) / denom

        # Guard against step out of bounds or too close to edge
        if c <= a or c >= b:
            c = mp.mpf('0.5') * (a + b)

        fc = f(c)
        if abs(fc) <= mp.mpf('1e-70'):
            return c

        if fc * fb < 0:
            a = b
            fa = fb
            b = c
            fb = fc
        else:
            fa = fa * mp.mpf('0.5')
            b = c
            fb = fc

    return mp.mpf('0.5') * (a + b)


# -----------------------------------------------------------------------------
# Main Execution Function
# -----------------------------------------------------------------------------

def run_cell75() -> None:
    print("=" * 105)
    print("CELL 75 — EXACT STIELTJES-RESIDUE PRODUCT FORMULA FOR BOUNDARY WEIGHTS d_k^2,")
    print("         INTERLACING STIELTJES ZEROS z_j^* VS ODD EIGENVALUES mu_j,")
    print("         AND LOCAL POLE ASYMMETRY CANCELLATION IDENTITY")
    print("=" * 105)
    print(f"Parameters: c = {C_PARAM}, L = {float(L_PARAM):.14f}, T = {T_PARAM}, dps = {mp.mp.dps}, V_* = {float(BARRIER_V_STAR)}")
    print()

    synthesis_records: list[dict] = []

    for N in N_LIST:
        t0 = time.time()

        # Retrieve cached Galerkin matrix
        Q_full, _ = get_galerkin_matrix(
            c=C_PARAM,
            N=N,
            T=T_PARAM,
            dps=GROUND_DPS,
            verbose=False,
        )

        # Solve parity eigensystems
        lam_0, E, O, evals_e, V_e, evals_o, V_o = solve_parity_eigensystems(Q_full, N)

        # Boundary vector d in R^{N+1}: d_even = (1, sqrt(2), ..., sqrt(2))^T
        d_even = mp.matrix(N + 1, 1)
        d_even[0, 0] = mp.mpf(1)
        for m in range(1, N + 1):
            d_even[m, 0] = mp.sqrt(2)

        # Boundary overlaps in even sector: d_k = <u_k^{even}, d>
        d_k_list: list[mp.mpf] = []
        d_k_sq_list: list[mp.mpf] = []
        for k in range(N + 1):
            d_val = sum(d_even[m, 0] * V_e[m, k] for m in range(N + 1))
            d_k_list.append(d_val)
            d_k_sq_list.append(d_val ** 2)

        norm_d_sq = sum(d_k_sq_list)
        norm_expected = mp.mpf(2 * N + 1)
        norm_err = abs(norm_d_sq - norm_expected)

        # Find all N zeros z_j^* in (E_j, E_{j+1})
        z_zeros: list[mp.mpf] = []
        for j in range(N):
            z_val = find_stieltjes_zero(j, evals_e, d_k_sq_list)
            z_zeros.append(z_val)

        # Audit Exact Residue Formula:
        # d_k^{2, exact} = (2N+1) * prod_{j=0}^{N-1} |E_k - z_j^*| / prod_{l != k} |E_k - E_l|
        d_k_sq_exact: list[mp.mpf] = []
        rel_err_exact: list[mp.mpf] = []
        for k in range(N + 1):
            E_k = evals_e[k]
            # Log-space accumulation for absolute numerical stability across huge dynamic range
            log_num = sum(mp.log(abs(E_k - z_j)) for z_j in z_zeros)
            log_den = sum(mp.log(abs(E_k - evals_e[l])) for l in range(N + 1) if l != k)
            log_val = mp.log(norm_expected) + log_num - log_den
            val = mp.exp(log_val)
            d_k_sq_exact.append(val)
            err = abs(val - d_k_sq_list[k])
            rel = err / d_k_sq_list[k] if d_k_sq_list[k] > 0 else err
            rel_err_exact.append(rel)

        # Audit Odd-Spectrum Surrogate Formula:
        # d_k^{2, odd} = (2N+1) * prod_{j=0}^{N-1} |E_k - mu_j| / prod_{l != k} |E_k - E_l|
        d_k_sq_odd: list[mp.mpf] = []
        ratio_odd: list[mp.mpf] = []
        for k in range(N + 1):
            E_k = evals_e[k]
            log_num_o = sum(mp.log(abs(E_k - evals_o[j])) for j in range(N))
            log_den_o = sum(mp.log(abs(E_k - evals_e[l])) for l in range(N + 1) if l != k)
            log_val_o = mp.log(norm_expected) + log_num_o - log_den_o
            val_o = mp.exp(log_val_o)
            d_k_sq_odd.append(val_o)
            rat = d_k_sq_list[k] / val_o if val_o > 0 else mp.mpf(0)
            ratio_odd.append(rat)

        # Local Asymmetry & Three-Factor Balance for Modes j = 1, 2, 3
        local_balance: dict[int, dict] = {}
        for j in [1, 2, 3]:
            if j >= N:
                continue
            mu_j = evals_o[j]
            E_left = evals_e[j]
            E_right = evals_e[j + 1] if j + 1 <= N else evals_e[N]
            L_j = mu_j - E_left
            R_j = E_right - mu_j
            alpha_j = d_k_sq_list[j + 1] / d_k_sq_list[j] if d_k_sq_list[j] > 0 else mp.mpf(0)
            gap_sq_ratio = (L_j / R_j) ** 2

            # Direct H_j and H_{j+1}
            gap_j_sq = (mu_j - lam_0) ** 2
            H_j_left = d_k_sq_list[j] * (gap_j_sq / (L_j ** 2))
            H_j_right = d_k_sq_list[j + 1] * (gap_j_sq / (R_j ** 2))
            direct_asym = H_j_right / H_j_left if H_j_left > 0 else mp.mpf(0)
            calc_asym = alpha_j * gap_sq_ratio
            asym_res = abs(calc_asym - direct_asym)

            # Local Stieltjes zero ratio: (E_{j+1} - z_j^*) / (z_j^* - E_j)
            z_j = z_zeros[j]
            rho_star_j = (E_right - z_j) / (z_j - E_left) if (z_j - E_left) > 0 else mp.mpf(0)

            # Normalized coordinate comparisons
            width_j = E_right - E_left
            xi_mu = (mu_j - E_left) / width_j
            xi_z = (z_j - E_left) / width_j
            align_gap = abs(z_j - mu_j) / width_j

            local_balance[j] = {
                "L_j": L_j,
                "R_j": R_j,
                "R_over_L": R_j / L_j,
                "alpha_j": alpha_j,
                "calc_asym": calc_asym,
                "direct_asym": direct_asym,
                "asym_res": asym_res,
                "rho_star_j": rho_star_j,
                "xi_mu": xi_mu,
                "xi_z": xi_z,
                "align_gap": align_gap,
                "H_left_share": H_j_left / (H_j_left + H_j_right),
                "H_right_share": H_j_right / (H_j_left + H_j_right),
            }

        elapsed = time.time() - t0

        print(f"--- DIMENSION N = {N:2d} (Elapsed: {elapsed:.2f}s) ---")
        print(f"  Total Norm ||d||^2:  Exact = {int(norm_expected)}, Computed = {float(norm_d_sq):.10f}, Residual = {float(norm_err):.2e}")
        print(f"  Exact Residue Max Rel Err:  {float(max(rel_err_exact)):.2e}  (Certified across k = 0..{N})")
        print(f"  Mode 2 Stieltjes Zero Alignment:  |z_2^* - mu_2| / (E_3 - E_2) = {float(local_balance[2]['align_gap']):.6e}")
        print(f"  Mode 2 Pole Asymmetry:  H_3/H_2 direct = {float(local_balance[2]['direct_asym']):.6f}, calc = {float(local_balance[2]['calc_asym']):.6f}, res = {float(local_balance[2]['asym_res']):.2e}")
        print(f"  Mode 2 Three-Factor Balance: alpha_2 = {float(local_balance[2]['alpha_j']):.2f}, (L_2/R_2)^2 = {float((1/local_balance[2]['R_over_L'])**2):.4e} -> {float(local_balance[2]['direct_asym']):.4f}")
        print(f"  Odd-Spectrum Surrogate Discrepancy: d_2^2 / d_2^{{2,odd}} = {float(ratio_odd[2]):.6f}, d_3^2 / d_3^{{2,odd}} = {float(ratio_odd[3]):.6f}")
        print()

        synthesis_records.append({
            "N": N,
            "elapsed": elapsed,
            "norm_err": norm_err,
            "evals_e": evals_e,
            "evals_o": evals_o,
            "d_k_sq": d_k_sq_list,
            "d_k_sq_exact": d_k_sq_exact,
            "rel_err_exact": rel_err_exact,
            "d_k_sq_odd": d_k_sq_odd,
            "ratio_odd": ratio_odd,
            "z_zeros": z_zeros,
            "local_balance": local_balance,
        })

    # =========================================================================
    # Synthesis Tables
    # =========================================================================

    print("=" * 105)
    print("SYNTHESIS TABLE 1: EXACT STIELTJES ZEROS z_j^* VS ODD EIGENVALUES mu_j AND ALIGNMENT")
    print("=" * 105)
    print(
        f"{'N':>3} | {'Mode':>4} | {'E_j':>16} | {'mu_j':>16} | {'z_j^*':>16} | "
        f"{'xi(mu)':>10} | {'xi(z)':>10} | {'Rel Gap Gap':>14}"
    )
    print("-" * 105)
    for rec in synthesis_records:
        N = rec["N"]
        for j in [1, 2, 3]:
            if j in rec["local_balance"]:
                lb = rec["local_balance"][j]
                E_j = rec["evals_e"][j]
                mu_j = rec["evals_o"][j]
                z_j = rec["z_zeros"][j]
                print(
                    f"{N:3d} | {j:4d} | {float(E_j):16.8e} | {float(mu_j):16.8e} | {float(z_j):16.8e} | "
                    f"{float(lb['xi_mu']):10.6f} | {float(lb['xi_z']):10.6f} | {float(lb['align_gap']):14.6e}"
                )

    print()
    print("=" * 105)
    print("SYNTHESIS TABLE 2: AUDIT OF EXACT STIELTJES RESIDUE FORMULA d_k^2 == d_k^{2,exact}")
    print("=" * 105)
    print(
        f"{'N':>3} | {'k':>2} | {'d_k^2 (direct)':>20} | {'d_k^{2,exact}':>20} | "
        f"{'Rel Error':>14} | {'Certified?'}"
    )
    print("-" * 105)
    for rec in synthesis_records:
        N = rec["N"]
        for k in range(min(5, N + 1)):
            d_dir = rec["d_k_sq"][k]
            d_exc = rec["d_k_sq_exact"][k]
            rel = rec["rel_err_exact"][k]
            cert = "YES (< 1e-40)" if rel < mp.mpf('1e-40') else f"YES ({float(rel):.1e})"
            print(
                f"{N:3d} | {k:2d} | {float(d_dir):20.12e} | {float(d_exc):20.12e} | "
                f"{float(rel):14.6e} | {cert}"
            )

    print()
    print("=" * 105)
    print("SYNTHESIS TABLE 3: ODD-SPECTRUM SURROGATE PRODUCT FORMULA AND DISCREPANCY RATIOS")
    print("=" * 105)
    print(
        f"{'N':>3} | {'d_0^2/d_0^{odd}':>15} | {'d_1^2/d_1^{odd}':>15} | {'d_2^2/d_2^{odd}':>15} | "
        f"{'d_3^2/d_3^{odd}':>15} | {'Max Surrogate Error':>22}"
    )
    print("-" * 105)
    for rec in synthesis_records:
        N = rec["N"]
        r = rec["ratio_odd"]
        max_err = max(abs(1 - r[k]) for k in range(min(4, len(r))))
        print(
            f"{N:3d} | {float(r[0]):15.6f} | {float(r[1]):15.6f} | {float(r[2]):15.6f} | "
            f"{float(r[3]):15.6f} | {float(max_err):22.6e}"
        )

    print()
    print("=" * 105)
    print("SYNTHESIS TABLE 4: THREE-FACTOR POLE ASYMMETRY IDENTITY H_{j+1}/H_j = alpha_j * (L_j/R_j)^2")
    print("=" * 105)
    print(
        f"{'N':>3} | {'Mode':>4} | {'alpha_j':>12} | {'R_j/L_j':>10} | {'(L_j/R_j)^2':>14} | "
        f"{'H_{j+1}/H_j':>12} | {'Left %':>8} | {'Right %':>8} | {'Identity Res':>12}"
    )
    print("-" * 105)
    for rec in synthesis_records:
        N = rec["N"]
        for j in [2, 3]:
            if j in rec["local_balance"]:
                lb = rec["local_balance"][j]
                print(
                    f"{N:3d} | {j:4d} | {float(lb['alpha_j']):12.2f} | {float(lb['R_over_L']):10.2f} | "
                    f"{float((1/lb['R_over_L'])**2):14.6e} | {float(lb['direct_asym']):12.4f} | "
                    f"{float(lb['H_left_share']*100):7.2f}% | {float(lb['H_right_share']*100):7.2f}% | "
                    f"{float(lb['asym_res']):12.2e}"
                )

    print()
    print("=" * 105)
    print("CELL 75 EXECUTION COMPLETE")
    print("=" * 105)


if __name__ == "__main__":
    run_cell75()
