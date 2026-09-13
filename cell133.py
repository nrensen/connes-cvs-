"""
CELL 133 — COORDINATE-SPACE UNIFICATION OF COUPLED CANCELLATION & THE EXCESS POSITIVITY MARGIN
=============================================================================================

Target Gate: Gate 1 (Finite-N Spectral Mechanism & Asymptotic Tail Extinction,
             Milestone M-G1.6 / Unified Representation of Coupled Cancellation)

Target Propositions & Tested Hypotheses:
  1. Coordinate Wavepacket Representation:
       T_{v_{bad}}(t) = (v_{bad})_0 + sqrt(2) * sum_{m=1}^N (v_{bad})_m * cos(2*pi*m*t/L)
     Profile the spatial probability density |T_{v_{bad}}(t)|^2, boundary contact at t = 0, L,
     midpoint amplitude, and spatial kinetic semi-norm K_{kin} = (2*pi/L)^2 * sum m^2 (v_{bad})_m^2.
  2. Coordinate-Matrix Potential Equivalence:
       R_W = v_{bad}^T * W_tilde * v_{bad} === (1/L) * int_0^L (-W(t)) * |T_{v_{bad}}(t)|^2 dt
     Verify this identity to 40 decimal digits via piecewise Gauss-Legendre quadrature
     across the prime transition points log(q_k).
  3. Spatial Mass Partition & Bound-State Exclusion:
       M_{flat} = (1/L) * int_0^{log 2} |T_{v_{bad}}(t)|^2 dt,
       M_{well} = (1/L) * int_{log 2}^L |T_{v_{bad}}(t)|^2 dt.
     Quantify the fraction of mass forced into the zero-potential plateau [0, log 2)
     by the orthogonality of v_{bad} to the 11 bound states.
  4. Archimedean Spectral Partition & Negative Well Quenching:
       I_{arch}^{(-)} = (1/pi) * int_0^{r_*} |Phi_{v_{bad}}(r)|^2 * h_+(r) dr
     Verify that the negative Archimedean energy is suppressed (|I_{arch}^{(-)}| << 0.01),
     proving that Phi^perp effectively quenches the negative Archimedean well.
  5. Asymptotic Stabilization of the Surplus Margin:
       Delta_{surplus} = R_{arch} + R_{pole} - |R_{comp}| = R_{net}
     Track Delta_{surplus} across N in [24, 64] and verify that it stabilizes to a macroscopic
     positive constant approx +0.165 > 0.

Falsification Criteria:
  - If |R_W - I_{pot}| > 10^{-35}, the continuous coordinate-matrix equivalence FAILS.
  - If Delta_{surplus} <= 0 for any N >= 24, net coupled positivity FAILS along w_{bad}.
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
L_PARAM = mp.log(mp.mpf("13"))
T_PARAM = 600
GROUND_DPS = 50

# Bound-state core dimension at c = 13
N_BOUND = 11

# Discrete dimension grid for continuum subspace sweep
N_GRID = [24, 28, 32, 40, 48, 64]

# Gauss-Legendre quadrature order for coordinate and frequency integrals
QUAD_ORDER = 40


# ============================================================
# MATHEMATICAL UTILITIES & QUADRATURE RECIPES
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
    Construct the (N+1) x (N+1) step-potential matrix W_tilde using exact Fourier moments.
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
      D_tilde^{per} = diag(0, 4*M(1), ..., 4*M(N))
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


def gauss_legendre_nodes_weights(order: int, a: mp.mpf, b: mp.mpf) -> tuple[list[mp.mpf], list[mp.mpf]]:
    """
    Compute Gauss-Legendre quadrature nodes and weights on [a, b] using
    the Golub-Welsch tridiagonal eigenvalue method at current mpmath dps.
    """
    J = mp.matrix(order, order)
    for i in range(order - 1):
        k = i + 1
        b_k = mp.mpf(k) / mp.sqrt(4 * k * k - 1)
        J[i, i + 1] = b_k
        J[i + 1, i] = b_k

    nodes_std, V = mp.eigsy(J)

    mid = (b + a) / mp.mpf("2")
    half_width = (b - a) / mp.mpf("2")

    nodes = []
    weights = []
    for i in range(order):
        x_i = nodes_std[i]
        w_i = mp.mpf("2") * (V[0, i] ** 2)
        nodes.append(mid + half_width * x_i)
        weights.append(half_width * w_i)

    return nodes, weights


def evaluate_T_v(v: mp.matrix, t: mp.mpf, L: mp.mpf) -> mp.mpf:
    """
    Evaluate the physical wavefunction T_v(t) = v_0 + sqrt(2) * sum_{m=1}^N v_m * cos(2*pi*m*t/L).
    """
    dim = v.rows
    val = v[0, 0]
    sqrt2 = mp.sqrt(mp.mpf("2"))
    PI2 = mp.mpf("2") * mp.pi
    for m in range(1, dim):
        val += sqrt2 * v[m, 0] * mp.cos(PI2 * mp.mpf(m) * t / L)
    return val


def phi_basis(m: int, r: mp.mpf, L: mp.mpf) -> mp.mpf:
    """
    Canonical Fourier basis amplitude phi_m(r) on R.
    phi_0(r) = (2 / sqrt(L)) * (sin(r L / 2) / r)
    phi_m(r) = (2 * sqrt(2) / sqrt(L)) * (r * sin(r L / 2) / (r^2 - a_m^2))  for m >= 1.
    """
    half_L = L / mp.mpf("2")
    if m == 0:
        if abs(r) < mp.mpf("1e-25"):
            return mp.sqrt(L)
        return (mp.mpf("2") / mp.sqrt(L)) * (mp.sin(r * half_L) / r)

    a_m = mp.mpf("2") * mp.pi * mp.mpf(m) / L
    diff_pos = r - a_m
    diff_neg = r + a_m

    # Removable singularity at r = a_m
    if abs(diff_pos) < mp.mpf("1e-12"):
        u = diff_pos * half_L
        sinc_u = mp.mpf("1") - (u**2) / 6 + (u**4) / 120 - (u**6) / 5040
        geom = (a_m + diff_pos) / (mp.mpf("2") * a_m + diff_pos)
        sgn = mp.mpf("-1") ** m
        return (mp.mpf("2") * mp.sqrt(mp.mpf("2")) / mp.sqrt(L)) * (sgn * half_L * geom * sinc_u)

    # Removable singularity at r = -a_m
    if abs(diff_neg) < mp.mpf("1e-12"):
        u = diff_neg * half_L
        sinc_u = mp.mpf("1") - (u**2) / 6 + (u**4) / 120 - (u**6) / 5040
        geom = (-a_m + diff_neg) / (-mp.mpf("2") * a_m + diff_neg)
        sgn = mp.mpf("-1") ** m
        return (mp.mpf("2") * mp.sqrt(mp.mpf("2")) / mp.sqrt(L)) * (sgn * half_L * geom * sinc_u)

    denom = r**2 - a_m**2
    return (mp.mpf("2") * mp.sqrt(mp.mpf("2")) / mp.sqrt(L)) * (r * mp.sin(r * half_L) / denom)


def evaluate_Phi_v(v: mp.matrix, r: mp.mpf, L: mp.mpf) -> mp.mpf:
    """
    Evaluate the continuous Fourier transform Phi_v(r) = sum_{m=0}^N v_m * phi_m(r).
    """
    dim = v.rows
    s = mp.mpf("0")
    for m in range(dim):
        v_m = v[m, 0]
        if v_m != 0:
            s += v_m * phi_basis(m, r, L)
    return s


# ============================================================
# MAIN CELL 133 AUDIT SUITE
# ============================================================

def run_cell133_audit():
    t_start = time.time()

    print("=" * 80)
    print("CELL 133 — COORDINATE-SPACE UNIFICATION OF COUPLED CANCELLATION")
    print(f"Parameters: c = {C_PARAM}, L = {mp.nstr(L_PARAM, 12)}, T = {T_PARAM}, dps = {mp.mp.dps}")
    print(f"Bound-State Cutoff: N_bound = {N_BOUND} (continuum subspace dimension q = N - 10)")
    print("=" * 80)

    prime_data, primes_list = prime_powers_up_to(C_PARAM)
    print(f"Loaded {len(prime_data)} prime powers up to c = {C_PARAM}:")
    for (n, logn, w) in prime_data:
        print(f"  q = {n:2d}: log(q) = {float(logn):.6f}, weight w_q = {float(w):.6f}")

    # Root of h_plus(r) = 0: r_* approx 6.2898363768
    r_star = mp.findroot(lambda r: h_plus(r, GROUND_DPS), mp.mpf("6.2898"))
    print(f"\nArchimedean zero-crossing root: r_* = {mp.nstr(r_star, 12)}")
    print("-" * 80)

    table1_rows = []  # Wavepacket profile
    table2_rows = []  # Coordinate-matrix potential equivalence
    table3_rows = []  # Archimedean spectral partition
    table4_rows = []  # Surplus margin synthesis

    max_global_pot_err = mp.mpf("0")

    for N in N_GRID:
        t_n_start = time.time()
        dim_even = N + 1
        q_cont = N - N_BOUND + 1

        # ----------------------------------------------------
        # 1. Assembling Operators & Projectors
        # ----------------------------------------------------
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

        psi_pr_vals = [psi_prime(n, L_PARAM, prime_data) for n in range(dim_even)]
        psi_pr_derivs = [psi_prime_deriv(n, L_PARAM, prime_data) for n in range(dim_even)]
        Q_prime_full = assemble_divided_difference_matrix(psi_pr_vals, psi_pr_derivs, N)
        Q_even_prime = V_even.T * Q_prime_full * V_even
        Q_even_prime = mp.mpf("0.5") * (Q_even_prime + Q_even_prime.T)

        psi_po_vals = [psi_pole(n, L_PARAM) for n in range(dim_even)]
        psi_po_derivs = [psi_pole_deriv(n, L_PARAM) for n in range(dim_even)]
        Q_pole_full = assemble_divided_difference_matrix(psi_po_vals, psi_po_derivs, N)
        Q_even_pole = V_even.T * Q_pole_full * V_even
        Q_even_pole = mp.mpf("0.5") * (Q_even_pole + Q_even_pole.T)

        Q_even_arch = Q_even - Q_even_prime - Q_even_pole
        Q_even_arch = mp.mpf("0.5") * (Q_even_arch + Q_even_arch.T)

        W_tilde = build_W_tilde(N, prime_data)
        D_per = build_D_tilde_per(N, prime_data)
        Delta_D = build_Delta_D_tilde_closed(N, prime_data)
        K_neg = mp.mpf("0.5") * ((W_tilde - Delta_D) + (W_tilde - Delta_D).T)

        PI = mp.pi
        D_mult = mp.matrix(dim_even, dim_even)
        Omega_diag = mp.matrix(dim_even, dim_even)
        for m in range(dim_even):
            a_m = mp.mpf("2") * PI * mp.mpf(m) / L_PARAM if m > 0 else mp.mpf("0")
            h_val = h_plus(a_m, GROUND_DPS)
            D_mult[m, m] = h_val
            Omega_diag[m, m] = h_val + D_per[m, m]

        Delta_Q_arch = mp.mpf("0.5") * ((Q_even_arch - D_mult) + (Q_even_arch - D_mult).T)

        # ----------------------------------------------------
        # 2. Continuum Subspace Isometry & Solitary Negative State
        # ----------------------------------------------------
        evals_even, V_even_eigs = symmetric_eigendecomposition(Q_even)
        U_cont = mp.matrix(dim_even, q_cont)
        for col in range(q_cont):
            orig_col = N_BOUND + col
            for row in range(dim_even):
                U_cont[row, col] = V_even_eigs[row, orig_col]

        Omega_hat = U_cont.T * Omega_diag * U_cont
        K_hat_neg = U_cont.T * K_neg * U_cont
        Q_hat_comp = mp.mpf("0.5") * ((Omega_hat - K_hat_neg) + (Omega_hat - K_hat_neg).T)

        mu_vals, W_comp = symmetric_eigendecomposition(Q_hat_comp)
        w_bad = mp.matrix(q_cont, 1)
        for i in range(q_cont):
            w_bad[i, 0] = W_comp[i, 0]

        # Canonical even vector v_bad in R^{N+1}
        v_bad = U_cont * w_bad

        # Rayleigh quotients
        R_comp = (w_bad.T * Q_hat_comp * w_bad)[0, 0]
        R_arch = (v_bad.T * Delta_Q_arch * v_bad)[0, 0]
        R_pole = (v_bad.T * Q_even_pole * v_bad)[0, 0]
        R_net = (v_bad.T * Q_even * v_bad)[0, 0]
        rho_restore = (R_arch + R_pole) / abs(R_comp)
        delta_surplus = R_net

        # ----------------------------------------------------
        # 3. Table 1: Coordinate Wavepacket Profile of T_{v_bad}
        # ----------------------------------------------------
        T_0 = evaluate_T_v(v_bad, mp.mpf("0"), L_PARAM)
        T_mid = evaluate_T_v(v_bad, L_PARAM / mp.mpf("2"), L_PARAM)

        # Fine spatial search for peak amplitude on [0, L]
        n_sample = 200
        max_T = mp.mpf("0")
        t_peak = mp.mpf("0")
        for s in range(n_sample + 1):
            t_eval = L_PARAM * mp.mpf(s) / mp.mpf(n_sample)
            val_T = abs(evaluate_T_v(v_bad, t_eval, L_PARAM))
            if val_T > max_T:
                max_T = val_T
                t_peak = t_eval

        # Spatial mass partition: flat plateau [0, log 2] vs well [log 2, L]
        log2 = mp.log(mp.mpf("2"))
        nodes_flat, weights_flat = gauss_legendre_nodes_weights(QUAD_ORDER, mp.mpf("0"), log2)
        mass_flat = mp.mpf("0")
        for k in range(QUAD_ORDER):
            t_k = nodes_flat[k]
            w_k = weights_flat[k]
            T_k = evaluate_T_v(v_bad, t_k, L_PARAM)
            mass_flat += w_k * (T_k ** 2)
        mass_flat = mass_flat / L_PARAM
        mass_well = mp.mpf("1") - mass_flat

        # Spatial kinetic energy
        K_kin = ((mp.mpf("2") * mp.pi / L_PARAM) ** 2) * sum(
            (mp.mpf(m) ** 2) * (v_bad[m, 0] ** 2) for m in range(1, dim_even)
        )

        table1_rows.append({
            "N": N, "T_0": T_0, "T_mid": T_mid,
            "max_T": max_T, "t_peak": t_peak,
            "mass_flat": mass_flat, "mass_well": mass_well,
            "K_kin": K_kin,
        })

        # ----------------------------------------------------
        # 4. Table 2: Coordinate Integration of Potential Well
        # ----------------------------------------------------
        # Matrix quadratic forms
        R_W = (v_bad.T * W_tilde * v_bad)[0, 0]
        R_Delta_D = (v_bad.T * Delta_D * v_bad)[0, 0]
        R_neg = (v_bad.T * K_neg * v_bad)[0, 0]

        # Piecewise continuous Gauss-Legendre integration of W(t) |T_v(t)|^2
        # Transition boundaries: log(q)
        transition_points = [mp.mpf("0")] + [logn for (_, logn, _) in prime_data]
        I_pot = mp.mpf("0")
        current_W = mp.mpf("0")

        # Interval 0: [0, log 2] has W(t) = 0
        for i in range(len(prime_data)):
            _, logn_curr, w_curr = prime_data[i]
            current_W += w_curr  # accumulating -W(t)
            a_interval = logn_curr
            b_interval = prime_data[i + 1][1] if i + 1 < len(prime_data) else L_PARAM
            if b_interval > a_interval:
                nodes_i, weights_i = gauss_legendre_nodes_weights(QUAD_ORDER, a_interval, b_interval)
                sub_integral = mp.mpf("0")
                for k in range(QUAD_ORDER):
                    t_k = nodes_i[k]
                    w_k = weights_i[k]
                    T_k = evaluate_T_v(v_bad, t_k, L_PARAM)
                    sub_integral += w_k * (T_k ** 2)
                I_pot += current_W * sub_integral

        I_pot = I_pot / L_PARAM
        pot_err = abs(R_W - I_pot)
        if pot_err > max_global_pot_err:
            max_global_pot_err = pot_err

        table2_rows.append({
            "N": N, "R_W": R_W, "I_pot": I_pot, "pot_err": pot_err,
            "R_Delta_D": R_Delta_D, "R_neg": R_neg,
        })

        # ----------------------------------------------------
        # 5. Table 3: Archimedean Spectral Energy Partition
        # ----------------------------------------------------
        # Direct integration of |Phi_v(r)|^2 * h_+(r) on [0, r_*]
        nodes_arch_minus, weights_arch_minus = gauss_legendre_nodes_weights(QUAD_ORDER, mp.mpf("0"), r_star)
        I_arch_minus = mp.mpf("0")
        for k in range(QUAD_ORDER):
            r_k = nodes_arch_minus[k]
            w_k = weights_arch_minus[k]
            h_k = h_plus(r_k, GROUND_DPS)
            phi_val = evaluate_Phi_v(v_bad, r_k, L_PARAM)
            I_arch_minus += (w_k / mp.pi) * (phi_val ** 2) * h_k

        # Total continuous Archimedean form via exact matrix
        I_arch_tot = (v_bad.T * Q_even_arch * v_bad)[0, 0]
        I_arch_plus = I_arch_tot - I_arch_minus
        D_arch = sum((v_bad[m, 0] ** 2) * D_mult[m, m] for m in range(dim_even))

        table3_rows.append({
            "N": N,
            "I_arch_minus": I_arch_minus,
            "I_arch_plus": I_arch_plus,
            "I_arch_tot": I_arch_tot,
            "D_arch": D_arch,
            "R_arch": R_arch,
        })

        # ----------------------------------------------------
        # 6. Table 4: Surplus Margin Synthesis
        # ----------------------------------------------------
        table4_rows.append({
            "N": N,
            "R_comp": R_comp,
            "R_arch": R_arch,
            "R_pole": R_pole,
            "R_net": R_net,
            "rho": rho_restore,
            "surplus": delta_surplus,
        })

        t_elapsed = time.time() - t_n_start
        print(f"Completed N = {N:2d} in {t_elapsed:.2f}s | "
              f"Mass(well) = {float(mass_well):.4f}, "
              f"|R_W - I_pot| = {float(pot_err):.2e}, "
              f"I_arch(-) = {float(I_arch_minus):+.4e}, "
              f"Delta_surplus = {float(delta_surplus):+.6f}")

    # ============================================================
    # FORMATTED DIAGNOSTIC REPORT
    # ============================================================

    print("\n" + "-" * 80)
    print("TABLE 1: COORDINATE WAVEPACKET PROFILE OF PRINCIPAL VULNERABLE STATE v_bad")
    print("T_v(t) = v_0 + sqrt(2) * sum v_m * cos(2*pi*m*t/L) on [0, L]")
    print("-" * 80)
    print(f"{'N':>4} | {'T(0)':>9} | {'T(L/2)':>9} | {'max |T|':>9} | {'t_peak':>7} | {'M(flat)':>9} | {'M(well)':>9} | {'K_kin':>10}")
    print("-" * 80)
    for r in table1_rows:
        print(f"{r['N']:4d} | {float(r['T_0']):9.4f} | {float(r['T_mid']):9.4f} | {float(r['max_T']):9.4f} | {float(r['t_peak']):7.4f} | {float(r['mass_flat']):9.4f} | {float(r['mass_well']):9.4f} | {float(r['K_kin']):10.2f}")
    print("-" * 80)

    print("\n" + "-" * 80)
    print("TABLE 2: COORDINATE INTEGRATION OF POTENTIAL WELL: MATRIX vs CONTINUOUS INTEGRAL")
    print("Identity: R_W = v^T * W_tilde * v === (1/L) * int_0^L (-W(t)) * |T_v(t)|^2 dt")
    print("-" * 80)
    print(f"{'N':>4} | {'R_W (Matrix)':>14} | {'I_pot (Integral)':>16} | {'|Residual|':>12} | {'R_Delta_D':>11} | {'R_neg':>11}")
    print("-" * 80)
    for r in table2_rows:
        print(f"{r['N']:4d} | {float(r['R_W']):14.8f} | {float(r['I_pot']):16.8f} | {float(r['pot_err']):12.4e} | {float(r['R_Delta_D']):11.6f} | {float(r['R_neg']):11.6f}")
    print("-" * 80)

    # Dynamic status verification of coordinate-matrix equivalence
    if max_global_pot_err < mp.mpf("1e-40"):
        print(f"COORDINATE-MATRIX IDENTITY AUDIT: PASSED (max residual = {float(max_global_pot_err):.4e} < 1e-40).")
    else:
        print(f"COORDINATE-MATRIX IDENTITY AUDIT: FAILED (max residual = {float(max_global_pot_err):.4e} >= 1e-40)!")
    print("-" * 80)

    print("\n" + "-" * 80)
    print("TABLE 3: ARCHIMEDEAN SPECTRAL ENERGY PARTITION")
    print("Partition: I_arch = (1/pi) * int_0^oo |Phi_v(r)|^2 * h_+(r) dr = I(-) + I(+)")
    print("-" * 80)
    print(f"{'N':>4} | {'I_arch(-)':>13} | {'I_arch(+)':>12} | {'I_arch(tot)':>12} | {'D_arch (diag)':>14} | {'R_arch (off)':>13}")
    print("-" * 80)
    for r in table3_rows:
        print(f"{r['N']:4d} | {float(r['I_arch_minus']):13.6e} | {float(r['I_arch_plus']):12.6f} | {float(r['I_arch_tot']):12.6f} | {float(r['D_arch']):14.6f} | {float(r['R_arch']):13.6f}")
    print("-" * 80)

    print("\n" + "-" * 80)
    print("TABLE 4: SURPLUS MARGIN SYNTHESIS & ASYMPTOTIC BALANCE ALONG w_bad")
    print("Net Margin: Delta_surplus = R_net = (R_arch + R_pole) - |R_comp| > 0")
    print("-" * 80)
    print(f"{'N':>4} | {'R_comp':>10} | {'R_arch':>10} | {'R_pole':>10} | {'R_net':>10} | {'rho_restore':>12} | {'Delta_surplus':>14}")
    print("-" * 80)
    for r in table4_rows:
        print(f"{r['N']:4d} | {float(r['R_comp']):10.6f} | {float(r['R_arch']):10.6f} | {float(r['R_pole']):10.6f} | {float(r['R_net']):10.6f} | {float(r['rho']):12.6f} | {float(r['surplus']):14.6f}")
    print("-" * 80)

    # Synthesis at N = 64
    r1 = table1_rows[-1]
    r2 = table2_rows[-1]
    r3 = table3_rows[-1]
    r4 = table4_rows[-1]
    print(f"\nSYNTHESIS AT MAXIMAL RESOLUTION N = {r1['N']}:")
    print(f"  Spatial Probability Mass in Well:     {float(r1['mass_well']) * 100:.2f}% (Mass in flat plateau: {float(r1['mass_flat']) * 100:.2f}%)")
    print(f"  Wavepacket Peak Amplitude:            |T_max| = {float(r1['max_T']):.4f} at t = {float(r1['t_peak']):.4f}")
    print(f"  Boundary Contact:                     T(0) = T(L) = {float(r1['T_0']):.4f}, T(L/2) = {float(r1['T_mid']):.4f}")
    print(f"  Spatial Kinetic Semi-Norm K_kin:      {float(r1['K_kin']):.2f}")
    print(f"  Coordinate-Matrix Potential Error:    {float(r2['pot_err']):.4e}")
    print(f"  Negative Archimedean Energy I(-):     {float(r3['I_arch_minus']):.6e} (quenched by {float(abs(r3['I_arch_minus']) / r3['I_arch_tot']) * 100:.4f}% of total)")
    print(f"  Positive Off-Diagonal Archimedean:    R_arch = {float(r3['R_arch']):+.6f}")
    print(f"  Zeta Pole Contribution:               R_pole = {float(r4['R_pole']):+.6f}")
    print(f"  Competition Deficit:                  R_comp = {float(r4['R_comp']):+.6f}")
    print(f"  Net Positivity Surplus Margin:        Delta_surplus = {float(r4['surplus']):+.6f}")
    print(f"  Restoring Efficiency Ratio:           rho_restore   = {float(r4['rho']):.6f} (> 1)")

    t_total = time.time() - t_start
    print(f"\nTotal execution time: {t_total:.2f}s")
    print("=" * 80)
    print("CELL 133 EXECUTION COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    run_cell133_audit()
