"""
CELL 48 — PHASE II: EXCITED STATES, NODAL ANATOMY, AND SPECTRAL RESONANCES WITH RIEMANN ZEROS

In Phase I (Cells 41–47), we established the complete theory of the Weil ground state v^{(0)}:
1. Universal eigenvalue scaling: lambda_min(N; c) ~ kappa * c^{-N} with universal kappa ~ 0.00238.
2. Semiclassical WKB barrier law: S_WKB ~ (pi N / 4) log c, matching suppression within 5.3%.
3. Super-polynomial resolvent decay and exact tri-partite zero-energy balance.

Cell 48 launches Phase II by investigating the excited spectrum of the Connes–van Suijlekom
Galerkin operator Q(c=13, N):
1. Full Galerkin Spectrum & Parity Sector Decomposition:
   Decompose Q into even (dimension N+1) and odd (dimension N) parity sectors.
   Track the combined lowest 10 eigenvalues E_0 <= E_1 <= ... <= E_9 across N in {8, 12, 16, 20}.
   Verify whether lambda_0^{even} -> 0 is an isolated zero mode, and evaluate the spectral gap
   Delta E = E_1 - E_0.

2. Spatial Profiles & Sturm–Liouville Nodal Ladder:
   Evaluate the spatial eigenfunctions T_k(t) on [0, L] for the lowest 6 states.
   Count interior zeros (nodes) in (0, L) to test the Sturm–Liouville nodal theorem:
   does state k possess exactly k interior nodes?

3. Boundary Confinement Across Excited States:
   Evaluate |T_{v_k}(0)| for excited even states across N to test whether all bound states
   develop Dirichlet boundary vanishing T_k(0) -> 0 (infinite-order confinement).

4. Fourier Resolvent Profiles & Resonance with Riemann Zeros:
   Evaluate the continuous Fourier amplitude Phi_k(r) for ground and excited states.
   Locate the real spectral zeros r^* in [0, 35] and test alignment with the low-lying
   non-trivial Riemann zeros:
       gamma_1 ~ 14.134725, gamma_2 ~ 21.022040, gamma_3 ~ 25.010858,
       gamma_4 ~ 30.424876, gamma_5 ~ 32.935062.

5. Tri-Partite Arithmetic Energy Decomposition for Excited States:
   Decompose Q(v_k) = Q_pole(v_k) + Q_prime(v_k) + Q_arch(v_k) for the lowest 4 states.
   Identify which arithmetic contribution drives the excited eigenvalues above zero.
"""

from __future__ import annotations

import mpmath as mp

from cell import (
    prime_power_terms,
)
from connes_cvs import build_galerkin_matrix


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

mp.mp.dps = 50

C_VAL = 13
L_VAL = mp.log(C_VAL)
N_LIST = [8, 12, 16, 20]
T_PARAM = 400
DPS_RUN = 50

# First 5 non-trivial Riemann zeros on the critical line
RIEMANN_ZEROS = [
    mp.mpf("14.13472514173469379045725198356247027078425711569924"),
    mp.mpf("21.02203963877155499262847959389690277733434052408000"),
    mp.mpf("25.01085758014568876321374348425422830881023772271871"),
    mp.mpf("30.42487612585951321031189753058409132018156002371544"),
    mp.mpf("32.93506158773918969066236896407490348881261560375788"),
]


# ---------------------------------------------------------------------------
# Sector Projectors and Eigensolver
# ---------------------------------------------------------------------------

def build_parity_projectors(N):
    dim = 2 * N + 1
    inv_sqrt2 = 1 / mp.sqrt(2)

    # Even projector: dim x (N + 1)
    V_even = mp.matrix(dim, N + 1)
    V_even[N, 0] = mp.mpf(1)
    for k in range(1, N + 1):
        V_even[N + k, k] = inv_sqrt2
        V_even[N - k, k] = inv_sqrt2

    # Odd projector: dim x N
    V_odd = mp.matrix(dim, N)
    for k in range(1, N + 1):
        col = k - 1
        V_odd[N + k, col] = inv_sqrt2
        V_odd[N - k, col] = -inv_sqrt2

    return V_even, V_odd


def diagonalize_sectors(Q, N):
    V_even, V_odd = build_parity_projectors(N)

    # Project
    Q_even = V_even.T * Q * V_even
    Q_odd = V_odd.T * Q * V_odd

    # Diagonalize
    eigs_ev, vecs_ev = mp.eigsy(Q_even)
    eigs_od, vecs_od = mp.eigsy(Q_odd)

    # Sort even sector
    idx_ev = sorted(range(len(eigs_ev)), key=lambda i: eigs_ev[i])
    even_states = []
    for i in idx_ev:
        lam = eigs_ev[i]
        # Coordinates in even basis: [v_0, v_1, ..., v_N]
        v_coord = [vecs_ev[row, i] for row in range(N + 1)]
        # Normalize
        nrm = mp.sqrt(sum(x ** 2 for x in v_coord))
        v_coord = [x / nrm for x in v_coord]
        # Standardize sign: make v_0 or central value positive
        if v_coord[0] < 0:
            v_coord = [-x for x in v_coord]
        # Lift to full 2N+1 basis
        v_full = mp.matrix(2 * N + 1, 1)
        for row in range(N + 1):
            col_vec = V_even[:, row]
            for r in range(2 * N + 1):
                v_full[r, 0] += v_coord[row] * col_vec[r, 0]
        even_states.append((lam, "even", v_coord, v_full))

    # Sort odd sector
    idx_od = sorted(range(len(eigs_od)), key=lambda i: eigs_od[i])
    odd_states = []
    for i in idx_od:
        lam = eigs_od[i]
        # Coordinates in odd basis: [w_1, ..., w_N]
        w_coord = [vecs_od[row, i] for row in range(N)]
        nrm = mp.sqrt(sum(x ** 2 for x in w_coord))
        w_coord = [x / nrm for x in w_coord]
        # Standardize sign: make first nonzero positive
        if w_coord[0] < 0:
            w_coord = [-x for x in w_coord]
        v_full = mp.matrix(2 * N + 1, 1)
        for row in range(N):
            col_vec = V_odd[:, row]
            for r in range(2 * N + 1):
                v_full[r, 0] += w_coord[row] * col_vec[r, 0]
        odd_states.append((lam, "odd", w_coord, v_full))

    # Merge and sort all states
    all_states = sorted(even_states + odd_states, key=lambda s: s[0])
    return even_states, odd_states, all_states


# ---------------------------------------------------------------------------
# Spatial Profiles T(t) and Nodal Counting
# ---------------------------------------------------------------------------

def T_eval_state(state, t, L):
    parity = state[1]
    kappa = 2 * mp.pi / L
    if parity == "even":
        v = state[2]
        val = v[0]
        for m in range(1, len(v)):
            val += mp.sqrt(2) * v[m] * mp.cos(kappa * m * t)
        return val
    else:
        w = state[2]
        val = mp.mpf(0)
        for m in range(1, len(w) + 1):
            val += mp.sqrt(2) * w[m - 1] * mp.sin(kappa * m * t)
        return val


def count_interior_nodes(state, L, num_grid=1000):
    """Count interior zeros of T(t) in (0, L) and locate them."""
    nodes = []
    dt = L / num_grid
    prev_t = dt * mp.mpf("0.5")
    prev_y = T_eval_state(state, prev_t, L)

    for i in range(1, num_grid):
        t_curr = dt * (i + mp.mpf("0.5"))
        if t_curr >= L:
            break
        y_curr = T_eval_state(state, t_curr, L)

        # Detect sign change
        if prev_y * y_curr < 0:
            # Refine root via bisection/findroot
            try:
                root = mp.findroot(
                    lambda t: T_eval_state(state, t, L),
                    (prev_t, t_curr),
                    solver="bisect",
                    tol=mp.mpf("1e-25"),
                )
                if mp.mpf(0) < root < L:
                    nodes.append(root)
            except Exception:
                nodes.append((prev_t + t_curr) / 2)

        prev_t = t_curr
        prev_y = y_curr

    return len(nodes), nodes


# ---------------------------------------------------------------------------
# Fourier Entire Amplitude Phi(r)
# ---------------------------------------------------------------------------

def Phi_eval_state(state, r, L):
    parity = state[1]
    kappa = 2 * mp.pi / L

    if parity == "even":
        v = state[2]
        if r == 0:
            return v[0] * mp.sqrt(L)
        sin_term = mp.sin(r * L / 2)
        sum_m = mp.mpf(0)
        for m in range(1, len(v)):
            am = kappa * m
            denom = r ** 2 - am ** 2
            if abs(denom) < mp.mpf("1e-20"):
                # Removable singularity at r = a_m
                term = mp.sqrt(2) * v[m] * ((-1) ** m * L / 4)
            else:
                term = mp.sqrt(2) * v[m] * r * sin_term / denom
            sum_m += term
        return (2 / mp.sqrt(L)) * (v[0] * sin_term / r + sum_m)

    else:
        w = state[2]
        if r == 0:
            return mp.mpf(0)
        sin_term = mp.sin(r * L / 2)
        sum_m = mp.mpf(0)
        for m in range(1, len(w) + 1):
            am = kappa * m
            denom = r ** 2 - am ** 2
            if abs(denom) < mp.mpf("1e-20"):
                term = mp.sqrt(2) * w[m - 1] * ((-1) ** m * L / 4)
            else:
                term = mp.sqrt(2) * w[m - 1] * am * sin_term / denom
            sum_m += term
        return (2 / mp.sqrt(L)) * sum_m


def find_spectral_roots(state, L, r_min=1.0, r_max=35.0, steps=1000):
    """Find real roots of Phi(r) = 0 in [r_min, r_max]."""
    roots = []
    dr = (r_max - r_min) / steps
    prev_r = mp.mpf(r_min)
    prev_val = Phi_eval_state(state, prev_r, L)

    for i in range(1, steps + 1):
        r_curr = mp.mpf(r_min) + i * dr
        val_curr = Phi_eval_state(state, r_curr, L)

        if prev_val * val_curr < 0:
            try:
                root = mp.findroot(
                    lambda r: Phi_eval_state(state, r, L),
                    (prev_r, r_curr),
                    solver="bisect",
                    tol=mp.mpf("1e-20"),
                )
                roots.append(root)
            except Exception:
                roots.append((prev_r + r_curr) / 2)

        prev_r = r_curr
        prev_val = val_curr

    return roots


# ---------------------------------------------------------------------------
# Prime and Pole Matrix Builders for Arithmetic Decomposition
# ---------------------------------------------------------------------------

def psi_prime_val(x, q, Lambda_q, L):
    a = 1 - mp.log(q) / L
    prefactor = -1 / mp.pi * Lambda_q / mp.sqrt(q)
    return prefactor * mp.sin(2 * mp.pi * x * a)


def psi_prime_deriv_val(x, q, Lambda_q, L):
    a = 1 - mp.log(q) / L
    prefactor = -1 / mp.pi * Lambda_q / mp.sqrt(q)
    return prefactor * 2 * mp.pi * a * mp.cos(2 * mp.pi * x * a)


def build_prime_matrix(N, c, L):
    size = 2 * N + 1
    Q_prime = mp.matrix(size, size)
    terms = prime_power_terms(c)

    for q, Lambda_q in terms:
        vals = [psi_prime_val(x, q, Lambda_q, L) for x in range(-N, N + 1)]
        ders = [psi_prime_deriv_val(x, q, Lambda_q, L) for x in range(-N, N + 1)]

        for i, m in enumerate(range(-N, N + 1)):
            for j, n in enumerate(range(-N, N + 1)):
                if m != n:
                    Q_prime[i, j] += (vals[i] - vals[j]) / mp.mpf(m - n)
                else:
                    Q_prime[i, j] += ders[i]
    return Q_prime


POLE_VALS_CACHE = {0: mp.mpf(0)}
POLE_DERS_CACHE = {}


def get_pole_val(x, L):
    if x not in POLE_VALS_CACHE:
        integrand = lambda y: 2 * mp.cosh(y / 2) * mp.sin(2 * mp.pi * abs(x) * (1 - y / L))
        v = (1 / mp.pi) * mp.quad(integrand, [0, L])
        POLE_VALS_CACHE[abs(x)] = v
        POLE_VALS_CACHE[-abs(x)] = -v
    return POLE_VALS_CACHE[x]


def get_pole_deriv(x, L):
    if x not in POLE_DERS_CACHE:
        integrand = lambda y: 2 * mp.cosh(y / 2) * (2 * mp.pi * (1 - y / L)) * mp.cos(2 * mp.pi * abs(x) * (1 - y / L))
        d = (1 / mp.pi) * mp.quad(integrand, [0, L])
        POLE_DERS_CACHE[abs(x)] = d
        POLE_DERS_CACHE[-abs(x)] = d
    return POLE_DERS_CACHE[x]


def build_pole_matrix(N, L):
    size = 2 * N + 1
    Q_pole = mp.matrix(size, size)

    for i, m in enumerate(range(-N, N + 1)):
        for j, n in enumerate(range(-N, N + 1)):
            if m != n:
                Q_pole[i, j] = (get_pole_val(m, L) - get_pole_val(n, L)) / mp.mpf(m - n)
            else:
                Q_pole[i, j] = get_pole_deriv(m, L)
    return Q_pole


# ---------------------------------------------------------------------------
# Main Analysis
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 80)
    print("CELL 48 — PHASE II: EXCITED STATES, NODAL ANATOMY, AND SPECTRAL RESONANCES")
    print("=" * 80)
    print(f"Cutoff c = {C_VAL}, L = {mp.nstr(L_VAL, 12)}")
    print(f"Dimensions N in {N_LIST}, dps = {mp.mp.dps}")

    # Build and diagonalize matrices across dimensions
    data_by_N = {}
    for N in N_LIST:
        print(f"\nConstructing and diagonalizing Galerkin matrix Q(c={C_VAL}, N={N})...")
        Q_mat = build_galerkin_matrix(c=C_VAL, N=N, T=T_PARAM, dps=DPS_RUN)
        ev_st, od_st, all_st = diagonalize_sectors(Q_mat, N)
        data_by_N[N] = {
            "Q": Q_mat,
            "even": ev_st,
            "odd": od_st,
            "all": all_st,
        }

    # =======================================================================
    # 1. Full Galerkin Spectrum & Level Spacing
    # =======================================================================
    print("\n" + "=" * 80)
    print("1. FULL GALERKIN SPECTRUM: LOWEST 8 EIGENVALUES ACROSS DIMENSIONS N")
    print("=" * 80)
    print(
        f"{'State':>5} "
        f"{'Parity':>7} "
        f"{'E(N=8)':>18} "
        f"{'E(N=12)':>18} "
        f"{'E(N=16)':>18} "
        f"{'E(N=20)':>18}"
    )
    print("-" * 80)

    num_display = 8
    for k in range(num_display):
        row_str = f"{'E_' + str(k):>5} "
        parity = data_by_N[20]["all"][k][1]
        row_str += f"{parity:>7} "
        for N in N_LIST:
            val = data_by_N[N]["all"][k][0]
            row_str += f"{mp.nstr(val, 10):>18} "
        print(row_str)

    # Spectral Gap
    E0_20 = data_by_N[20]["all"][0][0]
    E1_20 = data_by_N[20]["all"][1][0]
    gap_20 = E1_20 - E0_20
    print("-" * 80)
    print(f"Ground-state eigenvalue (N=20): E_0 = {mp.nstr(E0_20, 12)}")
    print(f"First excited state (N=20):     E_1 = {mp.nstr(E1_20, 12)}")
    print(f"Fundamental spectral gap:       Delta E = E_1 - E_0 = {mp.nstr(gap_20, 12)}")
    print(f"Ratio E_1 / E_0:                {mp.nstr(E1_20 / E0_20, 8)}")

    # =======================================================================
    # 2. Spatial Profiles & Sturm–Liouville Nodal Ladder
    # =======================================================================
    print("\n" + "=" * 80)
    print("2. SPATIAL PROFILES & STURM-LIOUVILLE NODAL LADDER (N = 20)")
    print("=" * 80)
    print(
        f"{'State':>5} "
        f"{'Parity':>7} "
        f"{'Eigenvalue E':>18} "
        f"{'Nodes':>6} "
        f"{'T(L/2)':>14} "
        f"{'|T(0)|':>14} "
        f"{'Interior Zeros in (0, L)':>20}"
    )
    print("-" * 90)

    for k in range(6):
        st = data_by_N[20]["all"][k]
        lam = st[0]
        parity = st[1]
        n_nodes, node_list = count_interior_nodes(st, L_VAL, num_grid=2000)
        t_mid = T_eval_state(st, L_VAL / 2, L_VAL)
        t_0 = T_eval_state(st, mp.mpf(0), L_VAL)

        node_str = ", ".join(f"{float(x):.3f}" for x in node_list[:4])
        if len(node_list) > 4:
            node_str += ", ..."

        print(
            f"{'E_' + str(k):>5} "
            f"{parity:>7} "
            f"{mp.nstr(lam, 10):>18} "
            f"{n_nodes:6d} "
            f"{mp.nstr(t_mid, 6):>14} "
            f"{mp.nstr(abs(t_0), 6):>14} "
            f"{node_str:>20}"
        )

    # =======================================================================
    # 3. Boundary Confinement Across Excited Even States
    # =======================================================================
    print("\n" + "=" * 80)
    print("3. BOUNDARY EXTINCTION ACROSS EVEN BOUND STATES: |T_{v_k}(0)|")
    print("=" * 80)
    print("(Note: All odd states satisfy T(0) = T(L) = 0 identically by parity)")
    print(
        f"{'State':>10} "
        f"{'N = 8':>16} "
        f"{'N = 12':>16} "
        f"{'N = 16':>16} "
        f"{'N = 20':>16}"
    )
    print("-" * 80)

    for k in range(4):
        state_name = f"Even #{k}"
        vals_by_N = []
        for N in N_LIST:
            st = data_by_N[N]["even"][k]
            t0 = abs(T_eval_state(st, mp.mpf(0), L_VAL))
            vals_by_N.append(t0)
        print(
            f"{state_name:>10} "
            f"{mp.nstr(vals_by_N[0], 8):>16} "
            f"{mp.nstr(vals_by_N[1], 8):>16} "
            f"{mp.nstr(vals_by_N[2], 8):>16} "
            f"{mp.nstr(vals_by_N[3], 8):>16}"
        )

    # =======================================================================
    # 4. Fourier Resolvent Profiles & Resonance with Riemann Zeros
    # =======================================================================
    print("\n" + "=" * 80)
    print("4. FOURIER AMPLITUDES Phi_k(r) AND SPECTRAL RESONANCE WITH RIEMANN ZEROS")
    print("=" * 80)
    print("Known Riemann zeros: gamma_1 ~ 14.1347, gamma_2 ~ 21.0220, gamma_3 ~ 25.0109")

    for k in range(4):
        st = data_by_N[20]["all"][k]
        parity = st[1]
        print(f"\n--- State E_{k} ({parity}, lambda = {mp.nstr(st[0], 8)}) ---")

        # Find roots of Phi(r) = 0
        roots = find_spectral_roots(st, L_VAL, r_min=1.0, r_max=35.0, steps=2000)
        root_str = ", ".join(f"{float(r):.4f}" for r in roots) if roots else "None detected"
        print(f"Roots of Phi_{k}(r) in [1, 35]: {root_str}")

        # Evaluate at the first 5 Riemann zeros
        print(
            f"{'Zero':>6} "
            f"{'gamma_k':>14} "
            f"{'Phi(gamma_k)':>18} "
            f"{'|Phi|^2':>16} "
            f"{'Nearest Root r*':>18} "
            f"{'|gamma_k - r*|':>16}"
        )
        print("-" * 88)
        for idx, g in enumerate(RIEMANN_ZEROS, 1):
            val = Phi_eval_state(st, g, L_VAL)
            val_sq = val ** 2
            nearest = None
            diff = None
            if roots:
                nearest = min(roots, key=lambda r: abs(r - g))
                diff = abs(nearest - g)
                diff_str = mp.nstr(diff, 6)
                near_str = mp.nstr(nearest, 8)
            else:
                near_str = "---"
                diff_str = "---"

            print(
                f"{'g_' + str(idx):>6} "
                f"{mp.nstr(g, 8):>14} "
                f"{mp.nstr(val, 8):>18} "
                f"{mp.nstr(val_sq, 8):>16} "
                f"{near_str:>18} "
                f"{diff_str:>16}"
            )

    # =======================================================================
    # 5. Tri-Partite Arithmetic Energy Decomposition for Excited States
    # =======================================================================
    print("\n" + "=" * 80)
    print("5. TRI-PARTITE ARITHMETIC ENERGY DECOMPOSITION FOR LOWEST 4 STATES (N = 20)")
    print("=" * 80)

    Q_20 = data_by_N[20]["Q"]
    Q_pr_20 = build_prime_matrix(20, C_VAL, L_VAL)
    Q_po_20 = build_pole_matrix(20, L_VAL)
    Q_ar_20 = Q_20 - Q_po_20 - Q_pr_20

    print(
        f"{'State':>5} "
        f"{'Parity':>7} "
        f"{'Q_pole':>16} "
        f"{'Q_prime':>16} "
        f"{'Q_arch':>16} "
        f"{'Q_total (Sum)':>18} "
        f"{'lambda_k':>18}"
    )
    print("-" * 88)

    for k in range(4):
        st = data_by_N[20]["all"][k]
        lam = st[0]
        parity = st[1]
        u_full = st[3]  # (2N+1) x 1

        pole_val = mp.fdot(u_full, Q_po_20 * u_full)
        prime_val = mp.fdot(u_full, Q_pr_20 * u_full)
        arch_val = mp.fdot(u_full, Q_ar_20 * u_full)
        total_val = pole_val + prime_val + arch_val

        print(
            f"{'E_' + str(k):>5} "
            f"{parity:>7} "
            f"{mp.nstr(pole_val, 8):>16} "
            f"{mp.nstr(prime_val, 8):>16} "
            f"{mp.nstr(arch_val, 8):>16} "
            f"{mp.nstr(total_val, 8):>18} "
            f"{mp.nstr(lam, 8):>18}"
        )

    print("\n" + "=" * 80)
    print("END OF CELL 48")
    print("=" * 80)
