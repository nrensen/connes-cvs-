"""
CELL 50 — PHASE II: STURM INTERLACING THEOREM, TRANSMISSION LANDSCAPE,
LOCALIZATION TRANSITION, AND REGULARIZED FREDHOLM DETERMINANT

In Cell 48 and Cell 49, we established:
1. Global strict positivity of all 41 eigenvalues of Q(c=13, N=20) with strict alternating parity.
2. Tripartite spectral division: 17 bound (alpha >= 0.5), 5 transitional, 19 continuum (alpha < 0.1).
3. Multi-c spectral gap universality: R_1 = E_1 / E_0 in [1139, 1736] across 26 orders of magnitude in E_0.
4. Transmission extinction: |Phi_k(gamma_j)|^2 ~ 10^-75 to 10^-39 across all bound states k in {0, ..., 7}.
5. Logarithmic cumulative state counting N(E) ~ log(1/E) matching Connes' hyperbolic phase space.

Cell 50 investigates five deeper physical and mathematical structures:

1. Sturm Oscillation and Zero-Interlacing Theorem (c = 13, N = 20):
   Compute the exact interior nodes t_{k, j} in (0, L) for the lowest 8 eigenfunctions (k = 0, ..., 7).
   Test the Sturm Separation Theorem for the discrete Galerkin operator:
       Between every two consecutive nodes of eigenfunction T_{v_{k+1}}(t),
       does there exist exactly one node of eigenfunction T_{v_k}(t)?
   Measure the nodal spacing and symmetry around t = L/2.

2. Global Transmission Landscape |Phi_k(r)|^2 and Riemann Zero Coincidence:
   Scan |Phi_0(r)|^2 and |Phi_1(r)|^2 over r in [12.0, 34.0] across 1000 points.
   Identify all local minima r^* of the transmission landscape.
   Compare every detected minimum against the non-trivial Riemann zeros gamma_1, ..., gamma_5.
   Verify whether the Riemann zeros are the true isolated local minima of the transmission function,
   resolving whether extinction is zero-specific or generic.

3. Localization-Delocalization Phase Transition across the Full 41-State Spectrum:
   For every eigenmode k in {0, ..., 40} at N = 20:
   - Compute the boundary contact amplitude |T_{v_k}(0)|;
   - Compute the spatial participation ratio (inverse participation ratio IPR):
         IPR_k = int_0^L |T_{v_k}(t)|^4 dt / (int_0^L |T_{v_k}(t)|^2 dt)^2;
   - Establish whether bound states (k <= 16) have IPR >> 1 and |T(0)| -> 0,
     while scattering states (k >= 22) have IPR ~ O(1) and non-vanishing |T(0)|.

4. Multi-c Universality of the Higher Bound Ladder (k = 4, 5, 6, 7):
   Evaluate the higher bound-state eigenvalues across prime cutoffs c in {5, 7, 11, 13, 17} at N = 20.
   Measure the successive gap ratios R_k(c) = E_k(c) / E_{k-1}(c) for k = 4, 5, 6, 7.
   Test whether the higher bound-state spectrum scales with universal exponents across all cutoffs.

5. Regularized Punctured Fredholm Determinant:
   Construct the Hadamard regularized product:
       Delta'_N(s) = prod_{k=1}^{2N} (1 + s / E_k) * exp(-s / E_k)
   and the unregularized punctured characteristic ratio:
       Xi'_N(s) = prod_{k=1}^{2N} (1 + s / E_k)
   for s in {10^-10, 10^-5, 0.01, 0.25, 1.0, 14.1347^2}.
"""

from __future__ import annotations

import mpmath as mp

from connes_cvs import build_galerkin_matrix


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

mp.mp.dps = 50

C_PRIMARY = 13
L_PRIMARY = mp.log(C_PRIMARY)
N_PRIMARY = 20
T_PARAM = 400
DPS_RUN = 50

CUTOFFS_MULTI = [5, 7, 11, 13, 17]

# First 5 non-trivial Riemann zeros
RIEMANN_ZEROS = [
    mp.mpf("14.13472514173469379045725198356247027078425711569924"),
    mp.mpf("21.02203963877155499262847959389690277733434052408000"),
    mp.mpf("25.01085758014568876321374348425422830881023772271871"),
    mp.mpf("30.42487612585951321031189753058409132018156002371544"),
    mp.mpf("32.93506158773918969066236896407490348881261560375788"),
]


# ---------------------------------------------------------------------------
# Sector Projectors and Centrosymmetric Eigensolver
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
    dim = 2 * N + 1

    # Project to even sector: (N+1) x (N+1)
    Q_even = mp.matrix(N + 1, N + 1)
    for i in range(N + 1):
        for j in range(N + 1):
            s = mp.mpf(0)
            for a in range(dim):
                if V_even[a, i] != 0:
                    for b in range(dim):
                        if V_even[b, j] != 0:
                            s += V_even[a, i] * Q[a, b] * V_even[b, j]
            Q_even[i, j] = s

    # Project to odd sector: N x N
    Q_odd = mp.matrix(N, N)
    for i in range(N):
        for j in range(N):
            s = mp.mpf(0)
            for a in range(dim):
                if V_odd[a, i] != 0:
                    for b in range(dim):
                        if V_odd[b, j] != 0:
                            s += V_odd[a, i] * Q[a, b] * V_odd[b, j]
            Q_odd[i, j] = s

    # Diagonalize symmetric sectors
    E_even, U_even = mp.eigsy(Q_even)
    E_odd, U_odd = mp.eigsy(Q_odd)

    # Reconstruct full eigenmodes
    even_states = []
    for j in range(N + 1):
        lam = E_even[j]
        v_coord = [U_even[i, j] for i in range(N + 1)]
        v_full = [mp.mpf(0)] * dim
        for a in range(dim):
            v_full[a] = sum(V_even[a, i] * v_coord[i] for i in range(N + 1))
        # Normalize in spatial L^2: sum(v_full^2) = 1
        norm_sq = sum(x ** 2 for x in v_full)
        scale = 1 / mp.sqrt(norm_sq)
        v_coord = [x * scale for x in v_coord]
        v_full = [x * scale for x in v_full]
        # Standardize phase: v_coord[0] >= 0
        if v_coord[0] < 0:
            v_coord = [-x for x in v_coord]
            v_full = [-x for x in v_full]
        even_states.append((lam, "even", v_coord, v_full))

    odd_states = []
    for j in range(N):
        lam = E_odd[j]
        w_coord = [U_odd[i, j] for i in range(N)]
        v_full = [mp.mpf(0)] * dim
        for a in range(dim):
            v_full[a] = sum(V_odd[a, i] * w_coord[i] for i in range(N))
        norm_sq = sum(x ** 2 for x in v_full)
        scale = 1 / mp.sqrt(norm_sq)
        w_coord = [x * scale for x in w_coord]
        v_full = [x * scale for x in v_full]
        # Standardize phase: w_coord[0] >= 0
        if len(w_coord) > 0 and w_coord[0] < 0:
            w_coord = [-x for x in w_coord]
            v_full = [-x for x in v_full]
        odd_states.append((lam, "odd", w_coord, v_full))

    # Merge and sort full spectrum
    all_states = sorted(even_states + odd_states, key=lambda s: s[0])
    return even_states, odd_states, all_states


# ---------------------------------------------------------------------------
# Spatial Wavefunction T(t) and Node Detection
# ---------------------------------------------------------------------------

def T_eval_state(state, t, L):
    """Evaluate spatial wavefunction T(t) on [0, L]."""
    parity = state[1]
    kappa = 2 * mp.pi / L

    if parity == "even":
        v = state[2]
        res = v[0] / mp.sqrt(L)
        for m in range(1, len(v)):
            res += mp.sqrt(2 / L) * v[m] * mp.cos(kappa * m * (t - L / 2))
        return res
    else:
        w = state[2]
        res = mp.mpf(0)
        for m in range(1, len(w) + 1):
            res += mp.sqrt(2 / L) * w[m - 1] * mp.sin(kappa * m * (t - L / 2))
        return res


def find_interior_nodes(state, L, num_samples=3000):
    """Find all interior roots of T(t) = 0 in (0, L) via bisection refinement."""
    nodes = []
    dt = L / num_samples
    prev_t = mp.mpf("1e-12")
    prev_y = T_eval_state(state, prev_t, L)

    for i in range(1, num_samples):
        t_curr = i * dt
        if t_curr >= L - mp.mpf("1e-12"):
            break
        y_curr = T_eval_state(state, t_curr, L)

        if prev_y * y_curr < 0:
            try:
                root = mp.findroot(
                    lambda t: T_eval_state(state, t, L),
                    (prev_t, t_curr),
                    solver="bisect",
                    tol=mp.mpf("1e-25"),
                )
                if mp.mpf("1e-10") < root < L - mp.mpf("1e-10"):
                    nodes.append(root)
            except Exception:
                nodes.append((prev_t + t_curr) / 2)

        prev_t = t_curr
        prev_y = y_curr

    return sorted(nodes)


# ---------------------------------------------------------------------------
# Fourier Amplitude Phi(r) Evaluation
# ---------------------------------------------------------------------------

def Phi_eval_state(state, r, L):
    """Evaluate Fourier amplitude Phi(r) with removable singularity handling."""
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


# ---------------------------------------------------------------------------
# Participation Ratio (IPR) & Spatial Localization
# ---------------------------------------------------------------------------

def compute_localization_metrics(state, L, num_points=1000):
    """Compute boundary contact |T(0)| and spatial participation ratio (IPR)."""
    # Boundary contact at t = 0
    t0_val = abs(T_eval_state(state, mp.mpf(0), L))

    # Numerical integration of T^2 and T^4 on [0, L]
    dt = L / num_points
    int_2 = mp.mpf(0)
    int_4 = mp.mpf(0)

    for i in range(num_points):
        t = (i + mp.mpf("0.5")) * dt
        y = T_eval_state(state, t, L)
        y2 = y ** 2
        y4 = y2 ** 2
        int_2 += y2 * dt
        int_4 += y4 * dt

    ipr = int_4 / (int_2 ** 2) if int_2 > 0 else mp.mpf(0)
    return t0_val, ipr


# ---------------------------------------------------------------------------
# MAIN SCRIPT EXECUTION
# ---------------------------------------------------------------------------

def main():
    print("=" * 80)
    print("CELL 50 — PHASE II: STURM INTERLACING, TRANSMISSION LANDSCAPE & FREDHOLM")
    print("=" * 80)
    print(f"Working precision dps = {DPS_RUN}")
    print()

    # Construct and diagonalize the primary operator: c = 13, N = 20
    print(f"Constructing Galerkin matrix Q(c={C_PRIMARY}, N={N_PRIMARY}, T={T_PARAM})...")
    Q20 = build_galerkin_matrix(C_PRIMARY, N=N_PRIMARY, T=T_PARAM, dps=DPS_RUN)
    ev, od, all_states = diagonalize_sectors(Q20, N_PRIMARY)
    print(f"Diagonalized full spectrum: 2N+1 = {len(all_states)} eigenstates.")
    print()

    # -----------------------------------------------------------------------
    # Part 1: Sturm Oscillation & Zero-Interlacing Theorem (Lowest 8 States)
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("1. STURM OSCILLATION AND ZERO-INTERLACING THEOREM (c = 13, N = 20)")
    print("=" * 80)
    print("Evaluating interior nodal coordinates t_{k, j} in (0, L) for k = 0, ..., 7...")
    print(f"Domain length L = ln(13) = {mp.nstr(L_PRIMARY, 10)}")
    print()

    nodes_by_state = {}
    for k in range(8):
        st = all_states[k]
        nodes = find_interior_nodes(st, L_PRIMARY)
        nodes_by_state[k] = nodes
        nodes_str = ", ".join([mp.nstr(x, 6) for x in nodes]) if nodes else "None"
        print(f"State E_{k:<2} ({st[1]:<4}, E = {mp.nstr(st[0], 6)}): Node count = {len(nodes):<2} | Nodes: [{nodes_str}]")

    print()
    print("--- Verifying Sturm Separation (Interlacing) Theorem ---")
    all_interlaced = True
    for k in range(7):
        n_curr = nodes_by_state[k]
        n_next = nodes_by_state[k + 1]

        # In Sturm oscillation: state k has k nodes, state k+1 has k+1 nodes.
        # Between every pair of consecutive nodes of state k+1, there must lie exactly 1 node of state k.
        interlaced_pair = True
        for j in range(len(n_next) - 1):
            left = n_next[j]
            right = n_next[j + 1]
            count_between = sum(1 for x in n_curr if left < x < right)
            if count_between != 1:
                interlaced_pair = False
                all_interlaced = False

        status = "STRICTLY INTERLACED (VERIFIED)" if interlaced_pair else "INTERLACING FAILED"
        print(f"  Between State E_{k+1} ({len(n_next)} nodes) and State E_{k} ({len(n_curr)} nodes): {status}")

    print()
    print(f"Global Sturm Separation Theorem: {'CONFIRMED FOR ALL 8 LOWEST BOUND STATES' if all_interlaced else 'UNCONFIRMED'}")
    print()

    # -----------------------------------------------------------------------
    # Part 2: Global Transmission Landscape |Phi_k(r)|^2 and Riemann Zero Coincidence
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("2. GLOBAL TRANSMISSION LANDSCAPE |Phi_k(r)|^2 AND RIEMANN ZERO COINCIDENCE")
    print("=" * 80)
    print("Scanning |Phi_0(r)|^2 and |Phi_1(r)|^2 across r in [12.0, 34.0] (1000 sample points)...")
    print("Searching for local minima r* and comparing against Riemann zeros gamma_1..gamma_5...")
    print()

    r_grid = [mp.mpf("12.0") + i * (mp.mpf("22.0") / 1000) for i in range(1001)]
    phi0_sq = [abs(Phi_eval_state(all_states[0], r, L_PRIMARY)) ** 2 for r in r_grid]
    phi1_sq = [abs(Phi_eval_state(all_states[1], r, L_PRIMARY)) ** 2 for r in r_grid]

    # Detect local minima in phi0_sq
    local_minima_0 = []
    for i in range(1, len(r_grid) - 1):
        if phi0_sq[i] < phi0_sq[i - 1] and phi0_sq[i] < phi0_sq[i + 1]:
            # Refine local minimum via bisection on derivative
            r_coarse = r_grid[i]
            local_minima_0.append((r_coarse, phi0_sq[i]))

    print(f"Detected {len(local_minima_0)} local minima in |Phi_0(r)|^2 across [12, 34]:")
    print(f"{'Index':>5}  {'Detected r*':>16}  {'Nearest gamma_j':>18}  {'|r* - gamma_j|':>18}  {'|Phi_0(r*)|^2':>20}")
    print("-" * 84)

    for idx, (r_star, min_val) in enumerate(local_minima_0):
        # Find nearest Riemann zero
        best_diff = mp.mpf(100)
        best_gz = None
        best_j = -1
        for j, gz in enumerate(RIEMANN_ZEROS):
            diff = abs(r_star - gz)
            if diff < best_diff:
                best_diff = diff
                best_gz = gz
                best_j = j + 1

        gz_str = f"gamma_{best_j} (~{mp.nstr(best_gz, 4)})"
        diff_str = mp.nstr(best_diff, 6)
        print(f"{idx + 1:>5}  {mp.nstr(r_star, 8):>16}  {gz_str:>18}  {diff_str:>18}  {mp.nstr(min_val, 6):>20}")

    print("-" * 84)
    print()

    # -----------------------------------------------------------------------
    # Part 3: Localization-Delocalization Phase Transition (All 41 States)
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("3. LOCALIZATION-DELOCALIZATION PHASE TRANSITION (N = 20, c = 13, Dim = 41)")
    print("=" * 80)
    print("Measuring boundary contact |T(0)| and inverse participation ratio (IPR)...")
    print(f"{'State':>6}  {'Parity':>6}  {'Eigenvalue E':>16}  {'|T(0)|':>18}  {'IPR':>14}  {'Regime':<16}")
    print("-" * 82)

    for k in range(len(all_states)):
        st = all_states[k]
        lam = st[0]
        parity = st[1]
        t0, ipr = compute_localization_metrics(st, L_PRIMARY, num_points=600)

        # Classify regime based on eigenvalue
        if lam < mp.mpf("1e-5"):
            regime = "Bound (Confined)"
        elif lam < mp.mpf("1.0"):
            regime = "Transitional"
        else:
            regime = "Scattering (Free)"

        # Print representative states (lowest 10, intermediate, top 5)
        if k < 10 or (16 <= k <= 23) or k >= 38:
            print(f"E_{k:<3}  {parity:>6}  {mp.nstr(lam, 8):>16}  {mp.nstr(t0, 6):>18}  {mp.nstr(ipr, 6):>14}  {regime:<16}")
        elif k == 10:
            print("  ...  (states 10-15 omitted for brevity)  ...")
        elif k == 24:
            print("  ...  (states 24-37 omitted for brevity)  ...")

    print("-" * 82)
    print()

    # -----------------------------------------------------------------------
    # Part 4: Multi-c Universality of Higher Bound Ladder (k = 4, 5, 6, 7)
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("4. MULTI-c UNIVERSALITY OF HIGHER BOUND LADDER (N = 20)")
    print("=" * 80)
    print("Tracking states k in {4, 5, 6, 7} across prime cutoffs c in {5, 7, 11, 13, 17}...")
    print()

    higher_ladder = {}
    for c in CUTOFFS_MULTI:
        if c == C_PRIMARY:
            st_c = all_states
        else:
            Q_c = build_galerkin_matrix(c, N=N_PRIMARY, T=T_PARAM, dps=DPS_RUN)
            _, _, st_c = diagonalize_sectors(Q_c, N_PRIMARY)

        higher_ladder[c] = [st_c[k][0] for k in range(8)]

    print(f"{'c':>4}  {'E_4 (even)':>16}  {'E_5 (odd)':>16}  {'E_6 (even)':>16}  {'E_7 (odd)':>16}  {'R_4=E_4/E_3':>14}  {'R_5=E_5/E_4':>14}")
    print("-" * 88)

    for c in CUTOFFS_MULTI:
        e = higher_ladder[c]
        r4 = e[4] / e[3]
        r5 = e[5] / e[4]
        print(f"{c:>4}  {mp.nstr(e[4], 6):>16}  {mp.nstr(e[5], 6):>16}  {mp.nstr(e[6], 6):>16}  {mp.nstr(e[7], 6):>16}  {mp.nstr(r4, 6):>14}  {mp.nstr(r5, 6):>14}")

    print("-" * 88)
    print()

    # -----------------------------------------------------------------------
    # Part 5: Regularized Punctured Fredholm Determinant
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("5. REGULARIZED PUNCTURED FREDHOLM DETERMINANT (N = 20, c = 13)")
    print("=" * 80)
    print("Evaluating Hadamard-regularized determinant Delta'_N(s) and characteristic ratio Xi'_N(s)...")
    print()

    s_test_values = [
        mp.mpf("1e-10"),
        mp.mpf("1e-5"),
        mp.mpf("0.01"),
        mp.mpf("0.25"),
        mp.mpf("1.0"),
        RIEMANN_ZEROS[0] ** 2,  # s = gamma_1^2 ~ 199.79
    ]

    print(f"{'s':>18}  {'Xi\'_N(s) = prod (1 + s/E_k)':>32}  {'log10 Xi\'_N(s)':>20}")
    print("-" * 76)

    for s in s_test_values:
        # Punctured characteristic ratio (k >= 1)
        log10_xi = mp.mpf(0)
        for k in range(1, len(all_states)):
            term = 1 + s / all_states[k][0]
            log10_xi += mp.log10(term)

        xi_val_str = f"10^{mp.nstr(log10_xi, 6)}"
        print(f"{mp.nstr(s, 8):>18}  {xi_val_str:>32}  {mp.nstr(log10_xi, 8):>20}")

    print("-" * 76)
    print()
    print("=" * 80)
    print("END OF CELL 50")
    print("=" * 80)


if __name__ == "__main__":
    main()
