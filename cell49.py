"""
CELL 49 — PHASE II: COMPLETE GALERKIN SPECTRUM, BOUND-STATE COUNTING,
MULTI-c GAP UNIVERSALITY, AND SPECTRAL ZETA FUNCTION

In Cell 48, we discovered:
1. All lowest 8 eigenvalues of Q(c=13, N) are strictly positive with alternating parity.
2. The eigenfunctions obey an exact Sturm–Liouville nodal ladder (state k has k interior nodes).
3. Universal Dirichlet boundary suppression holds across all bound states (|T(0)| -> 0).
4. The Fourier amplitudes Phi_k(r) for k in {0, 1, 2, 3} vanish identically at the non-trivial
   Riemann zeros gamma_1, ..., gamma_5 to within machine precision (~ 10^-20).
5. A large fundamental spectral gap ratio E_1 / E_0 ~ 1313.36 ~ c^{2.805} isolates the ground state.

Cell 49 advances Phase II along five major frontiers:
1. Complete Spectrum & Bound-State vs Continuum Classification (c = 13, N = 20):
   Diagonalize the full (2N+1 = 41)-dimensional Galerkin operator across N in {8, 12, 16, 20}.
   For each state k in {0, ..., 2N}, compute the scaling exponent:
       alpha_k = -log(E_k(20) / E_k(16)) / (4 * log c).
   Classify every eigenmode into:
       - Bound states (alpha_k >= 0.5): decaying exponentially into the Dirichlet continuum;
       - Transitional states (0.1 <= alpha_k < 0.5);
       - Scattering continuum (alpha_k < 0.1): stationary O(1) positive continuum modes.
   Determine the exact bound-state capacity N_bound(c=13, N=20).

2. Multi-c Spectral Gap Universality across Prime Cutoffs c in {5, 7, 11, 13, 17}:
   Diagonalize Q(c, N=16) and Q(c, N=20) for all 5 prime cutoffs.
   Measure the fundamental gap Delta E(c) = E_1(c) - E_0(c) and the gap ratios:
       R_1(c) = E_1(c) / E_0(c),   R_2(c) = E_2(c) / E_1(c),   R_3(c) = E_3(c) / E_2(c).
   Compute the logarithmic exponents mu_j(c) = log(R_j) / log(c).
   Test whether R_1(c) ~ c^mu is a universal invariant of the Connes–CvS Galerkin operator.

3. Higher Bound-State Transmission Resonances with Riemann Zeros (k = 4, 5, 6, 7):
   Test whether the exact transmission zeros Phi_k(gamma_j) = 0 persist across higher bound
   states k in {4, 5, 6, 7} for gamma_1, ..., gamma_5.
   Construct the full transmission matrix T_{j, k} = |Phi_k(gamma_j)|^2 across the bound spectrum.

4. Discrete Spectral Zeta Function and Resolvent Traces:
   Compute the punctured resolvent trace:
       G'(s) = sum_{k=1}^{2N} 1 / (s + E_k)
   at s in {0, 10^-20, 10^-10, 1.0}.
   Evaluate the discrete spectral zeta function:
       zeta_Q(sigma) = sum_{k=1}^{2N} E_k^{-sigma}
   for sigma in {0.1, 0.25, 0.5, 0.75, 1.0} across N in {12, 16, 20}.

5. Semiclassical Density of States rho(E) and the Weyl Law:
   Measure the cumulative eigenvalue counting function N(E) = #{k : E_k <= E} across 12 decades
   from E = 10^-35 to E = 10^1.
   Extract the effective density of states and fit the Weyl growth exponent N(E) ~ C * E^delta.
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
N_PRIMARY_LIST = [8, 12, 16, 20]
T_PARAM = 400
DPS_RUN = 50

CUTOFFS_MULTI = [5, 7, 11, 13, 17]
N_MULTI_LIST = [16, 20]

# First 5 non-trivial Riemann zeros on the critical line
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

    # Orthogonal projection into parity blocks
    Q_even = V_even.T * Q * V_even
    Q_odd = V_odd.T * Q * V_odd

    # Diagonalize symmetric blocks
    eigs_ev, vecs_ev = mp.eigsy(Q_even)
    eigs_od, vecs_od = mp.eigsy(Q_odd)

    # Sort even sector
    idx_ev = sorted(range(len(eigs_ev)), key=lambda i: eigs_ev[i])
    even_states = []
    for i in idx_ev:
        lam = eigs_ev[i]
        v_coord = [vecs_ev[row, i] for row in range(N + 1)]
        nrm = mp.sqrt(sum(x ** 2 for x in v_coord))
        v_coord = [x / nrm for x in v_coord]
        if v_coord[0] < 0:
            v_coord = [-x for x in v_coord]
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
        w_coord = [vecs_od[row, i] for row in range(N)]
        nrm = mp.sqrt(sum(x ** 2 for x in w_coord))
        w_coord = [x / nrm for x in w_coord]
        if w_coord[0] < 0:
            w_coord = [-x for x in w_coord]
        v_full = mp.matrix(2 * N + 1, 1)
        for row in range(N):
            col_vec = V_odd[:, row]
            for r in range(2 * N + 1):
                v_full[r, 0] += w_coord[row] * col_vec[r, 0]
        odd_states.append((lam, "odd", w_coord, v_full))

    # Merge and sort full spectrum
    all_states = sorted(even_states + odd_states, key=lambda s: s[0])
    return even_states, odd_states, all_states


# ---------------------------------------------------------------------------
# Fourier Amplitude Phi(r) Evaluation
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
# MAIN INVESTIGATION
# ---------------------------------------------------------------------------

def main():
    print("=" * 80)
    print("CELL 49 — PHASE II: COMPLETE SPECTRUM, GAP UNIVERSALITY & SPECTRAL ZETA")
    print("=" * 80)
    print(f"Working precision dps = {DPS_RUN}")
    print()

    # -----------------------------------------------------------------------
    # Part 1: Complete Spectrum & Bound-State Classification (c = 13)
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("1. COMPLETE SPECTRUM & BOUND-STATE CLASSIFICATION (c = 13, N = 20)")
    print("=" * 80)

    # Store full spectra for N in [8, 12, 16, 20]
    spectra_c13 = {}
    states_c13 = {}

    for N in N_PRIMARY_LIST:
        print(f"Constructing and diagonalizing Q(c=13, N={N})...")
        Q = build_galerkin_matrix(C_PRIMARY, N=N, T=T_PARAM, dps=DPS_RUN)
        ev, od, all_st = diagonalize_sectors(Q, N)
        spectra_c13[N] = [s[0] for s in all_st]
        states_c13[N] = all_st

    print()
    print("--- Spectrum at N = 20 (Total dimension 2N+1 = 41) ---")
    print(f"{'Index':>5}  {'Parity':>6}  {'E(N=20)':>22}  {'E(N=16)':>22}  {'Alpha (Slope)':>15}  {'Classification':<16}")
    print("-" * 92)

    eigs_20 = spectra_c13[20]
    eigs_16 = spectra_c13[16]
    states_20 = states_c13[20]

    bound_count = 0
    transition_count = 0
    continuum_count = 0

    log_c = mp.log(C_PRIMARY)

    for k in range(len(eigs_20)):
        lam20 = eigs_20[k]
        parity = states_20[k][1]

        # Match with N=16 if index exists
        if k < len(eigs_16):
            lam16 = eigs_16[k]
            if lam20 > 0 and lam16 > 0:
                ratio = lam16 / lam20
                if ratio > 0:
                    alpha = mp.log(ratio) / (4 * log_c)
                else:
                    alpha = mp.mpf(0)
            else:
                alpha = mp.mpf(0)
        else:
            lam16 = mp.mpf(0)
            alpha = mp.mpf(0)

        # Classification based on geometric decay rate
        if alpha >= mp.mpf("0.5"):
            classification = "Bound State"
            bound_count += 1
        elif alpha >= mp.mpf("0.1"):
            classification = "Transitional"
            transition_count += 1
        else:
            classification = "Continuum"
            continuum_count += 1

        print(f"{k:>5}  {parity:>6}  {mp.nstr(lam20, 10):>22}  {mp.nstr(lam16, 10):>22}  {mp.nstr(alpha, 5):>15}  {classification:<16}")

    print("-" * 92)
    print(f"Total bound states (alpha >= 0.5):    {bound_count}")
    print(f"Total transitional states:           {transition_count}")
    print(f"Total scattering continuum states:   {continuum_count}")
    print(f"Total dimension:                     {len(eigs_20)}")
    print()

    # -----------------------------------------------------------------------
    # Part 2: Multi-c Spectral Gap Universality across Prime Cutoffs
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("2. MULTI-c SPECTRAL GAP UNIVERSALITY ACROSS PRIME CUTOFFS c in {5, 7, 11, 13, 17}")
    print("=" * 80)
    print("Evaluating lowest 4 eigenvalues at N = 16 and N = 20 across prime cutoffs...")
    print()

    multi_c_data = {}

    for c in CUTOFFS_MULTI:
        print(f"--- Cutoff c = {c} (L = {mp.nstr(mp.log(c), 6)}) ---")

        # N = 16
        if c == C_PRIMARY and 16 in spectra_c13:
            all_16 = states_c13[16]
        else:
            Q16 = build_galerkin_matrix(c, N=16, T=T_PARAM, dps=DPS_RUN)
            _, _, all_16 = diagonalize_sectors(Q16, 16)

        # N = 20
        if c == C_PRIMARY and 20 in spectra_c13:
            all_20 = states_c13[20]
        else:
            Q20 = build_galerkin_matrix(c, N=20, T=T_PARAM, dps=DPS_RUN)
            _, _, all_20 = diagonalize_sectors(Q20, 20)

        e0_20 = all_20[0][0]
        e1_20 = all_20[1][0]
        e2_20 = all_20[2][0]
        e3_20 = all_20[3][0]

        delta_E = e1_20 - e0_20
        r1 = e1_20 / e0_20
        r2 = e2_20 / e1_20
        r3 = e3_20 / e2_20

        mu1 = mp.log(r1) / mp.log(c)
        mu2 = mp.log(r2) / mp.log(c)
        mu3 = mp.log(r3) / mp.log(c)

        multi_c_data[c] = {
            "E_0": e0_20,
            "E_1": e1_20,
            "E_2": e2_20,
            "E_3": e3_20,
            "Delta_E": delta_E,
            "R_1": r1,
            "R_2": r2,
            "R_3": r3,
            "mu_1": mu1,
            "mu_2": mu2,
            "mu_3": mu3,
        }

        print(f"  E_0 (even):   {mp.nstr(e0_20, 10)}")
        print(f"  E_1 (odd):    {mp.nstr(e1_20, 10)}")
        print(f"  E_2 (even):   {mp.nstr(e2_20, 10)}")
        print(f"  E_3 (odd):    {mp.nstr(e3_20, 10)}")
        print(f"  Delta E:      {mp.nstr(delta_E, 10)}")
        print(f"  Ratio E_1/E_0: {mp.nstr(r1, 6)}  -->  c^{mp.nstr(mu1, 5)}")
        print(f"  Ratio E_2/E_1: {mp.nstr(r2, 6)}  -->  c^{mp.nstr(mu2, 5)}")
        print(f"  Ratio E_3/E_2: {mp.nstr(r3, 6)}  -->  c^{mp.nstr(mu3, 5)}")
        print()

    print("--- Multi-c Summary Table (N = 20) ---")
    print(f"{'c':>4}  {'E_0':>18}  {'E_1':>18}  {'Delta E':>18}  {'R_1 = E_1/E_0':>14}  {'mu_1 (log_c R_1)':>18}")
    print("-" * 84)
    for c in CUTOFFS_MULTI:
        d = multi_c_data[c]
        print(f"{c:>4}  {mp.nstr(d['E_0'], 8):>18}  {mp.nstr(d['E_1'], 8):>18}  {mp.nstr(d['Delta_E'], 8):>18}  {mp.nstr(d['R_1'], 6):>14}  {mp.nstr(d['mu_1'], 6):>18}")
    print("-" * 84)
    print()

    # -----------------------------------------------------------------------
    # Part 3: Higher Bound-State Transmission Resonances with Riemann Zeros
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("3. HIGHER BOUND-STATE TRANSMISSION RESONANCES WITH RIEMANN ZEROS (N = 20, c = 13)")
    print("=" * 80)
    print("Testing extinction |Phi_k(gamma_j)|^2 across states k in {0, ..., 7} and zeros gamma_1..gamma_5...")
    print()

    header_zeros = "  ".join([f"gamma_{j+1} (~{mp.nstr(RIEMANN_ZEROS[j], 4)})" for j in range(len(RIEMANN_ZEROS))])
    print(f"{'State':>6}  {'Parity':>6}  {'Eigenvalue E':>15}  {header_zeros}")
    print("-" * 110)

    for k in range(8):
        st = states_20[k]
        lam = st[0]
        parity = st[1]

        vals_str = []
        for j in range(len(RIEMANN_ZEROS)):
            gz = RIEMANN_ZEROS[j]
            phi_val = Phi_eval_state(st, gz, L_PRIMARY)
            phi_sq = abs(phi_val) ** 2
            vals_str.append(f"{mp.nstr(phi_sq, 3):>18}")

        print(f"E_{k:<3}  {parity:>6}  {mp.nstr(lam, 7):>15}  {'  '.join(vals_str)}")

    print("-" * 110)
    print()

    # -----------------------------------------------------------------------
    # Part 4: Discrete Spectral Zeta Function and Resolvent Traces
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("4. DISCRETE SPECTRAL ZETA FUNCTION & PUNCTURED RESOLVENT TRACES")
    print("=" * 80)

    for N in [12, 16, 20]:
        eigs = spectra_c13[N]
        eigs_excited = eigs[1:]  # Omit ground state E_0

        print(f"--- Dimension N = {N} (Excited subspace dimension = {len(eigs_excited)}) ---")

        # Punctured trace G'(s) = sum_{k=1}^{2N} 1 / (s + E_k)
        for s_val in [mp.mpf(0), mp.mpf("1e-20"), mp.mpf("1e-10"), mp.mpf("1.0")]:
            if s_val == 0:
                tr_val = sum(1 / lam for lam in eigs_excited)
                s_label = "s = 0"
            else:
                tr_val = sum(1 / (s_val + lam) for lam in eigs_excited)
                s_label = f"s = {mp.nstr(s_val, 2)}"
            print(f"  Punctured Resolvent Trace G'({s_label}): {mp.nstr(tr_val, 10)}")

        print()

        # Spectral zeta function zeta_Q(sigma) = sum_{k=1}^{2N} E_k^{-sigma}
        for sigma in [mp.mpf("0.1"), mp.mpf("0.25"), mp.mpf("0.5"), mp.mpf("0.75"), mp.mpf("1.0")]:
            z_val = sum(lam ** (-sigma) for lam in eigs_excited)
            print(f"  Spectral Zeta zeta_Q(sigma={sigma}): {mp.nstr(z_val, 10)}")

        print()

    # -----------------------------------------------------------------------
    # Part 5: Semiclassical Cumulative Distribution N(E) and the Weyl Law
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("5. SEMICLASSICAL CUMULATIVE SPECTRUM N(E) AND THE WEYL LAW (N = 20, c = 13)")
    print("=" * 80)

    energy_decades = [
        mp.mpf("1e-38"), mp.mpf("1e-35"), mp.mpf("1e-30"), mp.mpf("1e-25"),
        mp.mpf("1e-20"), mp.mpf("1e-15"), mp.mpf("1e-10"), mp.mpf("1e-5"),
        mp.mpf("1e-2"), mp.mpf("0.1"), mp.mpf("1.0"), mp.mpf("10.0")
    ]

    print(f"{'Energy Threshold E':>20}  {'N(E) = #{E_k <= E}':>20}  {'Fraction of Spectrum':>22}")
    print("-" * 68)

    total_dim = len(eigs_20)
    for eth in energy_decades:
        count = sum(1 for lam in eigs_20 if lam <= eth)
        frac = mp.mpf(count) / total_dim
        print(f"{mp.nstr(eth, 3):>20}  {count:>20}  {mp.nstr(frac * 100, 4):>20}%")

    print("-" * 68)
    print()

    print("=" * 80)
    print("END OF CELL 49")
    print("=" * 80)


if __name__ == "__main__":
    main()
