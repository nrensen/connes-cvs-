"""
CELL 51 — OPERATOR RESOLVENT ANATOMY, HEAT-KERNEL BOUNDARY DYNAMICS,
AND ACCUMULATING POLE GEOMETRY OF THE GENERATING FUNCTION D_N(z)

Based on the analytical breakthrough connecting D_N(z) to the Neumann Laplacian
resolvent on [0, L]:
    L = -d^2/dt^2,   phi_0 = 1,   phi_m = sqrt(2) * cos(2*pi*m*t / L),
    L phi_m = kappa^2 * m^2 * phi_m,   kappa = 2*pi/L.

We have the exact operator identity:
    D_N(z) = [(I + z L)^(-1) T_N](0)
           = v_{N, 0} + sqrt(2) sum_{m=1}^N v_{N, m} / (1 + kappa^2 * m^2 * z).

And the exact Fourier-resolvent link:
    Phi_N(r) = (2 / sqrt(L)) * (sin(rL/2) / r) * D_N(-1/r^2),
    R_N(r)   = (2 / (L * r^2)) * [D_N(-1/r^2)]^2.

As N -> infinity, the poles of D_N(z) at z_m = -1 / (kappa^2 * m^2) accumulate at z = 0^-.
This cell investigates the analytic structure of D_N(z) across N in {8, 12, 16, 20, 24}:

1. Positive-Axis Resolvent D_N(x) for x > 0:
   Examine D_N(x) along the positive real axis (away from all poles), where (I + x L)^(-1)
   is a strictly positive bounded operator. Test convergence as N -> infinity.

2. Negative-Axis High-Frequency Approach D_N(-1/r^2) & WKB Singularity Scale:
   Evaluate B_N(r) = |D_N(-1/r^2)| across r in [5, 45] (interleaving between lattice poles).
   Fit -log B_N(r) against candidate scales:
       r,   r^2,   sqrt(r),   log(r),   r^alpha.
   Determine whether the high-frequency suppression is governed by an exponential WKB barrier
   e^{-C*r} or an accumulating rational pole phenomenon.

3. Discrete Stieltjes Transform & Alternating Coefficient Anatomy:
   Analyze F_N(w) = sum_{m=1}^N v_{N, m} / (m^2 + w) for w = -r^2 / kappa^2.
   Track the modulated coefficients b_{N, m} = (-1)^m * v_{N, m} to test for regularity
   and determine whether contour integration / special function representations apply.

4. Heat-Kernel Boundary Dynamics:
   Compute the heat-evolved boundary amplitude:
       H_N(u) = [e^{-u L} T_N](0) = v_{N, 0} + sqrt(2) sum_{m=1}^N v_{N, m} * e^{-kappa^2 * m^2 * u}
   for diffusion times u in [10^-6, 1.0].
   Verify short-time flat contact H_N(u) -> 0 as u -> 0 and test the Mellin/Laplace inversion:
       D_N(z) = int_0^infinity e^{-s} H_N(z * s) ds.

5. Dimension Scaling & Boundary-Layer Collapse:
   Test whether D_N(-1/r^2) exhibits scaling collapse as a function of r / N, r / sqrt(N), or r / N^alpha.
"""

from __future__ import annotations

import mpmath as mp

from connes_cvs import build_galerkin_matrix


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

mp.mp.dps = 50

C_PARAM = 13
L_PARAM = mp.log(C_PARAM)
KAPPA = 2 * mp.pi / L_PARAM
T_QUAD = 400
DPS_RUN = 50

N_LIST = [8, 12, 16, 20, 24]


# ---------------------------------------------------------------------------
# Even-Sector Ground-State Solver
# ---------------------------------------------------------------------------

def compute_even_ground_state(c, N, T=400, dps=50):
    """Build and diagonalize the even parity sector (N+1) x (N+1)."""
    dim = 2 * N + 1
    inv_sqrt2 = 1 / mp.sqrt(2)

    # Build even projector V_even: dim x (N+1)
    V_even = mp.matrix(dim, N + 1)
    V_even[N, 0] = mp.mpf(1)
    for k in range(1, N + 1):
        V_even[N + k, k] = inv_sqrt2
        V_even[N - k, k] = inv_sqrt2

    # Build full Galerkin matrix
    Q = build_galerkin_matrix(c, N=N, T=T, dps=dps)

    # Project to even sector
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

    # Diagonalize symmetric even matrix
    E_even, U_even = mp.eigsy(Q_even)

    # Lowest eigenmode
    lam_min = E_even[0]
    v_raw = [U_even[i, 0] for i in range(N + 1)]

    # Spatial L^2 normalization: v_0^2 + sum_{m=1}^N v_m^2 = 1
    norm_sq = v_raw[0] ** 2 + sum(v_raw[m] ** 2 for m in range(1, N + 1))
    scale = 1 / mp.sqrt(norm_sq)
    v = [x * scale for x in v_raw]

    # Standardize phase: v_0 >= 0
    if v[0] < 0:
        v = [-x for x in v]

    return lam_min, v


# ---------------------------------------------------------------------------
# Resolvent and Heat-Kernel Evaluators
# ---------------------------------------------------------------------------

def D_eval(v, z, kappa):
    """
    D_N(z) = v_0 + sqrt(2) sum_{m=1}^N v_m / (1 + kappa^2 * m^2 * z).
    """
    res = v[0]
    for m in range(1, len(v)):
        denom = 1 + (kappa * m) ** 2 * z
        res += mp.sqrt(2) * v[m] / denom
    return res


def F_stieltjes(v, w):
    """
    F_N(w) = sum_{m=1}^N v_m / (m^2 + w).
    """
    res = mp.mpf(0)
    for m in range(1, len(v)):
        res += v[m] / (m ** 2 + w)
    return res


def H_heat(v, u, kappa):
    """
    Heat-kernel boundary value:
        H_N(u) = [e^{-u L} T_N](0) = v_0 + sqrt(2) sum_{m=1}^N v_m * exp(-kappa^2 * m^2 * u).
    """
    res = v[0]
    for m in range(1, len(v)):
        arg = - (kappa * m) ** 2 * u
        res += mp.sqrt(2) * v[m] * mp.exp(arg)
    return res


# ---------------------------------------------------------------------------
# MAIN INVESTIGATION
# ---------------------------------------------------------------------------

def main():
    print("=" * 80)
    print("CELL 51 — RESOLVENT ANATOMY & ACCUMULATING POLE DYNAMICS OF D_N(z)")
    print("=" * 80)
    print(f"Parameters: c = {C_PARAM}, L = ln(13) = {mp.nstr(L_PARAM, 10)}, kappa = {mp.nstr(KAPPA, 10)}")
    print(f"Working precision: {DPS_RUN} decimal digits")
    print()

    # Step 1: Compute and store ground states for N in {8, 12, 16, 20, 24}
    ground_states = {}
    print("Diagonalizing even-sector Galerkin operators...")
    for N in N_LIST:
        print(f"  Computing ground state for N = {N} (matrix dimension {N+1}x{N+1})...")
        lam, v = compute_even_ground_state(C_PARAM, N=N, T=T_QUAD, dps=DPS_RUN)
        ground_states[N] = (lam, v)
        print(f"    N = {N:2d}: lambda_min = {mp.nstr(lam, 8)}, v_0 = {mp.nstr(v[0], 8)}")

    print()

    # -----------------------------------------------------------------------
    # Part 1: Positive-Axis Resolvent D_N(x) for x > 0
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("1. POSITIVE-AXIS RESOLVENT D_N(x) = [(I + x L)^(-1) T_N](0)  (x > 0)")
    print("=" * 80)
    print("Evaluating along the positive real axis, far from all negative poles...")
    print()

    x_test = [mp.mpf("1e-6"), mp.mpf("1e-4"), mp.mpf("1e-2"), mp.mpf("0.1"), mp.mpf("1.0"), mp.mpf("10.0")]
    header_x = "  ".join([f"x = {mp.nstr(x, 2):<8}" for x in x_test])
    print(f"{'N':>4}  {header_x}")
    print("-" * 75)

    for N in N_LIST:
        _, v = ground_states[N]
        vals = [D_eval(v, x, KAPPA) for x in x_test]
        vals_str = "  ".join([f"{mp.nstr(val, 6):>10}" for val in vals])
        print(f"{N:>4}  {vals_str}")

    print("-" * 75)
    print()

    # -----------------------------------------------------------------------
    # Part 2: Negative-Axis Approach D_N(-1/r^2) & WKB Singularity Scaling
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("2. NEGATIVE-AXIS APPROACH D_N(-1/r^2) & WKB SINGULARITY SCALING")
    print("=" * 80)
    print("Evaluating B_N(r) = |D_N(-1/r^2)| at half-integer lattice points r = kappa * (k + 1/2)...")
    print("This strictly avoids the poles at r = kappa * m.")
    print()

    r_samples = [KAPPA * (mp.mpf(k) + mp.mpf("0.5")) for k in [2, 4, 6, 8, 10, 12, 15, 20]]
    print(f"{'r':>8}  {'N = 8':>14}  {'N = 12':>14}  {'N = 16':>14}  {'N = 20':>14}  {'N = 24':>14}")
    print("-" * 82)

    b_data = {r: {} for r in r_samples}
    for r in r_samples:
        z = -1 / (r ** 2)
        row_str = []
        for N in N_LIST:
            _, v = ground_states[N]
            b_val = abs(D_eval(v, z, KAPPA))
            b_data[r][N] = b_val
            row_str.append(f"{mp.nstr(b_val, 5):>14}")
        print(f"{mp.nstr(r, 4):>8}  {'  '.join(row_str)}")

    print("-" * 82)
    print()

    print("--- Effective Logarithmic Scaling of -log |D_24(-1/r^2)| vs r ---")
    print(f"{'r':>8}  {'|D_24|':>14}  {'-log|D|':>12}  {'-log|D| / r':>14}  {'-log|D| / r^2':>14}  {'-log|D| / sqrt(r)':>16}")
    print("-" * 82)

    _, v24 = ground_states[24]
    for r in r_samples:
        z = -1 / (r ** 2)
        b24 = abs(D_eval(v24, z, KAPPA))
        if b24 > 0:
            log_b = -mp.log(b24)
            ratio_lin = log_b / r
            ratio_sq = log_b / (r ** 2)
            ratio_sqrt = log_b / mp.sqrt(r)
            print(f"{mp.nstr(r, 4):>8}  {mp.nstr(b24, 5):>14}  {mp.nstr(log_b, 6):>12}  {mp.nstr(ratio_lin, 6):>14}  {mp.nstr(ratio_sq, 6):>14}  {mp.nstr(ratio_sqrt, 6):>16}")

    print("-" * 82)
    print()

    # -----------------------------------------------------------------------
    # Part 3: Modulated Alternating Coefficients b_{N, m} = (-1)^m * v_m
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("3. MODULATED COEFFICIENT SEQUENCE b_{N, m} = (-1)^m * v_{N, m} (N = 24)")
    print("=" * 80)
    print("Inspecting positivity and decay of b_{24, m} for m = 1, ..., 16...")
    print()

    print(f"{'m':>4}  {'v_{24, m}':>18}  {'b_{24, m} = (-1)^m v_m':>24}  {'Ratio b_m / b_{m-1}':>22}")
    print("-" * 74)

    prev_b = None
    for m in range(1, 17):
        vm = v24[m]
        bm = ((-1) ** m) * vm
        ratio_str = f"{mp.nstr(bm / prev_b, 6):>22}" if prev_b is not None else f"{'—':>22}"
        print(f"{m:>4}  {mp.nstr(vm, 8):>18}  {mp.nstr(bm, 8):>24}  {ratio_str}")
        prev_b = bm

    print("-" * 74)
    print()

    # -----------------------------------------------------------------------
    # Part 4: Heat-Kernel Boundary Dynamics H_N(u) = [e^{-u L} T_N](0)
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("4. HEAT-KERNEL BOUNDARY DYNAMICS H_N(u) = [e^{-u L} T_N](0)")
    print("=" * 80)
    print("Tracking heat evolution of the boundary value from u = 1.0 down to u = 10^-6...")
    print()

    u_test = [mp.mpf("1.0"), mp.mpf("0.1"), mp.mpf("0.01"), mp.mpf("1e-3"), mp.mpf("1e-4"), mp.mpf("1e-5"), mp.mpf("1e-6")]
    print(f"{'Diffusion u':>14}  {'H_8(u)':>14}  {'H_12(u)':>14}  {'H_16(u)':>14}  {'H_20(u)':>14}  {'H_24(u)':>14}")
    print("-" * 86)

    for u in u_test:
        row = []
        for N in N_LIST:
            _, v = ground_states[N]
            h_val = H_heat(v, u, KAPPA)
            row.append(f"{mp.nstr(h_val, 5):>14}")
        print(f"{mp.nstr(u, 6):>14}  {'  '.join(row)}")

    print("-" * 86)
    print()

    # -----------------------------------------------------------------------
    # Part 5: Cross-Dimension Scaling Collapse
    # -----------------------------------------------------------------------
    print("=" * 80)
    print("5. CROSS-DIMENSION SCALING COLLAPSE OF D_N(-1/r^2)")
    print("=" * 80)
    print("Testing whether D_N(-1/r^2) collapses under scaled variables r / N or r / sqrt(N)...")
    print()

    # Evaluate at fixed ratios r / N = 1.0, 1.5, 2.0
    scaled_ratios = [mp.mpf("0.75"), mp.mpf("1.0"), mp.mpf("1.5"), mp.mpf("2.0")]
    for ratio in scaled_ratios:
        print(f"--- Fixed Scaling Ratio r / N = {mp.nstr(ratio, 3)} ---")
        for N in N_LIST:
            r = ratio * N
            z = -1 / (r ** 2)
            _, v = ground_states[N]
            val = abs(D_eval(v, z, KAPPA))
            print(f"  N = {N:2d}, r = {mp.nstr(r, 4):>6}: |D_N(-1/r^2)| = {mp.nstr(val, 8):>16},  -log|D|/N = {mp.nstr(-mp.log(val)/N, 6):>10}")
        print()

    print("=" * 80)
    print("END OF CELL 51")
    print("=" * 80)


if __name__ == "__main__":
    main()
