# ============================================================
# CELL 102 — DETERMINISTIC RESOLVENT COMPARISON & WASSERSTEIN AUDIT
# ============================================================
#
# Gate:       Gate 1 (Finite-N Spectral Mechanism)
#
# Target Proposition:
#
#   Theorems 1 & 2 in cell102.md:
#   For the Feshbach coupling block B_M, the resolvent discrepancy
#   Delta S_M = |S_M(E11) - S_M^iso(E11)| satisfies the exact
#   discrete summation-by-parts identity:
#
#       S_M - S_M^iso = a_bar * sum_{j=0}^{q-2} K_j (g_j - g_{j+1}),
#
#   where K_k = sum_{j=0}^k (r_j - 1) and g_j = 1/(mu_j - E11),
#   and is unconditionally bounded by the Kolmogorov-Smirnov distance:
#
#       Delta S_M <= B_KS(M) = ||B_M||_F^2 * D_KS(M) * (g_0 - g_{q-1}),
#
#   and by the Wasserstein-1 (optimal transport) distance:
#
#       Delta S_M <= B_W(M)  = (||B_M||_F^2 / delta_M^2) * W_1.
#
# Falsification Criterion:
#
#   If the exact identity fails beyond numerical precision (> 1e-60),
#   or if Delta S_M > B_KS(M) or Delta S_M > B_W(M) for any cutoff M,
#   the deterministic comparison theorems are refuted.
#
# Design:
#
#   Evaluates M in {24, 32, 48, 64} at N=192, c=13, T=600.
#   Vectorized projection V = B*U for high numerical efficiency.
#   Computes exact W_1 distance via piecewise integration of CDF deviation.
#   Audits bounding slack ratios sigma_KS = B_KS / Delta S_M and
#   sigma_W = B_W / Delta S_M.
#
# ============================================================

import time
import mpmath as mp

from cell import get_galerkin_matrix

# ============================================================
# PARAMETERS
# ============================================================

mp.mp.dps = 70

C = 13
T = 600
N = 192
DPS = 50

M_VALUES = [24, 32, 48, 64]

DIST_TOL = mp.mpf("1e-40")

# ============================================================
# MATRIX UTILITIES
# ============================================================

def full_to_canonical_matrix(Q, N):
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


def eigvals_sym(A):
    vals, _ = mp.eigsy(A)
    return [vals[i] for i in range(vals.rows)]


def eigsys_sym(A):
    vals, V = mp.eigsy(A)
    dim = A.rows
    idx = sorted(range(dim), key=lambda i: vals[i])

    evals = [vals[i] for i in idx]
    evecs = mp.matrix(dim, dim)

    for col_out, col_in in enumerate(idx):
        for row in range(dim):
            evecs[row, col_out] = V[row, col_in]

    return evals, evecs


def frobenius_norm_sq(A):
    total = mp.mpf(0)
    for i in range(A.rows):
        for j in range(A.cols):
            total += A[i, j] ** 2
    return total


def feshbach_blocks(H, M):
    d = H.rows
    p = M + 1
    q = d - p

    A = mp.matrix(p, p)
    B = mp.matrix(p, q)
    C_high = mp.matrix(q, q)

    for i in range(p):
        for j in range(p):
            A[i, j] = H[i, j]

    for i in range(p):
        for j in range(q):
            B[i, j] = H[i, p + j]

    for i in range(q):
        for j in range(q):
            C_high[i, j] = H[p + i, p + j]

    return A, B, C_high


# ============================================================
# LOAD MATRIX AND GROUND ENERGY
# ============================================================

print()
print("=" * 72)
print("CELL 102 — DETERMINISTIC RESOLVENT COMPARISON & WASSERSTEIN AUDIT")
print("=" * 72)
print()
print(f"C={C}, N={N}, T={T}, dps={DPS}, working_dps={mp.mp.dps}")
print(f"Sweep cutoffs M: {M_VALUES}")
print()

t0 = time.perf_counter()

Q, Q_meta = get_galerkin_matrix(
    c=C,
    N=N,
    T=T,
    dps=DPS,
    verbose=True,
)

print(f"Full Galerkin matrix loaded in {time.perf_counter() - t0:.3f} s")

t0 = time.perf_counter()
H = full_to_canonical_matrix(Q, N)
print(f"Canonical matrix converted (dim={H.rows}) in {time.perf_counter() - t0:.3f} s")

t0 = time.perf_counter()
full_eigs = eigvals_sym(H)
E11 = full_eigs[0]
print(f"Canonical ground energy E11 = {mp.nstr(E11, 18)}")
print(f"Full spectrum solve in {time.perf_counter() - t0:.3f} s")

# ============================================================
# COMPARISON THEOREM SWEEP
# ============================================================

audit_summary = []

for M in M_VALUES:
    print()
    print("=" * 72)
    print(f"AUDIT AT M = {M}  (dim P = {M+1}, dim Q = {N-M})")
    print("=" * 72)

    tM = time.perf_counter()
    A, B, C_high = feshbach_blocks(H, M)
    p = M + 1
    q = N - M

    B_frob_sq = frobenius_norm_sq(B)

    # High-sector eigensystem
    t_eig = time.perf_counter()
    evals_high, evecs_high = eigsys_sym(C_high)
    print(f"  high-sector eigensolve ({q}x{q}) in {time.perf_counter() - t_eig:.3f} s")

    mu_0 = evals_high[0]
    mu_max = evals_high[-1]
    delta_M = mu_0 - E11

    # Vectorized projection V = B * U (p x q)
    t_v = time.perf_counter()
    V = B * evecs_high
    print(f"  coupling projection V = B*U in {time.perf_counter() - t_v:.3f} s")

    coupling_a = []
    for j in range(q):
        a_j = mp.mpf(0)
        for i in range(p):
            a_j += V[i, j] ** 2
        coupling_a.append(a_j)

    total_mass = mp.fsum(coupling_a)
    a_bar = total_mass / q

    # Normalized coupling ratios
    r_vals = [a / a_bar for a in coupling_a]

    # Resolvent weights g_j = 1/(mu_j - E11)
    g_vals = [mp.mpf(1) / (evals_high[j] - E11) for j in range(q)]
    g_0 = g_vals[0]
    g_last = g_vals[-1]
    g_mean = mp.fsum(g_vals) / q

    # Actual trace S_M and isotropic trace S_M^iso
    S_M = mp.fsum(coupling_a[j] * g_vals[j] for j in range(q))
    S_M_iso = a_bar * mp.fsum(g_vals)
    delta_S = abs(S_M - S_M_iso)
    signed_delta_S = S_M - S_M_iso
    eps_M = signed_delta_S / S_M_iso

    # --------------------------------------------------------
    # Theorem 1: Exact Discrete Summation-by-Parts Identity
    # K_k = sum_{j=0}^k (r_j - 1)
    # diff = g_j - g_{j+1}
    # --------------------------------------------------------
    K_vals = []
    cumul_excess = mp.mpf(0)
    for j in range(q):
        cumul_excess += (r_vals[j] - mp.mpf(1))
        K_vals.append(cumul_excess)

    # Check vanishing boundary term: K_{q-1} == 0
    boundary_K = abs(K_vals[-1])

    # Summation-by-parts sum: a_bar * sum_{j=0}^{q-2} K_j (g_j - g_{j+1})
    sbp_sum = mp.mpf(0)
    for j in range(q - 1):
        diff_g = g_vals[j] - g_vals[j + 1]
        sbp_sum += K_vals[j] * diff_g
    sbp_total = a_bar * sbp_sum

    sbp_residual = abs(sbp_total - signed_delta_S)

    # --------------------------------------------------------
    # Theorem 2: Kolmogorov-Smirnov Discrepancy Bound
    # B_KS = ||B||_F^2 * D_KS * (g_0 - g_{q-1})
    # --------------------------------------------------------
    max_abs_K = max(abs(K) for K in K_vals)
    D_KS = max_abs_K / q

    B_KS = total_mass * D_KS * (g_0 - g_last)
    slack_KS = B_KS / delta_S if delta_S > 0 else mp.mpf(0)

    # Geometric ratio C_geom = (g_0 - g_{q-1}) / <g>_unif
    C_geom = (g_0 - g_last) / g_mean
    rel_bound_KS = C_geom * D_KS

    # --------------------------------------------------------
    # Theorem 3: Wasserstein-1 Optimal Transport Distance & Bound
    # W_1 = int_{mu_0}^{mu_max} |F_coup(x) - F_unif(x)| dx
    # Between mu_j and mu_{j+1}, F_coup is constant = (1/q) * sum_{k=0}^j r_k
    # and F_unif is constant = (j+1)/q.
    # Deviation |F_coup - F_unif| = |K_j| / q on interval [mu_j, mu_{j+1}].
    # --------------------------------------------------------
    W_1 = mp.mpf(0)
    for j in range(q - 1):
        interval_width = evals_high[j + 1] - evals_high[j]
        dev_j = abs(K_vals[j]) / q
        W_1 += dev_j * interval_width

    B_W = (total_mass / (delta_M ** 2)) * W_1
    slack_W = B_W / delta_S if delta_S > 0 else mp.mpf(0)

    # --------------------------------------------------------
    # Output section diagnostics
    # --------------------------------------------------------
    print()
    print("  THEOREM 1: EXACT SUMMATION-BY-PARTS IDENTITY")
    print(f"    actual S_M - S_M^iso        = {mp.nstr(signed_delta_S, 12)}")
    print(f"    SbP formula total           = {mp.nstr(sbp_total, 12)}")
    print(f"    identity absolute residual  = {mp.nstr(sbp_residual, 6)}")
    print(f"    boundary term |K_{{q-1}}|       = {mp.nstr(boundary_K, 6)}  (exact zero: {boundary_K < mp.mpf('1e-60')})")

    print()
    print("  THEOREM 2: KOLMOGOROV-SMIRNOV UPPER BOUND")
    print(f"    Kolmogorov-Smirnov D_KS     = {mp.nstr(D_KS, 8)}")
    print(f"    spectral inverse-gap diff   = {mp.nstr(g_0 - g_last, 8)}")
    print(f"    actual |S_M - S_M^iso|      = {mp.nstr(delta_S, 10)}")
    print(f"    KS bound B_KS               = {mp.nstr(B_KS, 10)}")
    print(f"    bound holds (delta <= B_KS) = {delta_S <= B_KS}")
    print(f"    bounding slack B_KS / delta = {mp.nstr(slack_KS, 6)}")
    print(f"    relative excess eps_M       = {mp.nstr(eps_M, 6)}")
    print(f"    relative bound C_geom*D_KS  = {mp.nstr(rel_bound_KS, 6)}")
    print(f"    geometric factor C_geom     = {mp.nstr(C_geom, 6)}")

    print()
    print("  THEOREM 3: WASSERSTEIN-1 OPTIMAL TRANSPORT BOUND")
    print(f"    Wasserstein-1 distance W_1  = {mp.nstr(W_1, 8)}")
    print(f"    spectral distance delta_M   = {mp.nstr(delta_M, 8)}")
    print(f"    Wasserstein bound B_W       = {mp.nstr(B_W, 10)}")
    print(f"    bound holds (delta <= B_W)  = {delta_S <= B_W}")
    print(f"    bounding slack B_W / delta  = {mp.nstr(slack_W, 6)}")

    record = {
        "M": M,
        "p": p,
        "q": q,
        "delta_M": delta_M,
        "B_frob_sq": total_mass,
        "S_M": S_M,
        "S_M_iso": S_M_iso,
        "delta_S": delta_S,
        "eps_M": eps_M,
        "D_KS": D_KS,
        "B_KS": B_KS,
        "slack_KS": slack_KS,
        "C_geom": C_geom,
        "W_1": W_1,
        "B_W": B_W,
        "slack_W": slack_W,
        "sbp_res": sbp_residual,
        "elapsed": time.perf_counter() - tM,
    }
    audit_summary.append(record)

    print(f"  M={M} completed in {record['elapsed']:.2f} s")

# ============================================================
# MULTI-CUTOFF SYNTHESIS TABLE
# ============================================================

print()
print("=" * 88)
print("CELL 102 — SYNTHESIS: DETERMINISTIC DISCREPANCY BOUNDS & SLACK RATIOS")
print("=" * 88)
print()

header = (
    f"{'M':>4} {'q':>4} {'delta_S':>11} {'B_KS':>11} {'slack_KS':>10} "
    f"{'W_1':>10} {'B_W':>11} {'slack_W':>10} {'C_geom':>9} {'eps_M':>9}"
)
print(header)
print("-" * len(header))

for r in audit_summary:
    print(
        f"{r['M']:4d} {r['q']:4d} "
        f"{mp.nstr(r['delta_S'], 5):>11} "
        f"{mp.nstr(r['B_KS'], 5):>11} "
        f"{mp.nstr(r['slack_KS'], 4):>10} "
        f"{mp.nstr(r['W_1'], 4):>10} "
        f"{mp.nstr(r['B_W'], 5):>11} "
        f"{mp.nstr(r['slack_W'], 4):>10} "
        f"{mp.nstr(r['C_geom'], 4):>9} "
        f"{mp.nstr(r['eps_M'], 4):>9}"
    )

print()
print("=" * 80)
print("CELL 102 EXECUTION COMPLETE")
print("=" * 80)
