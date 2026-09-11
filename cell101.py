# ============================================================
# CELL 101 — COUPLING ISOTROPIZATION PROFILE & DISCREPANCY SCALING
# ============================================================
#
# Gate:       Gate 1 (Finite-N Spectral Mechanism)
#
# Target Proposition:
#
#   Asymptotic Isotropization of the Coupling Spectral Measure:
#   For the Feshbach coupling block B_M and high-sector eigenvectors
#   C_M u_j = mu_j u_j (j = 0, ..., q_M - 1), the normalized modal
#   coupling ratios
#
#       r_j = ||B_M u_j||^2 / (||B_M||_F^2 / q_M),   mean(r) = 1,
#
#   converge toward the uniform spectral measure of C_M as M increases.
#   Consequently, the discrepancy metrics
#
#       D_KS(M) = max_j |F_coup(j) - (j+1)/q|,
#       D_TV(M) = (1 / 2q) sum_j |r_j - 1|,
#       CV(M)   = sqrt( (1/q) sum_j (r_j - 1)^2 ),
#       eps_M   = S_M(E11) / S_M^iso(E11) - 1,
#
#   decay with M, confirming that the matrix Feshbach problem
#   asymptotically reduces to a scalar spectral average over the
#   high-sector density of states.
#
# Falsification Criterion:
#
#   If the discrepancy metrics D_KS, D_TV, or CV fail to decrease
#   as M increases from 32 to 64, or if eps_M concentrates permanently
#   away from zero, the asymptotic spectral isotropization hypothesis
#   is refuted.
#
# Design:
#
#   Targeted evaluation at M in {32, 48, 64} (N=192, c=13, T=600).
#   Vectorized B * U evaluation for high numerical efficiency.
#   Omission of redundant direct Feshbach LU solves to keep execution
#   streamlined (~2-3 minutes per cutoff).
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

M_VALUES = [32, 48, 64]
N_DECILES = 10

DIST_TOL = mp.mpf("1e-40")

# ============================================================
# MATRIX UTILITIES
# ============================================================

def full_to_canonical_matrix(Q, N):
    """
    Convert the symmetric full-space Galerkin matrix Q into the
    canonical real-even orthonormal basis:
        e_0,
        (e_k + e_-k)/sqrt(2), k=1,...,N.
    """
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
    """
    P-sector = indices 0,...,M (dim p = M + 1)
    Q-sector = indices M+1,...,N (dim q = N - M).
    """
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
print("CELL 101 — COUPLING ISOTROPIZATION PROFILE & DISCREPANCY SCALING")
print("=" * 72)
print()
print(f"C={C}, N={N}, T={T}, dps={DPS}, working_dps={mp.mp.dps}")
print(f"Target cutoffs M: {M_VALUES}")
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
# TARGETED ISOTROPIZATION SWEEP
# ============================================================

sweep_summary = []

for M in M_VALUES:
    print()
    print("=" * 72)
    print(f"ANALYSIS AT M = {M}  (dim P = {M+1}, dim Q = {N-M})")
    print("=" * 72)

    tM = time.perf_counter()
    A, B, C_high = feshbach_blocks(H, M)
    p = M + 1
    q = N - M

    B_frob_sq = frobenius_norm_sq(B)

    # --------------------------------------------------------
    # High-sector eigensystem
    # --------------------------------------------------------
    t_eig = time.perf_counter()
    evals_high, evecs_high = eigsys_sym(C_high)
    print(f"  high-sector eigensolve ({q}x{q}) = {time.perf_counter() - t_eig:.3f} s")

    mu_min = evals_high[0]
    mu_max = evals_high[-1]
    spectral_span = mu_max - mu_min
    delta_M = mu_min - E11

    # --------------------------------------------------------
    # Vectorized coupling evaluation: V = B * U (p x q)
    # Column j of V is B * u_j
    # --------------------------------------------------------
    t_v = time.perf_counter()
    V = B * evecs_high
    print(f"  coupling projection V = B*U ({p}x{q}) = {time.perf_counter() - t_v:.3f} s")

    coupling_a = []
    for j in range(q):
        a_j = mp.mpf(0)
        for i in range(p):
            a_j += V[i, j] ** 2
        coupling_a.append(a_j)

    total_coupling_mass = mp.fsum(coupling_a)
    parseval_err = abs(total_coupling_mass - B_frob_sq)

    # Baseline isotropic coupling per mode
    isotropic_baseline = total_coupling_mass / q

    # Normalized coupling ratio r_j = a_j / (||B||_F^2 / q)
    r_vals = [a / isotropic_baseline for a in coupling_a]

    # --------------------------------------------------------
    # Resolvent measures & Isotropic ratio
    # --------------------------------------------------------
    inv_gaps = [mp.mpf(1) / (evals_high[j] - E11) for j in range(q)]
    trace_inv_gaps = mp.fsum(inv_gaps)

    # Isotropic resolvent weights w_j^iso = (1/(mu_j - E)) / sum_k 1/(mu_k - E)
    w_iso = [inv_gaps[j] / trace_inv_gaps for j in range(q)]

    # Actual trace bound S_M and isotropic trace bound S_M^iso
    S_M = mp.fsum(coupling_a[j] * inv_gaps[j] for j in range(q))
    S_M_iso = isotropic_baseline * trace_inv_gaps

    ratio_S = S_M / S_M_iso
    eps_M = ratio_S - mp.mpf(1)

    # Verify exact expectation identity: E_{w_iso}[r] == S_M / S_M_iso
    exp_r_under_w = mp.fsum(w_iso[j] * r_vals[j] for j in range(q))
    identity_check = abs(exp_r_under_w - ratio_S)

    # --------------------------------------------------------
    # Uniformity & Discrepancy Metrics
    # --------------------------------------------------------
    # 1. Total Variation distance from uniform distribution
    #    D_TV = (1 / 2q) * sum_j |r_j - 1|
    diff_abs = [abs(r - mp.mpf(1)) for r in r_vals]
    D_TV = mp.fsum(diff_abs) / (2 * q)

    # 2. Kolmogorov-Smirnov distance from uniform CDF
    #    F_coup(j) = (1/q) * sum_{k=0}^j r_k = sum_{k=0}^j a_k / total_mass
    #    F_unif(j) = (j + 1) / q
    cumul_r = mp.mpf(0)
    ks_max = mp.mpf(0)
    for j in range(q):
        cumul_r += r_vals[j]
        F_coup = cumul_r / q
        F_unif = mp.mpf(j + 1) / q
        dev = abs(F_coup - F_unif)
        if dev > ks_max:
            ks_max = dev
    D_KS = ks_max

    # 3. Variance & Coefficient of Variation
    #    Var(r) = (1/q) * sum_j (r_j - 1)^2
    var_r = mp.fsum((r - mp.mpf(1)) ** 2 for r in r_vals) / q
    cv_r = mp.sqrt(var_r)

    # 4. Extremal values
    r_min = min(r_vals)
    r_max = max(r_vals)
    r_0 = r_vals[0]

    # --------------------------------------------------------
    # Output section diagnostics
    # --------------------------------------------------------
    print()
    print("  BASIC PARAMETERS & INTEGRITY")
    print(f"    high spectrum span = [{mp.nstr(mu_min, 10)}, {mp.nstr(mu_max, 10)}]")
    print(f"    spectral distance delta_M = {mp.nstr(delta_M, 10)}")
    print(f"    total coupling mass ||B||_F^2 = {mp.nstr(B_frob_sq, 12)}")
    print(f"    Parseval absolute error = {mp.nstr(parseval_err, 6)}")
    print(f"    isotropic weight per mode = {mp.nstr(isotropic_baseline, 10)}")

    print()
    print("  ISOTROPIC TRACE RATIO & EXACT RESOLVENT IDENTITY")
    print(f"    S_M(E11)                  = {mp.nstr(S_M, 12)}")
    print(f"    S_M^iso(E11)              = {mp.nstr(S_M_iso, 12)}")
    print(f"    S_M / S_M^iso             = {mp.nstr(ratio_S, 10)}")
    print(f"    excess ratio eps_M        = {mp.nstr(eps_M, 10)} ({mp.nstr(100*eps_M, 4)}%)")
    print(f"    E_{{w_iso}}[r] check diff  = {mp.nstr(identity_check, 6)}")

    print()
    print("  SPECTRAL UNIFORMITY & DISCREPANCY METRICS")
    print(f"    Kolmogorov-Smirnov distance D_KS = {mp.nstr(D_KS, 8)}")
    print(f"    Total Variation distance D_TV    = {mp.nstr(D_TV, 8)}")
    print(f"    Coefficient of variation CV(r)   = {mp.nstr(cv_r, 8)}")
    print(f"    lowest mode ratio r_0            = {mp.nstr(r_0, 6)}  (mu_0 = {mp.nstr(mu_min, 6)})")
    print(f"    extremal range [r_min, r_max]    = [{mp.nstr(r_min, 6)}, {mp.nstr(r_max, 6)}]")

    # --------------------------------------------------------
    # Decile Profile Across Spectrum
    # --------------------------------------------------------
    print()
    print(f"  DECILE SPECTRAL PROFILE (q={q} modes into {N_DECILES} bins)")
    print(f"  {'bin':>4} {'modes':>12} {'mu range':>24} {'mean r_j':>12} {'mass frac':>12} {'iso target':>12}")

    bin_size = q // N_DECILES
    remainder = q - N_DECILES * bin_size
    start = 0

    decile_means = []
    for b in range(N_DECILES):
        size = bin_size + (1 if b < remainder else 0)
        end = start + size

        slice_r = r_vals[start:end]
        mean_r_b = mp.fsum(slice_r) / size
        mass_frac_b = mp.fsum(coupling_a[start:end]) / total_coupling_mass

        mu_start = evals_high[start]
        mu_end = evals_high[end - 1]
        mu_label = f"[{mp.nstr(mu_start, 5)}, {mp.nstr(mu_end, 5)}]"
        mode_label = f"{start:3d}--{end-1:3d}"

        decile_means.append(mean_r_b)

        print(
            f"  {b+1:4d} {mode_label:>12} {mu_label:>24} "
            f"{mp.nstr(mean_r_b, 6):>12} "
            f"{mp.nstr(100*mass_frac_b, 5) + '%':>12} "
            f"{'10.0%':>12}"
        )
        start = end

    record = {
        "M": M,
        "p": p,
        "q": q,
        "delta_M": delta_M,
        "B_frob_sq": B_frob_sq,
        "S_M": S_M,
        "S_M_iso": S_M_iso,
        "ratio_S": ratio_S,
        "eps_M": eps_M,
        "D_KS": D_KS,
        "D_TV": D_TV,
        "CV": cv_r,
        "r_0": r_0,
        "r_max": r_max,
        "r_min": r_min,
        "decile_means": decile_means,
        "elapsed": time.perf_counter() - tM,
    }
    sweep_summary.append(record)

    print(f"  cutoff M={M} completed in {record['elapsed']:.2f} s")

# ============================================================
# MULTI-M SYNTHESIS AND CONVERGENCE RATE ANALYSIS
# ============================================================

print()
print("=" * 80)
print("CELL 101 — MULTI-CUTOFF SYNTHESIS & ISOTROPIZATION RATE")
print("=" * 80)
print()

header = (
    f"{'M':>4} {'q':>4} {'delta_M':>10} {'S_M':>10} {'S_iso':>10} "
    f"{'ratio':>9} {'eps_M':>10} {'D_KS':>9} {'D_TV':>9} {'CV(r)':>9} {'r_0':>8}"
)
print(header)
print("-" * len(header))

for rec in sweep_summary:
    print(
        f"{rec['M']:4d} {rec['q']:4d} "
        f"{mp.nstr(rec['delta_M'], 6):>10} "
        f"{mp.nstr(rec['S_M'], 6):>10} "
        f"{mp.nstr(rec['S_M_iso'], 6):>10} "
        f"{mp.nstr(rec['ratio_S'], 6):>9} "
        f"{mp.nstr(rec['eps_M'], 5):>10} "
        f"{mp.nstr(rec['D_KS'], 5):>9} "
        f"{mp.nstr(rec['D_TV'], 5):>9} "
        f"{mp.nstr(rec['CV'], 5):>9} "
        f"{mp.nstr(rec['r_0'], 4):>8}"
    )

print()
print("EMPIRICAL SCALING EXPONENTS BETWEEN SUCCESSIVE CUTOFFS (rate gamma where metric ~ M^-gamma):")
print()

for i in range(len(sweep_summary) - 1):
    r1 = sweep_summary[i]
    r2 = sweep_summary[i + 1]
    M1, M2 = r1["M"], r2["M"]
    log_M_ratio = mp.log(mp.mpf(M2) / mp.mpf(M1))

    gamma_eps = -mp.log(r2["eps_M"] / r1["eps_M"]) / log_M_ratio
    gamma_ks = -mp.log(r2["D_KS"] / r1["D_KS"]) / log_M_ratio
    gamma_tv = -mp.log(r2["D_TV"] / r1["D_TV"]) / log_M_ratio
    gamma_cv = -mp.log(r2["CV"] / r1["CV"]) / log_M_ratio

    print(f"  Step M={M1} -> M={M2} (factor {M2/M1:.2f}):")
    print(f"    excess ratio eps_M rate gamma  = {mp.nstr(gamma_eps, 4)}")
    print(f"    Kolmogorov-Smirnov D_KS rate   = {mp.nstr(gamma_ks, 4)}")
    print(f"    Total Variation D_TV rate      = {mp.nstr(gamma_tv, 4)}")
    print(f"    Coefficient of variation rate  = {mp.nstr(gamma_cv, 4)}")
    print()

# Overall M=32 to M=64 rate
if len(sweep_summary) >= 3:
    r_first = sweep_summary[0]
    r_last = sweep_summary[-1]
    log_M_span = mp.log(mp.mpf(r_last["M"]) / mp.mpf(r_first["M"]))

    gamma_eps_full = -mp.log(r_last["eps_M"] / r_first["eps_M"]) / log_M_span
    gamma_ks_full = -mp.log(r_last["D_KS"] / r_first["D_KS"]) / log_M_span
    gamma_tv_full = -mp.log(r_last["D_TV"] / r_first["D_TV"]) / log_M_span
    gamma_cv_full = -mp.log(r_last["CV"] / r_first["CV"]) / log_M_span

    print(f"  Overall scaling M={r_first['M']} -> M={r_last['M']}:")
    print(f"    excess ratio eps_M rate gamma  = {mp.nstr(gamma_eps_full, 4)}")
    print(f"    Kolmogorov-Smirnov D_KS rate   = {mp.nstr(gamma_ks_full, 4)}")
    print(f"    Total Variation D_TV rate      = {mp.nstr(gamma_tv_full, 4)}")
    print(f"    Coefficient of variation rate  = {mp.nstr(gamma_cv_full, 4)}")

# ============================================================
# COMPLETION SENTINEL
# ============================================================

print()
print("=" * 80)
print("CELL 101 EXECUTION COMPLETE")
print("=" * 80)
