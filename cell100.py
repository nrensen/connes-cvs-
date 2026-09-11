# ============================================================
# CELL 100 — COUPLING-WEIGHTED SPECTRAL MEASURE
# ============================================================
#
# Gate:       Gate 1 (Finite-N Spectral Mechanism)
#
# Purpose:
#
#   Cell 99 showed that the Feshbach correction
#
#       R_M(E) = B_M (C_M - E I)^(-1) B_M^T
#
#   is not controlled by the nearest high-mode pole alone.
#   The natural scalar object is the coupling-weighted spectral
#   measure
#
#       nu_M = sum_j a_j delta_{mu_j},
#       a_j = ||B_M u_j||^2,
#
#   and its Stieltjes transform
#
#       S_M(E) = sum_j a_j / (mu_j - E)
#
#   (here E = E11 and C_M - E I is positive in the measured
#   regime).
#
#   Cell 100 therefore does NOT assume a particular mechanism.
#   It asks a narrower structural question:
#
#       Does the coupling-weighted spectral measure itself move
#       to higher energies as M increases, in a way that can
#       explain the decrease of the Feshbach correction?
#
# Diagnostics:
#
#   1. Coupling-mass quantiles in mu.
#   2. Coupling-weighted mean / variance of mu.
#   3. Low-edge cumulative coupling mass.
#   4. Resolvent-weighted mass distribution.
#   5. Resolvent-weighted mean energy and effective denominator.
#   6. Contribution by fixed energy windows.
#   7. Contribution by eigenvalue-index quartiles.
#   8. Exact trace S_M(E) and direct ||R_M(E)|| cross-check.
#   9. Effective spectral rank trace(R)/||R||.
#  10. Comparison with an isotropic coupling model having the
#      same ||B||_F^2 and the same high-sector spectrum.
#
# The important distinction is:
#
#   ||B||_F^2 = total coupling mass,
#
#   S_M(E) = coupling mass weighted by 1/(mu-E).
#
# A decrease in S_M can therefore come from redistribution of
# coupling mass, from spectral motion, or both. This cell measures
# the two effects separately.
#
# ============================================================

import mpmath as mp
import time

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

# Number of leading resolvent contributors to print.
N_TOP = 10

DIST_TOL = mp.mpf("1e-40")

# Fixed energy thresholds. These are deliberately descriptive
# rather than fitted. They cover the lower part of the observed
# high-sector spectrum at M >= 24.
ENERGY_WINDOWS = [
    (mp.mpf("0"), mp.mpf("1")),
    (mp.mpf("1"), mp.mpf("2")),
    (mp.mpf("2"), mp.mpf("3")),
    (mp.mpf("3"), mp.mpf("4")),
    (mp.mpf("4"), mp.mpf("5")),
    (mp.mpf("5"), mp.inf),
]

# ============================================================
# MATRIX CONVERSION
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


def symmetric_operator_norm(A):
    vals = eigvals_sym(A)
    return max(abs(vals[0]), abs(vals[-1]))


def frobenius_norm(A):
    total = mp.mpf(0)
    for i in range(A.rows):
        for j in range(A.cols):
            total += A[i, j] ** 2
    return mp.sqrt(total)


# ============================================================
# FESHBACH BLOCKS
# ============================================================

def feshbach_blocks(H, M):
    """
    P-sector = indices 0,...,M
    Q-sector = indices M+1,...,N.
    """
    d = H.rows

    if not (0 < M < d - 1):
        raise ValueError(
            f"M={M} incompatible with matrix dimension {d}"
        )

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


def feshbach_correction(B, C_high, lam):
    """
    Compute R(lambda) = B(C-lambda I)^(-1)B^T.
    """
    q = C_high.rows
    p = B.rows

    shifted = mp.matrix(q, q)

    for i in range(q):
        for j in range(q):
            shifted[i, j] = C_high[i, j]
        shifted[i, i] -= lam

    X = mp.matrix(q, p)

    for j in range(p):
        rhs = mp.matrix(q, 1)
        for i in range(q):
            rhs[i] = B[j, i]

        x = mp.lu_solve(shifted, rhs)

        for i in range(q):
            X[i, j] = x[i]

    R = B * X

    Rsym = mp.matrix(p, p)
    for i in range(p):
        for j in range(p):
            Rsym[i, j] = (R[i, j] + R[j, i]) / 2

    return Rsym


# ============================================================
# WEIGHTED-QUANTILE HELPER
# ============================================================

def weighted_quantile(sorted_pairs, total_mass, fraction):
    """
    sorted_pairs = [(mu, weight), ...] sorted by mu.
    Returns the first mu whose cumulative weight reaches
    fraction * total_mass.
    """
    target = fraction * total_mass
    cumulative = mp.mpf(0)

    for mu, weight in sorted_pairs:
        cumulative += weight
        if cumulative >= target:
            return mu

    return sorted_pairs[-1][0]


# ============================================================
# LOAD MATRIX
# ============================================================

print()
print("=" * 72)
print("CELL 100 — COUPLING-WEIGHTED SPECTRAL MEASURE")
print("=" * 72)
print()
print(f"C={C}, N={N}, T={T}, dps={DPS}")
print()

t0 = time.perf_counter()

Q, Q_meta = get_galerkin_matrix(
    c=C,
    N=N,
    T=T,
    dps=DPS,
    verbose=True,
)

print()
print(
    f"Full Galerkin matrix loaded in "
    f"{time.perf_counter() - t0:.3f} s"
)

# ============================================================
# CANONICAL MATRIX AND GROUND ENERGY
# ============================================================

t0 = time.perf_counter()

H = full_to_canonical_matrix(Q, N)

print()
print("Canonical matrix:")
print(f"  dimension = {H.rows}")
print(f"  conversion = {time.perf_counter() - t0:.3f} s")

t0 = time.perf_counter()

full_eigs = eigvals_sym(H)
E11 = full_eigs[0]

print()
print("FULL CANONICAL SPECTRUM")
print(f"  E11 = {mp.nstr(E11, 18)}")
print(f"  eigensolve = {time.perf_counter() - t0:.3f} s")

# ============================================================
# SWEEP
# ============================================================

print()
print("=" * 72)
print("COUPLING-WEIGHTED SPECTRAL MEASURE SWEEP")
print("=" * 72)

for M in M_VALUES:

    print()
    print("-" * 72)
    print(f"M = {M}  |  dim(P) = {M+1}  |  dim(Q) = {N-M}")
    print("-" * 72)

    tM = time.perf_counter()

    A, B, C_high = feshbach_blocks(H, M)
    p = M + 1
    q = N - M

    # --------------------------------------------------------
    # High-sector eigensystem.
    # --------------------------------------------------------

    t_eig = time.perf_counter()

    high_evals, high_evecs = eigsys_sym(C_high)

    print(
        f"  high-sector eigensolve = "
        f"{time.perf_counter() - t_eig:.3f} s"
    )

    # --------------------------------------------------------
    # Basic spectral / coupling data.
    # --------------------------------------------------------

    high_min = high_evals[0]
    high_max = high_evals[-1]

    delta_M = min(abs(E11 - mu) for mu in high_evals)

    B_frob = frobenius_norm(B)
    B_norm = symmetric_operator_norm(B.T * B) ** mp.mpf("0.5")

    coupling = []
    resolvent_contrib = []

    total_mass = mp.mpf(0)
    trace_sum = mp.mpf(0)

    for j in range(q):
        uj = mp.matrix(q, 1)

        for k in range(q):
            uj[k] = high_evecs[k, j]

        Bu = B * uj

        a = mp.mpf(0)
        for i in range(p):
            a += Bu[i] ** 2

        mu = high_evals[j]
        denom = abs(mu - E11)

        coupling.append(a)
        total_mass += a

        if denom > DIST_TOL:
            r = a / denom
        else:
            r = mp.inf

        resolvent_contrib.append(r)
        trace_sum += r

    # --------------------------------------------------------
    # Parseval.
    # --------------------------------------------------------

    parseval_err = abs(total_mass - B_frob ** 2)

    # --------------------------------------------------------
    # Coupling-weighted moments.
    # --------------------------------------------------------

    mean_mu = mp.fsum(
        coupling[j] * high_evals[j]
        for j in range(q)
    ) / total_mass

    second_mu = mp.fsum(
        coupling[j] * high_evals[j] ** 2
        for j in range(q)
    ) / total_mass

    variance_mu = max(
        mp.mpf(0),
        second_mu - mean_mu ** 2,
    )

    std_mu = mp.sqrt(variance_mu)

    coupling_pairs = [
        (high_evals[j], coupling[j])
        for j in range(q)
    ]

    # eigsy already returns eigenvalues sorted, so this is sorted.
    q10 = weighted_quantile(coupling_pairs, total_mass, mp.mpf("0.10"))
    q25 = weighted_quantile(coupling_pairs, total_mass, mp.mpf("0.25"))
    q50 = weighted_quantile(coupling_pairs, total_mass, mp.mpf("0.50"))
    q75 = weighted_quantile(coupling_pairs, total_mass, mp.mpf("0.75"))
    q90 = weighted_quantile(coupling_pairs, total_mass, mp.mpf("0.90"))

    # --------------------------------------------------------
    # Low-edge cumulative mass.
    # --------------------------------------------------------

    print()
    print("  SPECTRAL / COUPLING SUMMARY")
    print(f"    high spectrum = [{mp.nstr(high_min, 10)}, "
          f"{mp.nstr(high_max, 10)}]")
    print(f"    delta_M = {mp.nstr(delta_M, 10)}")
    print(f"    ||B||_op = {mp.nstr(B_norm, 10)}")
    print(f"    ||B||_F^2 = {mp.nstr(B_frob**2, 12)}")
    print(f"    Parseval error = {mp.nstr(parseval_err, 6)}")
    print()
    print("    coupling-weighted mu moments:")
    print(f"      mean = {mp.nstr(mean_mu, 12)}")
    print(f"      std  = {mp.nstr(std_mu, 12)}")
    print(
        f"      quantiles 10/25/50/75/90% = "
        f"{mp.nstr(q10, 8)} / {mp.nstr(q25, 8)} / "
        f"{mp.nstr(q50, 8)} / {mp.nstr(q75, 8)} / "
        f"{mp.nstr(q90, 8)}"
    )

    print()
    print("  COUPLING MASS BY FIXED ENERGY WINDOW")

    for lo, hi in ENERGY_WINDOWS:
        mass = mp.mpf(0)

        for j in range(q):
            mu = high_evals[j]
            if lo <= mu < hi:
                mass += coupling[j]

        frac = mass / total_mass if total_mass else mp.mpf(0)

        if hi == mp.inf:
            label = f"[{mp.nstr(lo, 4)}, inf)"
        else:
            label = f"[{mp.nstr(lo, 4)}, {mp.nstr(hi, 4)})"

        print(
            f"    {label:>18} : "
            f"{mp.nstr(mass, 10)} "
            f"({mp.nstr(100*frac, 6)}%)"
        )

    # --------------------------------------------------------
    # Resolvent-weighted measure.
    #
    # w_j = [a_j/(mu_j-E)] / S_M(E)
    #
    # This is the exact scalar spectral distribution appearing
    # in the trace bound.
    # --------------------------------------------------------

    resolvent_mean_mu = mp.fsum(
        resolvent_contrib[j] * high_evals[j]
        for j in range(q)
    ) / trace_sum

    resolvent_mean_inv_gap = (
        trace_sum / total_mass
    )

    effective_denominator = (
        total_mass / trace_sum
    )

    print()
    print("  RESOLVENT-WEIGHTED SPECTRAL MEASURE")
    print(
        f"    S_M(E11) = "
        f"{mp.nstr(trace_sum, 12)}"
    )
    print(
        f"    mean mu under resolvent weights = "
        f"{mp.nstr(resolvent_mean_mu, 12)}"
    )
    print(
        f"    mean 1/(mu-E11) under coupling weights = "
        f"{mp.nstr(resolvent_mean_inv_gap, 12)}"
    )
    print(
        f"    effective denominator "
        f"||B||_F^2 / S_M = "
        f"{mp.nstr(effective_denominator, 12)}"
    )

    # --------------------------------------------------------
    # Resolvent contribution by fixed energy window.
    # --------------------------------------------------------

    print()
    print("  RESOLVENT CONTRIBUTION BY FIXED ENERGY WINDOW")

    for lo, hi in ENERGY_WINDOWS:
        mass = mp.mpf(0)

        for j in range(q):
            mu = high_evals[j]
            if lo <= mu < hi:
                mass += resolvent_contrib[j]

        frac = mass / trace_sum if trace_sum else mp.mpf(0)

        if hi == mp.inf:
            label = f"[{mp.nstr(lo, 4)}, inf)"
        else:
            label = f"[{mp.nstr(lo, 4)}, {mp.nstr(hi, 4)})"

        print(
            f"    {label:>18} : "
            f"{mp.nstr(mass, 10)} "
            f"({mp.nstr(100*frac, 6)}%)"
        )

    # --------------------------------------------------------
    # Coupling contribution by eigenvalue-index quartile.
    # --------------------------------------------------------

    print()
    print("  INDEX-QUARTILE COUPLING MASS")

    n_bins = min(4, q)
    bin_size = q // n_bins
    remainder = q - n_bins * bin_size

    start = 0

    for b in range(n_bins):
        size = bin_size + (1 if b < remainder else 0)
        end = start + size

        mass = mp.fsum(coupling[start:end])
        frac = mass / total_mass

        print(
            f"    modes {start:4d}--{end-1:4d}: "
            f"{mp.nstr(mass, 10)} "
            f"({mp.nstr(100*frac, 6)}%)"
        )

        start = end

    # --------------------------------------------------------
    # Index-quartile resolvent contribution.
    # --------------------------------------------------------

    print()
    print("  INDEX-QUARTILE RESOLVENT CONTRIBUTION")

    start = 0

    for b in range(n_bins):
        size = bin_size + (1 if b < remainder else 0)
        end = start + size

        mass = mp.fsum(resolvent_contrib[start:end])
        frac = mass / trace_sum

        print(
            f"    modes {start:4d}--{end-1:4d}: "
            f"{mp.nstr(mass, 10)} "
            f"({mp.nstr(100*frac, 6)}%)"
        )

        start = end

    # --------------------------------------------------------
    # Leading contributors.
    # --------------------------------------------------------

    print()
    print(
        f"  TOP {N_TOP} RESOLVENT CONTRIBUTORS"
    )
    print(
        f"  {'rank':>4} {'j':>4} {'mu_j':>14} "
        f"{'||Bu_j||^2':>14} {'contrib':>14} "
        f"{'cumul_frac':>12}"
    )

    sorted_idx = sorted(
        range(q),
        key=lambda j: resolvent_contrib[j],
        reverse=True,
    )

    cumulative = mp.mpf(0)

    for rank, j in enumerate(sorted_idx[:N_TOP]):
        cumulative += resolvent_contrib[j]
        frac = cumulative / trace_sum

        print(
            f"  {rank+1:4d} {j:4d} "
            f"{mp.nstr(high_evals[j], 10):>14} "
            f"{mp.nstr(coupling[j], 8):>14} "
            f"{mp.nstr(resolvent_contrib[j], 8):>14} "
            f"{mp.nstr(frac, 6):>12}"
        )

    # --------------------------------------------------------
    # Isotropic comparison.
    #
    # Give every high eigenmode the same coupling weight
    # total_mass/q. This isolates the effect of the spectral
    # locations from the actual directional coupling profile.
    # --------------------------------------------------------

    isotropic_weight = total_mass / q

    isotropic_trace = mp.fsum(
        isotropic_weight / abs(high_evals[j] - E11)
        for j in range(q)
    )

    print()
    print("  ISOTROPIC COUPLING COMPARISON")
    print(
        f"    isotropic trace bound = "
        f"{mp.nstr(isotropic_trace, 12)}"
    )
    print(
        f"    actual coupling trace bound = "
        f"{mp.nstr(trace_sum, 12)}"
    )
    print(
        f"    actual / isotropic = "
        f"{mp.nstr(trace_sum / isotropic_trace, 10)}"
    )

    # --------------------------------------------------------
    # Direct Feshbach cross-check and effective rank.
    # --------------------------------------------------------

    if delta_M > DIST_TOL:
        R_direct = feshbach_correction(
            B,
            C_high,
            E11,
        )
        R_norm = symmetric_operator_norm(R_direct)
    else:
        R_norm = mp.inf

    envelope = B_norm ** 2 / delta_M

    effective_rank = (
        trace_sum / R_norm
        if R_norm not in (0, mp.inf)
        else mp.inf
    )

    print()
    print("  FESHBACH COMPARISON")
    print(
        f"    crude envelope ||B||^2/delta = "
        f"{mp.nstr(envelope, 12)}"
    )
    print(
        f"    trace bound S_M(E11) = "
        f"{mp.nstr(trace_sum, 12)}"
    )
    print(
        f"    actual ||R(E11)|| = "
        f"{mp.nstr(R_norm, 12)}"
    )
    print(
        f"    actual / trace bound = "
        f"{mp.nstr(R_norm / trace_sum, 10)}"
    )
    print(
        f"    effective spectral rank "
        f"trace(R)/||R|| = "
        f"{mp.nstr(effective_rank, 10)}"
    )

    print()
    print(
        f"  section elapsed = "
        f"{time.perf_counter() - tM:.3f} s"
    )

# ============================================================
# COMPLETION SENTINEL
# ============================================================

print()
print("=" * 80)
print("CELL 100 EXECUTION COMPLETE")
print("=" * 80)

