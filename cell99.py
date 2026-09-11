# ============================================================
# CELL 99 — SPECTRAL DECOMPOSITION OF THE FESHBACH CORRECTION
# ============================================================
#
# Gate:       Gate 1 (Finite-N Spectral Mechanism)
#
# Target:
#
#   Investigate WHY the Feshbach / Schur correction
#
#       R_M(E) = B_M (C_M - E I)^(-1) B_M^T
#
#   becomes small at large M even though:
#
#       ||B_M|| = O(1)                (stays ~0.88--1.15)
#       delta_M = dist(E, sigma(C))   (saturates ~0.89)
#
#   by computing the spectral decomposition
#
#       R_M(E) = sum_j  (B_M u_j)(B_M u_j)^T / (mu_j - E)
#
#   where C_M u_j = mu_j u_j.
#
#   Cell 98 established that the crude resolvent bound
#
#       ||R_M(E)|| <= ||B_M||^2 / delta_M
#
#   saturates at O(1) and cannot explain the observed
#   decay of ||R_M(E)|| from ~2.6 to ~0.32.
#
#   The discarded information is the directional coupling
#
#       ||B_M u_j||^2
#
#   for each eigenvector u_j of C_M.
#
# Hypothesis:
#
#   The Galerkin matrix is a Loewner (divided-difference)
#   matrix of a smooth function psi.  The coupling block
#   B_M has entries that decay with mode separation.
#   When dotted against oscillatory high-mode eigenvectors,
#   Riemann-Lebesgue cancellation suppresses ||B_M u_j||
#   preferentially for the low-lying eigenmodes of C_M
#   (those closest to E).
#
# Falsification criterion:
#
#   If ||B_M u_j||^2 ~ ||B_M||^2 / dim(C_M) uniformly
#   (isotropic coupling), then the Loewner directional
#   suppression hypothesis is falsified.
#
# Verification:
#
#   1. sum_j ||B_M u_j||^2 must equal ||B_M||_F^2
#      (Parseval identity).
#
#   2. sum_j ||B_M u_j||^2 / |mu_j - E| must reproduce
#      the trace of |R_M(E)| computed by cell 98.
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

# Low-sector cutoffs M to inspect.
# Focus on the regime where delta_M is already O(1).
M_VALUES = [
    24,
    32,
    48,
    64,
]

# Number of top contributing modes to report in detail.
N_TOP = 10

# Small safety tolerance.
DIST_TOL = mp.mpf("1e-40")


# ============================================================
# FULL -> CANONICAL MATRIX
# ============================================================

def full_to_canonical_matrix(Q, N):
    """
    Convert the symmetric full-space Galerkin matrix Q into the
    canonical real-even orthonormal basis.

    Canonical basis:
        e_0
        (e_k + e_-k)/sqrt(2), k=1,...,N.
    """

    H = mp.matrix(N + 1, N + 1)

    centre = N

    # Zero mode.
    H[0, 0] = Q[centre, centre]

    # Zero / positive-mode coupling.
    for k in range(1, N + 1):
        H[0, k] = (
            mp.sqrt(2)
            * Q[centre, centre + k]
        )
        H[k, 0] = H[0, k]

    # Positive-mode block.
    for j in range(1, N + 1):
        for k in range(j, N + 1):
            value = (
                Q[centre + j, centre + k]
                + Q[centre + j, centre - k]
            )

            H[j, k] = value
            H[k, j] = value

    return H


# ============================================================
# SYMMETRIC EIGENVALUES AND EIGENVECTORS
# ============================================================

def eigvals_sym(A):
    """
    Return the eigenvalues of a real symmetric mpmath matrix
    in ascending order.
    """
    vals, _ = mp.eigsy(A)
    return [vals[i] for i in range(vals.rows)]


def eigsys_sym(A):
    """
    Return sorted (eigenvalues, eigenvectors_as_columns) of
    a real symmetric mpmath matrix.  Eigenvalues in ascending
    order.
    """
    vals, V = mp.eigsy(A)
    dim = A.rows

    idx = sorted(
        range(dim),
        key=lambda i: vals[i],
    )

    evals = [vals[i] for i in idx]

    evecs = mp.matrix(dim, dim)
    for col_out, col_in in enumerate(idx):
        for row in range(dim):
            evecs[row, col_out] = V[row, col_in]

    return evals, evecs


# ============================================================
# MATRIX OPERATOR NORM
# ============================================================

def symmetric_operator_norm(A):
    """
    Operator norm of a real symmetric matrix.
    """
    vals = eigvals_sym(A)

    return max(
        abs(vals[0]),
        abs(vals[-1]),
    )


# ============================================================
# FROBENIUS NORM
# ============================================================

def frobenius_norm(A):
    """
    Frobenius norm.
    """
    total = mp.mpf(0)

    for i in range(A.rows):
        for j in range(A.cols):
            total += A[i, j] ** 2

    return mp.sqrt(total)


# ============================================================
# BLOCK DECOMPOSITION
# ============================================================

def feshbach_blocks(H, M):
    """
    Split the canonical matrix H at low-mode cutoff M.

    P-sector:
        indices 0,...,M

    Q-sector:
        indices M+1,...,N
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
    C = mp.matrix(q, q)

    for i in range(p):
        for j in range(p):
            A[i, j] = H[i, j]

    for i in range(p):
        for j in range(q):
            B[i, j] = H[i, p + j]

    for i in range(q):
        for j in range(q):
            C[i, j] = H[p + i, p + j]

    return A, B, C


# ============================================================
# FESHBACH CORRECTION (from cell 98, for cross-check)
# ============================================================

def feshbach_correction(A, B, C, lam):
    """
    Compute R(lambda) = B (C - lambda I)^(-1) B^T
    using column-by-column linear solves.
    """
    q = C.rows
    p = B.rows

    shifted = mp.matrix(q, q)

    for i in range(q):
        for j in range(q):
            shifted[i, j] = C[i, j]

        shifted[i, i] -= lam

    X = mp.matrix(q, p)

    for j in range(p):
        rhs = mp.matrix(q, 1)

        for i in range(q):
            rhs[i] = B[j, i]

        x = mp.lu_solve(
            shifted,
            rhs,
        )

        for i in range(q):
            X[i, j] = x[i]

    R = B * X

    # Symmetrise.
    Rsym = mp.matrix(p, p)

    for i in range(p):
        for j in range(p):
            Rsym[i, j] = (
                R[i, j] + R[j, i]
            ) / 2

    return Rsym


# ============================================================
# LOAD FULL GALERKIN MATRIX
# ============================================================

print()
print("=" * 72)
print("CELL 99 — SPECTRAL DECOMPOSITION OF FESHBACH CORRECTION")
print("=" * 72)
print()
print(
    f"C={C}, N={N}, T={T}, dps={DPS}"
)
print()

t0 = time.perf_counter()

Q, Q_meta = get_galerkin_matrix(
    c=C,
    N=N,
    T=T,
    dps=DPS,
    verbose=True,
)

load_elapsed = time.perf_counter() - t0

print()
print(
    f"Full Galerkin matrix loaded in "
    f"{load_elapsed:.3f} s"
)


# ============================================================
# CANONICAL MATRIX
# ============================================================

t0 = time.perf_counter()

H = full_to_canonical_matrix(
    Q,
    N,
)

canonical_elapsed = (
    time.perf_counter() - t0
)

print()
print(
    "Canonical matrix:"
)
print(
    f"  dimension = {H.rows}"
)
print(
    f"  conversion = {canonical_elapsed:.3f} s"
)


# ============================================================
# FULL LOW-ENERGY SPECTRUM
# ============================================================

t0 = time.perf_counter()

full_eigs = eigvals_sym(H)

eig_elapsed = (
    time.perf_counter() - t0
)

E11 = full_eigs[0]

print()
print("FULL CANONICAL SPECTRUM")
print(
    f"  E11 = {mp.nstr(E11, 18)}"
)
print(
    f"  eigensolve = {eig_elapsed:.3f} s"
)


# ============================================================
# FESHBACH SPECTRAL DECOMPOSITION SWEEP
# ============================================================

print()
print("=" * 72)
print("SPECTRAL DECOMPOSITION SWEEP")
print("=" * 72)

for M in M_VALUES:

    if M >= N:
        continue

    print()
    print("-" * 72)
    print(f"M = {M}  |  dim(P) = {M+1}  |  dim(Q) = {N-M}")
    print("-" * 72)

    tM = time.perf_counter()

    A, B, C_high = feshbach_blocks(
        H,
        M,
    )

    p = M + 1
    q = N - M

    # --------------------------------------------------------
    # High-sector eigensystem.
    # --------------------------------------------------------

    t_eig = time.perf_counter()

    high_evals, high_evecs = eigsys_sym(
        C_high
    )

    eig_secs = time.perf_counter() - t_eig

    high_min = high_evals[0]

    delta_M = min(
        abs(E11 - mu)
        for mu in high_evals
    )

    print(
        f"  high-sector eigensolve = {eig_secs:.3f} s"
    )
    print(
        f"  min(sigma(C)) = {mp.nstr(high_min, 12)}"
    )
    print(
        f"  dist(E11, sigma(C)) = {mp.nstr(delta_M, 12)}"
    )

    # --------------------------------------------------------
    # Coupling norms.
    # --------------------------------------------------------

    B_norm = symmetric_operator_norm(
        B.T * B
    ) ** mp.mpf("0.5")

    B_frob = frobenius_norm(B)

    print(
        f"  ||B||_op = {mp.nstr(B_norm, 12)}"
    )
    print(
        f"  ||B||_F  = {mp.nstr(B_frob, 12)}"
    )

    # --------------------------------------------------------
    # Compute ||B u_j||^2 for each high-mode eigenvector.
    # --------------------------------------------------------
    #
    # B is p x q, u_j is q x 1.
    # B u_j is p x 1.
    # ||B u_j||^2 = sum_i (B u_j)_i^2.
    #
    # --------------------------------------------------------

    coupling_strengths = []  # ||B u_j||^2
    contributions = []       # ||B u_j||^2 / |mu_j - E11|

    parseval_sum = mp.mpf(0)
    trace_sum = mp.mpf(0)

    for j in range(q):
        # Extract eigenvector column.
        uj = mp.matrix(q, 1)
        for k in range(q):
            uj[k] = high_evecs[k, j]

        # Compute B u_j.
        Bu_j = B * uj

        # ||B u_j||^2.
        norm_sq = mp.mpf(0)
        for i in range(p):
            norm_sq += Bu_j[i] ** 2

        coupling_strengths.append(norm_sq)
        parseval_sum += norm_sq

        # Contribution to resolvent trace.
        denom = abs(high_evals[j] - E11)

        if denom > DIST_TOL:
            contrib = norm_sq / denom
        else:
            contrib = mp.inf

        contributions.append(contrib)
        trace_sum += contrib

    # --------------------------------------------------------
    # Cross-check: actual correction via direct solve.
    # --------------------------------------------------------

    if delta_M > DIST_TOL:
        R_direct = feshbach_correction(
            A, B, C_high, E11,
        )

        R_direct_norm = symmetric_operator_norm(
            R_direct
        )
    else:
        R_direct_norm = mp.inf

    # --------------------------------------------------------
    # Parseval check: sum_j ||B u_j||^2 = ||B||_F^2.
    # --------------------------------------------------------

    parseval_err = abs(
        parseval_sum - B_frob ** 2
    )

    print()
    print("  PARSEVAL CHECK")
    print(
        f"    sum_j ||B u_j||^2"
        f" = {mp.nstr(parseval_sum, 15)}"
    )
    print(
        f"    ||B||_F^2"
        f"         = {mp.nstr(B_frob ** 2, 15)}"
    )
    print(
        f"    absolute error"
        f"    = {mp.nstr(parseval_err, 6)}"
    )

    # --------------------------------------------------------
    # Summary statistics.
    # --------------------------------------------------------

    # Isotropic baseline: if coupling were uniform,
    # each mode would get ||B||_F^2 / q.
    isotropic_baseline = B_frob ** 2 / q

    # Find the index of the mode closest to E11.
    idx_closest = min(
        range(q),
        key=lambda j: abs(high_evals[j] - E11),
    )

    print()
    print("  COUPLING PROFILE")
    print(
        f"    isotropic baseline"
        f" ||B||_F^2/q = {mp.nstr(isotropic_baseline, 8)}"
    )
    print(
        f"    closest mode to E11:"
        f" j={idx_closest},"
        f" mu_j = {mp.nstr(high_evals[idx_closest], 12)},"
        f" ||B u_j||^2"
        f" = {mp.nstr(coupling_strengths[idx_closest], 8)}"
    )
    print(
        f"    ratio to isotropic:"
        f" {mp.nstr(coupling_strengths[idx_closest] / isotropic_baseline, 6)}"
    )

    # --------------------------------------------------------
    # Top N_TOP contributors to the resolvent trace.
    # --------------------------------------------------------

    sorted_idx = sorted(
        range(q),
        key=lambda j: contributions[j],
        reverse=True,
    )

    print()
    print(
        f"  TOP {N_TOP} CONTRIBUTORS"
        f" to sum_j ||B u_j||^2 / |mu_j - E11|"
    )
    print(
        f"  {'rank':>4}  {'j':>4}"
        f"  {'mu_j':>16}"
        f"  {'||Bu_j||^2':>16}"
        f"  {'|mu_j-E|':>16}"
        f"  {'contrib':>16}"
        f"  {'cumul_frac':>12}"
    )

    cumul = mp.mpf(0)

    for rank, j in enumerate(
        sorted_idx[:N_TOP]
    ):
        cumul += contributions[j]

        frac = (
            cumul / trace_sum
            if trace_sum > 0
            else mp.mpf(0)
        )

        print(
            f"  {rank+1:4d}  {j:4d}"
            f"  {mp.nstr(high_evals[j], 12):>16}"
            f"  {mp.nstr(coupling_strengths[j], 8):>16}"
            f"  {mp.nstr(abs(high_evals[j] - E11), 8):>16}"
            f"  {mp.nstr(contributions[j], 8):>16}"
            f"  {mp.nstr(frac, 6):>12}"
        )

    # --------------------------------------------------------
    # Coupling decay profile: ||B u_j||^2 as a function of
    # eigenvalue index j (sorted by eigenvalue).
    #
    # Report in bins: bottom 5, next 5, ..., to show
    # whether coupling concentrates on low or high modes.
    # --------------------------------------------------------

    print()
    print("  COUPLING STRENGTH BY EIGENVALUE QUARTILE")

    n_bins = min(4, q)
    bin_size = q // n_bins
    remainder = q - n_bins * bin_size

    bin_start = 0
    for b in range(n_bins):
        bs = bin_size + (1 if b < remainder else 0)
        bin_end = bin_start + bs

        bin_sum = sum(
            coupling_strengths[j]
            for j in range(bin_start, bin_end)
        )

        bin_frac = (
            bin_sum / parseval_sum
            if parseval_sum > 0
            else mp.mpf(0)
        )

        mu_lo = high_evals[bin_start]
        mu_hi = high_evals[bin_end - 1]

        print(
            f"    modes {bin_start:4d}--{bin_end-1:4d}"
            f"  mu in [{mp.nstr(mu_lo, 6)}"
            f", {mp.nstr(mu_hi, 6)}]"
            f"  sum ||Bu_j||^2"
            f" = {mp.nstr(bin_sum, 8)}"
            f"  ({mp.nstr(100*bin_frac, 4)}%)"
        )

        bin_start = bin_end

    # --------------------------------------------------------
    # Overall comparison.
    # --------------------------------------------------------

    envelope = B_norm ** 2 / delta_M

    print()
    print("  COMPARISON")
    print(
        f"    crude envelope"
        f" ||B||^2/delta_M"
        f" = {mp.nstr(envelope, 8)}"
    )
    print(
        f"    trace bound"
        f" sum ||Bu_j||^2/|mu_j-E|"
        f" = {mp.nstr(trace_sum, 8)}"
    )
    print(
        f"    actual ||R(E11)||"
        f"                 = {mp.nstr(R_direct_norm, 8)}"
    )
    print(
        f"    ratio actual/envelope"
        f"        = {mp.nstr(R_direct_norm / envelope, 6)}"
    )
    print(
        f"    ratio actual/trace_bound"
        f"     = {mp.nstr(R_direct_norm / trace_sum, 6)}"
    )

    elapsed = time.perf_counter() - tM

    print(
        f"  section elapsed = {elapsed:.3f} s"
    )


# ============================================================
# COMPLETION SENTINEL
# ============================================================

print()
print("=" * 80)
print("CELL 99 EXECUTION COMPLETE")
print("=" * 80)
