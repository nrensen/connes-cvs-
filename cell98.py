# ============================================================
# CELL 98 — FESH BACH / HIGH-MODE SEPARATION DIAGNOSTIC
# ============================================================
#
# Target:
#
#   Test the analytical mechanism
#
#       H_N
#         = [ A_M       B_MN ]
#           [ B_MN^T   C_MN ]
#
#   with Feshbach correction
#
#       R_MN(lambda)
#         = B_MN (C_MN - lambda I)^(-1) B_MN^T.
#
#   The intended analytical route is:
#
#       high-sector separation
#             +
#       tunnelling suppression of B_MN
#             =>
#       small change in the low-energy effective operator.
#
# This cell does NOT assume a formula for the tunnelling scale.
# It measures the actual Schur correction first.
#
# The canonical coordinates are
#
#       v = (v_0, ..., v_N),
#
# corresponding to the symmetric full-space vector
#
#       u_0 = v_0,
#       u_{+k} = u_{-k} = v_k / sqrt(2).
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
#
# M is the last retained canonical low-mode index.
#
# Thus:
#
#   P_M dimension = M + 1
#   Q_M dimension = N - M
#
M_VALUES = [
    4,
    8,
    12,
    16,
    24,
    32,
    48,
    64,
]

# Small safety tolerance for deciding whether lambda is too
# close to the high-sector spectrum.
DIST_TOL = mp.mpf("1e-40")


# ============================================================
# FULL -> CANONICAL MATRIX
# ============================================================
#
# The cached Galerkin matrix Q is expressed in the full symmetric
# coordinates
#
#     -N, ..., -1, 0, 1, ..., N.
#
# We convert it to the orthonormal canonical basis
#
#     e_0,
#     (e_k + e_-k)/sqrt(2),   k >= 1.
#
# The resulting matrix is (N+1) x (N+1).
#
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
    #
    # <e_j^+, Q e_k^+>
    #
    # = Q_{j,k} + Q_{j,-k}
    #
    # by symmetry.
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
# SYMMETRIC EIGENVALUES
# ============================================================

def eigvals_sym(A):
    """
    Return the eigenvalues of a real symmetric mpmath matrix
    in ascending order.
    """
    vals, _ = mp.eigsy(A)
    return [vals[i] for i in range(vals.rows)]


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
# HIGH-SECTOR SPECTRAL DISTANCE
# ============================================================

def high_sector_spectrum(C):
    """
    Return sorted eigenvalues of the high-sector block.
    """
    return eigvals_sym(C)


def distance_to_spectrum(lam, eigenvalues):
    """
    dist(lambda, spectrum(C)).
    """
    return min(
        abs(lam - mu)
        for mu in eigenvalues
    )


# ============================================================
# FESHBACH CORRECTION
# ============================================================

def feshbach_correction(A, B, C, lam):
    """
    Compute

        R(lambda)
          = B (C - lambda I)^(-1) B^T

    using a linear solve rather than explicitly forming the inverse.
    """

    q = C.rows

    shifted = mp.matrix(q, q)

    for i in range(q):
        for j in range(q):
            shifted[i, j] = C[i, j]

        shifted[i, i] -= lam

    # Solve
    #
    #     (C - lambda I) X = B^T.
    #
    # X has q rows and p columns.

    X = mp.lu_solve(
        shifted,
        B.T,
    )

    R = B * X

    # Remove tiny numerical asymmetry.
    p = R.rows

    Rsym = mp.matrix(p, p)

    for i in range(p):
        for j in range(p):
            Rsym[i, j] = (
                R[i, j] + R[j, i]
            ) / 2

    return Rsym


# ============================================================
# EFFECTIVE HAMILTONIAN
# ============================================================

def effective_matrix(A, R):
    """
    Feshbach effective Hamiltonian

        H_eff(lambda) = A - R(lambda).
    """

    H_eff = mp.matrix(
        A.rows,
        A.cols,
    )

    for i in range(A.rows):
        for j in range(A.cols):
            H_eff[i, j] = (
                A[i, j] - R[i, j]
            )

    return H_eff


# ============================================================
# EIGENVALUE REPORT
# ============================================================

def first_two_eigenvalues(A):
    vals = eigvals_sym(A)

    return vals[0], vals[1]


# ============================================================
# LOAD FULL GALERKIN MATRIX
# ============================================================

print()
print("=" * 72)
print("CELL 98 — FESHBACH / HIGH-MODE SEPARATION")
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

# Structural symmetry check.
sym_error = mp.mpf(0)

for i in range(H.rows):
    for j in range(H.cols):
        sym_error = max(
            sym_error,
            abs(H[i, j] - H[j, i]),
        )

print(
    "  symmetry error = "
    f"{mp.nstr(sym_error, 8)}"
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
E12 = full_eigs[1]
G11 = E12 - E11

print()
print("FULL CANONICAL SPECTRUM")
print(
    f"  E11 = {mp.nstr(E11, 18)}"
)
print(
    f"  E12 = {mp.nstr(E12, 18)}"
)
print(
    f"  g11 = {mp.nstr(G11, 18)}"
)
print(
    f"  eigensolve = {eig_elapsed:.3f} s"
)

# ============================================================
# FESHBACH SWEEP
# ============================================================

print()
print("=" * 72)
print("FESHBACH SWEEP")
print("=" * 72)
print()

print(
    "Columns:"
)
print(
    "M | dim(P) | min high eig | "
    "dist(E11,sigma(C)) | "
    "||B|| | "
    "envelope(E11) | "
    "||R(E11)|| | "
    "envelope(E12) | "
    "||R(E12)||"
)

print("-" * 160)

results = []

for M in M_VALUES:

    if M >= N:
        continue

    tM = time.perf_counter()

    A, B, C_high = feshbach_blocks(
        H,
        M,
    )

    # High-sector spectrum.
    high_eigs = high_sector_spectrum(
        C_high
    )

    high_min = high_eigs[0]
    high_max = high_eigs[-1]

    dist11 = distance_to_spectrum(
        E11,
        high_eigs,
    )

    dist12 = distance_to_spectrum(
        E12,
        high_eigs,
    )

    B_norm = symmetric_operator_norm(
        B.T * B
    ) ** mp.mpf("0.5")

    B_frob = frobenius_norm(B)

    # --------------------------------------------------------
    # Rigorous resolvent envelope:
    #
    # ||B(C-lambda)^(-1)B^T||
    #
    # <= ||B||^2 / dist(lambda,sigma(C)).
    #
    # --------------------------------------------------------

    envelope11 = (
        B_norm ** 2
        / dist11
    )

    envelope12 = (
        B_norm ** 2
        / dist12
    )

    # --------------------------------------------------------
    # Actual Feshbach corrections.
    # --------------------------------------------------------

    if dist11 <= DIST_TOL:
        R11 = None
        R11_norm = mp.inf
        eff11 = None
    else:
        R11 = feshbach_correction(
            A,
            B,
            C_high,
            E11,
        )

        R11_norm = symmetric_operator_norm(
            R11
        )

        eff11 = effective_matrix(
            A,
            R11,
        )

    if dist12 <= DIST_TOL:
        R12 = None
        R12_norm = mp.inf
        eff12 = None
    else:
        R12 = feshbach_correction(
            A,
            B,
            C_high,
            E12,
        )

        R12_norm = symmetric_operator_norm(
            R12
        )

        eff12 = effective_matrix(
            A,
            R12,
        )

    # --------------------------------------------------------
    # Effective eigenvalues.
    #
    # These are diagnostic because the Feshbach equation is
    # nonlinear in lambda. We evaluate the effective operator
    # at the actual full-space eigenvalues.
    # --------------------------------------------------------

    if eff11 is not None:
        eff11_eigs = eigvals_sym(
            eff11
        )

        # Find eigenvalue closest to E11.
        eff11_closest = min(
            eff11_eigs,
            key=lambda x:
                abs(x - E11),
        )

        eff11_residual = (
            abs(eff11_closest - E11)
        )
    else:
        eff11_closest = mp.nan
        eff11_residual = mp.nan

    if eff12 is not None:
        eff12_eigs = eigvals_sym(
            eff12
        )

        eff12_closest = min(
            eff12_eigs,
            key=lambda x:
                abs(x - E12),
        )

        eff12_residual = (
            abs(eff12_closest - E12)
        )
    else:
        eff12_closest = mp.nan
        eff12_residual = mp.nan

    elapsed = (
        time.perf_counter() - tM
    )

    results.append(
        {
            "M": M,
            "dim_P": M + 1,
            "dim_Q": N - M,
            "high_min": high_min,
            "high_max": high_max,
            "dist11": dist11,
            "dist12": dist12,
            "B_norm": B_norm,
            "B_frob": B_frob,
            "envelope11": envelope11,
            "envelope12": envelope12,
            "R11_norm": R11_norm,
            "R12_norm": R12_norm,
            "eff11_closest": eff11_closest,
            "eff12_closest": eff12_closest,
            "eff11_residual": eff11_residual,
            "eff12_residual": eff12_residual,
            "seconds": elapsed,
        }
    )

    print(
        f"{M:2d} | "
        f"{M+1:6d} | "
        f"{mp.nstr(high_min, 9):>12} | "
        f"{mp.nstr(dist11, 9):>16} | "
        f"{mp.nstr(B_norm, 9):>12} | "
        f"{mp.nstr(envelope11, 9):>15} | "
        f"{mp.nstr(R11_norm, 9):>12} | "
        f"{mp.nstr(envelope12, 9):>15} | "
        f"{mp.nstr(R12_norm, 9):>12}"
    )

# ============================================================
# SUMMARY
# ============================================================

print()
print("=" * 72)
print("SUMMARY")
print("=" * 72)
print()

print(
    "The quantity to watch is not merely ||B||."
)
print(
    "The Feshbach control parameter is"
)
print(
    "    ||B||^2 / dist(lambda, sigma(C))"
)
print(
    "and the actual correction is"
)
print(
    "    ||B(C-lambda)^(-1)B^T||."
)
print()

print(
    "If the high-mode mechanism is working as hoped, "
    "we want to see:"
)
print()
print(
    "  1. dist(E11, sigma(C)) bounded away from zero;"
)
print(
    "  2. ||B|| decreasing as M increases;"
)
print(
    "  3. the actual Schur correction decreasing;"
)
print(
    "  4. the actual correction substantially below "
    "its resolvent envelope;"
)
print(
    "  5. effective eigenvalue residuals becoming small."
)
print()

print(
    "The decisive analytical target is then a bound of the form"
)
print()
print(
    "    ||B_MN||^2 / gamma_M"
)
print()
print(
    "    <= C * Delta_M * R_spec(M,L),"
)
print()
print(
    "with the right-hand side -> 0."
)
print()

print(
    "This cell deliberately does not assume that bound."
)
print(
    "It measures the two sides needed to formulate it."
)
