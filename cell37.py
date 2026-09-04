"""
CELL 37 — REORGANISE THE ARCHIMEDEAN TAIL COEFFICIENTS

Purpose
-------
Cell 36 established the exact asymptotic coefficients

    R_v(r) ~ A_0/r^2 + A_1/r^4 + A_2/r^6 + ...

for the reduced finite-N Fourier kernel

    K_fourier(v,r,L) = (1 - cos(r L)) R_v(r).

This cell takes the next analytical step.  Rather than fitting anything
numerically, it asks whether the exact coefficient A_k can be reorganised
into endpoint data of

    T_v(t) = v_0 + mp.sqrt(2) sum_{m=1}^N v_m cos(2 pi m t/L).

The key objects are the even endpoint derivatives

    D_k = T_v^(2k)(0)
        = mp.sqrt(2) (-1)^k kappa^(2k) M_{2k},  k >= 1,

where kappa = 2 pi/L and M_p = sum m^p v_m.

We compare the exact A_k from Cell 36 with several natural quadratic
endpoint invariants.  The aim is structural discovery, not an assumed
N -> infinity law.

No numerical quadrature is used.
"""

from __future__ import annotations

import mpmath as mp

from cell import get_ground_state, SURVEY_GROUND_STATE


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SURVEY = SURVEY_GROUND_STATE
c = SURVEY["c"]
L = mp.log(c)
T = SURVEY["T"]
dps = SURVEY["dps"]

N_VALUES = range(1, SURVEY["N"] + 1)
K_MAX = 5


# ---------------------------------------------------------------------------
# Basic exact-in-formula helpers
# ---------------------------------------------------------------------------

def moment(v, p: int):
    """M_p = sum_{m=1}^N m^p v_m."""
    return sum((mp.mpf(m) ** p) * v[m] for m in range(1, len(v)))


def endpoint_derivative(v, k: int, L):
    """Return T_v^(2k)(0), with k=0 meaning T_v(0)."""
    if k == 0:
        return v[0] + mp.sqrt(2) * moment(v, 0)
    kappa = 2 * mp.pi / L
    return mp.sqrt(2) * (-1) ** k * (kappa ** (2 * k)) * moment(v, 2 * k)


def H_k(m: int, n: int, k: int):
    """
    H_k(m,n) = (n^(2k+2)-m^(2k+2))/(n^2-m^2)
             = sum_{j=0}^k n^(2(k-j)) m^(2j).
    """
    return sum(
        (mp.mpf(n) ** (2 * (k - j))) *
        (mp.mpf(m) ** (2 * j))
        for j in range(k + 1)
    )


def exact_A_coefficients(v, L, K=K_MAX):
    """
    Exact asymptotic coefficients A_k of the reduced kernel R_v(r).

    A_0 is handled separately because C_0 contributes only its leading
    1/(L r^2) term.

    For k >= 1:
        A_k =
            4(k+1)/L * kappa^(2k) * sum m^(2k) v_m^2
          + 4 mp.sqrt(2)/L * kappa^(2k) * v_0 M_(2k)
          + 4/pi * kappa^(2k+1)
                * sum_{m<n} v_m v_n H_k(m,n).
    """
    kappa = 2 * mp.pi / L
    out = []

    v0 = v[0]
    M0 = moment(v, 0)

    A0 = (
        (2 / L) * v0 ** 2
        + (4 / L) * sum(v[m] ** 2 for m in range(1, len(v)))
        + (4 * mp.sqrt(2) / L) * v0 * M0
        + (8 / L) * sum(
            v[m] * v[n]
            for m in range(1, len(v))
            for n in range(m + 1, len(v))
        )
    )
    out.append(A0)

    for k in range(1, K + 1):
        diag = sum(
            (mp.mpf(m) ** (2 * k)) * v[m] ** 2
            for m in range(1, len(v))
        )
        mixed0 = v0 * moment(v, 2 * k)

        cross = sum(
            v[m] * v[n] * H_k(m, n, k)
            for m in range(1, len(v))
            for n in range(m + 1, len(v))
        )

        Ak = (
            (4 * (k + 1) / L) * (kappa ** (2 * k)) * diag
            + (4 * mp.sqrt(2) / L) * (kappa ** (2 * k)) * mixed0
            + (4 / mp.pi) * (kappa ** (2 * k + 1)) * cross
        )
        out.append(Ak)

    return out


# ---------------------------------------------------------------------------
# Alternative exact representations
# ---------------------------------------------------------------------------

def endpoint_quadratic_sum(v, k: int):
    """
    Q_k = sum_{m,n >= 1} v_m v_n H_k(m,n), written directly from
    the polynomial H_k.

    This is kept separate so that the relation between the cross term
    and full double sums can be inspected algebraically.
    """
    N = len(v) - 1
    total = mp.mpf("0")
    for m in range(1, N + 1):
        for n in range(1, N + 1):
            if m == n:
                # H_k(m,m) is the continuous value (k+1)m^(2k).
                H = (k + 1) * (mp.mpf(m) ** (2 * k))
            else:
                H = H_k(min(m, n), max(m, n), k)
            total += v[m] * v[n] * H
    return total


def endpoint_moment_quadratic(v, k: int):
    """
    E_k = (sum m^(2j) v_m) combinations suggested by the polynomial
    H_k.  For inspection only; the cell prints the exact residual between
    Q_k and products of moments.
    """
    return sum(
        moment(v, 2 * j) * moment(v, 2 * (k - j))
        for j in range(k + 1)
    )


# ---------------------------------------------------------------------------
# Symbolic identity probe
# ---------------------------------------------------------------------------

def print_polynomial_identity_probe():
    """
    Display the first few H_k polynomials.  These reveal immediately
    whether the cross term can be converted into ordinary endpoint
    moments.
    """
    print("\n" + "=" * 78)
    print("POLYNOMIAL STRUCTURE OF H_k")
    print("=" * 78)

    for k in range(K_MAX + 1):
        terms = []
        for j in range(k + 1):
            if j == 0:
                terms.append(f"n^{2*k}")
            elif j == k:
                terms.append(f"m^{2*k}")
            else:
                terms.append(f"n^{2*(k-j)} m^{2*j}")
        print(f"H_{k}(m,n) = " + " + ".join(terms))


# ---------------------------------------------------------------------------
# Numerical survey of exact algebraic identities
# ---------------------------------------------------------------------------

def survey_one_N(N: int):
    lam, v, meta = get_ground_state(c=c, N=N, T=T, dps=dps)

    A = exact_A_coefficients(v, L, K_MAX)

    print("\n" + "-" * 78)
    print(f"N = {N}")
    print("-" * 78)

    print(f"{'k':>3} {'A_k':>26} {'D_k=T^(2k)(0)':>26}")

    for k in range(K_MAX + 1):
        D = endpoint_derivative(v, k, L)
        print(f"{k:3d} {mp.nstr(A[k], 16):>26} {mp.nstr(D, 16):>26}")

    # The cross term in A_k is compared with a full double sum.
    # Since
    #   sum_{m<n} v_m v_n H_k(m,n)
    # = 1/2 [sum_{m,n} v_m v_n H_k(m,n)
    #         - (k+1) sum_m m^(2k) v_m^2],
    # the diagonal part cancels a substantial part of the first term
    # in A_k.  We expose the resulting simplification.
    print("\nExact cross-term reorganisation:")
    print(f"{'k':>3} {'cross':>24} {'Q_k':>24} {'moment-product':>24}")

    for k in range(1, K_MAX + 1):
        cross = sum(
            v[m] * v[n] * H_k(m, n, k)
            for m in range(1, N + 1)
            for n in range(m + 1, N + 1)
        )
        Q = endpoint_quadratic_sum(v, k)
        P = endpoint_moment_quadratic(v, k)
        print(
            f"{k:3d} "
            f"{mp.nstr(cross, 14):>24} "
            f"{mp.nstr(Q, 14):>24} "
            f"{mp.nstr(P, 14):>24}"
        )

    # A particularly useful residual: compare Q_k with the symmetric
    # moment-product sum.  If it vanishes, the entire coefficient hierarchy
    # has a much simpler endpoint-moment representation.
    print("\nQ_k - sum_j M_(2j) M_(2k-2j):")
    for k in range(1, K_MAX + 1):
        residual = endpoint_quadratic_sum(v, k) - endpoint_moment_quadratic(v, k)
        print(f"k={k}: {mp.nstr(residual, 20)}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    mp.mp.dps = dps

    print("=" * 78)
    print("CELL 37 — REORGANISE THE ARCHIMEDEAN TAIL COEFFICIENTS")
    print("=" * 78)
    print(f"c={c}, T={T}, dps={dps}, N_max={SURVEY['N']}")
    print(f"K_MAX={K_MAX}")

    print_polynomial_identity_probe()

    for N in N_VALUES:
        survey_one_N(N)

    print("\n" + "=" * 78)
    print("END OF CELL 37")
    print("=" * 78)
    print(
        "Interpretation: do not infer an N->infinity law from this cell. "
        "The objective is to discover an exact endpoint-moment "
        "reorganisation of the finite-N coefficients."
    )

