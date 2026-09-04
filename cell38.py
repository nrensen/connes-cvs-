"""
CELL 38 — CLOSED FORM FOR THE ARCHIMEDEAN TAIL HIERARCHY

Cell 37 established the exact identity

    Q_k = sum_{j=0}^k M_(2j) M_(2k-2j),

where M_p = sum_{m=1}^N m^p v_m.

This cell performs the algebraic cancellation in A_k completely.

For k >= 1 it proves the finite-N identity

    A_k = (4/L) kappa^(2k)
          [ Q_k + sqrt(2) v_0 M_(2k) ]

and then reorganises it using

    T_0 = T_v(0) = v_0 + sqrt(2) M_0

to obtain the much cleaner endpoint-jet formula

    A_k = (4/L) (-1)^k
          [ T_0 D_k + (1/2) sum_{j=1}^{k-1} D_j D_(k-j) ],

where

    D_k = T_v^(2k)(0),  k >= 1.

Equivalently, with D_0 := T_v(0),

    A_k = (2/L) (-1)^k
          sum_{j=0}^k D_j D_(k-j),     k >= 1.

Thus every higher asymptotic coefficient is exactly a quadratic
convolution of the even endpoint derivatives.

This is an exact finite-N algebraic result.  No fitting and no
quadrature are used.
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
# Moments and endpoint jets
# ---------------------------------------------------------------------------

def moment(v, p: int):
    """M_p = sum_{m=1}^N m^p v_m."""
    return sum((mp.mpf(m) ** p) * v[m] for m in range(1, len(v)))


def endpoint_jet(v, k: int):
    """
    D_k = T_v^(2k)(0).

    D_0 = T_v(0).
    For k >= 1:
        D_k = sqrt(2) (-1)^k kappa^(2k) M_(2k).
    """
    if k == 0:
        return v[0] + mp.sqrt(2) * moment(v, 0)

    kappa = 2 * mp.pi / L
    return (
        mp.sqrt(2)
        * (-1) ** k
        * (kappa ** (2 * k))
        * moment(v, 2 * k)
    )


# ---------------------------------------------------------------------------
# Exact coefficient from the Cell 36 formula
# ---------------------------------------------------------------------------

def H_k(m: int, n: int, k: int):
    """H_k(m,n) = sum_{j=0}^k n^(2(k-j)) m^(2j)."""
    return sum(
        (mp.mpf(n) ** (2 * (k - j)))
        * (mp.mpf(m) ** (2 * j))
        for j in range(k + 1)
    )


def exact_A_original(v, L, K=K_MAX):
    """
    Original exact finite-N coefficient formula used in Cell 36.

    This is retained as the independent starting point for the
    algebraic reorganisation.
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

        out.append(
            (4 * (k + 1) / L) * (kappa ** (2 * k)) * diag
            + (4 * mp.sqrt(2) / L) * (kappa ** (2 * k)) * mixed0
            + (4 / mp.pi) * (kappa ** (2 * k + 1)) * cross
        )

    return out


# ---------------------------------------------------------------------------
# Algebraically reduced coefficient
# ---------------------------------------------------------------------------

def exact_A_moment_form(v, L, k: int):
    """
    A_k after the diagonal/cross-term cancellation:

        A_k = (4/L) kappa^(2k)
              [sum_{j=0}^k M_(2j) M_(2k-2j)
               + sqrt(2) v0 M_(2k)].

    Valid for k >= 1.
    """
    assert k >= 1

    kappa = 2 * mp.pi / L
    Q = sum(
        moment(v, 2 * j) * moment(v, 2 * (k - j))
        for j in range(k + 1)
    )

    return (
        (4 / L)
        * (kappa ** (2 * k))
        * (Q + mp.sqrt(2) * v[0] * moment(v, 2 * k))
    )


def exact_A_jet_form(v, L, k: int):
    """
    Endpoint-jet form:

        A_k = (2/L) (-1)^k
              sum_{j=0}^k D_j D_(k-j),

    where D_j = T_v^(2j)(0), D_0 = T_v(0).
    """
    assert k >= 1

    return (
        (2 / L)
        * (-1) ** k
        * sum(
            endpoint_jet(v, j) * endpoint_jet(v, k - j)
            for j in range(k + 1)
        )
    )


def exact_A_jet_expanded(v, L, k: int):
    """
    Same identity written with the endpoint term isolated:

        A_k = (4/L) (-1)^k
              [D_0 D_k + 1/2 sum_{j=1}^{k-1} D_j D_(k-j)].
    """
    assert k >= 1

    D0 = endpoint_jet(v, 0)
    inner = sum(
        endpoint_jet(v, j) * endpoint_jet(v, k - j)
        for j in range(1, k)
    )

    return (4 / L) * (-1) ** k * (D0 * endpoint_jet(v, k) + inner / 2)


# ---------------------------------------------------------------------------
# Main survey
# ---------------------------------------------------------------------------

def survey_one_N(N: int):
    lam, v, meta = get_ground_state(c=c, N=N, T=T, dps=dps)

    A = exact_A_original(v, L, K_MAX)

    print("\n" + "=" * 78)
    print(f"N = {N}")
    print("=" * 78)

    print(
        f"{'k':>3} "
        f"{'A_original':>24} "
        f"{'A_moment':>24} "
        f"{'A_jet':>24}"
    )

    for k in range(1, K_MAX + 1):
        Am = exact_A_moment_form(v, L, k)
        Aj = exact_A_jet_form(v, L, k)

        print(
            f"{k:3d} "
            f"{mp.nstr(A[k], 16):>24} "
            f"{mp.nstr(Am, 16):>24} "
            f"{mp.nstr(Aj, 16):>24}"
        )

    print("\nIndependent residuals:")
    for k in range(1, K_MAX + 1):
        Am = exact_A_moment_form(v, L, k)
        Aj = exact_A_jet_form(v, L, k)
        Ae = exact_A_jet_expanded(v, L, k)

        print(
            f"k={k}: "
            f"original-moment = {mp.nstr(A[k] - Am, 8)}, "
            f"moment-jet = {mp.nstr(Am - Aj, 8)}, "
            f"jet-expanded = {mp.nstr(Aj - Ae, 8)}"
        )

    print("\nEndpoint jets:")
    for k in range(K_MAX + 1):
        print(
            f"D_{k} = "
            f"{mp.nstr(endpoint_jet(v, k), 20)}"
        )


if __name__ == "__main__":
    mp.mp.dps = dps

    print("=" * 78)
    print("CELL 38 — CLOSED FORM FOR THE ARCHIMEDEAN TAIL HIERARCHY")
    print("=" * 78)
    print(f"c={c}, L=log(c)={mp.nstr(L, 20)}, T={T}, dps={dps}")
    print(f"N_max={SURVEY['N']}, K_MAX={K_MAX}")

    print(
        "\nClaim to test exactly at finite N:\n"
        "  A_k = (2/L)(-1)^k sum_{j=0}^k D_j D_(k-j),  k >= 1."
    )

    for N in N_VALUES:
        survey_one_N(N)

    print("\n" + "=" * 78)
    print("END OF CELL 38")
    print("=" * 78)
    print(
        "If the residuals are at working-precision noise level, the "
        "entire higher tail hierarchy has been reduced to an endpoint-jet "
        "convolution."
    )

