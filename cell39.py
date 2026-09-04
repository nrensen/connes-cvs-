"""
CELL 39 — GENERATING FUNCTION FOR THE ARCHIMEDEAN TAIL HIERARCHY

Cell 38 established, for the reduced finite-N kernel

    K_fourier(v,r,L) = (1 - cos(r L)) R_v(r),

the exact asymptotic coefficients

    R_v(r) ~ sum_{k>=0} A_k / r^(2k+2),

with

    A_0 = (2/L) D_0^2,

    A_k = (2/L) (-1)^k sum_{j=0}^k D_j D_{k-j},  k >= 1,

where

    D_k = T_v^(2k)(0).

Cell 39 asks whether this hierarchy has a compact generating function.

The answer suggested by Cell 38 is

    D(z) = sum_{k>=0} D_k z^k
         = v_0 + sqrt(2) sum_{m=1}^N v_m/(1 + kappa^2 m^2 z),

where kappa = 2 pi/L.

Therefore

    A(z) = sum_{k>=0} A_k z^k
         = (2/L) D(-z)^2

and hence

    A(z)
      = (2/L) [
          v_0 + sqrt(2) sum_m v_m/(1 - kappa^2 m^2 z)
        ]^2.

This cell verifies those identities numerically at high precision and
then compares the resulting coefficient-generating function with the
exact reduced rational kernel through its asymptotic variable z = 1/r^2.

No quadrature, fitting, or large-r numerical extraction is used to discover
the coefficient formula.
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

# Number of exact tail coefficients to test.
K_MAX = 8

# A few non-asymptotic z values, kept safely inside the nearest pole.
Z_VALUES = ("0.000001", "0.00001", "0.0001")


# ---------------------------------------------------------------------------
# Endpoint derivatives and moments
# ---------------------------------------------------------------------------

def moment(v, p: int):
    """M_p = sum_{m=1}^N m^p v_m."""
    return sum(
        (mp.mpf(m) ** p) * v[m]
        for m in range(1, len(v))
    )


def endpoint_derivative(v, k: int):
    """
    D_k = T_v^(2k)(0), with D_0 = T_v(0).
    """
    if k == 0:
        return v[0] + mp.sqrt(2) * moment(v, 0)

    kappa = 2 * mp.pi / L
    return (
        mp.sqrt(2)
        * (-1) ** k
        * kappa ** (2 * k)
        * moment(v, 2 * k)
    )


# ---------------------------------------------------------------------------
# Exact coefficient formula from Cell 38
# ---------------------------------------------------------------------------

def exact_A_coefficients(v, K=K_MAX):
    """
    Exact finite-N tail coefficients from the Cell-38 endpoint-jet formula.
    """
    D = [endpoint_derivative(v, k) for k in range(K + 1)]

    A = []
    for k in range(K + 1):
        convolution = sum(
            D[j] * D[k - j]
            for j in range(k + 1)
        )
        A.append(
            (mp.mpf(2) / L)
            * (-1) ** k
            * convolution
        )

    return D, A


# ---------------------------------------------------------------------------
# Generating functions
# ---------------------------------------------------------------------------

def D_generating(v, z):
    """
    D(z) = sum_{k>=0} D_k z^k.

    For finite N this geometric sum is exact wherever the series converges:

        D(z) = v_0 + sqrt(2) sum_m v_m/(1 + kappa^2 m^2 z).
    """
    kappa = 2 * mp.pi / L
    return v[0] + mp.sqrt(2) * sum(
        v[m] / (1 + (kappa * m) ** 2 * z)
        for m in range(1, len(v))
    )


def A_generating_closed(v, z):
    """
    A(z) = (2/L) D(-z)^2.
    """
    return (mp.mpf(2) / L) * D_generating(v, -z) ** 2


def A_generating_from_coefficients(A, z):
    """Truncated power series sum_{k=0}^K A_k z^k."""
    return sum(
        A[k] * z ** k
        for k in range(len(A))
    )


# ---------------------------------------------------------------------------
# Exact reduced rational kernel
# ---------------------------------------------------------------------------

def reduced_K(v, r):
    """
    Exact R_v(r), obtained from K_fourier after removing the common factor
    (1 - cos(rL)) analytically.

    This reproduces the finite-N rational kernel used in Cells 32 and 36.
    """
    kappa = 2 * mp.pi / L

    def S_bar(m):
        a = kappa * m
        return a / (a * a - r * r)

    def C_bar(m):
        a = kappa * m
        return (r * r + a * a) / (
            L * (r * r - a * a) ** 2
        )

    v0 = v[0]

    value = 2 * v0 * v0 / (L * r * r)

    for m in range(1, len(v)):
        value += 2 * v[m] ** 2 * C_bar(m)
        value -= (
            2 * mp.sqrt(2) * v0 * v[m] / mp.pi
        ) * S_bar(m) / m
        value -= (
            v[m] ** 2 / mp.pi
        ) * S_bar(m) / m

    for m in range(1, len(v)):
        for n in range(m + 1, len(v)):
            value += (
                4 * v[m] * v[n] / mp.pi
                * (
                    m * S_bar(m) - n * S_bar(n)
                )
                / (n * n - m * m)
            )

    return value


# ---------------------------------------------------------------------------
# Direct coefficient extraction from the rational expression
# ---------------------------------------------------------------------------

def rational_A_coefficients(v, K=K_MAX):
    """
    Obtain A_k directly from the elementary large-r expansions of the
    reduced rational kernel.

    This is an independent coefficient calculation: it does NOT use
    endpoint derivatives or the Cell-38 formula.
    """
    kappa = 2 * mp.pi / L
    N = len(v) - 1

    A = []

    # k = 0
    A0 = (
        2 * v[0] ** 2 / L
        + 4 * sum(v[m] ** 2 for m in range(1, N + 1)) / L
        + 4 * mp.sqrt(2) * v[0] * moment(v, 0) / L
        + 8 * sum(
            v[m] * v[n]
            for m in range(1, N + 1)
            for n in range(m + 1, N + 1)
        ) / L
    )
    A.append(A0)

    for k in range(1, K + 1):
        diag = sum(
            (mp.mpf(m) ** (2 * k)) * v[m] ** 2
            for m in range(1, N + 1)
        )

        mixed0 = v[0] * moment(v, 2 * k)

        cross = sum(
            v[m] * v[n] * sum(
                (mp.mpf(n) ** (2 * (k - j)))
                * (mp.mpf(m) ** (2 * j))
                for j in range(k + 1)
            )
            for m in range(1, N + 1)
            for n in range(m + 1, N + 1)
        )

        A.append(
            (4 * (k + 1) / L)
            * kappa ** (2 * k)
            * diag
            + (4 * mp.sqrt(2) / L)
            * kappa ** (2 * k)
            * mixed0
            + (4 / mp.pi)
            * kappa ** (2 * k + 1)
            * cross
        )

    return A


# ---------------------------------------------------------------------------
# Main survey
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    mp.mp.dps = dps

    print("=" * 78)
    print("CELL 39 — GENERATING FUNCTION FOR THE ARCHIMEDEAN TAIL HIERARCHY")
    print("=" * 78)
    print(f"c={c}, L={mp.nstr(L, 20)}, T={T}, dps={dps}")
    print(f"N_max={SURVEY['N']}, K_MAX={K_MAX}")

    print("\nThe predicted exact generating functions are:")
    print("  D(z) = v0 + sqrt(2) sum_m v_m/(1 + kappa^2 m^2 z)")
    print("  A(z) = (2/L) D(-z)^2")

    for N in N_VALUES:
        lam, v, meta = get_ground_state(
            c=c, N=N, T=T, dps=dps
        )

        D, A_jet = exact_A_coefficients(v, K_MAX)
        A_direct = rational_A_coefficients(v, K_MAX)

        print("\n" + "-" * 78)
        print(f"N = {N}")
        print("-" * 78)

        # ---------------------------------------------------------------
        # 1. Direct coefficient identity
        # ---------------------------------------------------------------
        max_err = max(
            abs(A_jet[k] - A_direct[k])
            for k in range(K_MAX + 1)
        )

        print(
            "max |A_jet - A_direct| =",
            mp.nstr(max_err, 8)
        )

        print(
            f"{'k':>3} "
            f"{'A_direct':>25} "
            f"{'A_jet':>25} "
            f"{'abs diff':>14}"
        )

        for k in range(K_MAX + 1):
            diff = abs(A_direct[k] - A_jet[k])
            print(
                f"{k:3d} "
                f"{mp.nstr(A_direct[k], 14):>25} "
                f"{mp.nstr(A_jet[k], 14):>25} "
                f"{mp.nstr(diff, 6):>14}"
            )

        # ---------------------------------------------------------------
        # 2. Generating function D(z)
        #
        # Compare its closed rational expression against a truncated
        # derivative series.  This checks the geometric resummation.
        # ---------------------------------------------------------------
        print("\nD(z) geometric resummation check:")

        for z_text in Z_VALUES:
            z = mp.mpf(z_text)

            D_series = sum(
                D[k] * z ** k
                for k in range(K_MAX + 1)
            )
            D_closed = D_generating(v, z)

            print(
                f"z={z_text:>8}: "
                f"|D_series-D_closed|="
                f"{mp.nstr(abs(D_series - D_closed), 8)}"
            )

        # ---------------------------------------------------------------
        # 3. A(z) = (2/L) D(-z)^2
        #
        # Compare the coefficient series against the closed generating
        # function at small z.  The discrepancy should scale with the
        # first omitted power.
        # ---------------------------------------------------------------
        print("\nA(z) generating-function check:")

        for z_text in Z_VALUES:
            z = mp.mpf(z_text)

            A_series = A_generating_from_coefficients(A_jet, z)
            A_closed = A_generating_closed(v, z)

            print(
                f"z={z_text:>8}: "
                f"|A_series-A_closed|="
                f"{mp.nstr(abs(A_series - A_closed), 8)}"
            )

        # ---------------------------------------------------------------
        # 4. Compare the exact rational kernel with the asymptotic
        #    generating function through z = 1/r^2.
        #
        # Since
        #     R(r) ~ z A(z),  z=1/r^2,
        # the quantity R(r)/z should approach A(z).
        #
        # We deliberately use several r values and report the residual.
        # This is a validation of the generating interpretation, not a
        # coefficient fit.
        # ---------------------------------------------------------------
        print("\nRational-kernel / generating-function check:")

        for r_text in ("40", "80", "120"):
            r = mp.mpf(r_text)
            z = 1 / (r * r)

            R = reduced_K(v, r)
            predicted = z * A_generating_closed(v, z)

            rel = abs(R - predicted) / max(abs(R), mp.mpf("1e-100"))

            print(
                f"r={r_text:>3}: "
                f"R={mp.nstr(R, 14):>20} "
                f"z*A(z)={mp.nstr(predicted, 14):>20} "
                f"rel.err={mp.nstr(rel, 8)}"
            )

    print("\n" + "=" * 78)
    print("END OF CELL 39")
    print("=" * 78)
    print(
        "No N->infinity law has been assumed. "
        "The purpose is to establish the compact generating-function "
        "form of the exact finite-N tail hierarchy."
    )
