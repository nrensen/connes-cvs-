"""
CELL 36 — EXACT ARCHIMEDEAN TAIL COEFFICIENTS

Purpose
-------
Derive and test the complete algebraic asymptotic hierarchy of the
Fourier-side quadratic kernel used in the Archimedean calculation.

For large r, write

    K_fourier(v, r, L)
      = (1 - cos(r L)) R_v(r)

and

    R_v(r)
      ~ A_0(v)/r^2
       + A_1(v)/r^4
       + A_2(v)/r^6
       + ...

Cell 32 identified the leading coefficient A_0.
Cell 35 investigated endpoint derivatives and the associated spectral
moments.  This cell connects those two observations analytically.

The coefficients are derived directly from the exact finite Fourier
formula.  No numerical integration and no symbolic algebra package
are required.

For

    kappa = 2*pi/L,
    M_(2k) = sum_{m=1}^N m^(2k) v_m,

define, for k >= 0,

    H_k(m,n)
      = (n^(2k+2) - m^(2k+2))/(n^2-m^2)
      = sum_{j=0}^k n^(2(k-j)) m^(2j).

Then

    A_0
      = (2/L) v_0^2
        + (4/L) sum_m v_m^2
        + (4*sqrt(2)/L) v_0 M_0
        + (8/L) sum_{m<n} v_m v_n,

where

    M_0 = sum_m v_m.

Equivalently,

    A_0 = 2/L * T_v(0)^2.

For k >= 1,

    A_k
      = [4(k+1)/L] kappa^(2k) sum_m m^(2k) v_m^2
        + [4*sqrt(2)/L] kappa^(2k) v_0 M_(2k)
        + [4/pi] kappa^(2k+1)
            sum_{m<n} v_m v_n H_k(m,n).

The first two terms are directly connected to endpoint even
derivatives

    T_v^(2k)(0)
      = sqrt(2) (-1)^k kappa^(2k) M_(2k).

The principal numerical checks are:

1. exact A_0 agrees with 2*T(0)^2/L;
2. exact A_1 and A_2 agree with finite-r asymptotic extraction;
3. the remainder after subtracting successive terms decays with the
   predicted next power of r;
4. the exact coefficients are evaluated for the same ground-state
   survey used by Cell 35.

This cell deliberately makes no claim about an N -> infinity law.
"""

import mpmath as mp

from cell import get_ground_state, SURVEY_GROUND_STATE


# ============================================================================
# SURVEY CONFIGURATION
# ============================================================================

SURVEY = SURVEY_GROUND_STATE

c = SURVEY["c"]
T_ground = SURVEY["T"]
GROUND_DPS = SURVEY["dps"]

N_VALUES = range(1, SURVEY["N"] + 1)

# Number of asymptotic coefficients:
#
# A_0/r^2 + A_1/r^4 + A_2/r^6 + A_3/r^8 + ...
#
K_MAX = 3

# Large-r sample points used for numerical extraction.
# They should avoid obvious zeros of 1-cos(r L).
R_VALUES = (mp.mpf(80), mp.mpf(120), mp.mpf(180), mp.mpf(260))

mp.mp.dps = GROUND_DPS

L = mp.log(c)
kappa = 2 * mp.pi / L


# ============================================================================
# HEADER
# ============================================================================

print("=" * 78)
print("CELL 36 — EXACT ARCHIMEDEAN TAIL COEFFICIENTS")
print("=" * 78)
print()

print(f"c = {c}")
print(f"T_ground = {T_ground}")
print(f"ground dps = {GROUND_DPS}")
print(f"working dps = {mp.mp.dps}")
print(f"N range = {N_VALUES.start} ... {N_VALUES.stop - 1}")
print(f"K_MAX = {K_MAX}")
print()

print(f"L = {mp.nstr(L, 30)}")
print(f"kappa = 2*pi/L = {mp.nstr(kappa, 30)}")
print()


# ============================================================================
# ENDPOINT DATA
# ============================================================================

def endpoint_even_moments(v, K_MAX):
    """
    Return

        M_(2k) = sum_m m^(2k) v_m

    for k = 0,...,K_MAX.

    M_0 is included because it gives T_v(0).
    """
    N = len(v) - 1

    return [
        mp.fsum(
            (mp.mpf(m) ** (2 * k)) * v[m]
            for m in range(1, N + 1)
        )
        for k in range(K_MAX + 1)
    ]


def endpoint_even_derivatives(v, K_MAX):
    """
    Return

        T(0), T''(0), T^(4)(0), ..., T^(2K_MAX)(0).
    """
    moments = endpoint_even_moments(v, K_MAX)

    values = [
        v[0] + mp.sqrt(2) * moments[0]
    ]

    for k in range(1, K_MAX + 1):
        values.append(
            mp.sqrt(2)
            * (-1) ** k
            * kappa ** (2 * k)
            * moments[k]
        )

    return values


# ============================================================================
# EXACT ASYMPTOTIC COEFFICIENTS
# ============================================================================

def H_polynomial(k, m, n):
    """
    Return

        H_k(m,n)
          = (n^(2k+2)-m^(2k+2))/(n^2-m^2)

          = sum_{j=0}^k n^(2(k-j)) m^(2j).

    The polynomial form avoids subtraction of nearly equal quantities.
    """
    return mp.fsum(
        (mp.mpf(n) ** (2 * (k - j)))
        * (mp.mpf(m) ** (2 * j))
        for j in range(k + 1)
    )


def tail_coefficients(v, K_MAX):
    """
    Return [A_0, ..., A_KMAX] in

        R_v(r) ~ sum_k A_k / r^(2k+2).

    The formula is obtained by expanding the exact finite Fourier
    representation algebraically in powers of 1/r^2.
    """
    N = len(v) - 1
    moments = endpoint_even_moments(v, K_MAX)

    coefficients = []

    # ------------------------------------------------------------------
    # A_0
    #
    # This is written in the expanded canonical form so that its
    # equality with 2*T(0)^2/L is independently testable.
    # ------------------------------------------------------------------
    A0 = (
        2 * v[0] ** 2 / L
        + 4 * mp.fsum(
            v[m] ** 2
            for m in range(1, N + 1)
        ) / L
        + 4 * mp.sqrt(2) * v[0] * moments[0] / L
        + 8 * mp.fsum(
            v[m] * v[n]
            for m in range(1, N)
            for n in range(m + 1, N + 1)
        ) / L
    )

    coefficients.append(A0)

    # ------------------------------------------------------------------
    # A_k, k >= 1
    # ------------------------------------------------------------------
    for k in range(1, K_MAX + 1):
        moment = moments[k]

        diagonal = mp.fsum(
            (mp.mpf(m) ** (2 * k)) * v[m] ** 2
            for m in range(1, N + 1)
        )

        off_diagonal = mp.fsum(
            v[m]
            * v[n]
            * H_polynomial(k, m, n)
            for m in range(1, N)
            for n in range(m + 1, N + 1)
        )

        coefficient = (
            (4 * (k + 1) / L)
            * kappa ** (2 * k)
            * diagonal
            + (4 * mp.sqrt(2) / L)
            * kappa ** (2 * k)
            * v[0]
            * moment
            + (4 / mp.pi)
            * kappa ** (2 * k + 1)
            * off_diagonal
        )

        coefficients.append(coefficient)

    return coefficients


def leading_coefficient_from_endpoint(v):
    """
    A_0 = 2*T_v(0)^2/L.
    """
    T0 = v[0] + mp.sqrt(2) * mp.fsum(
        v[m]
        for m in range(1, len(v))
    )

    return 2 * T0 ** 2 / L


# ============================================================================
# EXACT REDUCED KERNEL R_v(r)
# ============================================================================

def reduced_S(m, r, L):
    """
    S_m(r) / (1-cos(rL)).

    For integer m >= 1,

        S_m(r)
          = (1-cos(rL)) * a_m/(a_m^2-r^2).

    The oscillatory factor is removed analytically.
    """
    a = kappa * m
    return a / (a * a - r * r)


def reduced_C(m, r, L):
    """
    C_m(r) / (1-cos(rL)).

    For integer m >= 0,

        C_m(r)
          = (1-cos(rL))/L
            * (r^2+a_m^2)/(r^2-a_m^2)^2.
    """
    a = kappa * m
    return (
        (r * r + a * a)
        / (L * (r * r - a * a) ** 2)
    )


def reduced_K(v, r, L):
    """
    Return

        R_v(r) = K_fourier(v,r,L)/(1-cos(rL))

    directly, with the common oscillatory factor removed.

    This is algebraically the same expression as K_fourier, but is
    substantially better suited to asymptotic coefficient extraction.
    """
    r = mp.mpf(r)
    L = mp.mpf(L)
    N = len(v) - 1

    total = mp.mpf(0)

    # Diagonal C terms.
    for m in range(0, N + 1):
        total += (
            2
            * v[m] ** 2
            * reduced_C(m, r, L)
        )

    # S_m terms.
    S = {
        m: reduced_S(m, r, L)
        for m in range(1, N + 1)
    }

    # Positive-mode diagonal S term.
    total -= (
        mp.fsum(
            v[m] ** 2 * S[m] / m
            for m in range(1, N + 1)
        )
        / mp.pi
    )

    # Zero-mode / positive-mode term.
    total -= (
        2
        * mp.sqrt(2)
        * v[0]
        * mp.fsum(
            v[m] * S[m] / m
            for m in range(1, N + 1)
        )
        / mp.pi
    )

    # Positive-mode off-diagonal terms.
    total += (
        4
        * mp.fsum(
            v[m]
            * v[n]
            * (
                m * S[m]
                - n * S[n]
            )
            / (n * n - m * m)
            for m in range(1, N)
            for n in range(m + 1, N + 1)
        )
        / mp.pi
    )

    return total


# ============================================================================
# COLLECT GROUND STATES
# ============================================================================

results = []

print("=" * 78)
print("GROUND-STATE ASYMPTOTIC DATA")
print("=" * 78)
print()

for N in N_VALUES:
    lam, v, meta = get_ground_state(
        c=c,
        N=N,
        T=T_ground,
        dps=GROUND_DPS,
        verbose=False,
    )

    derivatives = endpoint_even_derivatives(v, K_MAX)
    coefficients = tail_coefficients(v, K_MAX)
    A0_endpoint = leading_coefficient_from_endpoint(v)

    results.append({
        "N": N,
        "lambda": lam,
        "v": v,
        "derivatives": derivatives,
        "coefficients": coefficients,
        "A0_endpoint": A0_endpoint,
    })


# ============================================================================
# CHECK 1 — A_0 = 2*T(0)^2/L
# ============================================================================

print("=" * 78)
print("CHECK 1 — LEADING COEFFICIENT")
print("=" * 78)
print()

print(
    "Compare the expanded A_0 with"
)
print(
    "    A_0 = 2*T(0)^2/L."
)
print()

print(
    f"{'N':>3}"
    f"{'A0 expanded':>26}"
    f"{'A0 endpoint':>26}"
    f"{'absolute error':>20}"
)
print("-" * 78)

max_A0_error = mp.mpf(0)

for row in results:
    A0 = row["coefficients"][0]
    A0_endpoint = row["A0_endpoint"]
    error = abs(A0 - A0_endpoint)

    max_A0_error = max(max_A0_error, error)

    print(
        f"{row['N']:3d}"
        f"{mp.nstr(A0, 18):>26}"
        f"{mp.nstr(A0_endpoint, 18):>26}"
        f"{mp.nstr(error, 10):>20}"
    )

print()
print(
    "Maximum absolute A_0 identity error:"
)
print(mp.nstr(max_A0_error, 30))
print()


# ============================================================================
# CHECK 2 — EXACT COEFFICIENTS
# ============================================================================

print("=" * 78)
print("CHECK 2 — EXACT ASYMPTOTIC COEFFICIENTS")
print("=" * 78)
print()

for k in range(K_MAX + 1):
    print(f"A_{k} in R_v(r) ~ A_0/r^2 + A_1/r^4 + ...")
    print()

    print(
        f"{'N':>3}"
        f"{'A_' + str(k):>30}"
    )
    print("-" * 38)

    for row in results:
        print(
            f"{row['N']:3d}"
            f"{mp.nstr(row['coefficients'][k], 22):>30}"
        )

    print()


# ============================================================================
# CHECK 3 — ENDPOINT DERIVATIVE CONNECTION
# ============================================================================

print("=" * 78)
print("CHECK 3 — ENDPOINT DERIVATIVE / MOMENT CONNECTION")
print("=" * 78)
print()

print(
    "For k >= 1,"
)
print(
    "    T^(2k)(0) = sqrt(2)*(-1)^k*kappa^(2k)*M_(2k)."
)
print()

print(
    f"{'N':>3}"
    f"{'T(0)':>20}"
    f"{'T''(0)':>20}"
    f"{'T^(4)(0)':>20}"
    f"{'T^(6)(0)':>20}"
)
print("-" * 78)

for row in results:
    d = row["derivatives"]

    print(
        f"{row['N']:3d}"
        f"{mp.nstr(d[0], 12):>20}"
        f"{mp.nstr(d[1], 12):>20}"
        f"{mp.nstr(d[2], 12):>20}"
        f"{mp.nstr(d[3], 12):>20}"
    )

print()


# ============================================================================
# CHECK 4 — NUMERICAL EXTRACTION OF A_0, A_1, A_2
# ============================================================================

print("=" * 78)
print("CHECK 4 — LARGE-r COEFFICIENT EXTRACTION")
print("=" * 78)
print()

print(
    "At each r, form"
)
print(
    "    E_0(r) = r^2 R(r),"
)
print(
    "    E_1(r) = r^4 [R(r)-A_0/r^2],"
)
print(
    "    E_2(r) = r^6 [R(r)-A_0/r^2-A_1/r^4]."
)
print()
print(
    "These should approach A_0, A_1, A_2 respectively."
)
print()

for row in results:
    N = row["N"]
    v = row["v"]
    A = row["coefficients"]

    print("-" * 78)
    print(f"N = {N}")
    print()

    print(
        f"{'r':>10}"
        f"{'E0(r)':>24}"
        f"{'E1(r)':>24}"
        f"{'E2(r)':>24}"
    )
    print("-" * 78)

    for r in R_VALUES:
        R = reduced_K(v, r, L)

        E0 = r ** 2 * R
        E1 = r ** 4 * (
            R
            - A[0] / r ** 2
        )
        E2 = r ** 6 * (
            R
            - A[0] / r ** 2
            - A[1] / r ** 4
        )

        print(
            f"{mp.nstr(r, 8):>10}"
            f"{mp.nstr(E0, 16):>24}"
            f"{mp.nstr(E1, 16):>24}"
            f"{mp.nstr(E2, 16):>24}"
        )

    print()
    print(
        "Exact coefficients:"
    )
    print(
        f"  A_0 = {mp.nstr(A[0], 24)}"
    )
    print(
        f"  A_1 = {mp.nstr(A[1], 24)}"
    )
    print(
        f"  A_2 = {mp.nstr(A[2], 24)}"
    )
    print()


# ============================================================================
# CHECK 5 — RESIDUAL ORDER
# ============================================================================

print("=" * 78)
print("CHECK 5 — RESIDUAL ORDER")
print("=" * 78)
print()

print(
    "After subtracting A_0, A_1, A_2,"
)
print(
    "the remaining quantity should behave as O(r^-8)."
)
print()

for row in results:
    N = row["N"]
    v = row["v"]
    A = row["coefficients"]

    print(f"N = {N}")

    previous = None

    for r in R_VALUES:
        R = reduced_K(v, r, L)

        residual = (
            R
            - A[0] / r ** 2
            - A[1] / r ** 4
            - A[2] / r ** 6
        )

        scaled = r ** 8 * residual

        if previous is None:
            ratio_text = "-"
        else:
            ratio_text = mp.nstr(
                abs(scaled / previous),
                12,
            )

        print(
            f"  r = {mp.nstr(r, 8):>8}"
            f"   r^8 residual = {mp.nstr(scaled, 18):>22}"
            f"   successive ratio = {ratio_text:>14}"
        )

        previous = scaled

    print()


# ============================================================================
# SUMMARY
# ============================================================================

print("=" * 78)
print("SUMMARY")
print("=" * 78)
print()

print(
    "The exact finite Fourier representation admits an algebraic"
)
print(
    "large-r expansion in even inverse powers of r."
)
print()

print(
    "The leading coefficient is exactly"
)
print(
    "    A_0 = 2*T(0)^2/L."
)
print()

print(
    "For k >= 1, A_k depends on the spectral moment"
)
print(
    "    M_(2k) = sum_m m^(2k) v_m,"
)
print(
    "which is proportional to T^(2k)(0)."
)
print()

print(
    "Thus Cell 35's endpoint-derivative survey is directly probing"
)
print(
    "the coefficients of the Archimedean tail expansion."
)
print()

print(
    "The next analytical question is whether the observed suppression"
)
print(
    "of these endpoint moments with N is strong enough to imply a"
)
print(
    "quantitative N-dependent tail bound."
)
print()

print(
    "No N -> infinity asymptotic law is fitted or assumed here."
)
print()
