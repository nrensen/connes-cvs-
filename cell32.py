# ============================================================
# cell32.py — ANALYTICAL TAIL STRUCTURE OF K_fourier
#
# Purpose
# -------
# Reduce the exact finite-N K_fourier expression using
#
#     a_m = 2*pi*m/L
#
# and the identity
#
#     a_m L = 2*pi*m.
#
# This shows that every S_m(r) and C_m(r) contains the same
# oscillatory factor
#
#     1 - cos(r L).
#
# The cell then:
#
#   1. Constructs the exact reduced rational factor R_v(r).
#   2. Verifies
#
#        K_fourier(v,r,L)
#          = (1-cos(rL)) R_v(r)
#
#      numerically at high precision.
#   3. Computes the predicted leading coefficient
#
#        A = (2/L) *
#            (v_0 + sqrt(2) * sum_{m>=1} v_m)^2
#
#      so that
#
#        R_v(r) = A/r^2 + O(r^-4).
#
#   4. Tests the large-r behaviour of r^2 R_v(r).
#   5. Tests the residual
#
#        r^4 [R_v(r) - A/r^2]
#
#      to expose the next coefficient numerically.
#
# No Archimedean integration is performed here.
# ============================================================

import time

import mpmath as mp

from cell import (
    DEFAULT_DPS,
    FORENSIC_GROUND_STATE,
    K_fourier,
    compute_L,
    get_ground_state,
)


# ============================================================
# PARAMETERS
# ============================================================

DPS = 100

c = FORENSIC_GROUND_STATE["c"]
N = FORENSIC_GROUND_STATE["N"]
T_GROUND = FORENSIC_GROUND_STATE["T"]
GROUND_DPS = FORENSIC_GROUND_STATE["dps"]

# Large-r test points.
#
# These are deliberately not tied to zeros of 1-cos(rL).
# The reduced rational factor R_v(r) is independent of that
# oscillatory factor, so this avoids any loss of significance
# associated with dividing by a tiny 1-cos(rL).
R_TEST = [
    mp.mpf("100"),
    mp.mpf("300"),
    mp.mpf("1000"),
    mp.mpf("3000"),
    mp.mpf("10000"),
    mp.mpf("30000"),
    mp.mpf("100000"),
]


# ============================================================
# GROUND STATE
# ============================================================

def load_forensic_state():
    """
    Load the fixed forensic ground state at the requested
    working precision.

    The ground state itself is generated/cached at the fixed
    forensic generation precision, while all calculations in
    this cell occur at DPS working precision.
    """
    lambda_min, v_star, metadata = get_ground_state(
        c=c,
        N=N,
        T=T_GROUND,
        dps=GROUND_DPS,
    )

    L = compute_L(c)

    return lambda_min, v_star, metadata, L


# ============================================================
# EXACT REDUCED MODE FUNCTIONS
# ============================================================

def S_reduced(m, r, L):
    """
    Exact reduced form of S_m(r).

    For m >= 1,

        S_m(r)
          = (1-cos(rL)) *
            a_m / (a_m^2-r^2),

    where

        a_m = 2*pi*m/L.

    For m = 0, S_0 = 0.
    """
    m = int(m)
    r = mp.mpf(r)
    L = mp.mpf(L)

    if m == 0:
        return mp.mpf(0)

    a_m = 2 * mp.pi * m / L

    return (
        a_m / (a_m * a_m - r * r)
    )


def C_reduced(m, r, L):
    """
    Exact reduced rational factor of C_m(r).

    C_m(r)
      = (1-cos(rL)) * Cbar_m(r),

    with

      Cbar_m(r)
        = (a_m^2+r^2)
          / [L (a_m^2-r^2)^2].

    This formula also applies to m=0.
    """
    m = int(m)
    r = mp.mpf(r)
    L = mp.mpf(L)

    a_m = 2 * mp.pi * m / L

    return (
        (a_m * a_m + r * r)
        / (
            L
            * (a_m * a_m - r * r) ** 2
        )
    )


# ============================================================
# REDUCED KERNEL
# ============================================================

def K_reduced(v, r, L):
    """
    Exact rational factor R_v(r) defined by

        K_fourier(v,r,L)
          = (1-cos(rL)) * K_reduced(v,r,L).

    This is obtained directly by substituting the exact
    reduced S_m and C_m expressions into the current
    K_fourier formula.
    """
    r = mp.mpf(r)
    L = mp.mpf(L)

    N_local = len(v) - 1

    # --------------------------------------------------------
    # Diagonal terms
    # --------------------------------------------------------

    diag = mp.mpf(0)

    for m in range(0, N_local + 1):
        diag += (
            v[m]
            * v[m]
            * C_reduced(
                m,
                r,
                L,
            )
        )

    total = 2 * diag

    # --------------------------------------------------------
    # S_m terms
    # --------------------------------------------------------

    S = {
        m: S_reduced(
            m,
            r,
            L,
        )
        for m in range(1, N_local + 1)
    }

    off_diag = mp.mpf(0)

    for m in range(1, N_local + 1):
        off_diag += (
            v[m]
            * v[m]
            * S[m]
            / m
        )

    total -= off_diag / mp.pi

    # --------------------------------------------------------
    # Zero-positive terms
    # --------------------------------------------------------

    off_zero = mp.mpf(0)

    for m in range(1, N_local + 1):
        off_zero += (
            v[m]
            * S[m]
            / m
        )

    total -= (
        2
        * mp.sqrt(2)
        * v[0]
        * off_zero
        / mp.pi
    )

    # --------------------------------------------------------
    # Positive-positive terms
    # --------------------------------------------------------

    off = mp.mpf(0)

    for m in range(1, N_local):
        for n in range(m + 1, N_local + 1):
            off += (
                v[m]
                * v[n]
                * (
                    (
                        m * S[m]
                        - n * S[n]
                    )
                    / (n * n - m * m)
                )
            )

    total += 4 * off / mp.pi

    return total


# ============================================================
# LEADING ASYMPTOTIC COEFFICIENT
# ============================================================

def leading_coefficient(v, L):
    """
    Predicted coefficient A in

        K_reduced(v,r,L)
          = A/r^2 + O(r^-4).

    A = (2/L)
        * (v_0 + sqrt(2) * sum_{m>=1} v_m)^2.
    """
    N_local = len(v) - 1

    endpoint_value = v[0]

    for m in range(1, N_local + 1):
        endpoint_value += (
            mp.sqrt(2) * v[m]
        )

    A = (
        2
        / L
        * endpoint_value
        * endpoint_value
    )

    return A, endpoint_value


# ============================================================
# FORMATTING
# ============================================================

def fmt(x, digits=20):
    return mp.nstr(x, digits)


def relerr(a, b):
    if b == 0:
        return abs(a - b)
    return abs(a - b) / abs(b)


# ============================================================
# MAIN ANALYSIS
# ============================================================

def main():

    start = time.perf_counter()

    print("=" * 72)
    print("CELL 32 — ANALYTICAL TAIL STRUCTURE")
    print("=" * 72)

    print()
    print(f"c              = {c}")
    print(f"N              = {N}")
    print(f"T_ground       = {T_GROUND}")
    print(f"ground dps     = {GROUND_DPS}")
    print(f"working dps    = {DPS}")

    with mp.workdps(DPS):

        # ----------------------------------------------------
        # Load forensic ground state
        # ----------------------------------------------------

        lambda_min, v_star, metadata, L = (
            load_forensic_state()
        )

        print()
        print("GROUND STATE")
        print("-" * 72)
        print(f"lambda_min     = {fmt(lambda_min, 30)}")
        print(f"L              = {fmt(L, 30)}")

        norm = mp.sqrt(
            mp.fdot(
                v_star,
                v_star,
            )
        )

        print(f"||v||          = {fmt(norm, 30)}")
        print(
            f"norm error     = "
            f"{fmt(abs(norm - 1), 10)}"
        )

        # ----------------------------------------------------
        # Leading coefficient
        # ----------------------------------------------------

        A, endpoint_value = leading_coefficient(
            v_star,
            L,
        )

        print()
        print("LEADING COEFFICIENT")
        print("-" * 72)

        print(
            "T_v(0)         = "
            f"{fmt(endpoint_value, 40)}"
        )

        print(
            "A              = "
            f"{fmt(A, 40)}"
        )

        print()
        print(
            "Predicted form:"
        )
        print(
            "  K_fourier(r)"
            " = (1-cos(rL))"
            " * [ A/r^2 + O(r^-4) ]"
        )

        # ----------------------------------------------------
        # Exact factorisation test
        # ----------------------------------------------------

        print()
        print("EXACT FACTORISATION TEST")
        print("-" * 72)

        phase_points = [
            mp.mpf("17.3"),
            mp.mpf("37.1"),
            mp.mpf("103.7"),
            mp.mpf("317.2"),
        ]

        print(
            f"{'r':>12}"
            f"  {'K_fourier':>25}"
            f"  {'reconstructed':>25}"
            f"  {'relative error':>12}"
        )

        print("-" * 72)

        for r in phase_points:

            phase = 1 - mp.cos(r * L)

            K_direct = K_fourier(
                v_star,
                r,
                L,
            )

            R = K_reduced(
                v_star,
                r,
                L,
            )

            K_reconstructed = phase * R

            error = relerr(
                K_reconstructed,
                K_direct,
            )

            print(
                f"{fmt(r, 8):>12}"
                f"  {fmt(K_direct, 25):>25}"
                f"  {fmt(K_reconstructed, 25):>25}"
                f"  {fmt(error, 8):>12}"
            )

        # ----------------------------------------------------
        # Large-r leading asymptotic test
        # ----------------------------------------------------

        print()
        print("LARGE-r TEST")
        print("-" * 72)

        print(
            f"{'r':>12}"
            f"  {'r^2 R(r)':>30}"
            f"  {'difference from A':>20}"
            f"  {'relative':>12}"
        )

        print("-" * 72)

        for r in R_TEST:

            R = K_reduced(
                v_star,
                r,
                L,
            )

            scaled = r * r * R
            difference = scaled - A

            print(
                f"{fmt(r, 10):>12}"
                f"  {fmt(scaled, 30):>30}"
                f"  {fmt(difference, 20):>20}"
                f"  {fmt(relerr(scaled, A), 8):>12}"
            )

        # ----------------------------------------------------
        # Next asymptotic coefficient
        # ----------------------------------------------------
        #
        # If
        #
        #   R(r) = A/r^2 + B/r^4 + O(r^-6),
        #
        # then
        #
        #   r^4 [R(r)-A/r^2] -> B.
        #
        # We do not assume this merely from the numerical
        # result: the purpose here is to expose the candidate
        # coefficient so that we can derive it analytically
        # next.
        # ----------------------------------------------------

        print()
        print("NEXT-ORDER COEFFICIENT PROBE")
        print("-" * 72)

        print(
            f"{'r':>12}"
            f"  {'r^4 [R-A/r^2]':>30}"
            f"  {'successive change':>20}"
        )

        print("-" * 72)

        previous_B = None

        for r in R_TEST:

            R = K_reduced(
                v_star,
                r,
                L,
            )

            B_est = (
                r ** 4
                * (
                    R
                    - A / r ** 2
                )
            )

            if previous_B is None:
                change = mp.nan
            else:
                change = B_est - previous_B

            print(
                f"{fmt(r, 10):>12}"
                f"  {fmt(B_est, 30):>30}"
                f"  {fmt(change, 20):>20}"
            )

            previous_B = B_est

        # ----------------------------------------------------
        # Direct comparison at very large r
        # ----------------------------------------------------

        print()
        print("TAIL-SCALE SUMMARY")
        print("-" * 72)

        r = R_TEST[-1]

        R = K_reduced(
            v_star,
            r,
            L,
        )

        leading = A / r ** 2
        residual = R - leading

        print(
            f"r              = {fmt(r, 20)}"
        )
        print(
            f"R(r)           = {fmt(R, 40)}"
        )
        print(
            f"A/r^2          = {fmt(leading, 40)}"
        )
        print(
            f"residual        = {fmt(residual, 40)}"
        )
        print(
            f"residual/lead  = "
            f"{fmt(abs(residual / leading), 20)}"
        )

    elapsed = time.perf_counter() - start

    print()
    print("=" * 72)
    print(
        f"Elapsed time: {elapsed:.3f} s"
    )
    print("=" * 72)


if __name__ == "__main__":
    main()
