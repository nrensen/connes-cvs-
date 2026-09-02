# ============================================================
# cell33.py — N-DEPENDENCE OF THE FOURIER TAIL COEFFICIENTS
#
# Purpose
# -------
# Investigate how the large-r coefficients of
#
#     K_fourier(v,r,L)
#       = (1-cos(rL)) R_v(r)
#
# depend on the Galerkin dimension N.
#
# For finite N we have
#
#     R_v(r) = A/r^2 + B/r^4 + O(r^-6),
#
# where
#
#     A = (2/L) *
#         (v_0 + sqrt(2) * sum_{m=1}^N v_m)^2.
#
# This cell:
#
#   1. Builds the ground state for several N.
#   2. Computes T_v(0).
#   3. Computes A exactly from the finite vector.
#   4. Extracts B numerically from the exact reduced kernel.
#   5. Checks stability of the extracted B with r.
#   6. Reports ||v|| and lambda_min as sanity checks.
#
# No Archimedean tail integration is performed.
#
# IMPORTANT:
# ----------
# The ground states are all generated at the fixed forensic
# precision.  The current working precision is controlled
# independently by WORK_DPS.
# ============================================================

import time

import mpmath as mp

from cell import (
    FORENSIC_GROUND_STATE,
    get_ground_state,
    compute_L,
    K_fourier,
)


# ============================================================
# PARAMETERS
# ============================================================

WORK_DPS = 100

c = FORENSIC_GROUND_STATE["c"]
T_GROUND = FORENSIC_GROUND_STATE["T"]
GROUND_DPS = FORENSIC_GROUND_STATE["dps"]

# Start with modest N values.
#
# The purpose of this first sweep is to identify the trend,
# not to establish the N -> infinity limit.
N_VALUES = [
    2,
    4,
    6,
    8,
    10,
    12,
]

# Large-r values used to extract B.
#
# We need r large compared with all a_m = 2*pi*m/L,
# while avoiding unnecessarily extreme values.
R_VALUES = [
    mp.mpf("1000"),
    mp.mpf("3000"),
    mp.mpf("10000"),
    mp.mpf("30000"),
    mp.mpf("100000"),
]


# ============================================================
# REDUCED KERNEL
# ============================================================

def S_reduced(m, r, L):
    """
    Rational factor in

        S_m(r)
          = (1-cos(rL)) * S_reduced(m,r,L).

    For m >= 1.
    """
    if m == 0:
        return mp.mpf(0)

    a_m = 2 * mp.pi * m / L

    return (
        a_m
        / (a_m * a_m - r * r)
    )


def C_reduced(m, r, L):
    """
    Rational factor in

        C_m(r)
          = (1-cos(rL)) * C_reduced(m,r,L).
    """
    a_m = 2 * mp.pi * m / L

    return (
        (a_m * a_m + r * r)
        / (
            L
            * (a_m * a_m - r * r) ** 2
        )
    )


def K_reduced(v, r, L):
    """
    Exact rational factor R_v(r) defined by

        K_fourier(v,r,L)
          = (1-cos(rL)) * K_reduced(v,r,L).
    """
    N_local = len(v) - 1

    total = mp.mpf(0)

    # --------------------------------------------------------
    # C terms
    # --------------------------------------------------------

    for m in range(N_local + 1):
        total += (
            2
            * v[m] ** 2
            * C_reduced(
                m,
                r,
                L,
            )
        )

    # --------------------------------------------------------
    # S terms
    # --------------------------------------------------------

    S = {
        m: S_reduced(
            m,
            r,
            L,
        )
        for m in range(1, N_local + 1)
    }

    # Diagonal S contribution
    for m in range(1, N_local + 1):
        total -= (
            1
            / mp.pi
            * v[m] ** 2
            * S[m]
            / m
        )

    # v_0 / positive-mode contribution
    for m in range(1, N_local + 1):
        total -= (
            2
            * mp.sqrt(2)
            / mp.pi
            * v[0]
            * v[m]
            * S[m]
            / m
        )

    # Positive-positive contribution
    for m in range(1, N_local):
        for n in range(m + 1, N_local + 1):
            total += (
                4
                / mp.pi
                * v[m]
                * v[n]
                * (
                    m * S[m]
                    - n * S[n]
                )
                / (
                    n * n
                    - m * m
                )
            )

    return total


# ============================================================
# ASYMPTOTIC COEFFICIENT A
# ============================================================

def compute_A(v, L):
    """
    Compute the exact leading coefficient

        A = (2/L) T_v(0)^2,

    where

        T_v(0)
          = v_0 + sqrt(2) * sum_{m>=1} v_m.
    """
    endpoint = v[0]

    for m in range(1, len(v)):
        endpoint += mp.sqrt(2) * v[m]

    A = (
        2
        / L
        * endpoint ** 2
    )

    return A, endpoint


# ============================================================
# EXTRACT B
# ============================================================

def extract_B(v, L, A, r):
    """
    If

        R(r) = A/r^2 + B/r^4 + O(r^-6),

    then

        B = lim_{r->infinity}
            r^4 [R(r) - A/r^2].

    This function gives the finite-r estimator.
    """
    R = K_reduced(
        v,
        r,
        L,
    )

    return (
        r ** 4
        * (
            R
            - A / r ** 2
        )
    )


# ============================================================
# FORMATTERS
# ============================================================

def fmt(x, digits=16):
    return mp.nstr(x, digits)


# ============================================================
# MAIN
# ============================================================

def main():

    start = time.perf_counter()

    print("=" * 78)
    print("CELL 33 — N-DEPENDENCE OF FOURIER TAIL COEFFICIENTS")
    print("=" * 78)

    print()
    print(f"c              = {c}")
    print(f"T_ground       = {T_GROUND}")
    print(f"ground dps     = {GROUND_DPS}")
    print(f"working dps    = {WORK_DPS}")
    print(f"N values       = {N_VALUES}")

    with mp.workdps(WORK_DPS):

        L = compute_L(c)

        print()
        print(f"L = {fmt(L, 30)}")

        # ----------------------------------------------------
        # Summary table
        # ----------------------------------------------------

        print()
        print("GROUND-STATE / LEADING-COEFFICIENT SUMMARY")
        print("-" * 78)

        print(
            f"{'N':>4}"
            f"  {'lambda_min':>25}"
            f"  {'||v||-1':>12}"
            f"  {'T_v(0)':>25}"
            f"  {'A':>25}"
        )

        print("-" * 78)

        states = {}

        for N in N_VALUES:

            print(
                f"Building N={N} ...",
                flush=True,
            )

            lambda_min, v, metadata = (
                get_ground_state(
                    c=c,
                    N=N,
                    T=T_GROUND,
                    dps=GROUND_DPS,
                )
            )

            norm = mp.sqrt(
                mp.fdot(
                    v,
                    v,
                )
            )

            A, endpoint = compute_A(
                v,
                L,
            )

            states[N] = {
                "lambda": lambda_min,
                "v": v,
                "norm": norm,
                "endpoint": endpoint,
                "A": A,
            }

            print(
                f"{N:4d}"
                f"  {fmt(lambda_min, 18):>25}"
                f"  {fmt(norm - 1, 8):>12}"
                f"  {fmt(endpoint, 18):>25}"
                f"  {fmt(A, 18):>25}"
            )

        # ----------------------------------------------------
        # A scaling
        # ----------------------------------------------------

        print()
        print("LEADING COEFFICIENT SCALING")
        print("-" * 78)

        print(
            f"{'N':>4}"
            f"  {'T_v(0)':>25}"
            f"  {'|T_v(0)|':>20}"
            f"  {'A':>25}"
        )

        print("-" * 78)

        for N in N_VALUES:

            endpoint = states[N]["endpoint"]
            A = states[N]["A"]

            print(
                f"{N:4d}"
                f"  {fmt(endpoint, 22):>25}"
                f"  {fmt(abs(endpoint), 16):>20}"
                f"  {fmt(A, 22):>25}"
            )

        # ----------------------------------------------------
        # B extraction
        # ----------------------------------------------------

        print()
        print("NEXT-ORDER COEFFICIENT B — LARGE-r EXTRACTION")
        print("-" * 78)

        print(
            "For each N, the displayed values are"
        )
        print(
            "  B_est(r) = r^4 [ R(r) - A/r^2 ]."
        )
        print(
            "Convergence with r is a check on the O(r^-6)"
        )
        print(
            "remainder, not yet a derivation of B."
        )

        for N in N_VALUES:

            v = states[N]["v"]
            A = states[N]["A"]

            print()
            print(f"N = {N}")
            print("-" * 78)

            print(
                f"{'r':>12}"
                f"  {'B_est(r)':>30}"
            )

            print("-" * 78)

            for r in R_VALUES:

                B_est = extract_B(
                    v,
                    L,
                    A,
                    r,
                )

                print(
                    f"{fmt(r, 10):>12}"
                    f"  {fmt(B_est, 30):>30}"
                )

        # ----------------------------------------------------
        # Factorisation sanity check for each N
        # ----------------------------------------------------

        print()
        print("FACTORISATION SANITY CHECK")
        print("-" * 78)

        # Use one generic point per N.
        r_check = mp.mpf("137.3")

        print(
            f"{'N':>4}"
            f"  {'relative factorisation error':>32}"
        )

        print("-" * 78)

        for N in N_VALUES:

            v = states[N]["v"]

            phase = (
                1
                - mp.cos(r_check * L)
            )

            K_direct = K_fourier(
                v,
                r_check,
                L,
            )

            K_reconstructed = (
                phase
                * K_reduced(
                    v,
                    r_check,
                    L,
                )
            )

            if K_direct == 0:
                error = abs(
                    K_reconstructed
                    - K_direct
                )
            else:
                error = abs(
                    (
                        K_reconstructed
                        - K_direct
                    )
                    / K_direct
                )

            print(
                f"{N:4d}"
                f"  {fmt(error, 12):>32}"
            )

        # ----------------------------------------------------
        # Ratios between successive N values
        # ----------------------------------------------------

        print()
        print("SUCCESSIVE-N RATIOS")
        print("-" * 78)

        print(
            "These ratios are exploratory only; they are not"
        )
        print(
            "intended as fits to an asymptotic law."
        )

        print()
        print(
            f"{'N -> N+2':>12}"
            f"  {'|T_N+2(0)| / |T_N(0)|':>30}"
            f"  {'A_N+2 / A_N':>25}"
        )

        print("-" * 78)

        for i in range(len(N_VALUES) - 1):

            N0 = N_VALUES[i]
            N1 = N_VALUES[i + 1]

            t0 = abs(
                states[N0]["endpoint"]
            )

            t1 = abs(
                states[N1]["endpoint"]
            )

            A0 = states[N0]["A"]
            A1 = states[N1]["A"]

            if t0 == 0:
                t_ratio = mp.inf
            else:
                t_ratio = t1 / t0

            if A0 == 0:
                A_ratio = mp.inf
            else:
                A_ratio = A1 / A0

            print(
                f"{N0:5d} -> {N1:<5d}"
                f"  {fmt(t_ratio, 20):>30}"
                f"  {fmt(A_ratio, 20):>25}"
            )

    elapsed = time.perf_counter() - start

    print()
    print("=" * 78)
    print(
        f"Elapsed time: {elapsed:.3f} s"
    )
    print("=" * 78)


if __name__ == "__main__":
    main()
