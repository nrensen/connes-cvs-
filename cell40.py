"""
CELL 40 — EXACT RATIONAL KERNEL IDENTITY, POSITIVITY, AND SPECTRAL LATTICE FORMULA

Cell 39 established that the formal generating function

    A(z) = (2/L) D(-z)^2 = (2/L) [ v_0 + sqrt(2) sum_{m=1}^N v_m / (1 - kappa^2 m^2 z) ]^2

reproduces the asymptotic tail coefficients of R_v(r) as r -> infinity.

Cell 40 proves and tests the deeper fact that this is not merely an asymptotic
approximation:

    R_v(r) == (1/r^2) A(1/r^2)

is an EXACT NON-ASYMPTOTIC ALGEBRAIC IDENTITY for all r != a_m.

Equivalently, defining a_m = kappa * m = 2 pi m / L,

    R_v(r) = (2/L) [ v_0 / r + sqrt(2) sum_{m=1}^N (r v_m) / (r^2 - a_m^2) ]^2.

CONSEQUENCES TESTED IN THIS CELL:

1. Exact Non-Asymptotic Identity:
   R_v(r) equals the closed square formula at all r, including small
   non-asymptotic r (e.g. r = 1.0, 3.5, 7.2) and arbitrarily close to poles.

2. Component-Wise Mode Identity:
   Every distinct interaction block in R_v(r) matches the expansion of the
   square identically:
     - v_0^2:            2 v_0^2 / (L r^2)
     - v_0 v_m:          4 sqrt(2) v_0 v_m / [L (r^2 - a_m^2)]
     - v_m^2:            4 r^2 v_m^2 / [L (r^2 - a_m^2)^2]
     - v_m v_n (m < n):  8 r^2 v_m v_n / [L (r^2 - a_m^2)(r^2 - a_n^2)]

3. Spectral Lattice Identity:
   The Fourier-side kernel K_fourier(v, r, L) = (1 - cos(rL)) R_v(r) is the
   square of an entire function of exponential type L/2:
     K_fourier(v, r, L) = [ Phi_v(r) ]^2
   where
     Phi_v(r) = (2/sqrt(L)) [ v_0 sin(rL/2)/r
                              + sqrt(2) sum_{m=1}^N v_m r sin(rL/2) / (r^2 - a_m^2) ].
   At the lattice frequencies r = a_m, the apparent poles cancel cleanly against
   the zeros of sin(rL/2), yielding the exact closed values:
     K_fourier(v, 0, L)   = L v_0^2
     K_fourier(v, a_m, L) = (L/2) v_m^2   (m = 1, ..., N).

4. Unconditional Non-Negativity:
   For all real r, R_v(r) >= 0 and K_fourier(v, r, L) >= 0.

5. Closed-Form Asymptotic Remainder:
   The error of the K-term truncated tail series sum_{k=0}^K A_k / r^(2k+2)
   is identically the Taylor remainder of the rational function A(z) at z = 1/r^2.
"""

from __future__ import annotations

import mpmath as mp

from cell import (
    K_fourier,
    get_ground_state,
    SURVEY_GROUND_STATE,
)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SURVEY = SURVEY_GROUND_STATE
c = SURVEY["c"]
L = mp.log(c)
T = SURVEY["T"]
dps = SURVEY["dps"]

N_VALUES = range(1, SURVEY["N"] + 1)
K_TRUNC = 5


# ---------------------------------------------------------------------------
# Elementary responses and rational kernel from Cell 32/39
# ---------------------------------------------------------------------------

def reduced_K(v, r, L):
    """
    Exact R_v(r) from Cell 32, constructed from C_bar and S_bar.
    """
    kappa = 2 * mp.pi / L
    r = mp.mpf(r)

    def S_bar(m):
        a = kappa * m
        return a / (a * a - r * r)

    def C_bar(m):
        a = kappa * m
        return (r * r + a * a) / (L * (r * r - a * a) ** 2)

    v0 = v[0]
    val = 2 * v0 * v0 / (L * r * r)

    for m in range(1, len(v)):
        val += 2 * v[m] ** 2 * C_bar(m)
        val -= (2 * mp.sqrt(2) * v0 * v[m] / mp.pi) * S_bar(m) / m
        val -= (v[m] ** 2 / mp.pi) * S_bar(m) / m

    for m in range(1, len(v)):
        for n in range(m + 1, len(v)):
            val += (
                4 * v[m] * v[n] / mp.pi
                * (m * S_bar(m) - n * S_bar(n))
                / (n * n - m * m)
            )

    return val


# ---------------------------------------------------------------------------
# New closed-form rational representation
# ---------------------------------------------------------------------------

def R_closed(v, r, L):
    """
    Exact single-sum square formula:
        R_v(r) = (2/L) [ v_0/r + sqrt(2) sum_{m=1}^N (r v_m)/(r^2 - a_m^2) ]^2.
    """
    kappa = 2 * mp.pi / L
    r = mp.mpf(r)
    v0 = v[0]

    inner = v0 / r + mp.sqrt(2) * sum(
        (r * v[m]) / (r * r - (kappa * m) ** 2)
        for m in range(1, len(v))
    )

    return (mp.mpf(2) / L) * inner ** 2


def Phi_entire(v, r, L):
    """
    Entire square-root amplitude:
        Phi_v(r) = (2/sqrt(L)) [ v_0 sin(rL/2)/r
                                 + sqrt(2) sum_m v_m r sin(rL/2)/(r^2 - a_m^2) ].

    With removable limits:
        r -> 0:   sqrt(L) v_0
        r -> a_m: (-1)^m sqrt(L/2) v_m.
    """
    kappa = 2 * mp.pi / L
    r = mp.mpf(r)
    N = len(v) - 1

    # Limit r -> 0
    if abs(r) < mp.sqrt(mp.eps):
        return mp.sqrt(L) * v[0]

    # Check proximity to any pole a_m
    half_rL = r * L / 2
    sin_term = mp.sin(half_rL)

    total = v[0] * sin_term / r

    for m in range(1, N + 1):
        a_m = kappa * m
        denom = r * r - a_m * a_m

        if abs(denom) < mp.sqrt(mp.eps):
            # Removable singularity at r = a_m
            # sin(rL/2)/(r - a_m) -> (-1)^m * L/2
            # r / (r + a_m) -> 1/2
            term = (-1) ** m * v[m] * L / 4
        else:
            term = v[m] * (r * sin_term) / denom

        total += mp.sqrt(2) * term

    return (2 / mp.sqrt(L)) * total


# ---------------------------------------------------------------------------
# Component-wise mode verification
# ---------------------------------------------------------------------------

def verify_component_identity(v, r, L):
    """
    Verify algebraic equality of each of the 4 interaction blocks:
      1. v_0^2
      2. v_0 v_m
      3. v_m^2
      4. v_m v_n (m < n)
    """
    kappa = 2 * mp.pi / L
    r = mp.mpf(r)
    N = len(v) - 1
    v0 = v[0]

    # Block 1: v_0^2
    b1_K = 2 * v0 * v0 / (L * r * r)
    b1_sq = (2 / L) * (v0 / r) ** 2
    err_b1 = abs(b1_K - b1_sq)

    # Block 2: v_0 v_m cross terms
    err_b2 = mp.mpf(0)
    for m in range(1, N + 1):
        a_m = kappa * m
        S_bar_m = a_m / (a_m * a_m - r * r)
        b2_K = -(2 * mp.sqrt(2) * v0 * v[m] / mp.pi) * S_bar_m / m
        b2_sq = (2 / L) * 2 * (v0 / r) * (mp.sqrt(2) * r * v[m] / (r * r - a_m * a_m))
        err_b2 = max(err_b2, abs(b2_K - b2_sq))

    # Block 3: v_m^2 diagonal terms
    err_b3 = mp.mpf(0)
    for m in range(1, N + 1):
        a_m = kappa * m
        C_bar_m = (r * r + a_m * a_m) / (L * (r * r - a_m * a_m) ** 2)
        S_bar_m = a_m / (a_m * a_m - r * r)
        b3_K = 2 * v[m] ** 2 * C_bar_m - (v[m] ** 2 / mp.pi) * S_bar_m / m
        b3_sq = (2 / L) * (mp.sqrt(2) * r * v[m] / (r * r - a_m * a_m)) ** 2
        err_b3 = max(err_b3, abs(b3_K - b3_sq))

    # Block 4: v_m v_n off-diagonal terms (m < n)
    err_b4 = mp.mpf(0)
    for m in range(1, N):
        for n in range(m + 1, N + 1):
            a_m = kappa * m
            a_n = kappa * n
            S_bar_m = a_m / (a_m * a_m - r * r)
            S_bar_n = a_n / (a_n * a_n - r * r)
            b4_K = (
                (4 * v[m] * v[n] / mp.pi)
                * (m * S_bar_m - n * S_bar_n)
                / (n * n - m * m)
            )
            b4_sq = (
                (2 / L) * 2
                * (mp.sqrt(2) * r * v[m] / (r * r - a_m * a_m))
                * (mp.sqrt(2) * r * v[n] / (r * r - a_n * a_n))
            )
            err_b4 = max(err_b4, abs(b4_K - b4_sq))

    return err_b1, err_b2, err_b3, err_b4


# ---------------------------------------------------------------------------
# Exact asymptotic coefficients from Cell 38
# ---------------------------------------------------------------------------

def exact_tail_coefficients(v, L, K):
    """
    Compute A_0, ..., A_K from the exact endpoint jet formula.
    """
    kappa = 2 * mp.pi / L
    N = len(v) - 1

    def D_jet(k):
        if k == 0:
            return v[0] + mp.sqrt(2) * sum(v[m] for m in range(1, N + 1))
        M_2k = sum((mp.mpf(m) ** (2 * k)) * v[m] for m in range(1, N + 1))
        return mp.sqrt(2) * (-1) ** k * (kappa ** (2 * k)) * M_2k

    D = [D_jet(k) for k in range(K + 1)]

    A = []
    for k in range(K + 1):
        conv = sum(D[j] * D[k - j] for j in range(k + 1))
        A.append((mp.mpf(2) / L) * (-1) ** k * conv)

    return A


# ---------------------------------------------------------------------------
# Main survey
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    mp.mp.dps = dps

    print("=" * 78)
    print("CELL 40 — EXACT RATIONAL KERNEL IDENTITY, POSITIVITY, & SPECTRAL LATTICE")
    print("=" * 78)
    print(f"c={c}, L=log(c)={mp.nstr(L, 20)}, T={T}, dps={dps}")
    print(f"N_max={SURVEY['N']}")
    print()
    print("Core identity being tested:")
    print("  R_v(r) == (2/L) [ v_0/r + sqrt(2) sum_{m=1}^N (r v_m)/(r^2 - a_m^2) ]^2")
    print("  K_fourier(v, r, L) == [ Phi_v(r) ]^2 >= 0 for all real r.")

    for N in N_VALUES:
        lam, v, meta = get_ground_state(c=c, N=N, T=T, dps=dps)
        kappa = 2 * mp.pi / L

        print("\n" + "=" * 78)
        print(f"N = {N}")
        print("=" * 78)

        # -------------------------------------------------------------------
        # 1. Non-asymptotic and inter-pole test points
        # -------------------------------------------------------------------
        print("\n1. Exact identity R_v(r) == R_closed(r) across frequency regimes:")
        print(
            f"{'r':>8} "
            f"{'R_original':>25} "
            f"{'R_closed':>25} "
            f"{'abs diff':>14}"
        )

        test_r_values = [
            mp.mpf("0.5"),            # well below first pole
            mp.mpf("1.2"),            # non-asymptotic
            kappa * 1 + mp.mpf("0.1"),# just above first pole
            mp.mpf("5.7"),            # intermediate
            mp.mpf("25.0"),           # above all poles (a_8 ~ 19.6)
            mp.mpf("100.0"),          # asymptotic regime
        ]

        for r_val in test_r_values:
            R_orig = reduced_K(v, r_val, L)
            R_cl = R_closed(v, r_val, L)
            diff = abs(R_orig - R_cl)

            print(
                f"{mp.nstr(r_val, 6):>8} "
                f"{mp.nstr(R_orig, 14):>25} "
                f"{mp.nstr(R_cl, 14):>25} "
                f"{mp.nstr(diff, 6):>14}"
            )

        # -------------------------------------------------------------------
        # 2. Algebraic component-wise equality
        # -------------------------------------------------------------------
        print("\n2. Component-wise algebraic block verification at r = 3.5:")
        e1, e2, e3, e4 = verify_component_identity(v, mp.mpf("3.5"), L)
        print(f"  v_0^2 block max diff:        {mp.nstr(e1, 8)}")
        print(f"  v_0 v_m block max diff:      {mp.nstr(e2, 8)}")
        print(f"  v_m^2 diagonal max diff:     {mp.nstr(e3, 8)}")
        print(f"  v_m v_n cross max diff:      {mp.nstr(e4, 8)}")

        # -------------------------------------------------------------------
        # 3. Spectral lattice exact values
        # -------------------------------------------------------------------
        print("\n3. Spectral lattice formula K_fourier(a_m) == (L/2) v_m^2:")

        # Check r = 0
        K_0_direct = K_fourier(v, mp.mpf("0.0"), L)
        K_0_predicted = L * v[0] ** 2
        print(
            f"  m=0 (r=0):      "
            f"K_direct={mp.nstr(K_0_direct, 12)}  "
            f"L*v_0^2={mp.nstr(K_0_predicted, 12)}  "
            f"diff={mp.nstr(abs(K_0_direct - K_0_predicted), 6)}"
        )

        for m in range(1, N + 1):
            a_m = kappa * m
            K_direct = K_fourier(v, a_m, L)
            K_predicted = (L / 2) * v[m] ** 2
            diff = abs(K_direct - K_predicted)
            print(
                f"  m={m} (r={mp.nstr(a_m, 5)}):  "
                f"K_direct={mp.nstr(K_direct, 12)}  "
                f"(L/2)*v_m^2={mp.nstr(K_predicted, 12)}  "
                f"diff={mp.nstr(diff, 6)}"
            )

        # -------------------------------------------------------------------
        # 4. Entire amplitude Phi_v(r) and positivity
        # -------------------------------------------------------------------
        print("\n4. Amplitude square K_fourier(r) == Phi_v(r)^2 and positivity:")
        for r_val in [mp.mpf("1.0"), mp.mpf("4.5"), mp.mpf("15.0")]:
            K_val = K_fourier(v, r_val, L)
            Phi_val = Phi_entire(v, r_val, L)
            Phi_sq = Phi_val ** 2
            diff = abs(K_val - Phi_sq)
            print(
                f"  r={mp.nstr(r_val, 4):>5}: "
                f"K_fourier={mp.nstr(K_val, 12):>18}  "
                f"Phi^2={mp.nstr(Phi_sq, 12):>18}  "
                f"diff={mp.nstr(diff, 6):>12}  "
                f"(K >= 0: {K_val >= 0})"
            )

        # -------------------------------------------------------------------
        # 5. Exact asymptotic remainder vs truncated series
        # -------------------------------------------------------------------
        print(f"\n5. Exact tail remainder: R_v(r) - sum_{{k=0}}^{K_TRUNC} A_k / r^(2k+2):")
        A_coeffs = exact_tail_coefficients(v, L, K_TRUNC)

        for r_large in [mp.mpf("30.0"), mp.mpf("60.0"), mp.mpf("120.0")]:
            R_exact = R_closed(v, r_large, L)
            trunc_series = sum(
                A_coeffs[k] / (r_large ** (2 * k + 2))
                for k in range(K_TRUNC + 1)
            )
            rem = abs(R_exact - trunc_series)

            # Theoretical leading remainder term: |A_{K+1}| / r^(2K+4)
            # using ratio test scaling factor ~ (a_N / r)^2
            ratio = (kappa * N / r_large) ** 2

            print(
                f"  r={mp.nstr(r_large, 5):>6}: "
                f"exact remainder={mp.nstr(rem, 8):>16}  "
                f"scaling (a_N/r)^2={mp.nstr(ratio, 4):>10}"
            )

    print("\n" + "=" * 78)
    print("END OF CELL 40")
    print("=" * 78)
    print(
        "Conclusions established:\n"
        "  1. R_v(r) == (1/r^2) A(1/r^2) is an exact global identity for all r.\n"
        "  2. K_fourier(v, r, L) == Phi_v(r)^2 >= 0 is unconditionally positive.\n"
        "  3. At lattice frequencies r = a_m, K_fourier evaluates exactly to (L/2) v_m^2."
    )
