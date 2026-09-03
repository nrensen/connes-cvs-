"""
CELL 34 — SYSTEMATIC N-SCAN OF GROUND-STATE TAIL STRUCTURE

Purpose
-------
Extend Cell 33 from even N only to every N in a moderate range.

The primary question is whether the rapid suppression of

    T_v(0) = v_0 + sqrt(2) * sum_{m=1}^N v_m

is a systematic property of the Galerkin ground state as N increases.

We also record:
    * T_v(L)
    * low-order spectral moments of the coefficients
    * A_N = 2 T_v(0)^2 / L
    * numerical B extraction
    * lambda_min

No Archimedean integration is performed here.

The ground state is the expensive part. Existing cached states are reused.
"""

import time
import mpmath as mp

from cell import get_ground_state


# ============================================================================
# PARAMETERS
# ============================================================================

mp.mp.dps = 50

c = 13
T_ground = 400
GROUND_DPS = 50

# Every N, deliberately including odd N.
N_VALUES = list(range(1, 25))

# Large-r points used only for extracting the next coefficient B.
R_VALUES = [
    mp.mpf("1e3"),
    mp.mpf("3e3"),
    mp.mpf("1e4"),
    mp.mpf("3e4"),
    mp.mpf("1e5"),
]

L = mp.log(c)

print("=" * 78)
print("CELL 34 — SYSTEMATIC N-SCAN OF GROUND-STATE TAIL STRUCTURE")
print("=" * 78)
print()
print(f"c = {c}")
print(f"T_ground = {T_ground}")
print(f"ground dps = {GROUND_DPS}")
print(f"working dps = {mp.mp.dps}")
print(f"N values = {N_VALUES}")
print()
print(f"L = {mp.nstr(L, 30)}")
print()


# ============================================================================
# EXACT REDUCED FOURIER-SIDE COMPONENTS
# ============================================================================

def S_reduced(m, r, L):
    """
    S_m(r) / (1 - cos(rL)).
    """
    a = 2 * mp.pi * m / L
    return a / (a*a - r*r)


def C_reduced(m, r, L):
    """
    C_m(r) / (1 - cos(rL)).
    """
    a = 2 * mp.pi * m / L
    return (a*a + r*r) / (L * (a*a - r*r)**2)


def R_reduced(v, r, L):
    """
    Rational factor R_v(r) defined by

        K_fourier(v,r,L) = (1 - cos(rL)) R_v(r).

    This is the exact finite-N reduced expression.
    """
    N = len(v) - 1

    total = mp.mpf("0")

    # Diagonal C terms.
    for m in range(N + 1):
        total += 2 * v[m]**2 * C_reduced(m, r, L)

    # v0-vm terms.
    for m in range(1, N + 1):
        total -= (
            2 * mp.sqrt(2) * v[0] / mp.pi
            * v[m] * S_reduced(m, r, L) / m
        )

    # vm^2 S terms.
    for m in range(1, N + 1):
        total -= (
            v[m]**2 / mp.pi
            * S_reduced(m, r, L) / m
        )

    # Off-diagonal terms.
    for m in range(1, N + 1):
        for n in range(m + 1, N + 1):
            total += (
                4 / mp.pi
                * v[m] * v[n]
                * (
                    m * S_reduced(m, r, L)
                    - n * S_reduced(n, r, L)
                )
                / (n*n - m*m)
            )

    return total


# ============================================================================
# GROUND-STATE SCAN
# ============================================================================

results = []

print("=" * 78)
print("GROUND-STATE / ENDPOINT SUMMARY")
print("=" * 78)
print()
print(
    " N"
    "              lambda_min"
    "              |T(0)|"
    "              |T(L)|"
    "                    A"
)
print("-" * 78)

for N in N_VALUES:

    t0 = time.perf_counter()

    print(f"Building N={N} ...", flush=True)

    lam, v, meta = get_ground_state(
        c=c,
        N=N,
        T=T_ground,
        dps=GROUND_DPS,
        verbose=True,
    )

    elapsed = time.perf_counter() - t0

    # ------------------------------------------------------------------------
    # Endpoint values of T_v.
    # ------------------------------------------------------------------------

    T0 = v[0]

    for m in range(1, N + 1):
        T0 += mp.sqrt(2) * v[m]

    TL = v[0]

    for m in range(1, N + 1):
        TL += mp.sqrt(2) * ((-1) ** m) * v[m]

    T0_abs = abs(T0)
    TL_abs = abs(TL)

    # ------------------------------------------------------------------------
    # Leading tail coefficient.
    # ------------------------------------------------------------------------

    A = 2 * T0**2 / L

    # ------------------------------------------------------------------------
    # Low-order coefficient moments.
    #
    # These are useful because:
    #
    #   T''(0) = -sqrt(2) * (2*pi/L)^2 * sum(m^2 v_m)
    #
    # and similarly for higher derivatives.
    # ------------------------------------------------------------------------

    M2 = mp.fsum(
        (mp.mpf(m)**2) * v[m]
        for m in range(1, N + 1)
    )

    M4 = mp.fsum(
        (mp.mpf(m)**4) * v[m]
        for m in range(1, N + 1)
    )

    # ------------------------------------------------------------------------
    # B extraction.
    #
    # R(r) = A/r^2 + B/r^4 + O(r^-6)
    #
    # so
    #
    # B_est(r) = r^4 [R(r) - A/r^2].
    # ------------------------------------------------------------------------

    B_values = []

    for r in R_VALUES:
        R = R_reduced(v, r, L)
        B_est = r**4 * (R - A / r**2)
        B_values.append(B_est)

    B_last = B_values[-1]

    results.append({
        "N": N,
        "lambda": lam,
        "T0": T0,
        "T0_abs": T0_abs,
        "TL": TL,
        "TL_abs": TL_abs,
        "A": A,
        "M2": M2,
        "M4": M4,
        "B_est": B_last,
        "elapsed": elapsed,
    })

    print(
        f"{N:2d}"
        f"  {mp.nstr(lam, 18):>22}"
        f"  {mp.nstr(T0_abs, 18):>22}"
        f"  {mp.nstr(TL_abs, 18):>22}"
        f"  {mp.nstr(A, 18):>22}"
    )

print()


# ============================================================================
# SUCCESSIVE-N SCALING
# ============================================================================

print("=" * 78)
print("SUCCESSIVE-N SCALING")
print("=" * 78)
print()

print(
    " N -> N+1"
    "       |T0(N+1)|/|T0(N)|"
    "              A(N+1)/A(N)"
    "          |B(N+1)|/|B(N)|"
)
print("-" * 78)

for i in range(len(results) - 1):

    a = results[i]
    b = results[i + 1]

    ratio_T0 = b["T0_abs"] / a["T0_abs"]
    ratio_A = b["A"] / a["A"]

    if a["B_est"] != 0:
        ratio_B = abs(b["B_est"]) / abs(a["B_est"])
    else:
        ratio_B = mp.nan

    print(
        f"{a['N']:2d} -> {b['N']:2d}"
        f"        {mp.nstr(ratio_T0, 16):>18}"
        f"        {mp.nstr(ratio_A, 16):>18}"
        f"        {mp.nstr(ratio_B, 16):>18}"
    )

print()


# ============================================================================
# LOG-SCALE SUPPRESSION
# ============================================================================

print("=" * 78)
print("LOG-SCALE SUPPRESSION")
print("=" * 78)
print()
print(
    " N"
    "          log10 |T_v(0)|"
    "             log10 A"
    "          log10 |B_est|"
)
print("-" * 78)

for row in results:

    log_T0 = mp.log10(row["T0_abs"])
    log_A = mp.log10(abs(row["A"]))
    log_B = mp.log10(abs(row["B_est"]))

    print(
        f"{row['N']:2d}"
        f"        {mp.nstr(log_T0, 18):>20}"
        f"        {mp.nstr(log_A, 18):>20}"
        f"        {mp.nstr(log_B, 18):>20}"
    )

print()


# ============================================================================
# ENDPOINT / MOMENT STRUCTURE
# ============================================================================

print("=" * 78)
print("ENDPOINT / LOW-ORDER MOMENT STRUCTURE")
print("=" * 78)
print()

print(
    " N"
    "                 T_v(0)"
    "                 T_v(L)"
    "                  M2"
    "                  M4"
)
print("-" * 78)

for row in results:

    print(
        f"{row['N']:2d}"
        f"  {mp.nstr(row['T0'], 22):>24}"
        f"  {mp.nstr(row['TL'], 22):>24}"
        f"  {mp.nstr(row['M2'], 22):>24}"
        f"  {mp.nstr(row['M4'], 22):>24}"
    )

print()


# ============================================================================
# B CONVERGENCE CHECK FOR EACH N
# ============================================================================

print("=" * 78)
print("B CONVERGENCE CHECK")
print("=" * 78)
print()

for row in results:

    N = row["N"]

    # Recover v again from cache; this is cheap.
    _, v, _ = get_ground_state(
        c=c,
        N=N,
        T=T_ground,
        dps=GROUND_DPS,
        verbose=False,
    )

    A = row["A"]

    print(f"N = {N}")
    print("-" * 78)
    print(" r                 B_est(r)")

    for r in R_VALUES:

        R = R_reduced(v, r, L)
        B_est = r**4 * (R - A / r**2)

        print(
            f"{mp.nstr(r, 8):>10}"
            f"        {mp.nstr(B_est, 30)}"
        )

    print()


# ============================================================================
# FINAL INTERPRETATION
# ============================================================================

print("=" * 78)
print("INTERPRETATION GUIDE")
print("=" * 78)
print()
print(
    "This cell is exploratory. No asymptotic law in N is assumed."
)
print()
print(
    "The principal diagnostic is whether |T_v(0)| continues to fall"
)
print(
    "systematically as N increases, and whether the behaviour is"
)
print(
    "approximately geometric, faster than geometric, or irregular."
)
print()
print(
    "T_v(L), M2 and M4 are included to distinguish a generic small"
)
print(
    "coefficient from a more structured endpoint cancellation."
)
print()
print(
    "B_est(r) is only a numerical extraction of the r^-4 coefficient."
)
print(
    "Its convergence with r checks the asymptotic expansion but does"
)
print(
    "not constitute an independent derivation of B."
)
print()
print(
    "No Archimedean tail integral is evaluated by this cell."
)
print()
