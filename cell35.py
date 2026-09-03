"""
CELL 35 — ENDPOINT DERIVATIVE SUPPRESSION

Purpose
-------
Investigate whether the Galerkin ground states increasingly suppress
successive even derivatives of

    T_v(t) = v_0 + sqrt(2) * sum_{m=1}^N v_m cos(2*pi*m*t/L)

at the endpoint t = 0 as N increases.

The experiment is purely algebraic once the ground states are available.
It performs no Archimedean integration and no new ground-state generation.

For k >= 1,

    T_v^(2k)(0)
      = sqrt(2) * (-1)^k * (2*pi/L)^(2k)
        * sum_m m^(2k) v_m.

The k=0 quantity is

    T_v(0) = v_0 + sqrt(2) * sum_m v_m.

The corresponding endpoint at t=L is also reported for T_v itself.

The principal question is whether suppression extends systematically
from T_v(0) to T_v''(0), T_v^(4)(0), ... as N increases.
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

# Number of even derivatives to inspect:
#
# k = 0 -> T(0)
# k = 1 -> T''(0)
# k = 2 -> T^(4)(0)
# ...
#
# k = 4 gives T^(8)(0).
K_MAX = 4

mp.mp.dps = GROUND_DPS

L = mp.log(c)
alpha = 2 * mp.pi / L


# ============================================================================
# HEADER
# ============================================================================

print("=" * 78)
print("CELL 35 — ENDPOINT DERIVATIVE SUPPRESSION")
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
print(f"2*pi/L = {mp.nstr(alpha, 30)}")
print()


# ============================================================================
# DERIVATIVE CALCULATION
# ============================================================================

def endpoint_even_derivatives(v, K_MAX):
    """
    Return

        T(0), T''(0), T^(4)(0), ..., T^(2K_MAX)(0)

    for

        T(t) = v_0 + sqrt(2) sum_{m>=1} v_m cos(2*pi*m*t/L).

    The returned list has entries indexed by k.
    """
    N = len(v) - 1

    values = []

    # k = 0
    T0 = v[0] + mp.sqrt(2) * mp.fsum(
        v[m] for m in range(1, N + 1)
    )

    values.append(T0)

    # k >= 1
    for k in range(1, K_MAX + 1):

        moment = mp.fsum(
            (mp.mpf(m) ** (2 * k)) * v[m]
            for m in range(1, N + 1)
        )

        derivative = (
            mp.sqrt(2)
            * (-1) ** k
            * alpha ** (2 * k)
            * moment
        )

        values.append(derivative)

    return values


def endpoint_even_moments(v, K_MAX):
    """
    Return the dimensionless spectral moments

        sum_m m^(2k) v_m

    for k = 1,...,K_MAX.
    """
    N = len(v) - 1

    return [
        mp.fsum(
            (mp.mpf(m) ** (2 * k)) * v[m]
            for m in range(1, N + 1)
        )
        for k in range(1, K_MAX + 1)
    ]


def T_at_L(v):
    """
    T_v(L) = v_0 + sqrt(2) sum_m (-1)^m v_m.
    """
    N = len(v) - 1

    return v[0] + mp.sqrt(2) * mp.fsum(
        ((-1) ** m) * v[m]
        for m in range(1, N + 1)
    )


# ============================================================================
# COLLECT RESULTS
# ============================================================================

results = []

print("=" * 78)
print("GROUND-STATE ENDPOINT DATA")
print("=" * 78)
print()

print(
    " N"
    "             T(0)"
    "             T(L)"
    "             T''(0)"
    "             T^(4)(0)"
)
print("-" * 78)

for N in N_VALUES:

    lam, v, meta = get_ground_state(
        c=c,
        N=N,
        T=T_ground,
        dps=GROUND_DPS,
        verbose=False,
    )

    derivatives = endpoint_even_derivatives(v, K_MAX)

    T0 = derivatives[0]
    T2 = derivatives[1]
    T4 = derivatives[2]

    TL = T_at_L(v)

    results.append({
        "N": N,
        "lambda": lam,
        "v": v,
        "derivatives": derivatives,
        "TL": TL,
    })

    print(
        f"{N:2d}"
        f"  {mp.nstr(T0, 18):>22}"
        f"  {mp.nstr(TL, 18):>22}"
        f"  {mp.nstr(T2, 18):>22}"
        f"  {mp.nstr(T4, 18):>22}"
    )

print()


# ============================================================================
# ABSOLUTE MAGNITUDES ON LOG SCALE
# ============================================================================

print("=" * 78)
print("LOG10 MAGNITUDES OF EVEN ENDPOINT DERIVATIVES")
print("=" * 78)
print()

header = f"{'N':>3}"

for k in range(K_MAX + 1):
    if k == 0:
        label = "log10|T(0)|"
    else:
        label = f"log10|T^({2*k})(0)|"
    header += f"{label:>22}"

print(header)
print("-" * 78)

for row in results:

    line = f"{row['N']:3d}"

    for k in range(K_MAX + 1):
        value = abs(row["derivatives"][k])

        if value == 0:
            text = "-inf"
        else:
            text = mp.nstr(mp.log10(value), 16)

        line += f"{text:>22}"

    print(line)

print()


# ============================================================================
# DIMENSIONLESS DERIVATIVE SCALES
# ============================================================================

print("=" * 78)
print("DIMENSIONLESS EVEN-DERIVATIVE SCALES")
print("=" * 78)
print()

print(
    "These divide T^(2k)(0) by (2*pi/L)^(2k),"
)
print(
    "isolating the coefficient moment sum_m m^(2k) v_m."
)
print()

header = f"{'N':>3}"

for k in range(1, K_MAX + 1):
    header += f"{'|D_' + str(2*k) + '| / alpha^' + str(2*k):>24}"

print(header)
print("-" * 78)

for row in results:

    line = f"{row['N']:3d}"

    for k in range(1, K_MAX + 1):

        value = abs(row["derivatives"][k]) / alpha ** (2 * k)

        line += f"{mp.nstr(value, 18):>24}"

    print(line)

print()


# ============================================================================
# SUCCESSIVE-N RATIOS
# ============================================================================

print("=" * 78)
print("SUCCESSIVE-N RATIOS OF ABSOLUTE DERIVATIVE MAGNITUDES")
print("=" * 78)
print()

print(
    "These are |D_(N+1)| / |D_N|."
)
print()

header = f"{'N -> N+1':>9}"

for k in range(K_MAX + 1):
    if k == 0:
        label = "|T0| ratio"
    else:
        label = f"|T^({2*k})| ratio"
    header += f"{label:>20}"

print(header)
print("-" * 78)

for i in range(len(results) - 1):

    row_a = results[i]
    row_b = results[i + 1]

    line = f"{row_a['N']:2d} -> {row_b['N']:2d}"

    for k in range(K_MAX + 1):

        a = abs(row_a["derivatives"][k])
        b = abs(row_b["derivatives"][k])

        if a == 0:
            ratio = mp.nan
        else:
            ratio = b / a

        line += f"{mp.nstr(ratio, 14):>20}"

    print(line)

print()


# ============================================================================
# MOMENT SIGN STRUCTURE
# ============================================================================

print("=" * 78)
print("SIGNED SPECTRAL MOMENTS")
print("=" * 78)
print()

print(
    "M_(2k) = sum_m m^(2k) v_m"
)
print()

header = f"{'N':>3}"

for k in range(1, K_MAX + 1):
    header += f"{'M_' + str(2*k):>24}"

print(header)
print("-" * 78)

for row in results:

    moments = endpoint_even_moments(row["v"], K_MAX)

    line = f"{row['N']:3d}"

    for moment in moments:
        line += f"{mp.nstr(moment, 18):>24}"

    print(line)

print()


# ============================================================================
# CROSS-CHECK: DERIVATIVE / MOMENT IDENTITY
# ============================================================================

print("=" * 78)
print("DERIVATIVE / MOMENT CROSS-CHECK")
print("=" * 78)
print()

max_error = mp.mpf("0")

for row in results:

    moments = endpoint_even_moments(row["v"], K_MAX)

    for k in range(1, K_MAX + 1):

        expected = (
            mp.sqrt(2)
            * (-1) ** k
            * alpha ** (2 * k)
            * moments[k - 1]
        )

        actual = row["derivatives"][k]

        error = abs(actual - expected)

        if error > max_error:
            max_error = error

print(
    "Maximum absolute derivative/moment reconstruction error:"
)
print(
    mp.nstr(max_error, 30)
)
print()


# ============================================================================
# SUMMARY
# ============================================================================

print("=" * 78)
print("SUMMARY")
print("=" * 78)
print()

print(
    "This is a structural survey only."
)
print(
    "No asymptotic law in N is fitted or assumed."
)
print()
print(
    "The key question is whether suppression of T(0) is accompanied"
)
print(
    "by suppression of T''(0), T^(4)(0), T^(6)(0), ... ."
)
print()
print(
    "If successive derivatives are suppressed in an ordered fashion,"
)
print(
    "that would be evidence for increasing endpoint flatness of the"
)
print(
    "Galerkin ground state, rather than an isolated cancellation."
)
print()
print(
    "T(L) is retained as a control quantity: if T(0) and its even"
)
print(
    "derivatives collapse while T(L) remains O(1), the phenomenon is"
)
print(
    "strongly endpoint-specific."
)
print()
print(
    "No conclusion about the N -> infinity limit is drawn by this cell."
)
print()
