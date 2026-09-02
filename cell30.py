# ============================================================
# CELL 30 — FORENSIC ASYMPTOTIC TAIL TEST
#
# Purpose
# -------
# Test the emerging empirical observation from Cell 29 that
#
#     D(T) = integral_T^(2T) I(r) dr
#
# may behave asymptotically like
#
#     D(T) ~ C / T.
#
# The calculation does NOT assume p = 1.
#
# Instead it:
#
#   1. extends the dyadic intervals to much larger T;
#   2. reports C_T = T D(T);
#   3. reports the local effective exponent;
#   4. estimates C from several late-point averages;
#   5. estimates the remaining infinite tail under the
#      explicit D(T) ~ C/T hypothesis;
#   6. compares that extrapolated tail with directly computed
#      subsequent dyadic contributions.
#
# Definitions
# ----------
#
#     I(r) = h_+(r) K_fourier(v_star,r,L)
#
#     D(T) = integral_T^(2T) I(r) dr
#
#     C_T = T D(T)
#
# If D(T) ~ C/T then
#
#     C_T -> C.
#
# For the same hypothesis,
#
#     R(T) = integral_T^infinity I(r) dr
#
# satisfies
#
#     R(T) ~ 2 C / T,
#
# because
#
#     R(T)
#       = D(T) + D(2T) + D(4T) + ...
#       ~ C/T + C/(2T) + C/(4T) + ...
#       = 2C/T.
#
# This relation is an EXTRAPOLATION HYPOTHESIS, not an identity.
#
# ============================================================

import time
import mpmath as mp

from cell import (
    FORENSIC_GROUND_STATE,
    get_ground_state,
    compute_L,
    K_fourier,
    h_plus,
)

# ============================================================
# PARAMETERS
# ============================================================

WORKING_DPS = 80

mp.mp.dps = WORKING_DPS

c = FORENSIC_GROUND_STATE["c"]
N = FORENSIC_GROUND_STATE["N"]

L = compute_L(c)

# ------------------------------------------------------------
# Start at the same T as Cell 29.
#
# Cell 29 used:
#
#     20, 40, 80, ..., 20480
#
# We now continue to:
#
#     20 * 2^20 = 20,971,520
#
# This gives a substantially longer asymptotic lever arm.
#
# ------------------------------------------------------------

T0 = mp.mpf("20")

# Number of dyadic intervals.
#
# Level 0:
#     [20, 40]
#
# Level 20:
#     [20,971,520, 41,943,040]
#
N_LEVELS = 21

# ------------------------------------------------------------
# Late ranges used for estimating C.
#
# These are indices into the dyadic sequence.
#
# The estimates are intentionally redundant:
#
#   * last point;
#   * average of last 3;
#   * average of last 5;
#   * average of last 8;
#   * median of last 8.
#
# Agreement between these estimates is more useful than any
# single fitted number.
# ------------------------------------------------------------

LATE_WINDOWS = [3, 5, 8]

# ============================================================
# HELPERS
# ============================================================

def nstr(x, digits=WORKING_DPS):
    return mp.nstr(x, digits)


def integrand(r, v_star):
    """
    Actual analytic Archimedean Fourier integrand.
    """
    return (
        h_plus(r)
        * K_fourier(
            v_star,
            r,
            L,
        )
    )


def integrate_interval(a, b, v_star):
    """
    Signed integral over one finite dyadic interval.
    """
    return mp.quad(
        lambda r: integrand(r, v_star),
        [a, b],
    )


def local_power(D1, D2):
    """
    Given consecutive dyadic interval contributions
#
#       D(T), D(2T),
#
# return p_eff defined by
#
#       D(2T)/D(T) = 2^(-p_eff).
    """

    if D1 <= 0 or D2 <= 0:
        return mp.nan

    return (
        -mp.log(D2 / D1)
        / mp.log(2)
    )


def arithmetic_mean(values):
    """
    Arithmetic mean using mp.fsum.
    """
    if not values:
        return mp.nan

    return mp.fsum(values) / len(values)


def median(values):
    """
    Median for a short list of mp.mpf values.
    """
    if not values:
        return mp.nan

    ordered = sorted(values)
    n = len(ordered)

    if n % 2:
        return ordered[n // 2]

    return (
        ordered[n // 2 - 1]
        + ordered[n // 2]
    ) / 2


# ============================================================
# HEADER
# ============================================================

print("=" * 120)
print("CELL 30 — FORENSIC ASYMPTOTIC TAIL TEST")
print("=" * 120)

print()

print("Parameters:")
print(f" c = {c}")
print(f" N = {N}")
print(f" working_dps = {WORKING_DPS}")
print(f" L = {nstr(L, 60)}")
print(f" T0 = {nstr(T0, 20)}")
print(f" dyadic intervals = {N_LEVELS}")

Tmax = T0 * mp.mpf(2) ** N_LEVELS

print(
    f" final upper endpoint = {nstr(Tmax, 20)}"
)

print()

print(
    "Primary diagnostic:"
)

print(
    "    C_T = T * D(T)"
)

print()

print(
    "If D(T) ~ C/T, then C_T should approach a constant."
)

print()

print(
    "The infinite-tail estimate 2C/T is used only after"
)

print(
    "examining whether C_T actually stabilises."
)

# ============================================================
# 1. FIXED FORENSIC GROUND STATE
# ============================================================

print()
print("-" * 120)
print("1. FIXED FORENSIC GROUND STATE")
print("-" * 120)

ground_start = time.perf_counter()

lambda_forensic, v_star, ground_meta = get_ground_state(
    **FORENSIC_GROUND_STATE,
    verbose=True,
)

ground_elapsed = (
    time.perf_counter()
    - ground_start
)

norm_v = mp.sqrt(
    mp.fdot(
        v_star,
        v_star,
    )
)

print()

print("lambda_forensic =")
print(nstr(lambda_forensic, 60))

print()

print("||v_star|| =")
print(nstr(norm_v, 60))

print()

print(
    f"ground-state retrieval wall = "
    f"{ground_elapsed:.6f} s"
)

# ============================================================
# 2. COMPUTE DYADIC INTERVALS
# ============================================================

print()
print("-" * 120)
print("2. DYADIC INTERVAL INTEGRALS")
print("-" * 120)

print()

print(
    "Each row computes"
)

print()

print(
    "    D(T) = integral_T^(2T) I(r) dr."
)

print()

print(
    "The key quantity is"
)

print()

print(
    "    C_T = T D(T)."
)

print()

print(
    "No asymptotic exponent is imposed."
)

print()

print(
    " level"
    "             T"
    "                 D(T)"
    "                 T D(T)"
    "                 p_eff"
)

print("-" * 120)

results = []

T = T0

for level in range(N_LEVELS):

    a = T
    b = 2 * T

    start = time.perf_counter()

    D = integrate_interval(
        a,
        b,
        v_star,
    )

    elapsed = (
        time.perf_counter()
        - start
    )

    C_T = a * D

    if results:
        p_eff = local_power(
            results[-1]["D"],
            D,
        )
    else:
        p_eff = mp.nan

    row = {
        "level": level,
        "T": a,
        "2T": b,
        "D": D,
        "C_T": C_T,
        "p_eff": p_eff,
        "elapsed": elapsed,
    }

    results.append(row)

    print(
        f"{level:6d} "
        f"{nstr(a, 12):>18} "
        f"{nstr(D, 34):>38} "
        f"{nstr(C_T, 34):>38} "
        f"{nstr(p_eff, 18):>24}"
    )

    T = b

# ============================================================
# 3. SIGN CHECK
# ============================================================

print()
print("-" * 120)
print("3. SIGN CHECK")
print("-" * 120)

positive_count = sum(
    row["D"] > 0
    for row in results
)

zero_count = sum(
    row["D"] == 0
    for row in results
)

negative_count = sum(
    row["D"] < 0
    for row in results
)

print()

print(f"positive D(T) = {positive_count}")
print(f"zero D(T)     = {zero_count}")
print(f"negative D(T) = {negative_count}")

if negative_count:
    print()
    print(
        "WARNING: negative dyadic interval contributions detected."
    )

    for row in results:
        if row["D"] < 0:
            print(
                f"  level {row['level']}: "
                f"T={nstr(row['T'], 15)}, "
                f"D={nstr(row['D'], 30)}"
            )
else:
    print()
    print(
        "All computed dyadic interval contributions are positive."
    )

# ============================================================
# 4. STABILITY OF C_T
# ============================================================

print()
print("-" * 120)
print("4. STABILITY OF C_T = T D(T)")
print("-" * 120)

print()

print(
    "For a true D(T) ~ C/T asymptotic regime, C_T should"
)

print(
    "approach a constant."
)

print()

print(
    "The following table reports the fractional change"
)

print(
    "between successive C_T values:"
)

print()

print(
    " level"
    "             T"
    "                 C_T"
    "                 fractional change"
)

print("-" * 120)

for i, row in enumerate(results):

    if i == 0:
        fractional_change = mp.nan
    else:
        previous = results[i - 1]["C_T"]

        if previous == 0:
            fractional_change = mp.nan
        else:
            fractional_change = (
                row["C_T"] / previous
                - 1
            )

    print(
        f"{row['level']:6d} "
        f"{nstr(row['T'], 12):>18} "
        f"{nstr(row['C_T'], 34):>38} "
        f"{nstr(fractional_change, 20):>28}"
    )

# ============================================================
# 5. LATE C ESTIMATES
# ============================================================

print()
print("-" * 120)
print("5. LATE ESTIMATES OF C")
print("-" * 120)

print()

print(
    "If C_T has entered an asymptotic plateau, estimates of C"
)

print(
    "from different late windows should agree."
)

print()

for window in LATE_WINDOWS:

    if len(results) < window:
        continue

    late = [
        row["C_T"]
        for row in results[-window:]
    ]

    mean_C = arithmetic_mean(late)
    median_C = median(late)

    first = late[0]
    last = late[-1]

    if first != 0:
        drift = last / first - 1
    else:
        drift = mp.nan

    print(
        f"last {window:2d} points:"
    )

    print(
        f"    mean(C_T)   = {nstr(mean_C, 40)}"
    )

    print(
        f"    median(C_T) = {nstr(median_C, 40)}"
    )

    print(
        f"    endpoint drift = {nstr(drift, 24)}"
    )

    print()

# Use the median of the final 5 points as the principal
# empirical coefficient for the extrapolation below.
#
# Median is deliberately preferred over the final point so that
# one anomalous numerical interval cannot dominate the estimate.

C_est = median(
    [
        row["C_T"]
        for row in results[-5:]
    ]
)

print(
    "Principal empirical coefficient:"
)

print(
    f"    C_est = {nstr(C_est, 50)}"
)

# ============================================================
# 6. IMPLIED REMOTE TAIL
# ============================================================

print()
print("-" * 120)
print("6. IMPLIED INFINITE-TAIL ESTIMATE")
print("-" * 120)

print()

print(
    "Under the explicit hypothesis"
)

print()

print(
    "    D(T) ~ C/T,"
)

print()

print(
    "the remaining tail is"
)

print()

print(
    "    R(T) = integral_T^infinity I(r) dr"
)

print()

print(
    "    R(T) ~ 2C/T."
)

print()

print(
    "This is an EXTRAPOLATION, not a directly computed integral."
)

print()

print(
    "       T"
    "                    computed D(T)"
    "                    2 C_est / T"
)

print("-" * 120)

for row in results:

    R_est = (
        2 * C_est
        / row["T"]
    )

    print(
        f"{nstr(row['T'], 12):>14} "
        f"{nstr(row['D'], 32):>38} "
        f"{nstr(R_est, 32):>38}"
    )

# ============================================================
# 7. DIRECTLY COMPUTED REMOTE TAIL
# ============================================================

print()
print("-" * 120)
print("7. DIRECTLY COMPUTED FINITE REMOTE TAIL")
print("-" * 120)

print()

print(
    "For each cutoff T, sum the actually computed dyadic"
)

print(
    "intervals from T to the final endpoint."
)

print()

print(
    "This is"
)

print()

print(
    "    R(T,Tmax) = integral_T^Tmax I(r) dr."
)

print()

print(
    "It can be compared with the extrapolated 2C_est/T."
)

print()

print(
    "       T"
    "                    direct finite tail"
    "                    extrapolated infinite tail"
    "                    ratio"
)

print("-" * 120)

for i, row in enumerate(results):

    direct_tail = mp.fsum(
        item["D"]
        for item in results[i:]
    )

    extrapolated_tail = (
        2 * C_est
        / row["T"]
    )

    if extrapolated_tail != 0:
        ratio = (
            direct_tail
            / extrapolated_tail
        )
    else:
        ratio = mp.nan

    print(
        f"{nstr(row['T'], 12):>14} "
        f"{nstr(direct_tail, 34):>40} "
        f"{nstr(extrapolated_tail, 34):>40} "
        f"{nstr(ratio, 20):>24}"
    )

# ============================================================
# 8. HOW MUCH OF THE EXPECTED TAIL IS ALREADY RESOLVED?
# ============================================================

print()
print("-" * 120)
print("8. DIRECT / EXTRAPOLATED TAIL FRACTION")
print("-" * 120)

print()

print(
    "If the D(T) ~ C/T hypothesis is correct, then"
)

print(
    "the direct finite tail should approach the extrapolated"
)

print(
    "infinite tail from below as Tmax increases."
)

print()

print(
    "The ratio"
)

print(
    "    direct finite tail / extrapolated infinite tail"
)

print(
    "therefore gives a useful convergence diagnostic."
)

print()

for i, row in enumerate(results):

    direct_tail = mp.fsum(
        item["D"]
        for item in results[i:]
    )

    extrapolated_tail = (
        2 * C_est
        / row["T"]
    )

    if extrapolated_tail != 0:

        fraction = (
            direct_tail
            / extrapolated_tail
        )

        unresolved_fraction = (
            1 - fraction
        )

        print(
            f"T={nstr(row['T'], 12):>14} "
            f"resolved={nstr(fraction, 24):>28} "
            f"unresolved={nstr(unresolved_fraction, 24):>28}"
        )

# ============================================================
# 9. CONSISTENCY OF THE p_eff AND C_T DIAGNOSTICS
# ============================================================

print()
print("-" * 120)
print("9. ASYMPTOTIC CONSISTENCY CHECK")
print("-" * 120)

print()

print(
    "A genuine p=1 regime should show BOTH:"
)

print()

print(
    "    p_eff -> 1"
)

print()

print(
    "and"
)

print()

print(
    "    C_T = T D(T) -> constant."
)

print()

print(
    "The two diagnostics are independent views of the same"
)

print(
    "scaling behaviour."
)

print()

# Report the final five p_eff values for convenience.

print("Final p_eff values:")

for row in results[-5:]:

    print(
        f"    T={nstr(row['T'], 14):>16} "
        f"p_eff={nstr(row['p_eff'], 24)}"
    )

print()

print("Final C_T values:")

for row in results[-5:]:

    print(
        f"    T={nstr(row['T'], 14):>16} "
        f"C_T={nstr(row['C_T'], 32)}"
    )

# ============================================================
# 10. TIMINGS
# ============================================================

print()
print("-" * 120)
print("10. INTERVAL TIMINGS")
print("-" * 120)

print()

print(
    "       level"
    "                    T"
    "                    seconds"
)

print("-" * 120)

for row in results:

    print(
        f"{row['level']:>12d} "
        f"{nstr(row['T'], 14):>20} "
        f"{row['elapsed']:>28.6f}"
    )

# ============================================================
# 11. FINAL SUMMARY
# ============================================================

print()
print("=" * 120)
print("CELL 30 COMPLETE")
print("=" * 120)

print()

print(
    "Primary question:"
)

print(
    "Does C_T = T D(T) approach a constant?"
)

print()

print(
    "If yes, the data support"
)

print(
    "    D(T) ~ C/T"
)

print(
    "and consequently"
)

print(
    "    R(T) ~ 2C/T."
)

print()

print(
    "The latter remains an extrapolation until the direct"
)

print(
    "finite-tail sums demonstrate quantitative agreement."
)

print()

print(
    "Do NOT infer p=1 merely because the last few p_eff values"
)

print(
    "are close to one.  Require corresponding stabilisation"
)

print(
    "of C_T and agreement between the direct and extrapolated"
)

print(
    "remote tails."
)

print()

print(
    "The final C_est reported above is the median of the last"
)

print(
    "five C_T values and is used only as a robust diagnostic"
)

print(
    "coefficient for the tail extrapolation."
)
