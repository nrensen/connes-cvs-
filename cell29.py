# ============================================================
# CELL 29 — LOG-SCALE INTEGRATED TAIL SCALING
#
# Purpose
# -------
# Probe the asymptotic behaviour of the INTEGRATED tail without
# assuming a power law in advance.
#
# Define
#
#     I(r) = h_+(r) K_fourier(v_star, r, L)
#
# and, for geometrically expanding intervals,
#
#     D(T) = integral_T^(2T) I(r) dr.
#
# Cell 28 showed that the corresponding finite intervals remain
# positive while decaying slowly.  Cell 29 therefore asks:
#
#     How does D(T) itself scale?
#
# The principal diagnostic is the ratio
#
#     D(2T) / D(T).
#
# If eventually
#
#     D(T) ~ C T^(-p),
#
# then
#
#     D(2T) / D(T) -> 2^(-p).
#
# No value of p is assumed by the calculation.
#
# We additionally display T^p D(T) for several diagnostic
# exponents, but these columns are explicitly diagnostic only.
#
# A second diagnostic is the cumulative integral
#
#     C(T) = integral_T0^T I(r) dr.
#
# If C(T) approaches a finite value, the tail beyond T0 is
# integrable.  The difference between successive C(T) values is
# exactly the corresponding D(T) when the endpoints are powers
# of two.
#
# IMPORTANT:
#
# The forensic ground state is fixed.  No eigensolve is performed.
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
# Geometric tail intervals.
#
# Starting at 20 keeps this directly comparable with Cell 28.
#
# Each interval is [T, 2T], so successive intervals tile the
# positive-r axis:
#
#   [20,40], [40,80], [80,160], ...
#
# This is intentional.  It lets the cumulative sum be built
# directly from the logarithmic-scale increments.
# ------------------------------------------------------------

T0 = mp.mpf("20")
N_LEVELS = 10

# Diagnostic powers only.
#
# These are NOT assumptions about the asymptotic behaviour.
# p = 0.5, 1, 1.5 and 2 are useful reference cases.
#
# If one column becomes approximately constant over several
# successive levels, that provides an empirical indication of
# the corresponding power.
DIAGNOSTIC_POWERS = [
    mp.mpf("0"),
    mp.mpf("0.5"),
    mp.mpf("1"),
    mp.mpf("1.5"),
    mp.mpf("2"),
]

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
    Signed integral over one finite logarithmic interval.
    """
    return mp.quad(
        lambda r: integrand(r, v_star),
        [a, b],
    )


def local_power(D1, D2):
    """
    Empirical local power p defined by

        D(2T) / D(T) = 2^(-p).

    Hence

        p = -log_2(D(2T)/D(T)).

    This is only meaningful for positive D1 and D2.
    """
    if D1 <= 0 or D2 <= 0:
        return mp.nan

    return (
        -mp.log(D2 / D1)
        / mp.log(2)
    )


# ============================================================
# HEADER
# ============================================================

print("=" * 120)
print("CELL 29 — LOG-SCALE INTEGRATED TAIL SCALING")
print("=" * 120)

print()

print("Parameters:")
print(f" c = {c}")
print(f" N = {N}")
print(f" working_dps = {WORKING_DPS}")
print(f" L = {nstr(L, 60)}")
print(f" T0 = {nstr(T0, 20)}")
print(f" number of geometric intervals = {N_LEVELS}")

print()

print(
    "Principal quantity:"
)

print(
    " D(T) = integral_T^(2T) h_+(r) K_fourier(v_star,r,L) dr"
)

print()

print(
    "No asymptotic power law is assumed."
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
# 2. GEOMETRIC INTERVAL INTEGRALS
# ============================================================

print()
print("-" * 120)
print("2. GEOMETRIC INTERVAL INTEGRALS")
print("-" * 120)

print()

print(
    "Each row is one interval [T,2T]."
)

print()

print(
    "D(T) is the SIGNED integral."
)

print()

print(
    "The cumulative column is"
)

print(
    "C(T) = integral_20^(2T) I(r) dr."
)

print()

print(
    "     T"
    "          2T"
    "                    D(T)"
    "                    |D(T)|"
    "                    C(2T)"
)

print("-" * 120)

results = []

cumulative = mp.mpf("0")

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

    cumulative += D

    results.append(
        {
            "level": level,
            "T": a,
            "2T": b,
            "D": D,
            "cumulative": cumulative,
            "elapsed": elapsed,
        }
    )

    print(
        f"{nstr(a, 10):>10} "
        f"{nstr(b, 10):>10} "
        f"{nstr(D, 36):>42} "
        f"{nstr(abs(D), 28):>30} "
        f"{nstr(cumulative, 36):>42}"
    )

    T = b

# ============================================================
# 3. LOCAL EMPIRICAL POWER
# ============================================================

print()
print("-" * 120)
print("3. LOCAL EMPIRICAL POWER")
print("-" * 120)

print()

print(
    "For adjacent intervals:"
)

print()

print(
    "    p_eff(T) = -log_2(D(2T)/D(T))."
)

print()

print(
    "If D(T) ~ C T^(-p), then p_eff(T) -> p."
)

print()

print(
    "This is an OBSERVATION of the computed data, not an"
)

print(
    "assumption built into the integration."
)

print()

print(
    "       T"
    "                    D(T)"
    "                    D(2T)/D(T)"
    "                    p_eff"
)

print("-" * 120)

for i in range(len(results)):

    row = results[i]

    if i == 0:
        print(
            f"{nstr(row['T'], 10):>10} "
            f"{nstr(row['D'], 30):>34} "
            f"{'--':>30} "
            f"{'--':>24}"
        )
        continue

    previous = results[i - 1]

    ratio = (
        row["D"]
        / previous["D"]
    )

    p_eff = local_power(
        previous["D"],
        row["D"],
    )

    print(
        f"{nstr(row['T'], 10):>10} "
        f"{nstr(row['D'], 30):>34} "
        f"{nstr(ratio, 24):>30} "
        f"{nstr(p_eff, 24):>24}"
    )

# ============================================================
# 4. DIAGNOSTIC POWER SCALES
# ============================================================

print()
print("-" * 120)
print("4. DIAGNOSTIC POWER SCALES")
print("-" * 120)

print()

print(
    "The following columns show T^p D(T)."
)

print()

print(
    "They are diagnostic only; no asymptotic exponent is"
)

print(
    "assumed by the calculation."
)

print()

header = (
    f"{'T':>10}"
)

for p in DIAGNOSTIC_POWERS:
    header += (
        f"{'T^' + nstr(p, 4) + ' D(T)':>30}"
    )

print(header)

print("-" * 120)

for row in results:

    line = f"{nstr(row['T'], 10):>10}"

    for p in DIAGNOSTIC_POWERS:

        scaled = (
            row["T"] ** p
            * row["D"]
        )

        line += (
            f"{nstr(scaled, 24):>30}"
        )

    print(line)

# ============================================================
# 5. ABSOLUTE CONTRIBUTION PER LOG INTERVAL
# ============================================================

print()
print("-" * 120)
print("5. CONTRIBUTION PER LOGARITHMIC INTERVAL")
print("-" * 120)

print()

print(
    "Every interval has the same logarithmic width:"
)

print()

print(
    "log(2T) - log(T) = log(2)."
)

print()

print(
    "Therefore D(T)/log(2) is the average contribution per"
)

print(
    "unit of log(r) over that interval."
)

print()

print(
    "This is useful because a slowly decaying D(T) corresponds"
)

print(
    "directly to a slowly diminishing contribution on a"
)

print(
    "logarithmic r-axis."
)

print()

print(
    "       T"
    "                    D(T)/log(2)"
)

print("-" * 120)

log2 = mp.log(2)

for row in results:

    contribution_per_log = (
        row["D"]
        / log2
    )

    print(
        f"{nstr(row['T'], 10):>10} "
        f"{nstr(contribution_per_log, 36):>42}"
    )

# ============================================================
# 6. REMOTE TAIL SUMS
# ============================================================

print()
print("-" * 120)
print("6. REMOTE TAIL SUMS")
print("-" * 120)

print()

print(
    "For each cutoff T0, sum all computed intervals beginning"
)

print(
    "at that cutoff."
)

print()

print(
    "This gives a direct finite-window estimate of"
)

print(
    "integral_T0^(Tmax) I(r) dr."
)

print()

print(
    "      cutoff"
    "                    finite tail to Tmax"
)

print("-" * 120)

for i, row in enumerate(results):

    tail_sum = mp.fsum(
        item["D"]
        for item in results[i:]
    )

    print(
        f"{nstr(row['T'], 10):>12} "
        f"{nstr(tail_sum, 40):>46}"
    )

# ============================================================
# 7. POSITIVITY CHECK
# ============================================================

print()
print("-" * 120)
print("7. SIGN CHECK")
print("-" * 120)

print()

negative_intervals = [
    row
    for row in results
    if row["D"] < 0
]

zero_intervals = [
    row
    for row in results
    if row["D"] == 0
]

print(
    f"positive intervals = "
    f"{len(results) - len(negative_intervals) - len(zero_intervals)}"
)

print(
    f"zero intervals = "
    f"{len(zero_intervals)}"
)

print(
    f"negative intervals = "
    f"{len(negative_intervals)}"
)

if negative_intervals:
    print()

    print(
        "WARNING: negative interval contributions detected:"
    )

    for row in negative_intervals:
        print(
            f"  [{nstr(row['T'], 12)}, "
            f"{nstr(row['2T'], 12)}] "
            f"{nstr(row['D'], 30)}"
        )
else:
    print()

    print(
        "All computed logarithmic interval contributions are "
        "non-negative."
    )

# ============================================================
# 8. TIMINGS
# ============================================================

print()
print("-" * 120)
print("8. INTERVAL TIMINGS")
print("-" * 120)

print()

print(
    "       T"
    "                    seconds"
)

print("-" * 120)

for row in results:

    print(
        f"{nstr(row['T'], 10):>10} "
        f"{row['elapsed']:>28.6f}"
    )

# ============================================================
# 9. FINAL SUMMARY
# ============================================================

print()
print("=" * 120)
print("CELL 29 COMPLETE")
print("=" * 120)

print()

print(
    "The primary diagnostic is D(T) = integral_T^(2T) I(r) dr."
)

print()

print(
    "If the local empirical power p_eff settles toward a stable"
)

print(
    "positive value, this provides evidence for a corresponding"
)

print(
    "power-law decay of the INTEGRATED logarithmic interval."
)

print()

print(
    "If p_eff continues drifting substantially, no asymptotic"
)

print(
    "power law should yet be inferred."
)

print()

print(
    "The remote-tail sums show directly how much of the finite"
)

print(
    "integrated contribution remains beyond each cutoff."
)

print()

print(
    "In particular, positivity of D(T) means that the cumulative"
)

print(
    "tail cannot be reduced by cancellation between successive"
)

print(
    "geometric intervals."
)

print()

print(
    "Any apparent convergence must therefore come from genuine"
)

print(
    "decay of the integrated interval contribution."
)
