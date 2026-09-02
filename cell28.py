# ============================================================
# CELL 28 — DIRECT INTEGRATED ARCHIMEDEAN TAIL
#
# Purpose
# -------
# Integrate the actual analytic Archimedean integrand over
# successive large-r intervals for the FIXED forensic Galerkin
# ground state.
#
# The quantity being integrated is
#
#     I(r) = h_+(r) * K_fourier(v_star, r, L)
#
# This cell deliberately makes NO asymptotic power-law
# assumption.
#
# Instead it measures:
#
#     A(a,b) = integral_a^b I(r) dr
#
# over successive intervals and compares the cumulative
# integral against the interval increments.
#
# The central question is:
#
#     Does the apparent large-T ~1e-24 behaviour seen in
#     Cell 26 arise from a genuinely persistent tail, or
#     from an oscillatory integrand whose signed contributions
#     increasingly cancel?
#
# IMPORTANT:
#
# The forensic ground state is fixed:
#
#     c = 13
#     N = 8
#     Galerkin T = 400
#     generation dps = 150
#
# No ground-state regeneration is performed here.
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
# Integration intervals.
#
# The intervals deliberately become progressively wider.
#
# This lets us distinguish:
#
#   * local oscillation,
#   * interval cancellation,
#   * persistence of the cumulative tail.
#
# The first interval starts at r = 20 so that the output can
# also be compared against the region examined in Cell 27.
# ------------------------------------------------------------

INTERVALS = [
    (20, 40),
    (40, 60),
    (60, 80),
    (80, 100),
    (100, 150),
    (150, 200),
    (200, 300),
    (300, 400),
    (400, 500),
    (500, 600),
    (600, 800),
    (800, 1000),
    (1000, 1200),
    (1200, 1600),
    (1600, 2000),
    (2000, 3000),
    (3000, 4000),
    (4000, 5000),
    (5000, 6000),
    (6000, 8000),
    (8000, 10000),
]

# ============================================================
# HELPERS
# ============================================================

def nstr(x, digits=WORKING_DPS):
    return mp.nstr(x, digits)


def integrand(r, v_star):
    J = K_fourier(
        v_star,
        r,
        L,
    )

    h = h_plus(r)

    return h * J


def integrate_interval(a, b, v_star):
    """
    Compute the signed integral

        integral_a^b h_+(r) K_fourier(v_star,r,L) dr

    using mpmath's oscillatory-aware tanh-sinh/quadrature
    machinery.

    The interval is kept finite so that the result can be
    interpreted directly as a tail increment.
    """

    return mp.quad(
        lambda r: integrand(r, v_star),
        [mp.mpf(a), mp.mpf(b)],
    )


# ============================================================
# HEADER
# ============================================================

print("=" * 120)
print("CELL 28 — DIRECT INTEGRATED ARCHIMEDEAN TAIL")
print("=" * 120)

print()

print("Purpose:")
print(
    "Integrate the actual analytic Archimedean integrand "
    "over successive large-r intervals."
)

print()

print("Parameters:")
print(f" c = {c}")
print(f" N = {N}")
print(f" working_dps = {WORKING_DPS}")
print(f" L = {nstr(L, 60)}")

print()

print("No asymptotic power-law assumption is made.")

print()
print(
    "The reported interval integrals are SIGNED quantities."
)

# ============================================================
# 1. RETRIEVE FIXED FORENSIC GROUND STATE
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
# 2. INTEGRATED INTERVALS
# ============================================================

print()
print("-" * 120)
print("2. SUCCESSIVE SIGNED TAIL INTERVALS")
print("-" * 120)

print()

print(
    "For each interval [a,b], compute"
)

print()

print(
    " A(a,b) = integral_a^b I(r) dr"
)

print()

print(
    "where I(r) = h_+(r) K_fourier(v_star,r,L)."
)

print()

print(
    "Also report |A| and the cumulative integral from r=20."
)

print()

print(
    "       interval"
    "                    A(a,b)"
    "                    |A(a,b)|"
    "                    cumulative"
)

print("-" * 120)

interval_results = []

cumulative = mp.mpf("0")

for a, b in INTERVALS:

    start = time.perf_counter()

    A = integrate_interval(
        a,
        b,
        v_star,
    )

    elapsed = (
        time.perf_counter()
        - start
    )

    cumulative += A

    interval_results.append(
        (
            mp.mpf(a),
            mp.mpf(b),
            A,
            cumulative,
            elapsed,
        )
    )

    print(
        f"[{a:5d}, {b:5d}] "
        f"{nstr(A, 32):>38} "
        f"{nstr(abs(A), 24):>28} "
        f"{nstr(cumulative, 32):>38}"
    )

# ============================================================
# 3. CUMULATIVE TAIL DIFFERENCES
# ============================================================

print()
print("-" * 120)
print("3. CUMULATIVE TAIL DIFFERENCES")
print("-" * 120)

print()

print(
    "For each endpoint T, define"
)

print()

print(
    " C(T) = integral_20^T I(r) dr."
)

print()

print(
    "If C(T) stabilises while the pointwise integrand"
)

print(
    "continues to oscillate, this is direct evidence for"
)

print(
    "cancellation in the integrated tail."
)

print()

print(
    "     T"
    "                    C(T)"
    "                    |C(T)|"
    "                    last increment"
)

print("-" * 120)

previous_cumulative = mp.mpf("0")

for (
    a,
    b,
    A,
    cumulative,
    elapsed,
) in interval_results:

    increment = cumulative - previous_cumulative

    print(
        f"{int(b):6d} "
        f"{nstr(cumulative, 36):>40} "
        f"{nstr(abs(cumulative), 28):>30} "
        f"{nstr(increment, 28):>32}"
    )

    previous_cumulative = cumulative

# ============================================================
# 4. ABSOLUTE INTERVAL CONTRIBUTIONS
# ============================================================

print()
print("-" * 120)
print("4. MAGNITUDE OF SUCCESSIVE INTERVAL CONTRIBUTIONS")
print("-" * 120)

print()

print(
    "This section ignores the sign and reports |A(a,b)|."
)

print()

print(
    "It is diagnostic only: a decreasing |A| does not by"
)

print(
    "itself establish a particular asymptotic power."
)

print()

print(
    "       interval"
    "                    |A(a,b)|"
)

print("-" * 120)

for (
    a,
    b,
    A,
    cumulative,
    elapsed,
) in interval_results:

    print(
        f"[{int(a):5d}, {int(b):5d}] "
        f"{nstr(abs(A), 36):>42}"
    )

# ============================================================
# 5. TAIL-OF-TAIL COMPARISON
# ============================================================

print()
print("-" * 120)
print("5. NESTED TAIL COMPARISON")
print("-" * 120)

print()

print(
    "Compare cumulative integrals over increasingly remote"
)

print(
    "regions.  These are:"
)

print()

print(
    "   integral_400^T I(r) dr"
)

print(
    "   integral_1000^T I(r) dr"
)

print(
    "   integral_2000^T I(r) dr"
)

print(
    "   integral_4000^T I(r) dr"
)

print(
    "   integral_6000^T I(r) dr"
)

print()

print(
    "These quantities are reconstructed from the interval"
)

print(
    "integrals above, so no additional quadrature is needed."
)

print()

cutoffs = [
    400,
    1000,
    2000,
    4000,
    6000,
]

for cutoff in cutoffs:

    tail_sum = mp.mpf("0")

    for (
        a,
        b,
        A,
        cumulative,
        elapsed,
    ) in interval_results:

        if a >= cutoff:
            tail_sum += A

    print(
        f"T0 = {cutoff:5d}   "
        f"integral_{cutoff}^10000 I(r) dr = "
        f"{nstr(tail_sum, 40)}"
    )

# ============================================================
# 6. CHECK AGAINST POINTWISE SCALE
# ============================================================

print()
print("-" * 120)
print("6. INTERVAL INTEGRAL VERSUS POINTWISE SCALE")
print("-" * 120)

print()

print(
    "For each interval report the average signed integrand"
)

print(
    "A(a,b)/(b-a)."
)

print()

print(
    "This is NOT used as an asymptotic estimate."
)

print()

print(
    "       interval"
    "                    A/(b-a)"
)

print("-" * 120)

for (
    a,
    b,
    A,
    cumulative,
    elapsed,
) in interval_results:

    width = b - a

    average = A / width

    print(
        f"[{int(a):5d}, {int(b):5d}] "
        f"{nstr(average, 36):>42}"
    )

# ============================================================
# 7. TIMING
# ============================================================

print()
print("-" * 120)
print("7. INTERVAL TIMINGS")
print("-" * 120)

print()

print(
    "       interval"
    "                    seconds"
)

print("-" * 120)

for (
    a,
    b,
    A,
    cumulative,
    elapsed,
) in interval_results:

    print(
        f"[{int(a):5d}, {int(b):5d}] "
        f"{elapsed:>28.6f}"
    )

# ============================================================
# 8. FINAL SUMMARY
# ============================================================

print()
print("=" * 120)
print("CELL 28 COMPLETE")
print("=" * 120)

print()

print(
    "The cell has integrated the actual analytic"
)

print(
    "Archimedean integrand over successive finite intervals."
)

print()

print(
    "No asymptotic extrapolation has been performed."
)

print()

print(
    "The principal diagnostic is the behaviour of the"
)

print(
    "signed cumulative integral C(T) as T increases."
)

print()

print(
    "A stabilising C(T), despite continued pointwise"
)

print(
    "oscillation, would support an oscillatory-cancellation"
)

print(
    "interpretation of the Cell-26 tail."
)

print()

print(
    "A persistent drift in C(T) would instead indicate that"
)

print(
    "the integrated tail remains numerically significant."
)
