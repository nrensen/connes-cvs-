# ============================================================
# CELL 25 — FINITE-T ARCHIMEDEAN CONVERGENCE
#
# Purpose
# -------
# Re-evaluate the finite-T Archimedean functional using the
# validated analytic K_fourier reduction, first reproducing
# the key finite-T values from cell5_corrected and then
# extending the calculation to substantially larger T.
#
# The ground state is FIXED throughout the experiment.
#
# Primary cross-check configuration:
#
#     c = 13
#     N = 8
#     ground-state T = 60
#     generation dps = 80
#
# These are the parameters used by cell5_corrected.
#
# The T varied below is ONLY the upper limit of the
# Archimedean r-integral:
#
#     A_arch(T)
#       = 1/pi * int_0^T h_+(r) J_v(r) dr
#
# No Galerkin matrix is rebuilt as T changes.
#
# This cell therefore tests two things:
#
#   1. Does the efficient analytic calculation reproduce the
#      finite-T values obtained by cell5_corrected?
#
#   2. What happens to A_arch(T) as T is increased substantially
#      beyond the historical T=200 range?
#
# The ground state remains fixed while T is varied.
# ============================================================

import time

import mpmath as mp

from cell import (
    FORENSIC_GROUND_STATE,
    get_ground_state,
    compute_L,
    archimedean_integral,
)


# ============================================================
# PARAMETERS
# ============================================================

WORKING_DPS = 80

GROUND_STATE = {
    "c": 13,
    "N": 8,
    "T": 60,
    "dps": 80,
}

# Historical cell5_corrected values to reproduce.
#
# These are the explicit finite-T Archimedean values obtained
# from the older direct source-level integration.
#
# T = 20:
#   0.7125875193025273831377680272
#
# T = 40:
#   0.7421874519938487471064025458
#
# Keep these as strings so that they are not rounded before
# the comparison.
CELL5_REFERENCE = {
    20: mp.mpf(
        "0.7125875193025273831377680272"
    ),
    40: mp.mpf(
        "0.7421874519938487471064025458"
    ),
}

# Dense enough to see the initial convergence, followed by
# progressively larger cutoffs.
#
# The larger-T values are the main new experiment.
T_VALUES = (
    list(range(20, 101, 10))
    + [120, 150, 200, 250, 300, 400,
       600, 800, 1200, 1600, 2000]
)

mp.mp.dps = WORKING_DPS

c = GROUND_STATE["c"]
N = GROUND_STATE["N"]
L = compute_L(c)

DISPLAY_DIGITS = WORKING_DPS


def nstr(x):
    return mp.nstr(
        x,
        DISPLAY_DIGITS,
    )


def elapsed(start):
    return time.perf_counter() - start


# ============================================================
# HEADER
# ============================================================

print("=" * 78)
print("CELL 25 — FINITE-T ARCHIMEDEAN CONVERGENCE")
print("=" * 78)

print()
print("Ground-state parameters:")
print(f"  c              = {c}")
print(f"  N              = {N}")
print(f"  ground T       = {GROUND_STATE['T']}")
print(f"  generation dps = {GROUND_STATE['dps']}")
print(f"  working dps    = {WORKING_DPS}")

print()
print("L =")
print(nstr(L))

print()
print("T values:")
print(f"  {T_VALUES}")

print()
print("The ground state is held fixed.")
print("Only the upper limit of the Archimedean r-integral")
print("is varied.")
print("No Galerkin matrix is rebuilt as T changes.")


# ============================================================
# 1. GROUND STATE
# ============================================================

print()
print("-" * 78)
print("1. FIXED GROUND STATE")
print("-" * 78)

ground_wall_start = time.perf_counter()
ground_cpu_start = time.process_time()

lambda_ground, v_star, ground_meta = get_ground_state(
    **GROUND_STATE,
    verbose=True,
)

ground_wall_elapsed = elapsed(
    ground_wall_start
)

ground_cpu_elapsed = (
    time.process_time()
    - ground_cpu_start
)

print()
print("lambda_ground =")
print(nstr(lambda_ground))

print()
print("||v_star|| =")
print(
    nstr(
        mp.sqrt(
            mp.fdot(
                v_star,
                v_star,
            )
        )
    )
)

print()
print(
    f"ground-state retrieval wall = "
    f"{ground_wall_elapsed:.6f} s"
)

print(
    f"ground-state retrieval CPU  = "
    f"{ground_cpu_elapsed:.6f} s"
)


# ============================================================
# 2. FINITE-T SCAN
# ============================================================

print()
print("-" * 78)
print("2. FINITE-T ARCHIMEDEAN CONVERGENCE")
print("-" * 78)

print()
print(
    "Columns:"
)

print(
    "  T"
    "              A_arch(T)"
    "                          delta"
    "                 wall_s"
    "          cpu_s"
)

results = []

previous_value = None

total_wall_start = time.perf_counter()
total_cpu_start = time.process_time()

for T in T_VALUES:

    wall_start = time.perf_counter()
    cpu_start = time.process_time()

    value = archimedean_integral(
        T,
        v_star,
        L,
    )

    wall_elapsed = (
        time.perf_counter()
        - wall_start
    )

    cpu_elapsed = (
        time.process_time()
        - cpu_start
    )

    if previous_value is None:
        delta = mp.mpf("0")
    else:
        delta = value - previous_value

    results.append(
        (
            T,
            value,
            delta,
            wall_elapsed,
            cpu_elapsed,
        )
    )

    print(
        f"{T:6d}  "
        f"{nstr(value):>{DISPLAY_DIGITS + 4}}  "
        f"{nstr(delta):>{DISPLAY_DIGITS + 4}}  "
        f"{wall_elapsed:12.6f}  "
        f"{cpu_elapsed:12.6f}"
    )

    previous_value = value


total_wall_elapsed = (
    time.perf_counter()
    - total_wall_start
)

total_cpu_elapsed = (
    time.process_time()
    - total_cpu_start
)


# ============================================================
# 3. CELL5 CORRECTED CROSS-CHECK
# ============================================================

print()
print("-" * 78)
print("3. CELL5_CORRECTED CROSS-CHECK")
print("-" * 78)

print()
print(
    "The reference values below were obtained from the"
)
print(
    "older direct finite-T Archimedean calculation in"
)
print(
    "cell5_corrected."
)

print()
print(
    "T"
    "                  analytic A_arch(T)"
    "                         cell5 reference"
    "                         difference"
)

print("-" * 120)

crosscheck_results = {}

for T in sorted(CELL5_REFERENCE):

    analytic = next(
        value
        for (
            T_result,
            value,
            delta,
            wall,
            cpu,
        ) in results
        if T_result == T
    )

    reference = CELL5_REFERENCE[T]

    difference = analytic - reference

    crosscheck_results[T] = (
        analytic,
        reference,
        difference,
    )

    print(
        f"{T:6d}   "
        f"{nstr(analytic):>{DISPLAY_DIGITS + 4}}   "
        f"{nstr(reference):>{DISPLAY_DIGITS + 4}}   "
        f"{nstr(difference):>{DISPLAY_DIGITS + 4}}"
    )


# ============================================================
# 4. CONVERGENCE SUMMARY
# ============================================================

print()
print("-" * 78)
print("4. CONVERGENCE SUMMARY")
print("-" * 78)

print()
print(
    "T                  A_arch(T)"
)

for T, value, delta, wall, cpu in results:

    print(
        f"{T:6d}   "
        f"{nstr(value)}"
    )

print()
print(
    "Final T =",
    T_VALUES[-1],
)

print()
print("A_arch(final) =")
print(
    nstr(
        results[-1][1]
    )
)

print()
print(
    f"finite-T scan wall = "
    f"{total_wall_elapsed:.6f} s"
)

print(
    f"finite-T scan CPU  = "
    f"{total_cpu_elapsed:.6f} s"
)


# ============================================================
# 5. FINAL RESULT
# ============================================================

print()
print("=" * 78)
print("CELL 25 COMPLETE")
print("=" * 78)
