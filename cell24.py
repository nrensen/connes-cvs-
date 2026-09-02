# ============================================================
# CELL 24 — FINITE-T ARCHIMEDEAN CONVERGENCE MAP
#
# Purpose
# -------
# Investigate the finite-T behaviour of the Archimedean
# quadratic functional using the validated Cell-23 analytic
# reduction.
#
# IMPORTANT:
#
# The forensic ground state is held FIXED.
#
# The ground state is the canonical forensic ground state:
#
#     c = 13
#     N = 8
#     T = 400
#     generation_dps = 150
#
# The T varied in this cell is ONLY the upper limit of the
# Archimedean r-integral:
#
#     A_arch(T)
#       = 1/pi * int_0^T h_+(r) J_v(r) dr
#
# No Galerkin matrix is rebuilt as T changes.
#
# This isolates finite-T behaviour of the Archimedean functional
# for one fixed ground state.
#
# No numerical y-integration is performed.
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

WORKING_DPS = 40

# Initial finite-T survey.
#
# The relatively dense low-T grid is intended to show the
# shape of convergence rather than merely provide isolated
# endpoint values.
T_VALUES = (
    list(range(5, 105, 5))
    + [120, 150, 200, 250, 300, 350, 400]
)

mp.mp.dps = WORKING_DPS

c = FORENSIC_GROUND_STATE["c"]
N = FORENSIC_GROUND_STATE["N"]

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
print("CELL 24 — FINITE-T ARCHIMEDEAN CONVERGENCE MAP")
print("=" * 78)

print()
print("Parameters:")
print(f"  c              = {c}")
print(f"  N              = {N}")
print(f"  working_dps    = {WORKING_DPS}")

print()
print("T values:")
print(f"  {T_VALUES}")

print()
print("Forensic ground state:")
print(
    f"  c              = "
    f"{FORENSIC_GROUND_STATE['c']}"
)
print(
    f"  N              = "
    f"{FORENSIC_GROUND_STATE['N']}"
)
print(
    f"  T              = "
    f"{FORENSIC_GROUND_STATE['T']}"
)
print(
    f"  generation_dps = "
    f"{FORENSIC_GROUND_STATE['dps']}"
)

print()
print("L =")
print(nstr(L))


# ============================================================
# 1. GROUND STATE
# ============================================================

print()
print("-" * 78)
print("1. GROUND STATE")
print("-" * 78)

ground_wall_start = time.perf_counter()
ground_cpu_start = time.process_time()

lambda_forensic, v_star, ground_meta = get_ground_state(
    **FORENSIC_GROUND_STATE,
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
print("lambda_forensic =")
print(nstr(lambda_forensic))

print()
print(
    f"ground-state retrieval wall = "
    f"{ground_wall_elapsed:.6f} s"
)
print(
    f"ground-state retrieval CPU  = "
    f"{ground_cpu_elapsed:.6f} s"
)

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

# ============================================================
# 9. FINITE-T SCAN
# ============================================================

print()
print("-" * 78)
print("3. FINITE-T ARCHIMEDEAN CONVERGENCE")
print("-" * 78)

print()
print(
    "The ground state is held fixed while only the upper "
    "limit T of the Archimedean r-integral is varied."
)

print()
print(
    "No Galerkin matrix is rebuilt as T changes."
)

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

    value = archimedean_integral(T, v_star, L)

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
# 10. SUMMARY
# ============================================================

print()
print("-" * 78)
print("4. FINITE-T SUMMARY")
print("-" * 78)

print()
print("T                  A_arch(T)")

for T, value, delta, wall, cpu in results:

    print(
        f"{T:6d}   {nstr(value)}"
    )

print()
print(
    "Final T =",
    T_VALUES[-1],
)

print()
print(
    "A_arch(final) ="
)
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
# 11. FINAL RESULT
# ============================================================

print()
print("=" * 78)
print("CELL 24 COMPLETE")
print("=" * 78)
