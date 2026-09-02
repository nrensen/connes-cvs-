# ============================================================
# CELL 26 — LONG-RANGE FORENSIC ARCHIMEDEAN TAIL
#
# Purpose
# -------
# Investigate the behaviour of the Archimedean quadratic
# functional for the FIXED forensic Galerkin ground state as
# the upper limit T of the r-integral becomes large.
#
# IMPORTANT:
#
# The forensic ground state is held completely fixed.
#
#     c = 13
#     N = 8
#     Galerkin T = 400
#     generation dps = 150
#
# The T varied below is ONLY the upper limit of
#
#     A_arch(T)
#       = 1/pi * integral_0^T h_+(r) J_v(r) dr
#
# No Galerkin matrix is rebuilt as T changes.
#
# No numerical y-integration is performed.
#
# This is the first post-forensic investigation using the
# analytic Archimedean machinery.
#
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


# ------------------------------------------------------------
# T GRID
#
# Dense enough at low T to connect with Cell 24, then
# progressively coarser as we investigate the long tail.
#
# The final range is deliberately well beyond the historical
# T=400 forensic calculation.
# ------------------------------------------------------------

T_VALUES = (
    list(range(5, 105, 5))
    + [
        120,
        150,
        200,
        250,
        300,
        350,
        400,
        500,
        600,
        700,
        800,
        900,
        1000,
        1200,
        1400,
        1600,
        1800,
        2000,
        2250,
        2500,
        2750,
        3000,
        3500,
        4000,
        4500,
        5000,
        6000,
        7000,
        8000,
        10000,
    ]
)


mp.mp.dps = WORKING_DPS


# ============================================================
# PARAMETERS DERIVED FROM FORENSIC GROUND STATE
# ============================================================

c = FORENSIC_GROUND_STATE["c"]
N = FORENSIC_GROUND_STATE["N"]

L = compute_L(c)

DISPLAY_DIGITS = WORKING_DPS


# ============================================================
# FORMATTING / TIMING HELPERS
# ============================================================


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

print("=" * 90)
print("CELL 26 — LONG-RANGE FORENSIC ARCHIMEDEAN TAIL")
print("=" * 90)

print()
print("Purpose:")
print(
    "Investigate A_arch(T) for the fixed forensic ground state "
    "as T becomes large."
)

print()
print("Parameters:")
print(f" c = {c}")
print(f" N = {N}")
print(f" working_dps = {WORKING_DPS}")

print()
print("T values:")
print(f" {T_VALUES}")

print()
print("Forensic ground state:")
print(
    f" c = {FORENSIC_GROUND_STATE['c']}"
)
print(
    f" N = {FORENSIC_GROUND_STATE['N']}"
)
print(
    f" Galerkin T = {FORENSIC_GROUND_STATE['T']}"
)
print(
    f" generation_dps = "
    f"{FORENSIC_GROUND_STATE['dps']}"
)

print()
print("L =")
print(nstr(L))


# ============================================================
# 1. RETRIEVE FIXED FORENSIC GROUND STATE
# ============================================================

print()
print("-" * 90)
print("1. FORENSIC GROUND STATE")
print("-" * 90)

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
    f"ground-state retrieval CPU = "
    f"{ground_cpu_elapsed:.6f} s"
)


norm_v = mp.sqrt(
    mp.fdot(
        v_star,
        v_star,
    )
)


print()
print("||v_star|| =")
print(nstr(norm_v))


# ============================================================
# 2. LONG-RANGE FINITE-T SCAN
# ============================================================

print()
print("-" * 90)
print("2. LONG-RANGE FINITE-T ARCHIMEDEAN SCAN")
print("-" * 90)

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
    " T"
    " A_arch(T)"
    " delta"
    " |delta|"
    " wall_s"
    " cpu_s"
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
        delta = (
            value
            - previous_value
        )


    abs_delta = abs(delta)


    results.append(
        (
            T,
            value,
            delta,
            abs_delta,
            wall_elapsed,
            cpu_elapsed,
        )
    )


    print(
        f"{T:6d} "
        f"{nstr(value):>{DISPLAY_DIGITS + 4}} "
        f"{nstr(delta):>{DISPLAY_DIGITS + 4}} "
        f"{nstr(abs_delta):>{DISPLAY_DIGITS + 4}} "
        f"{wall_elapsed:12.6f} "
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
# 3. SELECTED LONG-RANGE VALUES
# ============================================================

print()
print("-" * 90)
print("3. SELECTED LONG-RANGE VALUES")
print("-" * 90)

print()
print(
    "This section repeats the tail in compact form, making "
    "the large-T behaviour easier to inspect."
)

print()
print(
    " T"
    " A_arch(T)"
    " delta"
    " |delta|"
)

print("-" * 90)


for (
    T,
    value,
    delta,
    abs_delta,
    wall,
    cpu,
) in results:

    if T >= 400:

        print(
            f"{T:6d} "
            f"{nstr(value):>{DISPLAY_DIGITS + 4}} "
            f"{nstr(delta):>{DISPLAY_DIGITS + 4}} "
            f"{nstr(abs_delta):>{DISPLAY_DIGITS + 4}}"
        )


# ============================================================
# 4. FINAL RESULT
# ============================================================

final_T = results[-1][0]
final_value = results[-1][1]
final_delta = results[-1][2]


print()
print("-" * 90)
print("4. FINAL RESULT")
print("-" * 90)

print()
print("Final T =")
print(final_T)

print()
print("A_arch(final) =")
print(nstr(final_value))

print()
print("delta at final T =")
print(nstr(final_delta))

print()
print("|delta| at final T =")
print(nstr(abs(final_delta)))


# ============================================================
# 5. TIMING SUMMARY
# ============================================================

print()
print("-" * 90)
print("5. TIMING SUMMARY")
print("-" * 90)

print()
print(
    f"finite-T scan wall = "
    f"{total_wall_elapsed:.6f} s"
)

print(
    f"finite-T scan CPU = "
    f"{total_cpu_elapsed:.6f} s"
)


# ============================================================
# 6. FINAL STATEMENT
# ============================================================

print()
print("=" * 90)
print("CELL 26 COMPLETE")
print("=" * 90)
