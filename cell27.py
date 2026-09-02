# ============================================================
# CELL 27 — ARCHIMEDEAN TAIL ANATOMY
#
# Purpose
# -------
# Investigate the large-r behaviour of the analytic
# Archimedean integrand for the FIXED forensic Galerkin
# ground state.
#
# This cell DOES NOT perform a long-range integral.
#
# Instead it examines:
#
#     J(r)       = K_fourier(v,r,L)
#     I(r)       = h_plus(r) * J(r)
#
# together with scaled versions of these quantities.
#
# The purpose is to understand the oscillatory ~1e-24 tail
# seen in Cell 26 before attempting any further large-T
# integration or asymptotic extrapolation.
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
# The r values below are NOT used to regenerate the ground
# state. They are only evaluation points for the analytic
# Archimedean integrand.
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


# ============================================================
# SAMPLE POINTS
# ============================================================

# Ordinary large-r points.
#
# These overlap with the region where Cell 26 showed the
# persistent ~1e-24 oscillation.

R_ORDINARY = [
    20,
    25,
    30,
    35,
    40,
    50,
    60,
    80,
    100,
    150,
    200,
    300,
    400,
    500,
    600,
    800,
    1000,
    1200,
    1600,
    2000,
    3000,
    4000,
    5000,
    6000,
    8000,
    10000,
]


# ------------------------------------------------------------
# Phase-related points.
#
# The Fourier-side expressions contain phases involving r L.
#
# We therefore sample points satisfying
#
#     r L = k*pi
#
# and
#
#     r L = (k + 1/2)*pi.
#
# These are useful for separating ordinary decay from
# oscillatory phase structure.
# ------------------------------------------------------------

K_PHASE = [
    20,
    40,
    80,
    120,
    160,
    200,
    300,
    400,
    600,
    800,
]


R_PHASE_INTEGER = [
    mp.mpf(k) * mp.pi / L
    for k in K_PHASE
]

R_PHASE_HALF = [
    (mp.mpf(k) + mp.mpf("0.5")) * mp.pi / L
    for k in K_PHASE
]


# ============================================================
# HELPERS
# ============================================================


def nstr(x, digits=WORKING_DPS):
    return mp.nstr(
        x,
        digits,
    )


def evaluate(r, v_star):
    J = K_fourier(
        v_star,
        r,
        L,
    )

    h = h_plus(r)

    I = h * J

    return J, h, I


# ============================================================
# HEADER
# ============================================================

print("=" * 110)
print("CELL 27 — ARCHIMEDEAN TAIL ANATOMY")
print("=" * 110)

print()
print("Purpose:")
print(
    "Inspect the large-r analytic Archimedean integrand "
    "for the fixed forensic ground state."
)

print()
print("Parameters:")
print(f" c = {c}")
print(f" N = {N}")
print(f" working_dps = {WORKING_DPS}")
print(f" L = {nstr(L, 60)}")

print()
print("No long-range integration is performed in this cell.")


# ============================================================
# 1. RETRIEVE FIXED FORENSIC GROUND STATE
# ============================================================

print()
print("-" * 110)
print("1. FIXED FORENSIC GROUND STATE")
print("-" * 110)

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
# 2. ORDINARY LARGE-r SAMPLE
# ============================================================

print()
print("-" * 110)
print("2. ORDINARY LARGE-r SAMPLE")
print("-" * 110)

print()
print(
    "J(r) = K_fourier(v,r,L)"
)

print(
    "I(r) = h_+(r) K_fourier(v,r,L)"
)

print()
print(
    "The scaled columns are diagnostic only; "
    "no asymptotic power is assumed."
)

print()
print(
    " r"
    "             J(r)"
    "             h_+(r)"
    "             I(r)"
    "           r^2 J"
    "           r^3 J"
    "           r^2 I"
    "           r^3 I"
)

print("-" * 110)


ordinary_results = []


for r in R_ORDINARY:

    start = time.perf_counter()

    J, h, I = evaluate(
        mp.mpf(r),
        v_star,
    )

    elapsed = (
        time.perf_counter()
        - start
    )

    r_mp = mp.mpf(r)

    row = (
        r_mp,
        J,
        h,
        I,
        r_mp**2 * J,
        r_mp**3 * J,
        r_mp**2 * I,
        r_mp**3 * I,
        elapsed,
    )

    ordinary_results.append(row)

    print(
        f"{r:6d} "
        f"{nstr(J, 28):>30} "
        f"{nstr(h, 20):>22} "
        f"{nstr(I, 28):>30} "
        f"{nstr(r_mp**2 * J, 20):>22} "
        f"{nstr(r_mp**3 * J, 20):>22} "
        f"{nstr(r_mp**2 * I, 20):>22} "
        f"{nstr(r_mp**3 * I, 20):>22}"
    )


# ============================================================
# 3. PHASE-LOCKED SAMPLE
# ============================================================

print()
print("-" * 110)
print("3. PHASE-LOCKED SAMPLE")
print("-" * 110)

print()
print(
    "First sample family: r L = k pi"
)

print()
print(
    " k"
    "          r"
    "             J(r)"
    "             I(r)"
    "             r^2 I"
)

print("-" * 110)


phase_integer_results = []


for k, r in zip(
    K_PHASE,
    R_PHASE_INTEGER,
):

    J, h, I = evaluate(
        r,
        v_star,
    )

    phase_integer_results.append(
        (
            k,
            r,
            J,
            I,
        )
    )

    print(
        f"{k:4d} "
        f"{nstr(r, 24):>26} "
        f"{nstr(J, 32):>36} "
        f"{nstr(I, 32):>36} "
        f"{nstr(r**2 * I, 24):>28}"
    )


print()
print(
    "Second sample family: r L = (k + 1/2) pi"
)

print()
print(
    " k"
    "          r"
    "             J(r)"
    "             I(r)"
    "             r^2 I"
)

print("-" * 110)


phase_half_results = []


for k, r in zip(
    K_PHASE,
    R_PHASE_HALF,
):

    J, h, I = evaluate(
        r,
        v_star,
    )

    phase_half_results.append(
        (
            k,
            r,
            J,
            I,
        )
    )

    print(
        f"{k:4d} "
        f"{nstr(r, 24):>26} "
        f"{nstr(J, 32):>36} "
        f"{nstr(I, 32):>36} "
        f"{nstr(r**2 * I, 24):>28}"
    )


# ============================================================
# 4. CONSECUTIVE LARGE-r RATIOS
# ============================================================

print()
print("-" * 110)
print("4. LARGE-r RATIOS")
print("-" * 110)

print()
print(
    "For consecutive ordinary sample points, inspect"
)

print()
print(
    " |J(r_i)| / |J(r_{i-1})|"
)

print(
    "and"

)

print(
    " |I(r_i)| / |I(r_{i-1})|"
)

print()
print(
    "These are only empirical diagnostics."
)

print()
print(
    " r"
    "          |J_i|/|J_prev|"
    "          |I_i|/|I_prev|"
)

print("-" * 110)


for i in range(1, len(ordinary_results)):

    r_prev, J_prev, _, I_prev, *_ = (
        ordinary_results[i - 1]
    )

    r, J, _, I, *_ = (
        ordinary_results[i]
    )

    if J_prev != 0:
        J_ratio = (
            abs(J)
            / abs(J_prev)
        )
    else:
        J_ratio = mp.nan

    if I_prev != 0:
        I_ratio = (
            abs(I)
            / abs(I_prev)
        )
    else:
        I_ratio = mp.nan

    print(
        f"{int(r):6d} "
        f"{nstr(J_ratio, 24):>28} "
        f"{nstr(I_ratio, 24):>28}"
    )


# ============================================================
# 5. DIRECT CHECK OF THE CELL-26 TAIL SCALE
# ============================================================

print()
print("-" * 110)
print("5. CELL-26 TAIL-SCALE CHECK")
print("-" * 110)

print()
print(
    "Cell 26 showed finite-T changes at roughly the 1e-24 "
    "level at large T."
)

print()
print(
    "This section reports the magnitude of the integrand "
    "itself at comparable r."
)

print()
print(
    " r"
    "             |J(r)|"
    "             |I(r)|"
)

print("-" * 110)


for (
    r,
    J,
    h,
    I,
    r2J,
    r3J,
    r2I,
    r3I,
    elapsed,
) in ordinary_results:

    if r >= 400:

        print(
            f"{int(r):6d} "
            f"{nstr(abs(J), 32):>36} "
            f"{nstr(abs(I), 32):>36}"
        )


# ============================================================
# 6. FINAL SUMMARY
# ============================================================

print()
print("=" * 110)
print("CELL 27 COMPLETE")
print("=" * 110)

print()
print(
    "No asymptotic extrapolation has been performed."
)

print(
    "No finite-T Archimedean integral has been extended."
)

print()
print(
    "The purpose of this cell is to expose the actual "
    "large-r structure before choosing the next experiment."
)
