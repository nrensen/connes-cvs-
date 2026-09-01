# ============================================================
# CELL 23 — OPTIMISED ANALYTIC ARCHIMEDEAN REDUCTION
#
# Purpose
# -------
# Cell 22 demonstrated that the y-integral in Cell 21 can be
# performed analytically, reducing the nested quadrature to a
# single r-integral.
#
# Cell 23 makes that same calculation more efficient and
# numerically better conditioned:
#
#   1. exploit (m,n) <-> (n,m) symmetry;
#   2. evaluate each Fourier-mode integral only once;
#   3. use sinc-style expressions for removable singularities.
#
# No mathematical change is intended relative to Cell 22.
#
# This is the candidate efficient implementation that will later
# be compared against:
#
#   Cell 22 @ 40 dps
#   Cell 21 @ 40 dps
#
# The latter remains the independent nested numerical control.
# ============================================================

import time

import mpmath as mp

from cell import (
    FORENSIC_GROUND_STATE,
    get_ground_state,
    compute_L,
    archimedean_integral,
    canonical_to_full,
    full_to_canonical,
)


# ============================================================
# PARAMETERS
# ============================================================

WORKING_DPS = 120
T = 60

mp.mp.dps = WORKING_DPS

c = FORENSIC_GROUND_STATE["c"]
N = FORENSIC_GROUND_STATE["N"]

L = compute_L(c)

DISPLAY_DIGITS = 120


def nstr(x):
    return mp.nstr(x, DISPLAY_DIGITS)


def elapsed(start):
    return time.perf_counter() - start


# ============================================================
# HEADER
# ============================================================

print("=" * 78)
print("CELL 23 — OPTIMISED ANALYTIC ARCHIMEDEAN REDUCTION")
print("=" * 78)

print()
print("Parameters:")
print(f"  c              = {c}")
print(f"  N              = {N}")
print(f"  T              = {T}")
print(f"  working_dps    = {WORKING_DPS}")

print()
print("Forensic ground state:")
print(f"  c              = {FORENSIC_GROUND_STATE['c']}")
print(f"  N              = {FORENSIC_GROUND_STATE['N']}")
print(f"  T              = {FORENSIC_GROUND_STATE['T']}")
print(f"  generation_dps = {FORENSIC_GROUND_STATE['dps']}")

print()
print("L =")
print(nstr(L))


# ============================================================
# 1. FORENSIC GROUND STATE
# ============================================================

print()
print("-" * 78)
print("1. FORENSIC GROUND STATE")
print("-" * 78)

ground_start = time.perf_counter()

lambda_forensic, v_star, ground_meta = get_ground_state(
    **FORENSIC_GROUND_STATE,
    verbose=True,
)

ground_elapsed = elapsed(ground_start)

u_star = canonical_to_full(v_star)
v = full_to_canonical(u_star)
u = canonical_to_full(v)

print()
print("lambda_forensic =")
print(nstr(lambda_forensic))

print()
print("||u_star|| =")
print(
    nstr(
        mp.sqrt(
            mp.fdot(
                u_star,
                u_star,
            )
        )
    )
)

print()
print(
    f"ground-state retrieval elapsed = "
    f"{ground_elapsed:.6f} s"
)


# ============================================================
# 2. CANONICAL / FULL REPRESENTATIONS
# ============================================================

print()
print("-" * 78)
print("2. CANONICAL / FULL REPRESENTATIONS")
print("-" * 78)

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
print("||u - u_star|| =")
print(
    nstr(
        mp.sqrt(
            mp.fdot(
                u - u_star,
                u - u_star,
            )
        )
    )
)


# ============================================================
# 3. EXPLICIT ARCHIMEDEAN FUNCTIONAL
#
# Only the outer r-integral remains numerical:
#
#   A_arch
#     = 1/pi int_0^T h_+(r) J_v(r) dr.
#
# ============================================================

print()
print("-" * 78)
print("3. OPTIMISED ANALYTIC ARCHIMEDEAN")
print("-" * 78)

print()
print(
    "Computing:"
)
print(
    "  J_v(r) from the symmetric finite Fourier sum"
)
print(
    "  A_arch = (1/pi) int_0^T h_+(r) J_v(r) dr"
)
print()
print(
    "No numerical y-integration is performed."
)

arch_start = time.perf_counter()

explicit_arch = archimedean_integral(T, v_star, L)

arch_elapsed = elapsed(
    arch_start
)

print()
print("Analytic Archimedean =")
print(nstr(explicit_arch))

print()
print(
    f"analytic Archimedean elapsed = "
    f"{arch_elapsed:.6f} s"
)


# ============================================================
# 4. FINAL RESULT
# ============================================================

print()
print("=" * 78)
print("CELL 23 COMPLETE")
print("=" * 78)

print()
print("Result:")
print(
    "  analytic Archimedean =",
    nstr(explicit_arch),
)

print()
print("Timing:")
print(
    f"  ground state         = "
    f"{ground_elapsed:.6f} s"
)
print(
    f"  analytic Archimedean = "
    f"{arch_elapsed:.6f} s"
)

print()
print("=" * 78)
print("CELL 23 COMPLETE")
print("=" * 78)
