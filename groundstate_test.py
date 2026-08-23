# ============================================================
# groundstate_test.py
#
# Test the persistent ground-state cache machinery in cell.py.
# ============================================================

import os
import sys
import shutil

import mpmath as mp

from cell import (
    GROUND_STATE_CACHE_DIR,
    canonical_to_full,
    full_to_canonical,
    get_ground_state,
    ground_state_cache_key,
)


# ============================================================
# PARAMETERS
# ============================================================

c = 13
N = 8
T = 40
dps = 80

mp.mp.dps = dps


print("=" * 70)
print("GROUND STATE CACHE TEST")
print("=" * 70)

print()
print("Parameters:")
print(f"c   = {c}")
print(f"N   = {N}")
print(f"T   = {T}")
print(f"dps = {dps}")

digest, record = ground_state_cache_key(
    c=c,
    N=N,
    T=T,
    dps=dps,
)

print()
print("Cache identity:")
print(f"key = {digest}")
print(f"record = {record}")


# ============================================================
# 1. FIRST REQUEST
# ============================================================

print()
print("-" * 70)
print("1. FIRST REQUEST")
print("-" * 70)

lambda_1, v_1, meta_1 = get_ground_state(
    c=c,
    N=N,
    T=T,
    dps=dps,
    cache=True,
    validate=True,
    verbose=True,
)

print()
print("lambda_1 =")
print(mp.nstr(lambda_1, 60))

print()
print("||v_1|| =")
print(mp.nstr(mp.sqrt(mp.fdot(v_1, v_1)), 60))

print()
print("metadata:")
print(meta_1)


# ============================================================
# 2. SECOND REQUEST — MUST BE CACHE HIT
# ============================================================

print()
print("-" * 70)
print("2. SECOND REQUEST")
print("-" * 70)

lambda_2, v_2, meta_2 = get_ground_state(
    c=c,
    N=N,
    T=T,
    dps=dps,
    cache=True,
    validate=True,
    verbose=True,
)

print()
print("lambda_2 =")
print(mp.nstr(lambda_2, 60))

print()
print("||v_2|| =")
print(mp.nstr(mp.sqrt(mp.fdot(v_2, v_2)), 60))


# ============================================================
# 3. REPEATABILITY
# ============================================================

print()
print("-" * 70)
print("3. REPEATABILITY")
print("-" * 70)

lambda_difference = abs(lambda_1 - lambda_2)

vector_difference = mp.sqrt(
    mp.fdot(
        v_1 - v_2,
        v_1 - v_2,
    )
)

print()
print("|lambda_1 - lambda_2| =")
print(mp.nstr(lambda_difference, 60))

print()
print("||v_1 - v_2|| =")
print(mp.nstr(vector_difference, 60))

print()
print("cache_hit on second request =")
print(meta_2["cache_hit"])


# ============================================================
# 4. CANONICAL <-> FULL ROUND TRIP
# ============================================================

print()
print("-" * 70)
print("4. CANONICAL / FULL ROUND TRIP")
print("-" * 70)

v_can = full_to_canonical(
    v_2,
    N,
)

v_roundtrip = canonical_to_full(
    v_can,
    N,
)

roundtrip_error = mp.sqrt(
    mp.fdot(
        v_roundtrip - v_2,
        v_roundtrip - v_2,
    )
)

print()
print("||v_canonical|| =")
print(mp.nstr(mp.sqrt(mp.fdot(v_can, v_can)), 60))

print()
print("round-trip error =")
print(mp.nstr(roundtrip_error, 60))


# ============================================================
# 5. CACHE FILE
# ============================================================

print()
print("-" * 70)
print("5. CACHE FILE")
print("-" * 70)

cache_path = meta_2["cache_path"]

print()
print(f"path = {cache_path}")

print()
print(f"exists = {cache_path.exists()}")

if cache_path.exists():
    print(f"size = {cache_path.stat().st_size} bytes")


# ============================================================
# 6. PARAMETER HASH SENSITIVITY
# ============================================================

print()
print("-" * 70)
print("6. PARAMETER HASH SENSITIVITY")
print("-" * 70)

tests = [
    ("N+1", dict(c=c, N=N + 1, T=T, dps=dps)),
    ("T+1", dict(c=c, N=N, T=T + 1, dps=dps)),
    ("dps+1", dict(c=c, N=N, T=T, dps=dps + 1)),
    ("c+1", dict(c=c + 1, N=N, T=T, dps=dps)),
]

for label, params in tests:
    changed_digest, _ = ground_state_cache_key(**params)

    print()
    print(f"{label}:")
    print(f"  key = {changed_digest}")
    print(f"  differs = {changed_digest != digest}")


# ============================================================
# 7. SUMMARY
# ============================================================

print()
print("-" * 70)
print("7. SUMMARY")
print("-" * 70)

print()
print("First request cache hit:")
print(meta_1["cache_hit"])

print()
print("Second request cache hit:")
print(meta_2["cache_hit"])

print()
print("Eigenvalue reproducibility:")
print(mp.nstr(lambda_difference, 10))

print()
print("Eigenvector reproducibility:")
print(mp.nstr(vector_difference, 10))

print()
print("Round-trip error:")
print(mp.nstr(roundtrip_error, 10))

print()
print("GROUND STATE CACHE TEST COMPLETE")
print("=" * 70)
