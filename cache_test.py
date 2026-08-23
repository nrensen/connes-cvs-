# ============================================================
# cache_test.py
#
# Test the generic JSON result cache and benchmark cache speed.
# ============================================================

import time

import mpmath as mp

from cell import (
    CELL_CACHE_DIR,
    cache_key,
    cache_load,
    cache_save,
    cache_get_or_compute,
    get_ground_state,
)


print("=" * 72)
print("GENERIC CACHE + GROUND STATE CACHE TEST")
print("=" * 72)


# ============================================================
# 1. GENERIC CACHE
# ============================================================

print()
print("-" * 72)
print("1. GENERIC JSON CACHE")
print("-" * 72)

parameters = {
    "purpose": "cache_test",
    "integer": 123,
    "nested": {
        "alpha": "1.234567890123456789",
        "beta": 42,
    },
    "list": [
        "a",
        "b",
        17,
    ],
}

results = {
    "answer": "cache machinery works",
    "value": "3.141592653589793238462643383279",
}

digest_1, record_1 = cache_key(
    "generic_test",
    parameters,
)

print()
print("digest:")
print(digest_1)

print()
print("canonical record:")
print(record_1)

save_start = time.perf_counter()

digest_2, path = cache_save(
    "generic_test",
    parameters,
    results,
)

save_time = (
    time.perf_counter()
    - save_start
)

print()
print(
    f"save time = {save_time:.6f} s"
)

print()
print(
    f"path = {path}"
)

load_start = time.perf_counter()

loaded_results, metadata = cache_load(
    "generic_test",
    parameters,
)

load_time = (
    time.perf_counter()
    - load_start
)

print()
print(
    f"load time = {load_time:.6f} s"
)

print()
print(
    f"loaded results = {loaded_results}"
)

assert digest_1 == digest_2
assert loaded_results == results


# ============================================================
# 2. PARAMETER IDENTITY
# ============================================================

print()
print("-" * 72)
print("2. PARAMETER IDENTITY")
print("-" * 72)

modified_parameters = dict(parameters)
modified_parameters["integer"] = 124

modified_digest, _ = cache_key(
    "generic_test",
    modified_parameters,
)

print()
print(
    "original digest:"
)
print(digest_1)

print()
print(
    "modified digest:"
)
print(modified_digest)

print()
print(
    "digest changed =",
    digest_1 != modified_digest,
)

assert digest_1 != modified_digest


# ============================================================
# 3. CACHE GET-OR-COMPUTE
# ============================================================

print()
print("-" * 72)
print("3. CACHE GET-OR-COMPUTE")
print("-" * 72)

counter = {
    "calls": 0,
}


def deliberately_slow_calculation():
    counter["calls"] += 1

    time.sleep(0.25)

    return {
        "computed_value": 987654321,
        "calls_at_compute": counter["calls"],
    }


slow_parameters = {
    "test": "deliberately_slow",
    "version": 1,
}

print()
print("First call:")

result_1, meta_1 = cache_get_or_compute(
    "generic_test",
    slow_parameters,
    deliberately_slow_calculation,
)

print()
print("Second call:")

result_2, meta_2 = cache_get_or_compute(
    "generic_test",
    slow_parameters,
    deliberately_slow_calculation,
)

print()
print(
    "compute function calls =",
    counter["calls"],
)

print()
print(
    "first cache_hit =",
    meta_1["cache_hit"],
)

print()
print(
    "second cache_hit =",
    meta_2["cache_hit"],
)

print()
print(
    f"first total  = "
    f"{meta_1['total_seconds']:.6f} s"
)

print()
print(
    f"second total = "
    f"{meta_2['total_seconds']:.6f} s"
)

assert counter["calls"] == 1
assert meta_1["cache_hit"] is False
assert meta_2["cache_hit"] is True
assert result_1 == result_2


# ============================================================
# 4. GROUND STATE — FAST CACHE
# ============================================================

print()
print("-" * 72)
print("4. GROUND STATE — FAST CACHE")
print("-" * 72)

c = 13
N = 8
T = 40
dps = 80

print()
print(
    f"parameters: c={c}, N={N}, T={T}, dps={dps}"
)

print()
print("Requesting ground state with fast validation:")

gs1_start = time.perf_counter()

lambda_1, v_1, gs_meta_1 = (
    get_ground_state(
        c=c,
        N=N,
        T=T,
        dps=dps,
        cache=True,
        validation="fast",
        verbose=True,
    )
)

gs1_elapsed = (
    time.perf_counter()
    - gs1_start
)

print()
print(
    f"wall time = {gs1_elapsed:.6f} s"
)

print()
print(
    "lambda_min ="
)
print(
    mp.nstr(lambda_1, 60)
)

print()
print(
    "||v|| ="
)
print(
    mp.nstr(
        mp.sqrt(mp.fdot(v_1, v_1)),
        60,
    )
)


# ============================================================
# 5. GROUND STATE — SECOND FAST REQUEST
# ============================================================

print()
print("-" * 72)
print("5. GROUND STATE — SECOND FAST REQUEST")
print("-" * 72)

gs2_start = time.perf_counter()

lambda_2, v_2, gs_meta_2 = (
    get_ground_state(
        c=c,
        N=N,
        T=T,
        dps=dps,
        cache=True,
        validation="fast",
        verbose=True,
    )
)

gs2_elapsed = (
    time.perf_counter()
    - gs2_start
)

print()
print(
    f"wall time = {gs2_elapsed:.6f} s"
)

lambda_difference = abs(
    lambda_1 - lambda_2
)

vector_difference = mp.sqrt(
    mp.fdot(
        v_1 - v_2,
        v_1 - v_2,
    )
)

print()
print(
    "|lambda_1 - lambda_2| ="
)
print(
    mp.nstr(
        lambda_difference,
        30,
    )
)

print()
print(
    "||v_1 - v_2|| ="
)
print(
    mp.nstr(
        vector_difference,
        30,
    )
)


# ============================================================
# 6. GROUND STATE — FULL VALIDATION
# ============================================================

print()
print("-" * 72)
print("6. GROUND STATE — FULL VALIDATION")
print("-" * 72)

print()
print(
    "This deliberately rebuilds Q."
)
print(
    "It is expected to be much slower than the fast cache hit."
)

full_start = time.perf_counter()

lambda_3, v_3, gs_meta_3 = (
    get_ground_state(
        c=c,
        N=N,
        T=T,
        dps=dps,
        cache=True,
        validation="full",
        verbose=True,
    )
)

full_elapsed = (
    time.perf_counter()
    - full_start
)

print()
print(
    f"wall time = {full_elapsed:.6f} s"
)

print()
print(
    "full-validation residual ="
)
print(
    mp.nstr(
        gs_meta_3["full_validation"][
            "residual_norm"
        ],
        30,
    )
)


# ============================================================
# 7. SPEEDUP
# ============================================================

print()
print("-" * 72)
print("7. SPEEDUP")
print("-" * 72)

if gs2_elapsed > 0:
    print()
    print(
        "cold/fast-request ratio =",
        f"{gs1_elapsed / gs2_elapsed:.3f}x",
    )

if gs2_elapsed > 0:
    print()
    print(
        "full-validation/fast-cache ratio =",
        f"{full_elapsed / gs2_elapsed:.3f}x",
    )


# ============================================================
# 8. SUMMARY
# ============================================================

print()
print("-" * 72)
print("8. SUMMARY")
print("-" * 72)

print()
print(
    "generic cache: PASS"
)

print(
    "parameter hashing: PASS"
)

print(
    "get-or-compute: PASS"
)

print(
    "ground-state first request cache_hit =",
    gs_meta_1["cache_hit"],
)

print(
    "ground-state second request cache_hit =",
    gs_meta_2["cache_hit"],
)

print(
    "ground-state full request cache_hit =",
    gs_meta_3["cache_hit"],
)

print()
print("=" * 72)
print("CACHE TEST COMPLETE")
print("=" * 72)
