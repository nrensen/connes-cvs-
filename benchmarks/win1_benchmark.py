#!/usr/bin/env python3
"""Historical-workload timer using the current checkout.

The frozen numbers below came from an April 2026 v0.1.0-era run. This script
now imports the current ``_compute_psi_pair`` and may include later changes.
Ratios printed against those constants are cross-run historical context, not a
same-environment A/B benchmark for the current checkout.

The historical baseline at c=13 N=50 dps=50 T=200 was recorded as:

    psi cache  : 122.112 s  (98.8%)
    TOTAL      : 123.607 s
    lambda_min : 4.8051148795216734485219293678933754790769597762326e-51

The script reports the current timing and printed-digit agreement with the
frozen lambda fingerprint.

Usage
-----
    python win1_benchmark.py [c] [N] [dps]

Defaults: c=13 N=50 dps=50 (same workload as the baseline).
"""
from __future__ import annotations

import os
import platform
import sys
import time

# Allow running directly with a bare python3 without PYTHONPATH set.
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

import mpmath as mp  # noqa: E402

from connes_cvs import (  # noqa: E402
    compute_ground_state,
    extract_zeros,
)
from connes_cvs.operator import (  # noqa: E402
    HAS_FLINT,
    _compute_psi_pair,
    prime_powers_up_to,
)

if HAS_FLINT:
    from flint import ctx as flint_ctx


BASELINE_LAMBDA_MIN_C13_N50_DPS50_T200 = (
    "4.8051148795216734485219293678933754790769597762326e-51"
)
BASELINE_T_CACHE_C13_N50_DPS50_T200 = 122.112  # seconds
BASELINE_T_TOTAL_C13_N50_DPS50_T200 = 123.607  # seconds


def bench(c: int, N: int, T: int, dps: int) -> dict:
    """Time each phase; return dict of timings + correctness fingerprint."""
    mp.mp.dps = dps
    if HAS_FLINT:
        flint_ctx.prec = int(dps * 3.5)

    # Phase 1: psi cache (with WIN 1 h_plus memoization active)
    t0 = time.perf_counter()
    L = mp.log(c)
    prime_data, _ = prime_powers_up_to(int(c))
    n_indices = list(range(-N, N + 1))
    psi_vals = {}
    psi_deriv_vals = {}
    for n_idx in n_indices:
        psi, psi_d = _compute_psi_pair(n_idx, L, T, dps, prime_data)
        psi_vals[n_idx] = psi
        psi_deriv_vals[n_idx] = psi_d
    t_cache = time.perf_counter() - t0

    # Phase 2: matrix assembly (same logic as build_galerkin_matrix)
    t0 = time.perf_counter()
    DIM = 2 * N + 1
    Q = mp.matrix(DIM, DIM)
    for i in range(DIM):
        m_idx = i - N
        for j in range(DIM):
            n_idx = j - N
            if m_idx == n_idx:
                Q[i, j] = psi_deriv_vals[n_idx]
            else:
                Q[i, j] = (psi_vals[m_idx] - psi_vals[n_idx]) / (m_idx - n_idx)
    for i in range(DIM):
        for j in range(i + 1, DIM):
            avg = (Q[i, j] + Q[j, i]) / 2
            Q[i, j] = avg
            Q[j, i] = avg
    t_assemble = time.perf_counter() - t0

    # Phase 3: eigensolver
    t0 = time.perf_counter()
    lam_min, eigvec = compute_ground_state(Q)
    t_eig = time.perf_counter() - t0

    # Phase 4: findroot
    t0 = time.perf_counter()
    zeros = extract_zeros(eigvec, c=c, n_zeros=1, dps=dps)
    t_findroot = time.perf_counter() - t0

    return {
        "c": c, "N": N, "T": T, "dps": dps, "DIM": DIM,
        "t_cache": t_cache, "t_assemble": t_assemble,
        "t_eig": t_eig, "t_findroot": t_findroot,
        "t_total": t_cache + t_assemble + t_eig + t_findroot,
        "lambda_min": lam_min,
        "psi_vals": psi_vals,
        "psi_deriv_vals": psi_deriv_vals,
        "zeros": zeros,
    }


def main() -> None:
    c = int(sys.argv[1]) if len(sys.argv) > 1 else 13
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    dps = int(sys.argv[3]) if len(sys.argv) > 3 else 50
    T = 200 if N <= 60 else 400

    print("=" * 72)
    print("Connes-van Suijlekom historical-workload timer (current checkout)")
    print("=" * 72)
    print(f"Python   : {platform.python_version()} ({sys.executable})")
    print(f"mpmath   : {mp.__version__}")
    print(f"HAS_FLINT: {HAS_FLINT}")
    print(f"Params   : c={c}, N={N}, T={T}, dps={dps}  (matrix DIM={2*N+1})")
    print("-" * 72)

    result = bench(c=c, N=N, T=T, dps=dps)

    def pct(value):
        return 100.0 * value / result["t_total"]
    print("Timings (seconds, wall, time.perf_counter()):")
    print(f"  1. psi cache (prime+pole+arch)  : {result['t_cache']:9.3f}  ({pct(result['t_cache']):5.1f}%)")
    print(f"  2. matrix assembly + symmetrize : {result['t_assemble']:9.3f}  ({pct(result['t_assemble']):5.1f}%)")
    print(f"  3. eigensolver (even sector)    : {result['t_eig']:9.3f}  ({pct(result['t_eig']):5.1f}%)")
    print(f"  4. findroot (first zero)        : {result['t_findroot']:9.3f}  ({pct(result['t_findroot']):5.1f}%)")
    print(f"  TOTAL                           : {result['t_total']:9.3f}  (100.0%)")
    print("-" * 72)

    # Baseline comparison (only meaningful at the exact baseline params)
    if c == 13 and N == 50 and dps == 50 and T == 200:
        baseline_cache = BASELINE_T_CACHE_C13_N50_DPS50_T200
        baseline_total = BASELINE_T_TOTAL_C13_N50_DPS50_T200
        speedup_cache = baseline_cache / result["t_cache"]
        speedup_total = baseline_total / result["t_total"]
        print("Ratio vs frozen April 2026 baseline (NOT a contemporary A/B):")
        print(f"  psi cache : {baseline_cache:8.3f}s -> {result['t_cache']:8.3f}s  ({speedup_cache:5.2f}x)")
        print(f"  total     : {baseline_total:8.3f}s -> {result['t_total']:8.3f}s  ({speedup_total:5.2f}x)")
        print("-" * 72)

        # Printed-digit comparison against the frozen baseline lambda.
        baseline_lm = mp.mpf(BASELINE_LAMBDA_MIN_C13_N50_DPS50_T200)
        new_lm = result["lambda_min"]
        abs_diff = abs(new_lm - baseline_lm)
        rel_diff = abs_diff / abs(baseline_lm) if baseline_lm != 0 else abs_diff
        print("Historical lambda comparison:")
        print(f"  baseline : {mp.nstr(baseline_lm, dps)}")
        print(f"  new      : {mp.nstr(new_lm, dps)}")
        print(f"  abs diff : {mp.nstr(abs_diff, 5)}")
        print(f"  rel diff : {mp.nstr(rel_diff, 5)}")
        tol = mp.mpf("1e-50")
        if rel_diff < tol:
            print("  PASS: relative difference < 1e-50; all 50 displayed digits agree")
        else:
            print("  FAIL: rel diff exceeds 1e-50")
        print("-" * 72)

    lam = result["lambda_min"]
    print("Correctness fingerprint:")
    print(f"  lambda_min = {mp.nstr(lam, dps)}")
    z0 = result["zeros"][0]
    gd, gt = z0["gamma_detected"], z0["gamma_true"]
    print(f"  gamma_1 true     = {mp.nstr(gt, 25)}")
    if gd is not None:
        print(f"  gamma_1 detected = {mp.nstr(gd, 25)}")
        print(f"  |gamma_1 error|  = {mp.nstr(z0['error'], 6)}")
    else:
        print(f"  gamma_1 detected = None (findroot failed; expected at N={N} dps={dps})")
    print("=" * 72)


if __name__ == "__main__":
    main()
