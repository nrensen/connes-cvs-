#!/usr/bin/env python3
"""Time one cell through the current public multiprocessing sweep path.

This is not directly comparable with a historical single-process v0.1.0 run.
A performance claim requires a matched before/after run with the same worker
count, versions, backend, hardware, workload and timing boundary.

Usage:
    python win1_pool_benchmark.py [c] [N] [T] [dps] [workers]

Defaults: c=13 N=80 T=400 dps=80 workers=8.
"""
from __future__ import annotations
import os
import platform
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

import mpmath as mp  # noqa: E402
from connes_cvs.sweep import run_sweep  # noqa: E402
from connes_cvs.operator import HAS_FLINT  # noqa: E402


def main() -> None:
    c   = int(sys.argv[1]) if len(sys.argv) > 1 else 13
    N   = int(sys.argv[2]) if len(sys.argv) > 2 else 80
    T   = int(sys.argv[3]) if len(sys.argv) > 3 else 400
    dps = int(sys.argv[4]) if len(sys.argv) > 4 else 80
    workers = int(sys.argv[5]) if len(sys.argv) > 5 else 8

    print("=" * 72)
    print("Connes-van Suijlekom POOL benchmark (production-style, multiprocessing)")
    print("=" * 72)
    print(f"Python   : {platform.python_version()}")
    print(f"mpmath   : {mp.__version__}")
    print(f"HAS_FLINT: {HAS_FLINT}")
    print(f"Params   : c={c} N={N} T={T} dps={dps}  workers={workers}  DIM={2*N+1}")
    print("-" * 72)
    print(f"START at {time.strftime('%H:%M:%S')}")
    sys.stdout.flush()

    t_start = time.perf_counter()
    result = run_sweep(
        cutoffs=[c], N=N, T=T, dps=dps, workers=workers, n_zeros=1
    )[c]
    t_end = time.perf_counter()

    print(f"END   at {time.strftime('%H:%M:%S')}")
    print("-" * 72)

    timing = result.get("timing", {})
    print("Phase timings (seconds, from run_sweep):")
    for k in (
        "psi_cache_s",
        "matrix_assembly_s",
        "diagonalize_s",
        "zeros_s",
        "total_with_zeros_s",
    ):
        v = timing.get(k)
        if v is not None:
            print(f"  {k:<14}: {v:9.3f}")
    print(f"  WALL TOTAL    : {(t_end - t_start):9.3f}")
    print("-" * 72)
    lam = result.get("lambda_min") or result.get("lambda_min_str")
    print(f"lambda_min : {lam}")
    print("=" * 72)


if __name__ == "__main__":
    main()
