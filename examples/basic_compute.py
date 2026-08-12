#!/usr/bin/env python3
"""Run a seconds-level smoke cell, or an extended c=13 validation cell.

Default smoke cell (measured in the release environment):

    c=13, N=8, T=60, dps=30
    lambda_min approximately 4.43043e-23
    |gamma_1 error| approximately 2.52738e-17

It completed in about 1.8 s with python-flint 0.8.0 and 4.8 s through the
mpmath fallback on the release machine. Runtime depends on versions and
hardware. This small cell checks the public API; it is not a paper benchmark.

``--extended`` runs the historical validation cell
``c=13, N=100, T=400, dps=80`` through the process-based runner. A v0.2.2
single-process measurement on an Apple M2 Max with python-flint 0.8.0 took
316 s for matrix construction and 14 s for the eigensolve. The extended mode
uses up to the runner's default eight-worker cap, prints progress, and locates
three finite-test roots after the cell completes.

Usage:
    python examples/basic_compute.py
    python examples/basic_compute.py --extended
"""

from __future__ import annotations

import argparse

import mpmath as mp

from connes_cvs import build_galerkin_matrix, compute_ground_state, extract_zeros


def _print_zeros(zeros: list[dict]) -> None:
    """Print located finite-test roots without float-only format codes."""
    for result in zeros:
        detected = result["gamma_detected"]
        detected_text = "N/A" if detected is None else mp.nstr(detected, 16)
        error_text = "N/A" if result["error"] is None else mp.nstr(result["error"], 7)
        print(
            f"  k={result['k']}: reference={mp.nstr(result['gamma_true'], 16)}, "
            f"detected={detected_text}, error={error_text}, "
            f"converged={result['converged']}"
        )


def run_smoke() -> None:
    """Run the small public-API smoke cell."""
    dps = 30
    mp.mp.dps = dps
    print("Smoke cell: c=13, N=8, T=60, dps=30")
    Q = build_galerkin_matrix(c=13, N=8, T=60, dps=dps)
    lam_min, eigvec = compute_ground_state(Q)
    zeros = extract_zeros(eigvec, c=13, n_zeros=1, dps=dps)
    print(f"  matrix size: {Q.rows} x {Q.cols}")
    print(f"  lambda_min: {mp.nstr(lam_min, 12)}")
    _print_zeros(zeros)
    print("This is a smoke test, not the published c=13 reference cell.")


def run_extended() -> None:
    """Run the longer validation cell through the guarded production runner."""
    from connes_cvs.runner import CellConfig, GalerkinCell

    dps = 80
    print("Extended cell: c=13, N=100, T=400, dps=80")
    cell = GalerkinCell(CellConfig(c=13, N=100, T=400, dps=dps))
    artifact = cell.run()
    if cell.eigvec_full is None:
        raise RuntimeError("runner completed without an eigenvector")
    zeros = extract_zeros(cell.eigvec_full, c=13, n_zeros=3, dps=dps)
    print(f"  lambda_min: {artifact['lambda_even']}")
    _print_zeros(zeros)


def main() -> None:
    """Parse mode and run under the multiprocessing-safe module guard."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--extended",
        action="store_true",
        help="run the multi-minute N=100 validation cell with progress",
    )
    args = parser.parse_args()
    if args.extended:
        run_extended()
    else:
        run_smoke()


if __name__ == "__main__":
    main()
