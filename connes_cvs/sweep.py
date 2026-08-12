"""Validated multi-cutoff orchestration for CvS Galerkin cells.

The sweep API delegates every cutoff to :mod:`connes_cvs.runner`, so serial
and multiprocessing runs share one arithmetic path, spawn-safe process setup,
lossless mpf transport, per-point heartbeat, and optional atomic checkpoints.

Multiprocessing entry points must be protected on spawn platforms::

    from connes_cvs.sweep import run_sweep

    if __name__ == "__main__":
        results = run_sweep([7, 8, 9], N=60, T=400, dps=80, workers=4)
        for cutoff, row in sorted(results.items()):
            print(cutoff, row["lambda_min"])

The numerical core mutates the process-global mpmath and Flint precision
contexts. Run independent cells in processes, not concurrent Python threads.
"""

from __future__ import annotations

import os
import time
from multiprocessing import cpu_count
from numbers import Integral
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

import mpmath as mp

from connes_cvs.operator import extract_zeros
from connes_cvs.runner import CellConfig, GalerkinCell


def _require_int(name: str, value: int, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{name} must be an integer, got {type(value).__name__}")
    value_int = int(value)
    if value_int < minimum:
        raise ValueError(f"{name} must be >= {minimum}, got {value_int}")
    return value_int


def _validated_cutoffs(cutoffs: Sequence[int]) -> List[int]:
    if isinstance(cutoffs, (str, bytes)) or not isinstance(cutoffs, Sequence):
        raise TypeError("cutoffs must be a nonempty sequence of integers")
    values = [_require_int(f"cutoffs[{i}]", value, 2) for i, value in enumerate(cutoffs)]
    if not values:
        raise ValueError("cutoffs must not be empty")
    if len(set(values)) != len(values):
        raise ValueError("cutoffs must not contain duplicates")
    return values


def _optional_directory(
    name: str, value: Optional[Union[str, os.PathLike]]
) -> Optional[Path]:
    if value is None:
        return None
    try:
        path = Path(value).expanduser().resolve()
    except TypeError as exc:
        raise TypeError(f"{name} must be path-like or None") from exc
    path.mkdir(parents=True, exist_ok=True)
    if not path.is_dir():
        raise ValueError(f"{name} is not a directory: {path}")
    return path


def run_sweep(
    cutoffs: Sequence[int],
    N: int = 100,
    T: int = 800,
    dps: int = 150,
    workers: Optional[int] = None,
    flint_bits: Optional[int] = None,
    n_zeros: int = 10,
    zero_tol: Optional[Union[mp.mpf, str]] = None,
    strict_zeros: bool = False,
    checkpoint_dir: Optional[Union[str, os.PathLike]] = None,
    artifact_dir: Optional[Union[str, os.PathLike]] = None,
    overwrite_artifacts: bool = False,
) -> Dict[int, Dict[str, Any]]:
    """Run validated, resumable Galerkin cells for several integer cutoffs.

    Cutoffs are processed sequentially; each cell may use up to ``workers``
    spawn processes for its independent basis-index quadratures. The returned
    mapping retains the v0.2 public fields while adding each runner artifact.

    ``checkpoint_dir`` stores one integrity-checked psi cache per cutoff.
    ``artifact_dir`` stores one atomic JSON runner artifact per cutoff. Those
    files are single-writer local inputs/outputs; a live lock prevents two
    processes from sharing the same target accidentally.
    """
    values = _validated_cutoffs(cutoffs)
    N = _require_int("N", N, 1)
    T = _require_int("T", T, 1)
    dps = _require_int("dps", dps, 15)
    if flint_bits is not None:
        flint_bits = _require_int(
            "flint_bits", flint_bits, (3322 * dps + 999) // 1000
        )
    n_zeros = _require_int("n_zeros", n_zeros, 1)
    if workers is None:
        workers = min(cpu_count(), 8)
    workers = _require_int("workers", workers, 1)
    if not isinstance(strict_zeros, bool):
        raise TypeError("strict_zeros must be a bool")
    if not isinstance(overwrite_artifacts, bool):
        raise TypeError("overwrite_artifacts must be a bool")
    checkpoint_root = _optional_directory("checkpoint_dir", checkpoint_dir)
    artifact_root = _optional_directory("artifact_dir", artifact_dir)

    results: Dict[int, Dict[str, Any]] = {}
    for cutoff in values:
        checkpoint_path = (
            checkpoint_root / f"c{cutoff}_N{N}_T{T}_dps{dps}.checkpoint.json"
            if checkpoint_root is not None
            else None
        )
        artifact_path = (
            artifact_root / f"c{cutoff}_N{N}_T{T}_dps{dps}.json"
            if artifact_root is not None
            else None
        )
        cell = GalerkinCell(
            CellConfig(
                c=cutoff,
                N=N,
                T=T,
                dps=dps,
                tag="run_sweep",
                flint_bits=flint_bits,
            ),
            workers=workers,
            ground_state="minimum",
            checkpoint_path=checkpoint_path,
            artifact_path=artifact_path,
            overwrite_artifact=overwrite_artifacts,
        )
        artifact = cell.run()
        if cell.eigvec_full is None or cell.lambda_even is None:
            raise RuntimeError("runner completed without an eigenpair")
        zero_started = time.perf_counter()
        zeros = extract_zeros(
            cell.eigvec_full,
            c=cutoff,
            n_zeros=n_zeros,
            dps=dps,
            tol=zero_tol,
            strict=strict_zeros,
        )
        zero_seconds = time.perf_counter() - zero_started
        timing = dict(artifact["timings_seconds"])
        timing["zeros_s"] = zero_seconds
        timing["total_with_zeros_s"] = timing["total_s"] + zero_seconds
        # Preserve the v0.2 public timing schema for downstream consumers while
        # retaining the more explicit v0.3 names from the runner artifact.
        timing["cache_sec"] = timing["psi_cache_s"]
        timing["matrix_sec"] = timing["matrix_assembly_s"]
        timing["diag_sec"] = timing["diagonalize_s"]
        timing["zeros_sec"] = timing["zeros_s"]
        timing["total_sec"] = timing["total_with_zeros_s"]
        results[cutoff] = {
            "cutoff": cutoff,
            "L": mp.log(cutoff),
            "N": N,
            "T": T,
            "dps": dps,
            "lambda_min": cell.lambda_even,
            "eigvec": cell.eigvec_full,
            "zeros": zeros,
            "gamma1_error": (
                zeros[0]["error"]
                if zeros and zeros[0].get("converged") is True
                else None
            ),
            "wall_time": timing["total_with_zeros_s"],
            "timing_scope": (
                "runner compute pipeline plus zero extraction; excludes runner "
                "artifact hashing, JSON serialization, and disk-write overhead"
            ),
            "timing": timing,
            "runner_artifact": artifact,
        }

    return results
