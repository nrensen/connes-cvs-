# `connes_cvs/` - Python package source

**Package:** [`connes-cvs` on PyPI](https://pypi.org/project/connes-cvs/) - an open-source arbitrary-precision implementation of the Connes-van Suijlekom Galerkin matrix.

**Package version:** `0.3.1`. The published papers pin the preserved historical release `0.2.2` (2026-04-19).

## Source files

| File | Content |
|---|---|
| `__init__.py` | Public top-level API: `build_galerkin_matrix`, `compute_ground_state`, `extract_zeros`, `arb_eigenpair_residual_bound`, and `__version__`. |
| `operator.py` | Core arithmetic: prime-power data, the prime/pole/archimedean pieces, Galerkin assembly, even-sector eigensolve and seeded finite-test-function root location. `extract_zeros(..., c=...)` avoids the float64-`L` pitfall. |
| `kernels.py` | Stable low-level Fourier-kernel helpers. |
| `runner.py` | Explicit-precision production cell runner (`CellConfig`, `GalerkinCell`, `run_cell`) with spawn-safe parallel psi-cache computation, progress, lossless mpmath transport, integrity-checked atomic checkpoints and JSON artifacts. Checkpoint hashes detect corruption, not malicious modification: resume only trusted local files. Locks are single-writer guards, not authentication. This module is imported explicitly rather than re-exported at package top level. |
| `sweep.py` | Validated multi-cutoff orchestration (`run_sweep`) over the runner. Each cutoff is sequential; each cell may use a process pool. |
| `validation.py` | `arb_eigenpair_residual_bound`, a rigorously scoped residual bound for an exact supplied finite real-symmetric mpmath matrix. It is not an infinite-operator or truncation-error certificate. |
| `py.typed` | PEP 561 marker. |

The package intentionally has no automatic precision recommender or a-priori “sufficient precision” certificate. Choose `dps` and optional `flint_bits` explicitly, then validate through a committed reference or a cross-precision recomputation.

## Version history

| Version | Date | Distribution status | Headline |
|---|---|---|---|
| 0.1.0 | 2026-04-13 | Repository tag only; **never uploaded to PyPI** | Initial public repository release and Paper Zenodo Version 1 lineage. |
| 0.2.0 | 2026-04-14 | Live on PyPI | Memoized/fused archimedean kernel; historical A/B record reports 2.06x on the psi-cache phase and agreement in all 80 printed lambda digits. |
| 0.2.1 | 2026-04-19 | **Yanked on PyPI** | Internal `__version__` drift; superseded the same day. |
| 0.2.2 | 2026-04-19 | Preserved on PyPI | Version-string fix. This is the package version pinned by the published papers. |
| 0.3.0 | 2026-08-12 | Live on PyPI | Full-precision `extract_zeros(c=...)`, float-`L` warning, validated runner/sweep, finite-matrix Arb residual bound, hardening and real regression gates. |
| 0.3.1 | 2026-08-14 | Current package version | Documentation and metadata only; no code or numerical change. Corrects the `build_galerkin_matrix` basis docstring, replaces paper version DOIs with concept DOIs so the frozen PyPI description cannot go stale, and ships the current `ERRATA.md`. |

## Install

```bash
pip install connes-cvs
pip install connes-cvs==0.3.1          # this package version
pip install connes-cvs==0.2.0          # retained historical PyPI release
pip install connes-cvs==0.2.2          # version pinned by the papers
```

To test a source checkout, clone the repository and install it in an isolated environment:

```bash
pip install -e '.[dev]'
```

## Release-owner commands

Build, upload, tag and push operations are owner-run only after review. Preparing a source diff does not authorize any of those public actions.

## Cross-references

- [Top-level README](https://github.com/akivag613/connes-cvs-/blob/main/README.md)
- [Regression and identity tests](https://github.com/akivag613/connes-cvs-/tree/main/tests)
- [Committed reference values](https://github.com/akivag613/connes-cvs-/blob/main/tests/reference_values.json)
- [Runnable smoke and extended example](https://github.com/akivag613/connes-cvs-/blob/main/examples/basic_compute.py)
- [Historical benchmarks](https://github.com/akivag613/connes-cvs-/tree/main/benchmarks)
- [Riemann-zeros paper reference data](https://github.com/akivag613/connes-cvs-/blob/main/data/results_15pt_T800.json)
- [Changelog](https://github.com/akivag613/connes-cvs-/blob/main/CHANGELOG.md)

## Numerical contract

- The slow suite recomputes the same-precision A/B cell (`c=13`, `N=80`, `T=400`, `dps=80`) and requires every one of the 80 stored printed digits. It separately recomputes the extended validation cell (`c=13`, `N=100`, `T=400`, `dps=80`) and requires at least 20 lambda digits plus 18 gamma-error digits against its higher-precision reference.
- The slow suite also recomputes the exact paper-reference runner cell (`c=13`, `N=100`, `T=800`, `dps=150`). The run performs a fatal check against the 39-digit lambda reference (at least 22 significant digits) and checks the extracted first-zero error against the published row.
- A small `c=13`, `N=8`, `T=60`, `dps=30` test compares the classic and runner paths on raw mpmath tuples for every matrix entry, the selected eigenvalue and every eigenvector component. A second small test compares pooled and serial runner output exactly.
- These scopes do not establish blanket raw-value equality for every cutoff, dependency version or platform. Published `c=100` artifacts used `flint_bits=4*dps`; state that explicitly when reproducing them.
- The numerical core mutates process-global mpmath/Flint precision contexts and caches. Do not run different-precision cells concurrently in Python threads; use `runner`/`sweep`, which isolate arithmetic in processes.
- Never change version metadata independently. `pyproject.toml`, `connes_cvs.__version__`, `CHANGELOG.md`, and `CITATION.cff` must agree before an owner-run release.
- Do not yank or delete a published release without a separate, evidence-backed release decision. The yanked 0.2.1 remains historical metadata; 0.2.0 and 0.2.2 remain live.
