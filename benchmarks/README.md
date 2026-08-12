# `benchmarks/` - Historical records and current-path timers

This directory separates two different kinds of evidence:

1. the frozen April 2026 v0.1.0-to-v0.2.0 A/B summary; and
2. scripts that execute the **current checkout** for profiling or comparison with frozen historical numbers.

Do not treat a current run divided by a historical time as a same-environment A/B result.

## Distributed record

| File | Scope |
|---|---|
| `AB_VERIFIED_2026-04-14.md` | Historical same-environment A/B summary at c=13, N=80, T=400, dps=80, 12 workers. It records 2.06x on the psi-cache phase, 1.83x end-to-end, and agreement in all 80 printed lambda digits. Raw console logs are not distributed, so the summary is the surviving evidence. |
| `baseline_benchmark.py` | Single-process phase timer for the current checkout using the historical full-index assembly shape. Despite its filename, it cannot reconstruct the pre-optimization v0.1.0 source. |
| `win1_benchmark.py` | Runs the current checkout at the historical small workload and prints a clearly labelled comparison with frozen April 2026 numbers. It is archaeological context, not a clean contemporary A/B. |
| `win1_pool_benchmark.py` | Times one current public `run_sweep` cell with a configurable process count. It does not claim comparability with a v0.1.0 single-process run. |

Example current-path timing:

```bash
python benchmarks/win1_pool_benchmark.py 13 80 400 80 8
```

## Interpreting the historical terminology

The April summary described two decimal strings as identical when they agree in all 80 printed digits. This repository reserves raw identity language for direct `_mpf_` tuple comparisons. The historical table has therefore been relabelled as an 80-printed-digit agreement; no underlying value or timing was changed.

Performance claims require a matched before/after harness with identical code scope, workload, precision, backend, worker count, hardware and timing boundary. The distributed evidence does not support a universal python-flint multiplier, so none is claimed.

## Cross-references

- [Package source](../connes_cvs/)
- [Regression and identity tests](../tests/)
- [Top-level performance section](../README.md#performance)
