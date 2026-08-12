# WIN 1 A/B Verification - Production Scale

**Date:** 2026-04-14
**Branch:** `feature/v0.2.0-matrix-microopt` (WIN 1 still uncommitted)
**Method:** Stash WIN 1, run baseline; pop stash, run WIN 1; compare in the same recorded environment.

**Evidence boundary:** this is a frozen historical summary. Its raw console logs were local development artifacts and are not distributed. “80 digits agree” below means the two stored decimal renderings agree in all 80 printed significant digits; it does not assert equality of raw mpmath tuples.

## Test config

- `c=13, N=80, T=400, dps=80`
- 12-way multiprocessing Pool (production code path via `sweep._run_single_cutoff`)
- Python 3.12.11, mpmath 1.4.1, python-flint enabled (HAS_FLINT=True)
- Matrix DIM = 2N+1 = 161

## Results

| Phase | Baseline (no WIN 1) | WIN 1 | Speedup |
|---|---|---|---|
| cache_sec | 57.554 | 27.941 | **2.06×** |
| matrix_sec | 0.112 | 0.119 | ~unchanged |
| diag_sec | 6.113 | 6.185 | ~unchanged |
| zeros_sec | 1.156 | 1.149 | ~unchanged |
| **WALL TOTAL** | **64.939** | **35.398** | **1.83×** |
| **λ_min** | 2.5282661401965756026025862533001704392434144948201908268289070778008968019858182e-59 | 2.5282661401965756026025862533001704392434144948201908268289070778008968019858182e-59 | **all 80 printed digits agree** |

## Implication

- Cache phase 2.06× confirmed at production scale.
- Total wall 1.83× (matrix/diag/zeros are unchanged ~7.5s constant overhead).
- The stored multiprocessing-path renderings agree in all 80 printed digits.

## Extrapolation against historical production records

| Historical record | c | N | T | dps | Old | Projected w/ WIN 1 (cache × 0.486 + other) |
|---|---|---|---|---|---|---|
| results_A_c13.json | 13 | 100 | 400 | 150 | 199s | ~108s |
| results_U_T800_c13.json | 13 | 100 | 800 | 150 | 215s | ~117s |
| results_U_T800_c67.json | 67 | 100 | 800 | 200 | 543s | ~280s |
| Full 15-cutoff sweep (sum) | - | - | - | - | ~5400s | **~2900s (saves ~40 min)** |

## Files

- `AB_baseline_no_win1.txt` - local-only baseline output (not distributed)
- `AB_with_win1.txt` - local-only WIN 1 output (not distributed)
- `win1_pool_benchmark.py` - current public sweep-path timer
- `win1_benchmark.py` - current checkout versus frozen small-workload numbers
- `baseline_benchmark.py` - current-checkout single-process phase timer

## Gold-standard verification at the published reference workload

The c=13 N=80 T=400 dps=80 A/B above proves the speedup at moderate scale. To prove the new code reproduces the *published headline number* in Paper 1 Table 18 (`λ_min(c=13) = 2.865 × 10⁻⁵⁹` at N=100 T=800 dps=150), v0.2.0 was run end-to-end at the exact paper-reference workload:

- **Computed:** `λ_min = 2.86545361493028029516151514986747977533798676783773219101029565377637421791530494377666704141009139776092287892559119370499915456183497252629004918875e-59`
- **Paper Table 18:** `2.865 × 10⁻⁵⁹` ✓ exact match to all reported precision
- **Wall time:** 127.3 s
- **Historical baseline (v0.1.0 at the same workload):** 214.8 s
- **Speedup at exact paper workload:** **1.69×** total wall

Raw output was recorded in the gitignored development log `AB_published_reference_check.txt`; it is not part of the public repository or source distribution.

## Verdict

The historical record supports the stated speedup in its recorded environment and agreement with the paper's reported lambda precision. Re-measure before making a performance claim about a later release.
