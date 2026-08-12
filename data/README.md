# `data/` - Canonical numerical data (public)

**Contents.**

| Path | Content |
|---|---|
| `results_15pt_T800.json` | The 15-cutoff Paper 1 summary: `lambda_even` and `gamma_1_abs_error` for each `c` in `{13, 14, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67}` at `N = 100`, `T = 800`, `dps = 150` or `200`. |
| `c100/` | The `c = 100` verification dataset: N-sweep results, gamma-extraction tables, and the Richardson/Aitken extrapolation inputs used by `examples/c100_aitken_check.py`. See `c100/README.md`. |

**Purpose.** The public-facing numerical record of Paper 1 and of the `c = 100` verification. Linked from [`../README.md`](../README.md) and cross-validated by the examples and tests. Papers 2 and 3 carry their own reproducibility artifacts inside their [paper folders on GitHub](https://github.com/akivag613/connes-cvs-/tree/main/papers), not here.

## Schema of `results_15pt_T800.json`

```json
{
  "description": "CvS operator reference results at T=800, N=100",
  "source": "https://github.com/akivag613/connes-cvs-",
  "cutoffs": 15,
  "span_oom": 113.1,
  "results": [
    {
      "cutoff": 13,
      "L": 2.564949,
      "lambda_even": "2.865453614930280295161515e-59",
      "gamma_1_abs_error": "2.0054614440503604741098249893e-55",
      "gamma_1_over_lambda": 6998.8,
      "log10_abs_error": -54.7,
      "T": 800,
      "N": 100,
      "dps": 150
    }
  ]
}
```

(`results` holds one such row per cutoff; `cutoffs` is the row count.) The precision-critical values (`lambda_even`, `gamma_1_abs_error`) are stored as decimal strings, not floats; they span 113 orders of magnitude.

## Packaging note (wheel vs sdist)

The PyPI **wheel** contains only the `connes_cvs/` package - installing `connes-cvs` from a wheel does not install `data/`. The reference JSON ships in the **source distribution** (sdist) and in the GitHub repository; clone the repository (or download the sdist) to obtain the data files.

## Reproducing the 15-cutoff dataset

```python
from connes_cvs.sweep import run_sweep

LOW_DPS = [13, 14, 17, 19, 23, 29, 31, 37]
HIGH_DPS = [41, 43, 47, 53, 59, 61, 67]

if __name__ == "__main__":
    results = run_sweep(LOW_DPS, N=100, T=800, dps=150, workers=8)
    results.update(
        run_sweep(HIGH_DPS, N=100, T=800, dps=200, workers=8)
    )
```

`dps = 200` is used for `c >= 41`. The guard is required because `run_sweep` starts a multiprocessing pool on spawn platforms. See [`../README.md`](../README.md) for serialization and comparison details.

## Discipline

- This directory is **public** (git-tracked). Never include in-progress or revision-pending findings here.
- The `results_15pt_T800.json` file is byte-for-byte unchanged across the published Paper 1 versions that carry it (concept DOI [10.5281/zenodo.19546514](https://doi.org/10.5281/zenodo.19546514), currently resolving to Version 3.3, version DOI [10.5281/zenodo.20931069](https://doi.org/10.5281/zenodo.20931069)). If a correction is ever needed, deposit a new Zenodo version and update this file in lockstep.
