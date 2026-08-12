# `examples/` - Minimal examples for the `connes-cvs` package

## What's here

| File | Purpose |
|---|---|
| `basic_compute.py` | Defaults to the seconds-level `c=13, N=8, T=60, dps=30` API smoke cell and one seeded finite-test root. `--extended` runs the `N=100, T=400, dps=80` validation cell through the process-based runner, with progress and the default eight-worker cap, then locates three roots. Heavy work is protected by `main()`, so importing the module is safe. |
| `c100_aitken_check.py` | Loads the published `c = 100` N-sweep data from [`../data/c100/`](../data/c100/) and reproduces the Aitken extrapolation check in under a second. |
| `make_fig9_c100_aitken.py`, `make_fig10_c100_gamma_digits.py` | Generators for the two `c = 100` figures of the Paper 1 reproducibility package (Aitken extrapolation; canonical `floor(-log10(error))` counts for `gamma_1..gamma_10`). They require Matplotlib in addition to the package dependencies. Papers 2 and 3 keep their own figure generators inside their [paper folders on GitHub](https://github.com/akivag613/connes-cvs-/tree/main/papers). |

## Running

```bash
pip install connes-cvs
python examples/basic_compute.py

# Multi-minute validation cell, with progress and process isolation:
python examples/basic_compute.py --extended

# Figure generators additionally require matplotlib:
python -m pip install matplotlib
```

Expected output for the default `c = 13`, `N = 8`, `T = 60`, `dps = 30` smoke cell:

- `lambda_min` approximately `4.43043e-23`
- `gamma_1` detected as `14.13472514173469...` with `|gamma_1 error|` approximately `2.52738e-17`

Measured release-environment runtime was about 1.8 seconds with python-flint 0.8.0 and 4.8 seconds through the mpmath fallback; runtime varies by versions and hardware. This is only an API smoke test. The `--extended` mode is the `T=400` validation cell whose committed slow gate requires at least 20 lambda digits and 18 gamma-error digits. The paper's c=13 reference datum (`|gamma_1 error|` approximately `2.005e-55`) is the separate `N=100`, `T=800`, `dps=150` cell.

## Cross-references (public)

- Top-level README + headline result: [`../README.md`](../README.md)
- Package source: [`../connes_cvs/`](../connes_cvs/)
- Canonical paper data: [`../data/results_15pt_T800.json`](../data/results_15pt_T800.json)
- Tests: [`../tests/`](../tests/)
- Benchmarks: [`../benchmarks/`](../benchmarks/)
