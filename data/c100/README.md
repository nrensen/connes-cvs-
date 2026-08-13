# `data/c100/` - Verification data for the $c = 100$ out-of-sample test (public)

**Contents.** Twelve small JSON files containing the public $c=100$ finite-cell results and analysis inputs. The four-point $N$-sweep supports two Aitken-$\Delta^2$ estimates that approach the Connes 2026 §6.4 heuristic continuum prediction $\approx -530.38$ as the anchor moves from $(100,150,200)$ to $(150,200,250)$ (gaps 6.39 and 3.32 OOM). The consecutive first-difference ratios 0.8373 and 0.8355 are consistent with a local geometric model over these four points; they do not rule out other convergence laws.

A minimal verification script is at [`../../examples/c100_aitken_check.py`](../../examples/c100_aitken_check.py); it loads the JSONs in this directory and reproduces both Aitken triples plus the Connes prediction in under a second.

## Files

| File | Content |
|---|---|
| `results_c100_N100_T800_dps500_v020.json` | $c=100$, $N=100$, $T=800$, $\mathrm{dps}=500$. $\lambda_{\min}^{\mathrm{even}} \approx 1.22 \times 10^{-191}$. |
| `results_c100_N150_T800_dps500_v020.json` | $c=100$, $N=150$, $T=800$, $\mathrm{dps}=500$. $\lambda_{\min}^{\mathrm{even}} \approx 6.42 \times 10^{-248}$. |
| `results_c100_N200_T800_dps500_v020.json` | $c=100$, $N=200$, $T=800$, $\mathrm{dps}=500$. $\lambda_{\min}^{\mathrm{even}} \approx 4.87 \times 10^{-295}$. |
| `results_c100_N250_T800_dps500_v020.json` | $c=100$, $N=250$, $T=800$, $\mathrm{dps}=500$. $\lambda_{\min}^{\mathrm{even}} \approx 2.08 \times 10^{-334}$. |
| `results_c100_N150_T800_dps1000_v020.json` | $c=100$, $N=150$, $T=800$, $\mathrm{dps}=1000$. Precision retest of the $N=150$ row above. $\lambda_{\min}^{\mathrm{even}}$ agrees with the $\mathrm{dps}=500$ value to 25 leading significant digits. |
| `results_c67_N150_T800_dps500_v020.json` | $c=67$, $N=150$, $T=800$, $\mathrm{dps}=500$. Corroborative measurement at the deepest in-sample cutoff. $\lambda_{\min}^{\mathrm{even}} \approx 5.33 \times 10^{-219}$, a 46-OOM drop versus the same-$c$, $N=100$ value reported in [`../results_15pt_T800.json`](../results_15pt_T800.json). |
| `richardson_n_extrapolation.json` | Aitken-$\Delta^2$ extrapolation of the three-point sequence $\{\log_{10}\lvert\lambda_N\rvert\}$ at $c=100$, $N\in\{100,150,200\}$. Reports `aitken: -536.965`. Note: the four-point analysis (incorporating $N=250$) is computed in [`../../examples/c100_aitken_check.py`](../../examples/c100_aitken_check.py); the deeper-anchored triple gives Aitken $\approx -533.70$. |
| `c100_N150_dps1000_gamma_extraction.json` | $\gamma_1$ through $\gamma_{10}$ extraction from the $c=100$, $N=150$, $\mathrm{dps}=1000$ smallest-positive eigenvector. Matching-digit counts range 219–242. Per-$\gamma_k$: detected value, true `mp.zetazero(k).imag` reference, error magnitude, log10 error. Independent verification path documented in the embedded `verification_protocol` field. |
| `c100_N150_gamma_extraction_retight.json` | $N=150$, $\mathrm{dps}=500$ root-extraction retest with tolerance $10^{-100}$; stores per-zero error and log10 error for the retight curve in `make_fig10_c100_gamma_digits.py`. |
| `c100_N250_dps500_gamma_extraction.json` | $\gamma_1$ through $\gamma_{10}$ extraction from the $c=100$, $N=250$, $\mathrm{dps}=500$ smallest-positive eigenvector with tight findroot tolerance $10^{-380}$. Matching-digit counts range 307–329 (deeper than the $N=150$, $\mathrm{dps}=1000$ extraction). Both detected and reference values stored to 400 significant digits; the `verification_protocol` field describes how to independently confirm `matching_digits = floor(-log10(error))` using `mp.zetazero(k)` at $\mathrm{dps}=400$. |
| `c100_N120_dps560_gamma_extraction.json` | $\gamma_1$ through $\gamma_{10}$ from the $c=100$, $N=120$, $T=800$, $\mathrm{dps}=560$ smallest-positive eigenvector (tight findroot tolerance $10^{-520}$). $\gamma_1$ matches `mp.zetazero(1)` to 210 digits; per-$k$ counts range 186–210. Stores **detected** $\gamma_k$ and the true reference to 500 significant digits (verified dps-stable to $>515$ digits), enabling an eigenvalue-vs-eigenvalue diff. Off-sweep cross-validation cell at an independent grid point (not part of the four-point Aitken sweep); reproduction recipe in the embedded `verification_protocol`. |
| `c100_N160_dps560_gamma_extraction.json` | $\gamma_1$ through $\gamma_{10}$ from the $c=100$, $N=160$, $T=800$, $\mathrm{dps}=560$ smallest-positive eigenvector (tight findroot tolerance $10^{-520}$). $\gamma_1$ matches `mp.zetazero(1)` to 253 digits; per-$k$ counts range 230–253. Stores **detected** $\gamma_k$ and the true reference to 500 significant digits, enabling an eigenvalue-vs-eigenvalue diff. Off-sweep cross-validation cell at an independent grid point; reproduction recipe in the embedded `verification_protocol`. |

## Schema (representative)

The five `results_c100_*_v020.json` files share a common shape:

```json
{
  "tag":          "...",
  "c":            100,
  "N":            150,
  "T":            800,
  "dps":          500,
  "flint_prec":   2000,
  "engine":       "v020",
  "lambda_even":  "<decimal string as recorded by the run>",
  "t_cache_s":    <float>,
  "t_mat_s":      <float>,
  "t_diag_s":     <float>,
  "t_total_s":    <float>,
  "n_workers":    12,
  "version":      "v0.2.3-local"
}
```

`lambda_even` is always stored as a decimal string rather than a Python float, so magnitudes down to $10^{-334}$ are not underflowed. String lengths vary by artifact: most production rows retain hundreds of digits, while the $N=150$, $\mathrm{dps}=500$ row intentionally records a shorter 48-significant-digit value. Do not infer unrecorded trailing digits.

`richardson_n_extrapolation.json` carries the three-point $N$-sweep array (legacy; pre-$N{=}250$), fitted models (exponential, power, $1/N$, stretched-exponential), and the Aitken-$\Delta^2$ acceleration scalar for that three-point case.

`c100_N*_gamma_extraction.json` files list per-$k$ detected $\gamma_k$ (decimal-string), true `mp.zetazero(k).imag` reference, absolute error, $\log_{10}$ error, and floor matching-digit count. The $N=250$ extraction stores all $\gamma$ values to 400 significant digits, sufficient for independent verification past the 329-digit headline.

## Provenance

The production result files were produced on a 12-worker Apple M-series workstation using the v0.2.0 mathematical core with a local cell runner. Each result JSON records `flint_prec=4*dps`; this differs from the package's historical default and must be supplied explicitly for an exact reproduction attempt. Recorded wall times are:

- $c=100$, $N=100$, $\mathrm{dps}=500$ - 13.9 min
- $c=100$, $N=150$, $\mathrm{dps}=500$ - ~21 min
- $c=100$, $N=200$, $\mathrm{dps}=500$ - 28.4 min
- $c=100$, $N=250$, $\mathrm{dps}=500$ - ~38 min
- $c=100$, $N=150$, $\mathrm{dps}=1000$ - ~111 min
- $c=67$, $N=150$, $\mathrm{dps}=500$ - ~26 min

The two $N=150$ rows at $\mathrm{dps}\in\{500,1000\}$ agree for 25 leading significant digits and differ beginning at digit 26. This supports the recorded leading digits at that $N=150$ cell only; it is not a precision certificate for $N=250$.

## Exact package recipe and branch scope

```python
from connes_cvs.runner import CellConfig, GalerkinCell

if __name__ == "__main__":
    dps = 500
    cell = GalerkinCell(
        CellConfig(c=100, N=250, T=800, dps=dps, flint_bits=4*dps),
        ground_state="smallest_positive",
    )
    artifact = cell.run()
    print(artifact["lambda_even"])
```

At $c=100$, $T=800$, the raw finite-$T$ matrix has negative-sign eigenvalues that disappear at larger archimedean cutoff. The published rows therefore follow the empirically distinguished smallest-positive even-sector branch; `compute_ground_state`, which returns the raw minimum, is not the selector for these cells. Full gamma regeneration uses `cell.eigvec_full` after `cell.run()` and `extract_zeros(..., c=100, tol=...)`. It is expensive and should use the tolerance and exact precision configuration stated in each JSON protocol.

## Cross-references (public)

- Top-level headline + Aitken match: [`../../README.md`](../../README.md)
- Verification script: [`../../examples/c100_aitken_check.py`](../../examples/c100_aitken_check.py)
- 15-cutoff sweep data: [`../results_15pt_T800.json`](../results_15pt_T800.json)
- Package source: [`../../connes_cvs/`](../../connes_cvs/)
