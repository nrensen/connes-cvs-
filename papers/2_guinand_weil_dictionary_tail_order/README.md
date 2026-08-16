[**← `connes-cvs`**](../../README.md) · [**Papers**](../README.md) &nbsp;|&nbsp; [Riemann zeros](../1_high_precision_riemann_zeros/) · **Guinand-Weil dictionary** · [von Mangoldt measure](../3_matrix_von_mangoldt_measure/)

<div align="center">

# A finite Guinand-Weil dictionary and archimedean tail order<br>for the truncated Weil quadratic form

**The finite Guinand-Weil dictionary - _the structure_ · Akiva Groskin, 2026**

[![arXiv](https://img.shields.io/badge/arXiv-2607.02828-b31b1b.svg)](https://arxiv.org/abs/2607.02828)
[![Zenodo DOI](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.21124802-1682D4.svg?logo=zenodo&logoColor=white)](https://doi.org/10.5281/zenodo.21124802)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-blue.svg)](LICENSE-PAPER-CC-BY-4.0.txt)

</div>

> An exact finite **Guinand-Weil zero-source dictionary** for the truncated Weil form: every value
> of the quadratic form is an exact sum over the nontrivial zeros of ζ. Plus a finite-cutoff
> **archimedean tail-order theorem** with a two-sided certification rule. No claim regarding the
> Riemann Hypothesis.

Part of the [`connes-cvs` series](../../README.md#papers): [**Riemann zeros** - the numerics](../1_high_precision_riemann_zeros/) (the [`connes-cvs`](../../README.md) package) · **Guinand-Weil dictionary - the structure** · [**von Mangoldt measure** - the arithmetic](../3_matrix_von_mangoldt_measure/). Published on arXiv, [arXiv:2607.02828](https://arxiv.org/abs/2607.02828) (math.NT, math.SP), and archived on Zenodo, concept DOI [10.5281/zenodo.21124802](https://doi.org/10.5281/zenodo.21124802) (resolves to the latest version).

**Corrections.** This paper has no errata: **no mathematical statement has been corrected
since publication.** The PDF here is byte-identical to the manuscript in the Zenodo deposit
that the concept DOI above resolves to. One published ancillary value was corrected on
2026-08-14: the derivative envelope ladder in
`artifacts/arch_tail_exact_vs_asymptotic.json` reported `h_+'(t)` at `t = 50`, `100` and
`1000` from a series evaluation that `mp.nsum` had accelerated incorrectly at `dps = 25`, most
visibly at `t = 1000`, where the published `0.0002384087979` should read `0.001000000083`.
The formula in the script was correct and is now evaluated directly from the trigamma
function; Lemma 3.1's envelope holds at every `t` under the corrected values, and
`B_quadrature` (renamed from `B_exact` on 2026-08-16 because it is a numerical evaluation with a remainder bound, not an interval-exact value), `B_asym` and the certification-floor solve are unchanged. **No statement in the
manuscript quotes these values, so this correction does not change the manuscript.** The defect
was identified by B. W. A. Silva.

Separately, and also on 2026-08-14, two references were updated to the journal versions in which
they have since appeared: Connes and van Suijlekom to *Communications in Mathematical Physics*
**406** (2025), article 312, and Connes, Consani and Moscovici to a chapter in *Applications of
Noncommutative Geometry to Gauge Theories, Field Theories, and Quantum Space-Time*, EMS Series of
Lectures in Mathematics, EMS Press (2026), pages 39-76. That is the only reason the PDF differs
from the one deposited before that date; **no mathematical statement changes.** The full change
history is in the repository [CHANGELOG](../../CHANGELOG.md) and in the Zenodo record.

## What this paper proves

For a real even finite Connes-van Suijlekom / Connes-Consani-Moscovici Galerkin coefficient
vector at cutoff `c` and band `N`:

1. an **exact finite Guinand-Weil zero-source dictionary** (Theorem 2.5): every value of the
   truncated Weil quadratic form is an exact sum over the nontrivial zeros of the Riemann zeta
   function, via a band-limited Guinand-Weil test function, with an entry-identification lemma
   pinning the source assembly to the CCM closed forms at equation level (Lemma 2.1);
2. an exact finite **source quotient** of dimension `2N+1`, stated as an iff (Corollary 2.4),
   and a positive-dimensional non-collapsing **pole-neutral** subfamily (Corollary 2.7);
3. a finite-`T` **archimedean tail-order theorem** (Theorem 3.2): past the Galerkin band the
   omitted archimedean tail is a strictly positive definite, strictly totally positive
   Cauchy-Stieltjes increment;
4. a **two-sided certification rule** with an explicit closed-form budget and asymptotic
   `B_T = (2N+1) rho (log(T/2pi)+1)/(pi^2 T)(1+o(1))` (Corollary 3.3), plus a self-contained
   `h_+` envelope lemma (Lemma 3.1).

The manuscript includes a worked example verified against the first 512 nontrivial zeros of
zeta and an eigenvalue-flow demonstration of the certification rule. It makes no claim
regarding the Riemann Hypothesis, Weil positivity, or a prime-location bound.

## Layout

```text
finite_guinand_weil_dictionary_tail_order.pdf   the paper (15 pp)
README.md                                       this file
VERIFICATION.md                                 what each guard checks
requirements.txt                                dependencies
LICENSE-PAPER-CC-BY-4.0.txt · LICENSE           licenses (paper CC-BY-4.0; scripts MIT)
SHA256SUMS                                      checksums for every file
source/      main.tex, main.bib, main.bbl, plainurl.bst   (LaTeX source)
figures/     fig_dictionary.pdf, fig_tailorder.pdf + make_figures.py
scripts/     verification guards (exact symbolic/integer + Arb interval)
artifacts/   guard outputs (JSON) and the 9000-bit certificate log
audit/       CLAIM_TRACE_AUDIT.md, NOVELTY_BOUNDARY_AUDIT.md
```

The `scripts/` guards, grouped by result:

- **Dictionary / source quotient** - `audit_exact_series_identity.py` (single-frequency
  identity), `audit_kernel_span_rank.py` (finite Volterra-kernel span), `audit_full_matrix_source_quotient.py`
  (factorization through the `2N+1` quotient), `audit_pole_neutral_survival.py` (pole-square
  factorization + dimension formula), `verify_finite_dictionary.py` (single-frequency source
  identity + pole normalization by three routes), `verify_dictionary_threeroute.py` (three-route
  dictionary confirmation over the first 512 zeros), `verify_zero_side.py` (original `c=13,N=4`
  zero-side check, retained for continuity).
- **Archimedean tail order** - `verify_arch_tail_order.py` (tail-order algebra + strict total
  positivity), `audit_arch_tail_dt_bridge.py` (rank-two Cauchy density = finite-`T` derivative),
  `arch_tail_budget.py` (Arb interval budget at `c=100, N=200, T=800`), `arch_tail_stress_ladder.py`
  (interval stress ladder across `T, N`, precision), `arch_tail_exact_asymptotic.py` (exact `B_T`
  vs the closed asymptotic; the `B_T = 1e-59` solve).
- **Cutoff-free inertia certificate** - `arb_ldlt_certify.py` (rigorous Arb interval `LDL^T`
  inertia certificate; generator of the 9000-bit `c=100, N=200` certificate).

## Reproduce

Install the pinned environment in `requirements.txt`. The exact integer guard
`audit_exact_series_identity.py` uses only the standard library. Four symbolic
guards (`audit_kernel_span_rank.py`, `audit_full_matrix_source_quotient.py`,
`audit_pole_neutral_survival.py`, and `audit_arch_tail_dt_bridge.py`) require
`sympy`; the numerical dictionary/asymptotic checks require `mpmath`. The four
Arb interval scripts (`verify_arch_tail_order.py`, `arch_tail_budget.py`,
`arch_tail_stress_ladder.py`, and `arb_ldlt_certify.py`) additionally require
`python-flint`. The figure generator requires both `mpmath` and `matplotlib`.

```bash
PYTHONDONTWRITEBYTECODE=1 python3 scripts/audit_exact_series_identity.py
PYTHONDONTWRITEBYTECODE=1 python3 scripts/audit_kernel_span_rank.py
PYTHONDONTWRITEBYTECODE=1 python3 scripts/audit_full_matrix_source_quotient.py
PYTHONDONTWRITEBYTECODE=1 python3 scripts/audit_pole_neutral_survival.py
PYTHONDONTWRITEBYTECODE=1 python3 scripts/verify_finite_dictionary.py
PYTHONDONTWRITEBYTECODE=1 python3 scripts/audit_arch_tail_dt_bridge.py
PYTHONDONTWRITEBYTECODE=1 python3 scripts/verify_zero_side.py
PYTHONDONTWRITEBYTECODE=1 python3 scripts/verify_dictionary_threeroute.py 512
PYTHONDONTWRITEBYTECODE=1 python3 scripts/arch_tail_exact_asymptotic.py
# with python-flint:
PYTHONDONTWRITEBYTECODE=1 python3 scripts/verify_arch_tail_order.py
PYTHONDONTWRITEBYTECODE=1 python3 scripts/arch_tail_budget.py
PYTHONDONTWRITEBYTECODE=1 python3 scripts/arch_tail_stress_ladder.py
PYTHONDONTWRITEBYTECODE=1 python3 scripts/arb_ldlt_certify.py --selftest --c 13 --N 8 --prec 300
```

The archived 2026-07-02 headline-certificate run took about 15 minutes on one
Apple M2 Max core under Python 3.12.11 and `python-flint` 0.8.0. This is a
historical measurement, not a performance guarantee; runtime varies by hardware
and dependency build. Regenerate the certificate with:

```bash
python3 scripts/arb_ldlt_certify.py --selftest --c 100 --N 200 --prec 9000 \
        --json-out artifacts/c100_N200_arb_ldlt_prec9000_provenance.json
```

See `VERIFICATION.md` for the guard-to-theorem map, and `audit/CLAIM_TRACE_AUDIT.md` for the
claim-to-artifact trace.

For archive integrity, run `shasum -a 256 -c SHA256SUMS` before regenerating an
artifact in place. In particular, the explicit `--json-out` command above writes
the current date and measured build/LDL timings into its provenance JSON, so a
fresh successful run is expected to differ from the archived JSON checksum in
those volatile provenance fields.

## Build the paper

The compiled PDF is provided at the top level. To rebuild, first regenerate the
figures if desired, then copy the figure PDFs next to `source/main.tex`:

```bash
python3 figures/make_figures.py  # optional; rewrites figures/fig_*.pdf
cp figures/fig_*.pdf source/
cd source && pdflatex main && bibtex main && pdflatex main && pdflatex main
```

## License

Manuscript and figures: [CC BY 4.0](LICENSE-PAPER-CC-BY-4.0.txt). Verification scripts:
[MIT](LICENSE). Checksums for every file are in `SHA256SUMS`.
