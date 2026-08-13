[**← `connes-cvs`**](../../README.md) · [**Papers**](../README.md) &nbsp;|&nbsp; **Riemann zeros** · [Guinand-Weil dictionary](../2_guinand_weil_dictionary_tail_order/) · [von Mangoldt measure](../3_matrix_von_mangoldt_measure/)

<div align="center">

# High-Precision Approximation of Riemann Zeros<br>via the Truncated Weil Form

**Riemann zeros at high precision - _the numerics_ · Akiva Groskin, 2026**

[![arXiv](https://img.shields.io/badge/arXiv-2605.20224-b31b1b.svg)](https://arxiv.org/abs/2605.20224)
[![Zenodo DOI](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.19546514-1682D4.svg?logo=zenodo&logoColor=white)](https://doi.org/10.5281/zenodo.19546514)
[![PyPI](https://img.shields.io/pypi/v/connes-cvs.svg?label=connes-cvs&color=4c1)](https://pypi.org/project/connes-cvs/)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-blue.svg)](https://doi.org/10.5281/zenodo.19546514)

</div>

> Builds and diagonalizes the finite Connes–van Suijlekom **Galerkin matrix** at high precision,
> extracting Riemann zeros to **hundreds of matching digits** and giving an independent
> out-of-sample test of the Connes 2026 §6.4 continuum asymptotic. Empirical measurements only;
> no claim regarding the Riemann Hypothesis.

Part of the [`connes-cvs` series](../../README.md#papers): **Riemann zeros - the numerics** · [**Guinand-Weil dictionary** - the structure](../2_guinand_weil_dictionary_tail_order/) · [**von Mangoldt measure** - the arithmetic](../3_matrix_von_mangoldt_measure/). Published on arXiv, [arXiv:2605.20224](https://arxiv.org/abs/2605.20224) (math.NT); archived on Zenodo, concept DOI [10.5281/zenodo.19546514](https://doi.org/10.5281/zenodo.19546514) (resolves to the latest version).

## About this folder

`high_precision_approximation_of_riemann_zeros.pdf` is the published manuscript (as on arXiv).
Unlike the two companion notes - self-contained theorem notes whose code lives beside them - **this paper's
reproducibility package is the `connes-cvs` package itself**, which is the whole repository. So
this folder holds the paper; the code and data live at the repository root:

| Artifact | Location |
| :--- | :--- |
| The implementation (builds + diagonalizes the CvS Galerkin matrix, extracts zeros) | [`../../connes_cvs/`](../../connes_cvs/) - installable as `pip install connes-cvs` |
| The `c = 100` datasets (`gamma` extractions, `N`-sweep) | [`../../data/c100/`](../../data/c100/) |
| Runnable examples (incl. the Aitken-Δ² check) | [`../../examples/`](../../examples/) |
| Regression tests against the committed `c = 13` references | [`../../tests/`](../../tests/) |
| Performance A/B benchmarks | [`../../benchmarks/`](../../benchmarks/) |
| Errata (four entries; see below) | [`../../ERRATA.md`](../../ERRATA.md) |

The [repository root README](../../README.md) documents the headline result, installation,
quick start, the `c = 100` verification, and how it works.

## What the paper does

Builds and diagonalizes the finite Connes-van Suijlekom Galerkin matrix of the truncated Weil
quadratic form at high precision. It computes the smallest-positive eigenvalue across a
15-cutoff sweep and at `c = 100`, extracts the first ten Riemann zeros to hundreds of matching
digits, and gives an independent out-of-sample numerical test of the Connes 2026 §6.4 continuum
asymptotic. It reports empirical measurements only; it makes no claim regarding the Riemann
Hypothesis. See [`../../ERRATA.md`](../../ERRATA.md) for the four recorded corrections, none of
which changes a measured value: the 2026-06-26 finite-cutoff sign correction; the withdrawal of
the Section 8.2 Paley-Wiener mechanism (Table 14 is measured at `T = 400`, not `T = 800`); the
`c = 67` matching-digit count, 167 rather than 168 under the paper's own definition; and two
summary digit-increment ranges in Section 6.5, 93-117 and 181-203.

All four are incorporated in the revised manuscript text deposited on Zenodo as **Version 3.6**
([10.5281/zenodo.21918387](https://doi.org/10.5281/zenodo.21918387)). **The PDF in this folder is
the currently announced arXiv v2**, which carries the 2026-06-26 correction but predates the three
later ones; its replacement has been submitted to arXiv and is awaiting announcement. Read the
arXiv v2 PDF against `ERRATA.md`, or use the Zenodo Version 3.6 manuscript, which needs no errata
overlay.

## License

Manuscript: CC BY 4.0. The `connes-cvs` package is [MIT](../../LICENSE).
