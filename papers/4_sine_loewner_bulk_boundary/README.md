[**← `connes-cvs`**](../../README.md) · [**Papers**](../README.md)

# Bulk trace laws and algebraic boundary interaction for the sine Loewner operator

Akiva Groskin. First published 2026-09-08.

[![Zenodo DOI](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.22662221-1682D4.svg?logo=zenodo&logoColor=white)](https://doi.org/10.5281/zenodo.22662221)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-blue.svg)](source/LICENSES.md)

The complete paper has **88 pages, three figures and full proofs**. It proves uniform
trace laws for the sine Loewner operator and its finite matrices, determines half-line
boundary flow and algebraic interaction, and certifies a boundary mode with nonzero
even–odd splitting. The results include infinite-dimensional spectral theorems and
make no claim of proving the Riemann Hypothesis.

Read the [paper](paper.pdf). The [Zenodo concept DOI](https://doi.org/10.5281/zenodo.22662221)
resolves to the current release; the [first-publication version DOI](https://doi.org/10.5281/zenodo.22662222)
identifies this release. The PDF and the 64 files under `source/` match that release.
See [CITATION.bib](CITATION.bib) for the paper citation.

## Contents

```text
paper.pdf                complete manuscript
source/main.tex          build entry point
source/references.tex    bibliography
source/sections/         all statements and proofs
source/figures/          TikZ inputs
source/anc/              verification programs, helper and exact witnesses
source/SOURCE_README.md  full build instructions
source/LICENSES.md       manuscript, figure, witness and software licenses
SHA256SUMS               integrity checks for this GitHub folder
```

The complete source layout is retained so its relative paths work directly. The
source ZIP is available from Zenodo; it is not duplicated in this folder.

## Build

With a standard TeX Live installation, run from this paper folder:

```bash
cd source
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error main.tex
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error main.tex
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error main.tex
```

The bibliography is already in `references.tex`; the figures use TikZ. No BibTeX,
external image conversion or shell escape is required. See
[source/SOURCE_README.md](source/SOURCE_README.md) for the recorded build environment.

## Reproduce the checks

Use Python 3.12 and the [pinned requirements](source/anc/requirements.txt). Read
[source/anc/README.md](source/anc/README.md), then follow the commands in the
[trace guide](source/anc/trace-README.md) and
[boundary guide](source/anc/boundary-README.md). The standard replay has **23 programs**:

- 18 floating-point diagnostics;
- two exact rational scalar checks;
- three outward-rounded interval certificate programs for existence, exclusion and residue.

The diagnostics test finite samples and do not replace proofs. The interval
certificates use the exact witnesses and the analytic operator-error estimates
proved in the paper. All proofs are in the manuscript; no external mathematical
supplement or private input is needed.

These programs are independent of the `connes-cvs` Python package and its PyPI
release. The repository's package CI does not run the 23 paper checks. Run them
sequentially with `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1` and
`VECLIB_MAXIMUM_THREADS=1`, as directed in the guides.

Before rebuilding or replaying, verify the released bytes from this folder:

```bash
shasum -a 256 -c SHA256SUMS
```

Keep regenerated outputs separate and choose fresh output names.

## License

Manuscript, figures and exact mathematical witness data: CC BY 4.0. Original
Python programs: MIT. See [source/LICENSES.md](source/LICENSES.md). Third-party
libraries retain their own licenses.
