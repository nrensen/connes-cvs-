[**← `connes-cvs`**](../../README.md) · [**Papers**](../README.md) &nbsp;|&nbsp; [Riemann zeros](../1_high_precision_riemann_zeros/) · [Guinand-Weil dictionary](../2_guinand_weil_dictionary_tail_order/) · **von Mangoldt measure**

# A matrix-valued von Mangoldt measure in the finite Connes-van Suijlekom path (version 2, corrected)

Akiva Groskin, 2026. Manuscript and full reproducibility package, version 2.

Archived on Zenodo, concept DOI
[10.5281/zenodo.21242028](https://doi.org/10.5281/zenodo.21242028), which always
resolves to the current version. The manuscript is the corrected version 2 of
2026-07-27; later versions have revised documentation in the reproducibility
archive only, leaving the manuscript unchanged. For theorem-sensitive use, cite
the version DOI of the specific version you consulted, which each Zenodo record
displays. Companion to the Riemann-zeros paper
([arXiv:2605.20224](https://arxiv.org/abs/2605.20224)) and the Guinand-Weil dictionary
([arXiv:2607.02828](https://arxiv.org/abs/2607.02828)). This paper is not on arXiv.

**Corrections.** Version 2 corrected version 1 substantively, and **version 1 remains publicly
citable** at [10.5281/zenodo.21242029](https://doi.org/10.5281/zenodo.21242029). Anyone holding
that v1 DOI should read [What version 2 corrects](#what-version-2-corrects-relative-to-the-version-1-deposit)
below before relying on it; the same list ships in the README of the deposited reproducibility
archive. The guard script `scripts/check_negative_controls.py` verifies that each corrected
version-1 statement fails on its counterexample and each version-2 statement holds. The PDF
here is the corrected version 2 manuscript.

On 2026-08-14 two references were updated to the journal versions in which they have since
appeared: Connes and van Suijlekom to *Communications in Mathematical Physics* **406** (2025),
article 312, and Connes, Consani and Moscovici to a chapter in *Applications of Noncommutative
Geometry to Gauge Theories, Field Theories, and Quantum Space-Time*, EMS Series of Lectures in
Mathematics, EMS Press (2026), pages 39-76. That is the only reason the PDF differs from the one
deposited before that date; the text is otherwise the corrected version 2 of 2026-07-27, and
**no theorem, proof, or numerical result changes.**

Fix a Galerkin level `N` in the finite Connes-van Suijlekom truncation of the Weil
quadratic form (no archimedean cutoff), and vary the prime cutoff `u = log c`.
Differentiating the finite matrix path `u -> Q_N(u)` across a prime-power threshold
`u = log q` returns the von Mangoldt weight exactly, in every matrix entry: the
first-derivative jump is `-2 Lambda(q)/(sqrt(q) log q)` times the all-ones rank-one
matrix. This identity is an elementary structural derivative of the path's
defining prime sum; the contribution is the isolation and naming of the
resulting finite matrix-valued measure and the exact finite geometry around it.
The paper proves the event is arithmetically rigid, develops the finite
source-to-jet dictionary (confluent Vandermonde, sharp `2N+1` window, universal
recurrence), a sharp finite vanishing-moment ceiling at the prime edge (a
prime-edge uncertainty principle in the band-limited sense), a
coincidence-resolvent generating identity stated by spectral projections, and a
rank-one Weyl-function increment; the Krein-string reading is retained only as
an explicitly labeled analogy. It proves no positivity, no Riemann Hypothesis,
and no prime-counting, next-prime, or factoring statement.

## What version 2 corrects (relative to the version 1 deposit)

- The coincidence-averaging corollary now carries its pairwise-uncorrelatedness
  hypothesis explicitly; a new remark gives the correlated-noise counterexample,
  the generalized least-squares statement for general covariance, and the
  symmetric-measurement variance.
- The directional-residual and singular-ratio diagnostics are no longer called
  equivalent; the text states what each measures.
- The coincidence-resolvent pole statement is stated by spectral projection,
  allowing removable points and repeated eigenvalues.
- The reciprocal Weyl identity is stated meromorphically with its pointwise
  domain made exact; the nonreciprocal Sherman-Morrison form is the primary
  identity.
- The Krein-string boundary-mass reading is downgraded from a theorem to an
  explicitly labeled analogy, with the missing Stieltjes/positivity hypotheses
  stated.
- Scope calibration: the Selberg-class transfer is conditioned on an actually
  constructed finite path, conditional statements are marked conditional (the
  spectral-barrier statement is now a conditional remark stating its
  positive-definiteness hypothesis), the vanishing-moment result carries its
  calibrated name (with the uncertainty-principle reading kept as a gloss),
  the verification table maps each row to a named guard script and JSON
  artifact and marks the variance reduction as an IID-model statement, and
  point-of-use attributions are added (Connes-van Suijlekom Props. 4.1-4.2;
  Andrade's pole and prime notes; Andrews' quantitative convergence law;
  Kac-Krein and Gesztesy-Simon for the string and Jacobi classification
  facts).

The central rank-one prime-power jump, the distributional matrix-valued measure,
the dictionary, uncertainty, and deflation results are unchanged.

## Layout

```text
matrix_valued_von_Mangoldt_measure_finite_CvS_path.pdf      the paper (18 pp)
README.md                                                this file
VERIFICATION.md                                          what each check verifies
requirements.txt                                         optional dependencies
LICENSE-PAPER-CC-BY-4.0.txt                              license (CC BY 4.0)
SHA256SUMS                                               checksums for every file
source/      main.tex, main.bib, main.bbl, plainurl.bst  (LaTeX source)
figures/     fig_*.pdf + make_figures.py                 (the three figures + generator)
scripts/     check_*.py                                  (13 reproducibility guards + 1 negative-control guard)
artifacts/   *.json                                      (14 check outputs, status PASS)
```

## Reproduce

The thirteen version 1 reproducibility guards are separate checks (exact symbolic,
exact modular over prime fields to `N=1000`, and floating-point checks including
the canonical scale `N=200`). Eight use only the Python 3 standard library; five
use `sympy` and/or `numpy` (see `requirements.txt`). The fourteenth script,
`check_negative_controls.py` (standard library only, new in version 2), verifies
that each corrected version-1 statement fails on its counterexample and that each
corrected version-2 statement holds (seven controls, NC1-NC7, including a
symmetry-correlated-noise control and an independent-inversion Weyl check). Each
script reports its status and exits non-zero on failure; some write JSON files,
while the byte-preserved guards described in `VERIFICATION.md` print their
status for explicit redirection. See that file for the check-to-theorem map.
Because the thirteen version 1 guards are kept
byte-identical to the version 1 deposit, the internal label "L1 Krein
boundary-mass identity" inside `check_elevations.py` and its artifact predates
the version 2 downgrade of the Krein-string reading; that legacy label refers
only to the verified reciprocal-increment algebra `1/W_+ - 1/W_- = -a_q`, not to
a Krein-string boundary-mass interpretation, which version 2 records as an
analogy only.

```bash
for s in scripts/check_*.py; do python3 "$s"; done
```

## Integrity and regenerated artifacts

Run the archive-integrity check **before** running the guards:

```bash
shasum -a 256 -c SHA256SUMS
```

Eleven guard scripts write their results directly into `artifacts/`. The three
byte-preserved version 1 scripts `check_coincidence_readout.py`,
`check_dirichlet_readout.py`, and `check_universal_jet.py` print JSON to stdout;
regenerate their archived files with explicit redirection:

```bash
python3 scripts/check_coincidence_readout.py > artifacts/coincidence_readout_audit.json
python3 scripts/check_dirichlet_readout.py > artifacts/dirichlet_readout_audit.json
python3 scripts/check_universal_jet.py > artifacts/universal_jet_audit.json
```

Six long-running guards
(`check_event_jet_determinant.py`, `check_event_jet_largeN.py`,
`check_event_jet_recurrence.py`, `check_event_prony_reconstruction.py`,
`check_source_quotient_and_transport.py`, and `check_spectral_barrier_jump.py`)
also record the measured `runtime_seconds`. That field is intentionally different
on every machine and run, so rerunning those guards makes the corresponding six
artifact hashes differ even when every mathematical result and `"status": "PASS"`
field reproduces exactly. To preserve an archive-integrity check, run the guards
in a disposable copy or check `SHA256SUMS` first and keep regenerated artifacts as
separate run evidence.

## Build the paper

```bash
python3 figures/make_figures.py  # optional; rewrites figures/fig_*.pdf
cp figures/fig_*.pdf source/
cd source && pdflatex main && bibtex main && pdflatex main && pdflatex main
```

The copy step is required because the preserved version 2 TeX source searches
beside `main.tex`; the figure generator writes to the sibling `figures/` directory.

## License

Manuscript and figures: CC BY 4.0 (see
[`LICENSE-PAPER-CC-BY-4.0.txt`](LICENSE-PAPER-CC-BY-4.0.txt)). In this GitHub
repository, the verification scripts and other software files are distributed
under the repository's [MIT License](../../LICENSE). The archived Zenodo v2 bundle
permits the verification scripts to be used freely for reproduction.
