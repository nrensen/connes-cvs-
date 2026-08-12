# Errata

Corrections to the accompanying paper, [arXiv:2605.20224](https://arxiv.org/abs/2605.20224), *High-Precision Approximation of Riemann Zeros via the Truncated Weil Form*. This note records them for transparency. The 2026-06-26 correction is incorporated into the current published version (arXiv:2605.20224v2, Zenodo Version 3.3); the 2026-08-12 corrections below are recorded here first and will be folded into the next version of the paper. Entries are newest first. No entry changes a quantitative result.

## 2026-08-12 - the Section 8.2 Paley-Wiener mechanism is withdrawn, and Table 14 is measured at T = 400

**Cutoff provenance.** The caption of Table 14 (`tab:sobolev-scaling`, the Sobolev regularity exponent s(c)) does not state the archimedean integration cutoff. Its rows are measured at **T = 400**, inherited from the per-cutoff N-convergence ladders, not at the T = 800 used for the 15-cutoff sweep. Section 8.2 nevertheless derives A = 55 * 2*pi / 800 = 0.432 as though the slope were a T = 800 measurement. That substitution is incorrect, and the caption should state T = 400.

**Mechanism withdrawn.** Section 8.2 proposes s = sigma_eff * T / (2*pi), which predicts s proportional to T at fixed c, and Section 11 preregisters exactly that test. The test has now been run at c = 23 on the N in {40, 60, 80} grid at dps = 150:

| T | s |
|---:|---|
| 400 | 46.140 |
| 800 | 46.031 |
| 1600 | 45.934 |

Across a fourfold increase in T the exponent moves by 0.206, where proportionality would require roughly a fourfold increase. The Paley-Wiener reading of Section 8.2 is therefore **withdrawn**, and the constant A = 0.432 has no meaning as derived.

The empirical scaling s(c) ~ 55 log c - 128 is a fit to the measured exponents, not a consequence of the withdrawn mechanism, and is unaffected; so are the measured values in Table 14, once correctly labelled T = 400. The three-cutoff test was first carried out by M. Osman ([github.com/Osman209/prime-number-studies](https://github.com/Osman209/prime-number-studies), issue #2) and has since been reproduced independently here.

## 2026-08-12 - matching-digit count at c = 67

Section 6 defines the matching-digit count as the floor of -log10 of the absolute error. Applied to the c = 67, N = 100, dps = 200 datum, whose error is 1.478e-168, that definition gives **167** matching digits; the text prints 168, which labels the error by its decade instead. The measured error value is unchanged and is reported correctly throughout, and every other matching-digit count in the paper follows the stated definition exactly: 307 to 329 at c = 100, N = 250, dps = 500, and 219 to 242 at N = 150, dps = 1000. Only the c = 67 corroborative comparison figure is affected; no quantitative result changes.

## 2026-06-26 - negative-sign eigenvalue blocks are a finite-cutoff artifact

The paper reported small blocks of negative-sign even-sector eigenvalues in two places, and interpreted them as features of the finite-N truncation:

- at c=100 (abstract, §2.4, §6.6, the N-sweep table): a block of dps-stable negative-sign eigenvalues, taken to mean the matrix-level smallest eigenvalue is negative, attributed tentatively to condition-driven sign loss; and
- for L(s, χ₃) at c=23, 29 (§8.10, and the Future Directions section): negative even-sector eigenvalues, interpreted as a character-dependent positivity breakdown.

Both are artifacts of the finite archimedean integration cutoff T, not features of the operator. The negatives are stable under increasing working precision but not under increasing T:

- c=100 (T=800): absent at T=1200, and a cutoff-free evaluation of the archimedean entries leaves the even sector with no negative eigenvalues.
- χ₃ (T=400): re-running the exact original computation with only T varied, the c=23 negative (−6.46×10⁻²³) becomes positive by T=800, and the c=29 negative (−5.82×10⁻¹⁷) by T=1200; both even sectors are then non-negative.

dps-stability was mistaken for correctness; the correct diagnostic of a deep-spectrum value is agreement between two values of T. The "structural character dependence" reading of the χ₃ case, and the suggestion that the c=100 and χ₃ blocks arise from different mechanisms, are withdrawn: both are the same archimedean-truncation artifact.

No quantitative result changes. The reported smallest-positive branch is the genuine smallest eigenvalue, and the recovery of γ₁ through γ₁₀ to 307–329 matching digits at c=100, N=250 stands, as do the Aitken extrapolation and all convergence data.

Professor A. Connes prompted this investigation through his questions about the c=100 spectrum. The cutoff sensitivity was then independently identified by B. W. A. Silva ([zenodo.org/records/20650146](https://zenodo.org/records/20650146)) and is consistent with the naturally even, positive ground state reported by R. Andrews ([zenodo.org/records/20427500](https://zenodo.org/records/20427500)).
