# Verification package (version 2)

Reproducibility checks for *A matrix-valued von Mangoldt measure in the finite
Connes-van Suijlekom path*, version 2 (A. Groskin). The proofs in the paper are
algebraic; these scripts are separate reproducibility checks for signs, indexing,
determinant factors, and finite-field behavior. They are not proof substitutes.

There are fourteen scripts: the thirteen reproducibility guards of version 1
(byte-identical to the version 1 deposit) plus one negative-control guard new in
version 2. Nine use only the Python 3 standard library (no third-party
dependencies); five (`check_canonical_scale`, `check_elevations`,
`check_universal_jet`, `check_coincidence_readout`, and
`check_dirichlet_readout`) use `sympy` and/or `numpy`. Eleven scripts write a
JSON artifact under `artifacts/`; the three byte-preserved version 1 scripts
named below print JSON to stdout and need explicit redirection. Every script
reports a PASS status and exits non-zero on failure. All run in seconds to a
few minutes each.

## Run

```
python3 scripts/check_rank_one_jump.py
python3 scripts/check_uncertainty_ceiling.py
python3 scripts/check_canonical_scale.py
python3 scripts/check_elevations.py
python3 scripts/check_universal_jet.py > artifacts/universal_jet_audit.json
python3 scripts/check_coincidence_readout.py > artifacts/coincidence_readout_audit.json
python3 scripts/check_dirichlet_readout.py > artifacts/dirichlet_readout_audit.json
python3 scripts/check_event_jet_largeN.py
python3 scripts/check_event_jet_determinant.py
python3 scripts/check_event_jet_recurrence.py
python3 scripts/check_event_prony_reconstruction.py
python3 scripts/check_source_quotient_and_transport.py
python3 scripts/check_spectral_barrier_jump.py
python3 scripts/check_negative_controls.py
```

## Integrity order

On an untouched copy, first verify the release bytes:

```bash
shasum -a 256 -c SHA256SUMS
```

Then run the guards. Six long-running guards include the current
`runtime_seconds` in their JSON artifact, so their regenerated artifact hashes
will normally differ from the archived hashes even when all substantive fields
reproduce exactly. Run the guards in a disposable copy if both the original hash
gate and regenerated artifacts are needed. The six scripts are listed explicitly
in the README's “Integrity and regenerated artifacts” section.

## What each check verifies

Theorem numbers refer to the compiled version 2 PDF.

| Check | Manuscript result |
|---|---|
| `check_rank_one_jump.py` | Lemma 2.3 (edge derivative) and Theorem 3.1 (the `-2 Lambda(q)/(sqrt q log q) 11^T` first-derivative jump), by evaluating `A_N` and finite-differencing the assembled prime path for `q in {3,4,5,7,8,9,25}`, `N<=6`. |
| `check_uncertainty_ceiling.py` | Theorem 6.1 (finite vanishing-moment ceiling at the prime edge): even-moment Vandermonde invertible, the centered finite-difference stencil uniquely attains `e=2N` / visibility order `4N+1`, and a family realizes `1,5,...,4N+1`. Exact rational arithmetic, `N<=10`. |
| `check_canonical_scale.py` | Theorem 3.1 at the program's canonical scale: at `N=200` (dimension 401) every prime power `q<=100` has the rank-one first jump `-2 Lambda(q)/(sqrt q log q) 11^T`, floating point. |
| `check_elevations.py` | The divided-difference identity for `A_N`, Proposition 3.3 (second-order event law, `+4 Lambda/(sqrt q (log q)^2) 11^T`, PSD), and the Theorem 7.5 rank-one Weyl-function increment `1/W_+ - 1/W_- = -a_q` (z-independent), evaluated at generic sample points, all inside the exact pointwise domain `W_+(z) W_-(z) != 0` stated in version 2. The guard's internal label "L1 Krein boundary-mass identity" is the version 1 legacy name for this same reciprocal-increment algebra; version 2 records the Krein-string reading as an analogy only. |
| `check_universal_jet.py` | Proposition 4.2 (closed form of the universal jet `B_{r,N}(u0)`), exact symbolic (sympy), entrywise, orders `r<=5`, incl. the `r=1,2` specializations. |
| `check_coincidence_readout.py` | Corollary 3.4 (coincidence-averaged weight readout): exact clean recovery + rank-one certificate, and the `(2N+1)^2` variance reduction of the matched average under an IID model, which satisfies the pairwise-uncorrelated hypothesis that Corollary 3.4 now states explicitly (Monte Carlo, fixed seed). |
| `check_dirichlet_readout.py` | Corollary 8.1 (residue-class readout across the Dirichlet family): character-orthogonality reconstruction of `Lambda(q) 1[q==a mod m]` from the per-character first jumps, cyclic moduli `m<=13`. |
| `check_event_jet_largeN.py` | Theorem 5.2 edge-jet rank and Lemma 5.3 transport rank, exact modular, `N<=200`. |
| `check_event_jet_determinant.py` | Theorem 5.2 determinant `(-1)^{N(N-1)/2} 2^N prod k^6 prod (j^2-i^2)^4` and the Lemma 5.3 `tau`-transport determinant, exact modular over four prime fields, `N<=200`. |
| `check_event_jet_recurrence.py` | Theorem 5.5(iii) recurrence `S prod (S-k^2)^2`, exact modular over four prime fields, `N<=1000`. |
| `check_event_prony_reconstruction.py` | Theorem 5.5 window / blind line / recurrence residues, exact modular, `N<=120` over four fields (`N=160,200` over one). |
| `check_source_quotient_and_transport.py` | Lemma 5.1 source quotient + Lemma 5.3 transport, exact modular, `N<=200`. |
| `check_spectral_barrier_jump.py` | Remark 7.4 (conditional: hypothesizes a positive-definite cell, which is not verified for the actual path at canonical parameters): elementary-symmetric barrier deceleration, exact integer, 50 positive-definite cases. |
| `check_negative_controls.py` | New in version 2. Seven negative controls (NC1-NC7) verifying that each corrected version-1 statement fails on its counterexample and each version-2 statement holds: correlated errors kill the variance reduction while IID and mixed pairwise-uncorrelated models satisfying Corollary 3.4 restore it (NC1); the singular-ratio and directional-residual diagnostics separate on a wrong-direction rank-one matrix (NC2); an eigenspace orthogonal to the all-ones vector gives a removable point of the Theorem 7.1 resolvent, and a repeated eigenvalue gives one grouped spectral-projection residue, both computed from the matrices by exact `Fraction` resolvent solves rather than from asserted spectral weights (NC3, NC4); a Weyl-function zero off both spectra breaks the reciprocal identity exactly where version 2 says it does, with `W_+` computed independently by exact inversion (not from the Sherman-Morrison formula under test) and the zero certified off both spectra by determinant evaluation (NC5); an indefinite velocity matrix shows the Krein-string reading needs Stieltjes/positivity hypotheses beyond the rank-one algebra, hence its downgrade to an analogy (NC6); symmetry-correlated noise breaks the full `(2N+1)^2` reduction while the Remark 3.5(ii) variance `sigma^2 (2d-1)/d^3` holds exactly, in exact `Fraction` arithmetic for `d = 3, 5, 7, 9` (NC7). Exact `Fraction` arithmetic except the seeded Monte Carlo. |

## Figures

`make_figures.py` regenerates the three manuscript figures
(`fig_event_signal.pdf`, `fig_reconstruction.pdf`, `fig_uncertainty.pdf`); it
needs `matplotlib`, `mpmath`, `numpy` (see `requirements.txt`).

## License

Manuscript: CC-BY-4.0 (see `LICENSE-PAPER-CC-BY-4.0.txt`). In this GitHub
repository, the verification scripts and other software files are under the
repository's MIT License (`../../LICENSE`); the archived Zenodo v2 bundle permits
the scripts to be used freely for reproduction.
