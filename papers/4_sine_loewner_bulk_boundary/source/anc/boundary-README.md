# Boundary certificates and diagnostics

Run from this directory with Python 3.12 and the versions in requirements.txt.
The scripts are self-contained; no private research repository or separate paper is needed to execute them.

## Three rigorous certificate chains

Set the numerical thread limits once for this shell:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
```

Run sequentially:

    python3 verify_boundary_gap_certificate.py
    python3 verify_negative_window_certificate.py
    python3 verify_boundary_residue_certificate.py

The first check uses the exact frequency witness, full-Hilbert exterior
projection bound, and physical synthesis norm to enclose an eigenvalue in
(0.405,0.428), at 192 and 256 bits. The second uses the exact model trial,
global midpoint/tail error and full residual quadratic to exclude
[-0.4765,-0.3565], at 128 and 192 bits. Combined with the analytic signed
alternation theorem, this proves uniqueness in [0.3565,0.4765].

The third verifies the global residue-row and derivative bounds, the exact
trial functional, and the physical spectral-complement transfer at 192 and
256 bits. Its isolation input is supplied by the first two checks and the
analytic alternation theorem. The three passing programs certify finite
inequalities within the proofs; their dependency on the full operator
arguments is explicit in the manuscript.

No quadrature, floating-point eigensolver or saved approximate eigenvalue
is accepted as a certificate. The optional --output argument writes a new
JSON record and refuses an existing output. The negative-window and residue
scripts also accept --bits followed by the requested precision values.

The negative-window helper parametrix_defect_step_bound.py consumes exact
rational meshes. The witnesses, rational meshes and rational model parameters
are read exactly. Transcendental constants and cell primitives are enclosed
with outward-rounded Arb arithmetic.

| Witness | SHA-256 |
|---|---|
| boundary_gap_witness.json | 61946f85902ba95c0b695fb7c43e8e814b84c081918e7602c145b739dbfea57b |
| step_parametrix_witness.json | 2b5d1dca6fbcdf3ae73a2a8f7bb0f75e28cd3a514395e4ee176e3696e87d2d1c |

The generation-time status EXACT_POINT_WITNESS_UNVERIFIED in
step_parametrix_witness.json describes the saved candidate before verification.
It is retained as provenance, not used as a certificate. The supplied verifier
reads the exact rational coefficients and independently establishes the stated
inequalities; its successful output is PASS. Reproduction requires these
witnesses and the supplied verifiers, not the original witness generator.

## Four floating-point diagnostics

Run these independently of the rigorous checks:

    python3 verify_phase_dependence.py
    python3 verify_parity_limits.py
    BOUNDARY_DIAGNOSTIC_OUTPUT=harmonic_homogenization.json python3 verify_harmonic_homogenization.py
    BOUNDARY_DIAGNOSTIC_OUTPUT=closed_interaction.json python3 verify_closed_interaction.py

- verify_phase_dependence.py reproduces the plotted selections with two
  quadrature families. It compares Gauss–Legendre orders 256/384 with
  Clenshaw–Curtis order 384. Agreement below 1e-10 is a diagnostic, not an
  enclosure or proof of an exact eigenvalue count.
- verify_parity_limits.py compares parity powers and determinant ratios
  at two quadrature orders for three scales. It is the parity-only part
  of verify_boundary_extensions.py; unrelated finite-matrix and
  multifrequency tests are omitted.
- verify_closed_interaction.py independently evaluates the real-space
  and frequency-space scalar integrals at 40/65 decimal digits and checks
  the three-dimensional inverse identities. It writes its output only
  if the destination is new. It does not certify the tail amplitude of
  a physical boundary eigenfunction.

The predicted asymptotic limits and splitting coefficient are proved
analytically. Numerical agreement at finitely many scales does not
establish those limits, their uniformity, or an improved error rate.

The harmonic-homogenization check compares the physical oscillatory kernel with independently integrated polynomial Hilbert models and the logit Fourier integral with its beta expression. It writes a protected JSON output; choose a fresh BOUNDARY_DIAGNOSTIC_OUTPUT for each replay. It is a diagnostic and does not certify a boundary tail amplitude.
