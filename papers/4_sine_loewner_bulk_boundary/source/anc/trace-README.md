# Trace and finite-matrix diagnostics

Run from this directory with Python 3.12 and the versions in requirements.txt. The floating-point checks use NumPy 2.4.4, SciPy 1.17.1 and, when needed, mpmath 1.3.0. The exact scalar checkers use only the Python standard library.

Set OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1 and VECLIB_MAXIMUM_THREADS=1 before running numerical scripts. Run sequentially. These are bounded standard runs; the optional spectral-gap --full sweep is not part of the standard checks. The scripts use POSIX date and Unix resource facilities; native Windows is not supported.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python3 verify_critical_scaling.py
python3 verify_extended_traces.py
python3 verify_general_trace.py
python3 verify_trace_transfer.py
python3 verify_determinant.py
python3 verify_parity_and_covariance.py
python3 verify_boundary_extensions.py
python3 verify_finite_and_signed_bounds.py
python3 verify_supercritical_counts.py
python3 verify_edge_rates.py
python3 verify_second_order_constants.py
python3 verify_quartic_law.py
python3 verify_structural_bounds.py
python3 verify_spectral_gap.py
python3 verify_trace_constant.py --output trace_constant.json
python3 verify_trace_constant_independent.py --output trace_constant_independent.json
```

Each program must exit 0 and report PASS. Preserve its raw output and environment when reporting a failure. Exact scalar outputs refuse overwrite; choose a fresh output filename for a replay. Keep verify_parity_and_covariance.py in this layout: it records the adjacent source hashes, mapped to the actual compiled unified modules.

| Script | Scope |
|---|---|
| `verify_critical_scaling.py` | Finite entries, exact critical sampling, reflection, contraction, kernel-error estimates, exact finite trace, two-sign certificates and small-parameter coefficients. |
| `verify_extended_traces.py` | Every retained higher-even-moment octave fixture, the nonpolynomial tests and positive-prolate comparison, using two Gauss–Legendre resolutions. Also the Schatten coefficients against their closed forms at p=1,2,3,4, their monotonicity and convexity in p, the sharp Lipschitz bound on the coefficient functional, the parity-resolved second moment against its stated deviation bound, and the exact finite reflection identity. |
| `verify_general_trace.py` | Jump-integral coefficients against the independent Mellin integral, rational even-moment coefficients through degree 14, finite star and Hilbert-corner reductions, and continuum traces at two quadrature orders. It independently reconstructs the phase kernel, checks the cyclic fourth-moment identities, extracts Ω₂ and Ω₃, and checks their explicit endpoint-strip bounds. General trace residuals are reported as diagnostics; no finite-sample bound on the theorem's remainder is asserted. |
| `verify_trace_transfer.py` | Sampled versus cell-averaged Fourier factorizations; the uniform nuclear and smooth finite trace bounds with exact linear correction; model-band lower bounds; the trace-norm coefficient; refined continuum trace norms, weighted gap mass and signed-rank observations. The finite samples do not certify asymptotic limits. |
| `verify_determinant.py` | The logarithmic jump integral versus its closed coefficient, the small-parameter coefficient limit, the arcsine form of the squared determinant coefficient against both the arctangent form and a direct integral, refined phase-dependent Fredholm logarithms and independent matrix determinants. |
| `verify_parity_and_covariance.py` | Independent parity convolution and leading oscillation with explicit error, signed ordered eigenvalue monotonicity at two resolutions, scalar Schatten and determinant coefficients, and the actual finite compression losses and eigenvalue-list comparisons behind Lipschitz transfer. |
| `verify_boundary_extensions.py` | Noncommuting covariance comparison, twelve exact finite/continuum pairs, high-precision seam coefficients, three commensurate frequency families, and polynomial/holomorphic parity limits at two quadrature orders. |
| `verify_finite_and_signed_bounds.py` | Direct finite divided differences versus the independent exact Fourier sum; the global minimum defect and refined scalar Lipschitz bounds; seven Lipschitz tests; signed concentration identities and explicit errors; translated-phase kernel symmetries and antinode pairing; the adjacent-band quadratic identity. Two modest quadrature resolutions are compared. It neither estimates an unknown asymptotic constant nor tests the quantified theorem by sampling. |
| `verify_trace_constant.py` | Exact rational checks of the scalar inequalities for the full-fiber gauge bound and admissible trace constant 143; conditional on the analytic operator estimates. |
| `verify_trace_constant_independent.py` | Independent multinomial derivative coefficients, rational Machin enclosure for pi, and safe trace/band/gap constants. |
| `verify_supercritical_counts.py` | Independent modulated-prolate reconstruction, reflection reduction, the concentration comparison; prints the four signed fixed-threshold parity counts on a finite grid and checks the historical concentrated-subspace transfer estimate (a diagnostic retained for comparison; that auxiliary proof is omitted from this edition). |
| `verify_edge_rates.py` | Four signed edge defects versus half-band prolate leakage on sampled finite cases, and the proved upper bound `2 mu (1-mu)/(2 mu - 1)` of the dual comparison at every cell where `mu > 1/2`. The sampled cells have `d*eps <= 8.02`, below the `d*eps >= 25` hypothesis of the factor-3 statement; the ratio window is empirical and does not replace the uniform analytic exponential proof. |
| `verify_second_order_constants.py` | Exact continuum quadratic formulas, defect traces and cross terms at integer and off-integer scales; exact finite identity and its asymptotic law. |
| `verify_quartic_law.py` | The algebraic fourth-moment reduction (an implementation-consistency identity) and the proved logarithm-centered bounds 7.33, 15.676 and 32, with the observed/bound ratios printed; the sampled remainders lie far below these constants, which are not claimed to be sharp. It refines the fourth moment itself, samples off-integer phases, and independently reconstructs the boundary operator from its exterior integral, with an explicit omitted-tail budget. |
| `verify_structural_bounds.py` | SciPy Wasserstein distance versus the finite bound; four Lipschitz test functions; refined continuum absolute traces and sampled extreme-edge monotonicity; trial-packet Fourier convolution and Parseval against Toeplitz losses and independent tail quadrature. Finite observations do not prove infinite inertia, strict monotonicity at every parameter or tiny exponential defects. |
| `verify_spectral_gap.py` | All five census columns and empty guard bands through tau=128 in its default run; the proved (not sharp) window-count bounds and fourth-moment increments. |

Floating-point diagnostics test finite samples and formula normalizations. They are not interval certificates, proofs of asymptotic convergence, or evidence that a universal error constant is optimal. Some combined checks cover identities whose complete proofs are in the unified paper. Their scope is preserved so no previous diagnostic evidence is discarded.

`verify_spectral_gap.py` is also included byte-for-byte as the local spectrum helper imported by `verify_extended_traces.py`. Its standalone census check is included once in the list above; its role as an imported helper is not counted as an additional independent validation.
