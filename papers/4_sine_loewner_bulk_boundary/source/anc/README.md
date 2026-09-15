# Reproducibility programs for the unified sine Loewner paper

All mathematical statements and proofs are in the paper. This directory supplies exact rational witnesses and reproducibility programs; it is not a mathematical supplement.

Use Python 3.12 and the pinned versions in requirements.txt. Read trace-README.md for sixteen trace/finite checks and boundary-README.md for seven boundary checks. Run the 23 standard checks sequentially, with OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1 and VECLIB_MAXIMUM_THREADS=1. Optional expanded sweeps are outside the standard replay.

The checks comprise eighteen floating-point diagnostics, two exact rational scalar checks and three outward-rounded interval certificate programs. Their distinct evidential scopes are documented beside the commands. Passing finite samples does not prove a universal theorem. The three interval certificates rely on the analytic inequalities and infinite-dimensional reductions proved in the paper.

Two programs record hashes of the compiled manuscript modules: verify_trace_constant.py and verify_parity_and_covariance.py. These anchors point into sections/, not to uncompiled copies of old sources.

No private path or input is required. Choose fresh output names for repeated runs; the exact and boundary output writers refuse overwrite.
