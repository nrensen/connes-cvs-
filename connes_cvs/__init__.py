"""
connes-cvs: arbitrary-precision construction and validation of the
Connes-van Suijlekom Galerkin matrix.

This package constructs and diagonalizes the truncated Weil operator Q(c)
from Connes & van Suijlekom (arXiv:2511.23257), whose ground-state eigenvalue
measures proximity to the Riemann Hypothesis.

Basic usage::

    from connes_cvs import build_galerkin_matrix, compute_ground_state, extract_zeros

    Q = build_galerkin_matrix(c=13, N=100, T=400, dps=80)
    lam_min, eigvec = compute_ground_state(Q)
    zeros = extract_zeros(eigvec, c=13, n_zeros=1, dps=80)

Version 0.3.0:

- ``extract_zeros(..., c=13)`` computes L = log(c) internally at full
  precision (avoids the float64 L pitfall); passing a Python float as L
  now emits a UserWarning.
- ``connes_cvs.runner`` - an explicit-precision production cell runner
  (``CellConfig``, ``run_cell``) with parallel psi-cache computation,
  atomic checkpoints, progress reporting, and narrowly scoped identity
  tests against the classic path.
- ``connes_cvs.validation.arb_eigenpair_residual_bound`` - a rigorously
  scoped Arb residual bound for an exact supplied finite symmetric matrix.
"""

__version__ = "0.3.0"

from connes_cvs.operator import build_galerkin_matrix, compute_ground_state, extract_zeros
from connes_cvs.validation import arb_eigenpair_residual_bound

__all__ = [
    "__version__",
    "build_galerkin_matrix",
    "compute_ground_state",
    "extract_zeros",
    "arb_eigenpair_residual_bound",
]
