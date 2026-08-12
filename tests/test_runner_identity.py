"""
Bit-identity of the v0.3.0 runner against the classic path.

The runner (``connes_cvs.runner``) must be arithmetically
indistinguishable from the classic

    build_galerkin_matrix(c, N, T, dps) -> compute_ground_state(Q)

pipeline: same kernels, same assembly conventions, same flint precision
convention, lossless worker transport. This test computes a tiny cell
(c=13, N=8, T=60, dps=30) through both paths and asserts exact equality:

- every Galerkin matrix entry, compared both as ``mp.nstr`` strings and
  on the raw mpf bit patterns (``_mpf_``);
- the ground-state eigenvalue;
- every eigenvector component.

Runs in the fast suite (seconds with python-flint; still tractable
without it, since both paths use the same digamma backend either way).
"""

from __future__ import annotations

import pytest


@pytest.mark.timeout(600)
def test_runner_matches_classic_path_entrywise():
    import mpmath as mp
    from connes_cvs.operator import build_galerkin_matrix, compute_ground_state
    from connes_cvs.runner import CellConfig, GalerkinCell

    c, N, T, dps = 13, 8, 60, 30

    Q_classic = build_galerkin_matrix(c=c, N=N, T=T, dps=dps)
    lam_classic, vec_classic = compute_ground_state(Q_classic)

    # Serial runner path (workers=1): no pool, identical arithmetic.
    cell = GalerkinCell(
        CellConfig(c=c, N=N, T=T, dps=dps),
        workers=1,
        progress_callback=False,
    )
    cell.run()

    DIM = 2 * N + 1
    assert cell.Q.rows == Q_classic.rows == DIM

    mp.mp.dps = dps
    mismatches = []
    for i in range(DIM):
        for j in range(DIM):
            a, b = Q_classic[i, j], cell.Q[i, j]
            if mp.nstr(a, dps + 5) != mp.nstr(b, dps + 5) or a._mpf_ != b._mpf_:
                mismatches.append((i, j, mp.nstr(a, 20), mp.nstr(b, 20)))
    assert not mismatches, (
        f"{len(mismatches)} matrix entries differ between classic and "
        f"runner paths; first: {mismatches[0]}"
    )

    assert mp.nstr(lam_classic, dps + 5) == mp.nstr(cell.lambda_even, dps + 5)
    assert lam_classic._mpf_ == cell.lambda_even._mpf_, (
        "lambda_even differs at the bit level between classic and runner paths"
    )

    for i in range(DIM):
        assert vec_classic[i, 0]._mpf_ == cell.eigvec_full[i, 0]._mpf_, (
            f"eigenvector component {i} differs between classic and runner paths"
        )


@pytest.mark.timeout(600)
def test_runner_parallel_matches_serial():
    """The pooled path (workers=2) is bit-identical to the serial path."""
    from connes_cvs.runner import CellConfig, GalerkinCell

    c, N, T, dps = 13, 6, 60, 30

    serial = GalerkinCell(
        CellConfig(c=c, N=N, T=T, dps=dps), workers=1, progress_callback=False
    )
    serial.run()
    pooled = GalerkinCell(
        CellConfig(c=c, N=N, T=T, dps=dps), workers=2, progress_callback=False
    )
    pooled.run()

    assert serial.lambda_even._mpf_ == pooled.lambda_even._mpf_
    DIM = 2 * N + 1
    for i in range(DIM):
        for j in range(DIM):
            assert serial.Q[i, j]._mpf_ == pooled.Q[i, j]._mpf_
