"""
Unit test suite for the Galerkin matrix caching infrastructure in cell.py.

Verifies:
1. Cache miss generates matrix and stores JSON cache entry.
2. Cache hit loads from disk with matching dimensions and verified structure.
3. Numerical identity: output of get_galerkin_matrix exactly matches build_galerkin_matrix.
4. Sub-dimension reuse: requesting N <= N_cached reuses basis points without quadratures.
5. Incremental evaluation: requesting N > N_cached computes only the missing basis points.
"""

from __future__ import annotations

import time
from pathlib import Path

import mpmath as mp
import pytest

import cell
from connes_cvs import build_galerkin_matrix
from cell import get_galerkin_matrix


@pytest.fixture(autouse=True)
def isolate_cell_cache(tmp_path, monkeypatch):
    """Ensure all test cache reads/writes occur in an isolated temporary directory."""
    isolated_dir = tmp_path / ".cell_cache"
    monkeypatch.setattr(cell, "CELL_CACHE_DIR", isolated_dir)
    return isolated_dir


def test_galerkin_matrix_cache_miss_then_hit(isolate_cell_cache):
    """Verify cache miss generates matrix, and subsequent call is a clean cache hit."""
    mp.mp.dps = 25
    c_val = 13
    N_val = 2
    T_val = 50
    dps_val = 25

    # First call: cache miss
    Q_first, meta_first = get_galerkin_matrix(
        c=c_val,
        N=N_val,
        T=T_val,
        dps=dps_val,
        verbose=False,
    )
    assert not meta_first["cache_hit"]
    assert meta_first["computed_points"] == N_val + 1
    assert meta_first["reused_points"] == 0
    assert Q_first.rows == 2 * N_val + 1
    assert Q_first.cols == 2 * N_val + 1

    # Second call: cache hit
    Q_second, meta_second = get_galerkin_matrix(
        c=c_val,
        N=N_val,
        T=T_val,
        dps=dps_val,
        verbose=False,
    )
    assert meta_second["cache_hit"]
    assert meta_second["total_seconds"] < 0.5

    # Numerical identity between miss and hit
    dim = 2 * N_val + 1
    for i in range(dim):
        for j in range(dim):
            assert Q_first[i, j] == Q_second[i, j]


def test_galerkin_matrix_numerical_identity_with_build_galerkin_matrix(isolate_cell_cache):
    """Verify that cached get_galerkin_matrix produces identical entries to build_galerkin_matrix."""
    mp.mp.dps = 25
    c_val = 13
    N_val = 2
    T_val = 50
    dps_val = 25

    Q_direct = build_galerkin_matrix(c=c_val, N=N_val, T=T_val, dps=dps_val)
    Q_cached, _ = get_galerkin_matrix(c=c_val, N=N_val, T=T_val, dps=dps_val, verbose=False)

    dim = 2 * N_val + 1
    max_diff = mp.mpf(0)
    for i in range(dim):
        for j in range(dim):
            diff = abs(Q_direct[i, j] - Q_cached[i, j])
            if diff > max_diff:
                max_diff = diff

    assert max_diff < mp.mpf("1e-20")


def test_subdimension_reuse(isolate_cell_cache):
    """Verify that an existing N=3 cache entry enables N=1 to be assembled without quadratures."""
    mp.mp.dps = 25
    c_val = 13
    T_val = 50
    dps_val = 25

    # Generate N=3 (computes points 0, 1, 2, 3)
    Q_3, meta_3 = get_galerkin_matrix(
        c=c_val,
        N=3,
        T=T_val,
        dps=dps_val,
        verbose=False,
    )
    assert meta_3["computed_points"] == 4

    # Request N=1 (should reuse points 0, 1 from the N=3 entry)
    Q_1, meta_1 = get_galerkin_matrix(
        c=c_val,
        N=1,
        T=T_val,
        dps=dps_val,
        verbose=False,
    )
    assert meta_1["reused_points"] == 2
    assert meta_1["computed_points"] == 0
    assert Q_1.rows == 3
    assert Q_1.cols == 3

    # Direct check of submatrix extraction
    # Indices in Q_3: row/col 0..6 map to m=-3..3.
    # For N=1, indices map to m=-1..1, which correspond to indices 2, 3, 4 of Q_3.
    for i in range(3):
        for j in range(3):
            assert Q_1[i, j] == Q_3[i + 2, j + 2]


def test_incremental_evaluation(isolate_cell_cache):
    """Verify that requesting higher N reuses available points and evaluates only missing points."""
    mp.mp.dps = 25
    c_val = 13
    T_val = 50
    dps_val = 25

    # Step 1: N=1 computes points 0, 1 (2 points)
    _, meta_1 = get_galerkin_matrix(
        c=c_val,
        N=1,
        T=T_val,
        dps=dps_val,
        verbose=False,
    )
    assert meta_1["computed_points"] == 2
    assert meta_1["reused_points"] == 0

    # Step 2: N=2 should reuse points 0, 1 and only compute point 2 (1 point)
    _, meta_2 = get_galerkin_matrix(
        c=c_val,
        N=2,
        T=T_val,
        dps=dps_val,
        verbose=False,
    )
    assert meta_2["reused_points"] == 2
    assert meta_2["computed_points"] == 1
