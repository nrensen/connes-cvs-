#!/usr/bin/env python3
"""Deterministic normalization and finite-example replay for the sine Loewner paper.

This script builds the exact finite matrix from equation (1.1), independently
reconstructs it from its modulated prolate factorization, checks
reflection reduction, the contraction bound, and the comparison A^2 <= C, then
prints all four fixed-threshold counts.  It uses ordinary double precision and
is a falsifier/reproduction aid, not a proof or an interval certificate.
"""

from __future__ import annotations

import numpy as np


def parity_bases(n: int) -> tuple[np.ndarray, np.ndarray]:
    d = 2 * n + 1
    even = np.zeros((d, n + 1))
    odd = np.zeros((d, n))
    even[n, 0] = 1.0
    for j in range(1, n + 1):
        even[n - j, j] = 2.0**-0.5
        even[n + j, j] = 2.0**-0.5
        odd[n + j, j - 1] = 2.0**-0.5
        odd[n - j, j - 1] = -(2.0**-0.5)
    return even, odd


def source_matrix(n: int, insertion_width: float) -> np.ndarray:
    modes = np.arange(-n, n + 1, dtype=float)
    difference = modes[:, None] - modes[None, :]
    sine = np.sin(2.0 * np.pi * modes * insertion_width)
    matrix = np.divide(
        sine[:, None] - sine[None, :],
        np.pi * difference,
        out=np.zeros_like(difference),
        where=difference != 0.0,
    )
    np.fill_diagonal(
        matrix,
        2.0 * insertion_width * np.cos(2.0 * np.pi * modes * insertion_width),
    )
    return matrix


def centered_concentration(n: int, interval_length: float) -> np.ndarray:
    modes = np.arange(-n, n + 1, dtype=float)
    difference = modes[:, None] - modes[None, :]
    matrix = np.divide(
        np.sin(np.pi * interval_length * difference),
        np.pi * difference,
        out=np.zeros_like(difference),
        where=difference != 0.0,
    )
    np.fill_diagonal(matrix, interval_length)
    return matrix


def counts(values: np.ndarray, threshold: float) -> tuple[int, int]:
    return int(np.sum(values > threshold)), int(np.sum(values < -threshold))


def modulated_prolate_reconstruction(
    n: int, insertion_width: float
) -> np.ndarray:
    modes = np.arange(-n, n + 1, dtype=float)
    prolate = centered_concentration(n, insertion_width)
    phase = np.exp(1j * np.pi * modes * insertion_width)
    modulated = phase[:, None] * prolate * phase[None, :]
    return 2.0 * np.real(modulated)


def run_case(n: int, insertion_width: float, threshold: float) -> None:
    d = 2 * n + 1
    source = source_matrix(n, insertion_width)
    reconstructed = modulated_prolate_reconstruction(n, insertion_width)
    reflection = np.fliplr(np.eye(d))
    even, odd = parity_bases(n)
    even_values = np.linalg.eigvalsh(even.T @ source @ even)
    odd_values = np.linalg.eigvalsh(odd.T @ source @ odd)
    all_values = np.linalg.eigvalsh(source)

    endpoint_concentration = centered_concentration(n, 2.0 * insertion_width)
    comparison_floor = np.linalg.eigvalsh(endpoint_concentration - source @ source)[0]
    factorization_error = np.linalg.norm(source - reconstructed, ord=2)
    reflection_error = np.linalg.norm(source @ reflection - reflection @ source, ord=2)

    # Transfer lemma budget: rho = (1-theta)/16 <= 1/32 gives d_rho = sqrt(1-rho) - sqrt(2 rho) >= 1/2,
    # hence 4 rho / d_rho^2 <= 16 rho = 1 - theta, as the lemma on concentrated subspaces requires.
    leakage_budget = (1.0 - threshold) / 16.0
    base = centered_concentration(n, insertion_width)
    base_even = np.linalg.eigvalsh(even.T @ base @ even)
    base_odd = np.linalg.eigvalsh(odd.T @ base @ odd)
    base_counts = (
        int(np.sum(base_even > 1.0 - leakage_budget)),
        int(np.sum(base_odd > 1.0 - leakage_budget)),
    )

    observed = (*counts(even_values, threshold), *counts(odd_values, threshold))
    target = d * insertion_width / 2.0
    print(
        f"N={n:3d} eps={insertion_width:.6f} d*eps/2={target:.6f} "
        f"counts(e+,e-,o+,o-)={observed} base(e,o)={base_counts} "
        f"||A||={np.max(np.abs(all_values)):.12f} "
        f"factor_error={factorization_error:.3e} "
        f"reflection_error={reflection_error:.3e} "
        f"lambda_min(C_endpoint-A^2)={comparison_floor:.3e}"
    )

    # The base even space transfers to (even,+) and (odd,-); the base odd
    # space transfers to (odd,+) and (even,-).
    assert observed[0] >= base_counts[0]
    assert observed[3] >= base_counts[0]
    assert observed[2] >= base_counts[1]
    assert observed[1] >= base_counts[1]
    assert factorization_error < 5.0e-12
    assert reflection_error < 5.0e-12
    assert np.max(np.abs(all_values)) <= 1.0 + 5.0e-12
    assert comparison_floor > -5.0e-12


def main() -> None:
    threshold = 0.8
    for n, insertion_width in (
        (100, 0.05),
        (100, 0.10),
        (100, 0.20),
        (200, 0.05),
        (200, 0.10),
    ):
        run_case(n, insertion_width, threshold)
    print("DIAGNOSTIC PASS (not an interval certificate)")


if __name__ == "__main__":
    main()
