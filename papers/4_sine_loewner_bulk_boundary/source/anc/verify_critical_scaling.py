#!/usr/bin/env python3
"""Deterministic falsifier for the critical sine-Loewner scaling paper.

The checks use ordinary double precision.  They test normalization and finite
examples; they are not interval certificates and do not replace the proofs.
"""

from __future__ import annotations

import numpy as np


def source_matrix(n: int, epsilon: float) -> np.ndarray:
    modes = np.arange(-n, n + 1, dtype=float)
    difference = modes[:, None] - modes[None, :]
    sine = np.sin(2.0 * np.pi * epsilon * modes)
    matrix = np.divide(
        sine[:, None] - sine[None, :],
        np.pi * difference,
        out=np.zeros_like(difference),
        where=difference != 0.0,
    )
    np.fill_diagonal(
        matrix,
        2.0 * epsilon * np.cos(2.0 * np.pi * epsilon * modes),
    )
    return matrix


def critical_kernel(tau: float, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    difference = x - y
    numerator = np.sin(2.0 * np.pi * tau * x) - np.sin(2.0 * np.pi * tau * y)
    kernel = np.divide(
        numerator,
        np.pi * difference,
        out=np.zeros(np.broadcast_shapes(x.shape, y.shape), dtype=float),
        where=difference != 0.0,
    )
    diagonal = 2.0 * tau * np.cos(2.0 * np.pi * tau * x)
    return np.where(difference == 0.0, diagonal, kernel)


def parity_bases(dimension: int) -> tuple[np.ndarray, np.ndarray]:
    if dimension % 2 == 0:
        half = dimension // 2
        even = np.zeros((dimension, half))
        odd = np.zeros((dimension, half))
        for j in range(half):
            even[j, j] = 2.0**-0.5
            even[dimension - 1 - j, j] = 2.0**-0.5
            odd[j, j] = 2.0**-0.5
            odd[dimension - 1 - j, j] = -(2.0**-0.5)
        return even, odd

    n = (dimension - 1) // 2
    even = np.zeros((dimension, n + 1))
    odd = np.zeros((dimension, n))
    even[n, 0] = 1.0
    for j in range(1, n + 1):
        even[n - j, j] = 2.0**-0.5
        even[n + j, j] = 2.0**-0.5
        odd[n - j, j - 1] = 2.0**-0.5
        odd[n + j, j - 1] = -(2.0**-0.5)
    return even, odd


def check_exact_sampling_and_symmetry() -> None:
    for n, tau in ((16, 0.25), (32, 0.75), (48, 1.25)):
        epsilon = tau / n
        matrix = source_matrix(n, epsilon)
        modes = np.arange(-n, n + 1, dtype=float)
        sample = critical_kernel(tau, modes[:, None] / n, modes[None, :] / n) / n
        reflection = np.fliplr(np.eye(2 * n + 1))
        sampling_error = np.linalg.norm(matrix - sample, ord=np.inf)
        reflection_error = np.linalg.norm(
            matrix @ reflection - reflection @ matrix,
            ord=np.inf,
        )
        spectral_norm = np.max(np.abs(np.linalg.eigvalsh(matrix)))
        assert sampling_error < 2.0e-13
        assert reflection_error < 2.0e-13
        assert spectral_norm <= 1.0 + 2.0e-13
        print(
            f"sampling N={n:2d} tau={tau:.2f} "
            f"error={sampling_error:.3e} reflection={reflection_error:.3e} "
            f"norm={spectral_norm:.12f}"
        )


def check_step_kernel_bound() -> None:
    for n, tau in ((24, 0.20), (32, 0.75), (48, 1.25)):
        d = 2 * n + 1
        h = 2.0 / d
        modes = np.arange(-n, n + 1, dtype=float)
        cell_centers = h * modes
        offsets = np.array((-0.499, 0.0, 0.499)) * h
        largest_error = 0.0
        for offset_x in offsets:
            x = (cell_centers + offset_x)[:, None]
            for offset_y in offsets:
                y = (cell_centers + offset_y)[None, :]
                embedded = (d / (2.0 * n)) * critical_kernel(
                    tau,
                    (modes / n)[:, None],
                    (modes / n)[None, :],
                )
                exact = critical_kernel(tau, x, y)
                largest_error = max(largest_error, float(np.max(np.abs(embedded - exact))))
        pointwise_bound = (tau + 4.0 * np.pi * tau**2) / n
        assert largest_error <= pointwise_bound * (1.0 + 2.0e-12)
        print(
            f"step-bound N={n:2d} tau={tau:.2f} "
            f"observed={largest_error:.6e} bound={pointwise_bound:.6e}"
        )


def signed_bounds(tau: float) -> tuple[float, float]:
    witness_width = min(1.0, 1.0 / (4.0 * tau))
    positive_even = 16.0 * witness_width * tau / np.pi**2
    negative_odd_magnitude = (
        128.0 * witness_width**3 * tau**3 / (3.0 * np.pi**2)
    )
    return positive_even, negative_odd_magnitude


def check_trace_and_finite_sign_certificate() -> None:
    for n, tau in ((160, 0.10), (160, 0.25), (160, 0.50)):
        matrix = source_matrix(n, tau / n)
        exact_trace = (2.0 * tau / n) * (
            np.sin((2 * n + 1) * np.pi * tau / n)
            / np.sin(np.pi * tau / n)
        )
        trace_error = abs(float(np.trace(matrix)) - exact_trace)

        even, odd = parity_bases(2 * n + 1)
        even_top = float(np.max(np.linalg.eigvalsh(even.T @ matrix @ even)))
        odd_bottom = float(np.min(np.linalg.eigvalsh(odd.T @ matrix @ odd)))
        positive_bound, negative_bound = signed_bounds(tau)
        finite_error = (2.0 * tau + 8.0 * np.pi * tau**2) / n

        assert trace_error < 2.0e-13
        assert finite_error < min(positive_bound, negative_bound)
        assert even_top >= positive_bound - finite_error - 2.0e-13
        assert odd_bottom <= -negative_bound + finite_error + 2.0e-13
        assert even_top > 0.0
        assert odd_bottom < 0.0
        print(
            f"signed N={n:3d} tau={tau:.2f} trace_error={trace_error:.3e} "
            f"even_top={even_top:.10f} odd_bottom={odd_bottom:.10f}"
        )


def continuum_parity_edges(tau: float, order: int = 192) -> tuple[float, float]:
    nodes, weights = np.polynomial.legendre.leggauss(order)
    weighted = critical_kernel(tau, nodes[:, None], nodes[None, :])
    weighted *= np.sqrt(weights[:, None] * weights[None, :])
    even, odd = parity_bases(order)
    even_values = np.linalg.eigvalsh(even.T @ weighted @ even)
    odd_values = np.linalg.eigvalsh(odd.T @ weighted @ odd)
    return float(np.max(np.abs(even_values))), float(np.min(odd_values))


def check_small_tau_coefficients() -> None:
    expected_even = 16.0 * np.pi**2 / 9.0
    expected_odd = 8.0 * np.pi**2 / 9.0
    even_errors: list[float] = []
    odd_errors: list[float] = []
    for tau in (0.08, 0.05, 0.03):
        even_edge, odd_edge = continuum_parity_edges(tau)
        even_coefficient = (4.0 * tau - even_edge) / tau**3
        odd_coefficient = -odd_edge / tau**3
        even_error = abs(even_coefficient - expected_even)
        odd_error = abs(odd_coefficient - expected_odd)
        even_errors.append(even_error)
        odd_errors.append(odd_error)
        print(
            f"small-tau tau={tau:.2f} even_coeff={even_coefficient:.10f} "
            f"odd_coeff={odd_coefficient:.10f}"
        )
    assert even_errors[-1] < even_errors[0]
    assert odd_errors[-1] < odd_errors[0]
    assert even_errors[-1] < 0.08
    assert odd_errors[-1] < 0.04


def main() -> None:
    check_exact_sampling_and_symmetry()
    check_step_kernel_bound()
    check_trace_and_finite_sign_certificate()
    check_small_tau_coefficients()
    print("CRITICAL SCALING FALSIFIER PASS (ordinary precision; not a certificate)")


if __name__ == "__main__":
    main()
