#!/usr/bin/env python3
"""Small independent diagnostics for the strengthened structural bounds.

Uses scipy.stats.wasserstein_distance, two continuum quadrature resolutions,
and Fourier-coefficient packet constructions. Ordinary double precision only:
finite checks cannot prove injectivity, infinite signed inertia, global strict
monotonicity, or a uniform asymptotic theorem. No tiny edge defect is certified.
No files are written by this script; capture stdout in the validation log.
"""
from __future__ import annotations

import os

# Set before importing numerical libraries: these modest checks use one thread.
for _name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ[_name] = "1"

import hashlib
import math
from pathlib import Path
import platform
import subprocess

import numpy as np
import scipy
from scipy.special import roots_legendre
from scipy.stats import wasserstein_distance

PI = math.pi
ROUND_TOL = 2e-11
FUNCTIONS = (
    ("abs", np.abs, 1.0),
    ("linear", lambda x: x, 1.0),
    ("cos_zero", lambda x: 1.0 - np.cos(PI * x / 2.0), PI / 2.0),
    ("square", lambda x: x * x, 2.0),
)


def source_matrix(n: int, epsilon: float) -> np.ndarray:
    """Direct divided-difference source, independent of the positive split."""
    modes = np.arange(-n, n + 1, dtype=float)
    difference = modes[:, None] - modes[None, :]
    sine = np.sin(2.0 * PI * epsilon * modes)
    result = np.divide(sine[:, None] - sine[None, :], PI * difference,
                       out=np.zeros_like(difference), where=difference != 0.0)
    np.fill_diagonal(result, 2.0 * epsilon * np.cos(2.0 * PI * epsilon * modes))
    return result


def check_trace_functions(values: np.ndarray, budget: float,
                          published_bound: float) -> float:
    defect = budget - float(values @ values)
    signed_trace = float(np.sum(values))
    tolerance = ROUND_TOL * max(1.0, budget)
    assert defect >= -tolerance, (budget, defect)
    assert float(np.sum(np.abs(values))) <= budget + tolerance
    smallest_slack = math.inf
    for name, function, lipschitz in FUNCTIONS:
        target = 0.5 * budget * float(function(1.0) + function(-1.0))
        error = abs(float(np.sum(function(values))) - target)
        generic = lipschitz * (2.0 * defect + abs(signed_trace))
        assert error <= generic + tolerance, (name, error, generic)
        assert error <= lipschitz * published_bound + tolerance, (
            name, error, lipschitz * published_bound)
        smallest_slack = min(smallest_slack,
                             (lipschitz * published_bound - error) / lipschitz)
    return smallest_slack


def finite_checks() -> None:
    print("FAMILY finite: SciPy W1, nuclear budget, defect, four Lipschitz tests", flush=True)
    cases = 0
    largest_ratio = 0.0
    for n in (3, 15, 64):
        for epsilon in (1e-6, 0.01, 0.1, 0.3, 0.499):
            values = np.linalg.eigvalsh(source_matrix(n, epsilon))
            d = 2 * n + 1
            budget = 2.0 * d * epsilon
            assert np.max(np.abs(values)) <= 1.0 + ROUND_TOL
            defect = budget - float(values @ values)
            defect_bound = 3.0 + 2.0 * math.log(max(1.0, 2.0 * d * epsilon))
            assert defect <= defect_bound + ROUND_TOL * max(1.0, budget)
            trace_bound = 7.0 + 4.0 * math.log(max(1.0, 2.0 * d * epsilon))
            slack = check_trace_functions(values, budget, trace_bound)
            observed = wasserstein_distance(
                values, np.array([-1.0, 0.0, 1.0]),
                u_weights=np.full(d, 1.0 / d),
                v_weights=np.array([epsilon, 1.0 - 2.0 * epsilon, epsilon]))
            bound = min(4.0 * epsilon, trace_bound / d)
            assert observed <= bound + 5e-13, (n, epsilon, observed, bound)
            largest_ratio = max(largest_ratio, observed / bound)
            cases += 1
            print(f"  N={n:2d} eps={epsilon:.6g} W1={observed:.9g} "
                  f"bound={bound:.9g} D={defect:.7g} "
                  f"trace_slack/L={slack:.6g}", flush=True)
    print(f"FINITE PASS: {cases} cases; maximum W1/bound={largest_ratio:.6g}", flush=True)


def continuum_spectra(tau: float, order: int) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = roots_legendre(order)
    # Build the full weighted kernel and then compress by orthonormal parity bases.
    x = nodes[:, None]
    y = nodes[None, :]
    kernel = 2.0 * tau * np.cos(PI * tau * (x + y)) * np.sinc(tau * (x - y))
    matrix = np.sqrt(weights[:, None] * weights[None, :]) * kernel
    half = order // 2
    e = np.zeros((order, half))
    o = np.zeros_like(e)
    index = np.arange(half)
    e[index, index] = e[order - 1 - index, index] = 1.0 / math.sqrt(2.0)
    o[index, index] = 1.0 / math.sqrt(2.0)
    o[order - 1 - index, index] = -1.0 / math.sqrt(2.0)
    return np.linalg.eigvalsh(e.T @ matrix @ e), np.linalg.eigvalsh(o.T @ matrix @ o)


def continuum_checks() -> None:
    print("FAMILY continuum: first absolute moment, signed counts, extreme monotonicity", flush=True)
    print("  Counts below use |lambda| > 1e-9. Finite quadrature cannot prove infinite inertia or injectivity.", flush=True)
    previous = None
    for tau in (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0):
        coarse = continuum_spectra(tau, 128)
        fine = continuum_spectra(tau, 192)
        extrema = np.array([fine[0][-1], fine[0][0], fine[1][-1], fine[1][0]])
        coarse_extrema = np.array([coarse[0][-1], coarse[0][0], coarse[1][-1], coarse[1][0]])
        assert np.max(np.abs(extrema - coarse_extrema)) < 2e-11
        values = np.concatenate(fine)
        values_coarse = np.concatenate(coarse)
        first = float(np.sum(np.abs(values)))
        first_coarse = float(np.sum(np.abs(values_coarse)))
        assert abs(first - first_coarse) < 2e-10
        assert first <= 4.0 * tau + ROUND_TOL * max(1.0, 4.0 * tau)
        exact_trace = 2.0 * math.sin(2.0 * PI * tau) / PI
        assert abs(float(np.sum(values)) - exact_trace) < 2e-11
        counts = tuple((int(np.count_nonzero(v > 1e-9)),
                        int(np.count_nonzero(v < -1e-9))) for v in fine)
        coarse_counts = tuple((int(np.count_nonzero(v > 1e-9)),
                               int(np.count_nonzero(v < -1e-9))) for v in coarse)
        assert counts == coarse_counts, (tau, counts, coarse_counts)
        assert all(positive >= 1 and negative >= 1 for positive, negative in counts)
        assert np.all(extrema[[0, 2]] > 0.0) and np.all(extrema[[1, 3]] < 0.0)
        assert np.max(np.abs(values)) <= 1.0 + ROUND_TOL
        if previous is not None:
            increments = np.array([extrema[0] - previous[0], previous[1] - extrema[1],
                                   extrema[2] - previous[2], previous[3] - extrema[3]])
            assert np.min(increments) > 1e-10, (tau, increments)
        previous = extrema
        if tau >= 2.0:
            bound = 4.0 / PI**2 * math.log(tau) + 2.712
            check_trace_functions(values, 4.0 * tau, bound)
        print(f"  tau={tau:4.2f} tr|K|={first:.10f} <= {4*tau:.3f} "
              f"counts(e+/e-,o+/o-)={counts} "
              f"edges={np.array2string(extrema, precision=9)}", flush=True)
    print("CONTINUUM PASS: sampled inequalities and increasing/decreasing extreme edges at two resolutions", flush=True)


def band_energy(coefficients: np.ndarray, epsilon: float) -> float:
    """Exact Fourier band integrals, evaluated as a Toeplitz quadratic form."""
    modes = np.arange(len(coefficients))
    band = epsilon * np.sinc(epsilon * (modes[:, None] - modes[None, :]))
    return float(np.real(np.vdot(coefficients, band @ coefficients)))


def tail_energy(coefficients: np.ndarray, epsilon: float, order: int) -> float:
    """Independent direct quadrature of the trigonometric polynomial outside the band."""
    nodes, weights = roots_legendre(order)
    lo, hi = epsilon / 2.0, 0.5
    x = (hi + lo) / 2.0 + (hi - lo) * nodes / 2.0
    modes = np.arange(len(coefficients)) - (len(coefficients) - 1) // 2
    evaluated = np.exp(2j * PI * x[:, None] * modes[None, :]) @ coefficients
    # Both packets have even modulus squared, so double the positive half.
    return float((hi - lo) * np.dot(weights, np.abs(evaluated)**2))


def packet_checks() -> None:
    print("FAMILY packets: Fourier convolution, Parseval norm, Toeplitz band loss, independent tail quadrature", flush=True)
    print("  Moderate k=2,3,4 test the explicit packet inequalities; the t>=20 exponential envelope is not a numerical tiny-defect certificate.", flush=True)
    for m in (5, 9, 17):
        epsilon = 2.0 / m
        h = 2 * m + 1
        a = h * epsilon
        q_coefficients = np.full(h, 1.0 / h)
        for k in (2, 3, 4):
            n = m * k + 1
            assert math.ceil(2.0 / epsilon) == m
            assert (n - 1) // m == k
            power = np.array([1.0])
            for _ in range(k):
                power = np.convolve(power, q_coefficients)
            odd = np.convolve(power, np.array([-1.0 / (2j), 0.0, 1.0 / (2j)]))
            delta = 1.0 / (PI * h * math.sqrt(k))
            lower_norms = (4.0 * delta / 3.0, 64.0 * delta**3 / 9.0)
            bounds = (3.0 * PI / 4.0 * a * math.sqrt(k) / (2*k - 1) * a**(-2*k),
                      9.0 * PI**5 / 64.0 * a**3 * k**1.5 / (2*k - 3) * a**(-2*k))
            losses = []
            for parity, coefficients, lower_norm, bound in zip(
                    ("even", "odd"), (power, odd), lower_norms, bounds):
                norm = float(np.real(np.vdot(coefficients, coefficients)))  # Parseval
                assert norm >= lower_norm * (1.0 - 2e-13)
                toeplitz_loss = 1.0 - band_energy(coefficients, epsilon) / norm
                loss = tail_energy(coefficients, epsilon, 256) / norm
                refined = tail_energy(coefficients, epsilon, 384) / norm
                assert abs(loss - refined) < 2e-12
                assert abs(toeplitz_loss - refined) < 2e-11
                assert refined >= 0.0 and refined <= bound + 2e-13, (
                    m, k, parity, refined, bound)
                losses.append((refined, bound))
            print(f"  m={m:2d} k={k} N={n:2d} eps={epsilon:.6g} "
                  f"even_loss={losses[0][0]:.8g} <= {losses[0][1]:.8g} "
                  f"odd_loss={losses[1][0]:.8g} <= {losses[1][1]:.8g}", flush=True)
    print("PACKET PASS: 18 normalized loss bounds and denominator bounds", flush=True)


def main() -> None:
    stamp = subprocess.check_output(["date", "+%Y-%m-%d %H:%M:%S %Z"],
                                    env={**os.environ, "TZ": "Asia/Jerusalem"}, text=True).strip()
    digest = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(f"Started {stamp}; Python {platform.python_version()}; NumPy {np.__version__}; SciPy {scipy.__version__}", flush=True)
    print(f"Source SHA256 {digest}; all numerical thread limits=1", flush=True)
    finite_checks()
    continuum_checks()
    packet_checks()
    print("STRUCTURAL BOUNDS DIAGNOSTIC PASS (ordinary precision; finite observations, not proof certificates)", flush=True)


if __name__ == "__main__":
    main()
