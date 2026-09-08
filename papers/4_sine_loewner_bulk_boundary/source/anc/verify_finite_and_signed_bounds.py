#!/usr/bin/env python3
"""Small independent diagnostics for the paper's refined assertions.

Ordinary-double error detection, not interval certification or asymptotic proof.
No external data or ancillary modules are read, and no files are written.
In particular, no numerical constant is asserted for an O_delta remainder.
"""
from __future__ import annotations

import os

for _key in (
    "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_key] = "1"

import math
import sys

import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigvalsh
from scipy.special import roots_legendre, sici


PI = math.pi
NS = (1, 3, 10, 80)
EPSILONS = (0.003, 0.05, 0.2, 0.49)
TAUS = (0.1, 0.7, 2.0, 3.125, 8.375)
ORDERS = (96, 144)
PHASES = (0.0, PI / 7, PI / 2, -PI / 5, 7 * PI / 3)
TESTS = {
    "absolute": (np.abs, 1.0),
    "positive_part": (lambda x: np.maximum(x, 0), 1.0),
    "truncated_absolute": (lambda x: np.minimum(np.abs(x), 0.35), 1.0),
    "shifted_cusp": (lambda x: np.abs(x - 0.37) - 0.37, 1.0),
    "clipped_odd": (lambda x: np.clip(x, -0.2, 0.2), 1.0),
    "square": (lambda x: x * x, 2.0),
    "exponential": (np.expm1, math.e),
}
MAXIMA: dict[str, float] = {}


def close(name, actual, expected, tolerance=2e-10):
    a, b = np.asarray(actual), np.asarray(expected)
    error = float(np.max(np.abs(a - b)))
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise AssertionError(f"{name}: non-finite value")
    MAXIMA[name] = max(MAXIMA.get(name, 0.0), error)
    if error > tolerance:
        raise AssertionError(f"{name}: error {error:.12g} > {tolerance:.12g}")


def bound(name, value, upper, tolerance=2e-10):
    if not math.isfinite(value) or not math.isfinite(upper) or upper <= 0:
        raise AssertionError(f"{name}: invalid value or upper bound")
    MAXIMA[name + " ratio"] = max(MAXIMA.get(name + " ratio", 0.0), value / upper)
    if value > upper + tolerance:
        raise AssertionError(f"{name}: {value:.12g} > {upper:.12g}")


def finite_matrix(n, epsilon):
    """Direct divided differences, not the cosine/sinc factorization."""
    nodes = np.arange(-n, n + 1, dtype=float)
    difference = nodes[:, None] - nodes[None, :]
    values = np.sin(2 * PI * epsilon * nodes)
    numerator = values[:, None] - values[None, :]
    result = numerator / (PI * np.where(difference == 0, 1, difference))
    np.fill_diagonal(result, 2 * epsilon * np.cos(2 * PI * epsilon * nodes))
    return result


def finite_sum(n, epsilon):
    k = np.arange(-2 * n, 2 * n + 1, dtype=float)
    s = np.sin(PI * epsilon * k) ** 2 / (PI * np.where(k == 0, 1, k)) ** 2
    s[k == 0] = epsilon ** 2
    width = 2 * n + 1 - np.abs(k)
    return float(np.sum(s * (2 * width + 2 * np.sin(2 * PI * epsilon * width)
                             / math.sin(2 * PI * epsilon))))


def global_defect_bound(x):
    return min(2 * x, 3 + 2 * max(0, math.log(2 * x)),
               8 / PI ** 2 * (3 + max(0, math.log(PI * x))))


def check_finite():
    for n in NS:
        d = 2 * n + 1
        for epsilon in EPSILONS:
            a = finite_matrix(n, epsilon)
            eigenvalues = eigvalsh(a)
            quadratic = float(np.sum(a * a))
            close("finite exact sum", quadratic, finite_sum(n, epsilon))
            exact_trace = 2 * epsilon * math.sin(PI * d * epsilon) / math.sin(PI * epsilon)
            close("finite diagonal trace", float(np.trace(a)), exact_trace)
            deficit = 2 * d * epsilon - quadratic
            if deficit < -2e-10:
                raise AssertionError("negative finite deficit")
            b = global_defect_bound(d * epsilon)
            bound("global finite deficit", deficit, b)
            bound("finite contraction", float(np.max(np.abs(eigenvalues))), 1.0)
            for function, lip in TESTS.values():
                residual = abs(float(np.sum(function(eigenvalues)))
                               - d * epsilon * float(function(1.0) + function(-1.0)))
                bound("finite scalar Lipschitz", residual,
                      lip * (2 * max(0, deficit) + abs(exact_trace)))
                bound("finite sharpened global Lipschitz", residual,
                      lip * (2 * b + abs(exact_trace)))


def continuum_kernel(nodes, weights, tau, phase=0.0):
    """Independent direct divided differences, with exact diagonal."""
    difference = nodes[:, None] - nodes[None, :]
    values = np.sin(2 * PI * tau * nodes + phase)
    result = (values[:, None] - values[None, :]) / (
        PI * np.where(difference == 0, 1, difference))
    np.fill_diagonal(result, 2 * tau * np.cos(2 * PI * tau * nodes + phase))
    return result * np.sqrt(weights[:, None] * weights[None, :])


def prolate_kernel(nodes, weights, tau):
    difference = nodes[:, None] - nodes[None, :]
    return tau * np.sinc(tau * difference) * np.sqrt(weights[:, None] * weights[None, :])


def exact_prolate_defect(a):
    si, ci = sici(4 * a)
    j = 0.5 * math.cos(4 * a) - 2 * a * (PI / 2 - si)
    return (math.log(4 * a) + 1 + np.euler_gamma - ci - 2 * j) / PI ** 2


def check_continuum():
    phase_bound = 4 * quad(lambda u: math.exp(PI * u * u / 2), 0, 1,
                           epsabs=1e-12, epsrel=1e-12)[0]  # 4 b(1/2)
    previous = {}
    for order in ORDERS:
        nodes, weights = roots_legendre(order)
        for tau in TAUS:
            k = continuum_kernel(nodes, weights, tau)
            c = prolate_kernel(nodes, weights, tau)
            c2 = prolate_kernel(nodes, weights, 2 * tau)
            cosine, sine = np.cos(PI * tau * nodes), np.sin(PI * tau * nodes)
            plus = 2 * cosine[:, None] * c * cosine[None, :]
            minus = 2 * sine[:, None] * c * sine[None, :]
            close("positive split sum", plus + minus, c2)
            close("positive split difference", plus - minus, k)
            dr = float(np.trace(c2) - np.sum(c2 * c2))
            dj = 4 * tau - float(np.sum(k * k))
            cross = float(np.trace(k) - np.sum(k * c2))
            dp = float(np.trace(plus) - np.sum(plus * plus))
            dm = float(np.trace(minus) - np.sum(minus * minus))
            product = float(np.sum(plus * minus))
            close("signed defect plus identity", dp, (dr + dj + 2 * cross) / 4)
            close("signed defect minus identity", dm, (dr + dj - 2 * cross) / 4)
            close("signed product identity", product, (dj - dr) / 4)
            bound("signed cross constant", abs(cross - 4 * math.log(2) / PI ** 2
                                                * math.cos(2 * PI * tau)),
                  3 / (PI ** 3 * tau))
            common = 3 * math.log(tau) + math.log(2 * PI ** 3) + 3 * (1 + np.euler_gamma)
            for value, wave in ((dp, math.cos(PI * tau)), (dm, math.sin(PI * tau))):
                expected = (common + 16 * math.log(2) * wave ** 4) / (4 * PI ** 2)
                bound("signed defect constants", abs(value - expected), 21 / (8 * PI ** 3 * tau))
            expected_product = (math.log(2 * PI * tau) + 1 + np.euler_gamma
                                + 2 * math.log(2) * math.cos(4 * PI * tau)) / (4 * PI ** 2)
            bound("signed product constant", abs(product - expected_product), 9 / (8 * PI ** 3 * tau))

            modulation = np.exp(2j * PI * tau * nodes)
            adjacent = modulation[:, None] * c * modulation.conj()[None, :]
            adjacent_trace = np.sum(c * adjacent.T)
            close("adjacent trace reality", adjacent_trace.imag, 0.0)
            adjacent_trace = float(adjacent_trace.real)
            close("adjacent quadratic identity", adjacent_trace,
                  exact_prolate_defect(PI * tau) - exact_prolate_defect(2 * PI * tau) / 2)
            bound("adjacent quadratic constant", abs(adjacent_trace -
                  (math.log(2 * PI * tau) + 1 + np.euler_gamma) / (2 * PI ** 2)),
                  5 / (4 * PI ** 3 * tau))

            base_eigenvalues = eigvalsh(k)
            observables = [dr, dj, cross, dp, dm, product, adjacent_trace]
            for function, lip in TESTS.values():
                trace = float(np.sum(function(base_eigenvalues)))
                observables.append(trace)
                if tau >= 2:
                    residual = abs(trace - 2 * tau * float(function(1.0) + function(-1.0)))
                    bound("continuum sharpened Lipschitz", residual,
                          lip * (4 / PI ** 2 * math.log(tau) + 2.712))

            for phase in PHASES:
                phased = continuum_kernel(nodes, weights, tau, phase)
                close("phase reflection", phased[::-1, ::-1],
                      continuum_kernel(nodes, weights, tau, -phase))
                close("phase sign reversal", -phased,
                      continuum_kernel(nodes, weights, tau, phase + PI))
                close("phase periodicity", phased,
                      continuum_kernel(nodes, weights, tau, phase + 2 * PI))
                close("phase trace", float(np.trace(phased)),
                      2 / PI * math.sin(2 * PI * tau) * math.cos(phase))
                eigenvalues = eigvalsh(phased)
                for function, lip in TESTS.values():
                    trace = float(np.sum(function(eigenvalues)))
                    observables.append(trace)
                    bound("uniform phase Lipschitz", abs(trace -
                          float(np.sum(function(base_eigenvalues)))), phase_bound * lip)
                if phase == PI / 2:
                    close("antinode eigenvalue pairing", eigenvalues, -eigenvalues[::-1])
                    for odd in (lambda x: x, lambda x: x ** 3,
                                lambda x: np.clip(x, -0.2, 0.2)):
                        close("antinode odd trace", float(np.sum(odd(eigenvalues))), 0.0)
            if tau in previous:
                close("quadrature refinement", observables, previous[tau], tolerance=2e-8)
            previous[tau] = observables


def main():
    print("Ordinary-double diagnostics; no interval or asymptotic certification.")
    print(f"Finite N={NS}; epsilon={EPSILONS}")
    print(f"Continuum tau={TAUS}; Gauss-Legendre orders={ORDERS}; phases={PHASES}")
    print(f"Lipschitz tests={tuple(TESTS)}; maximum matrix dimension={max(2 * max(NS) + 1, max(ORDERS))}")
    check_finite()
    check_continuum()
    for name in sorted(MAXIMA):
        print(f"{name}: {MAXIMA[name]:.12g}")
    print("PASS")


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, ValueError, FloatingPointError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        sys.exit(1)
