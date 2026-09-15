#!/usr/bin/env python3
"""Small independent diagnostics for the general trace and phase reductions.

Only NumPy/SciPy and the Python standard library are required. No external
datasets or other ancillary modules are read. All numerical work uses one
thread and ordinary precision. Finite checks are not proofs of uniform
remainders, trace-class estimates, or asymptotic counting laws.
"""
from __future__ import annotations

import os

for _thread_key in (
    "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_key] = "1"

from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.linalg import eigvalsh, toeplitz
from scipy.special import digamma, polygamma, roots_legendre


TESTS = {
    "x2": (lambda x: np.asarray(x) ** 2, 2.0),
    "x4": (lambda x: np.asarray(x) ** 4, 4.0),
    "x6": (lambda x: np.asarray(x) ** 6, 6.0),
    "exp_minus_one": (np.expm1, math.e),
    "cosine": (lambda x: 2 * np.sin(np.pi * np.asarray(x) / 4) ** 2, math.pi / 2),
    "abs_cubed_C2": (lambda x: np.abs(x) ** 3, 3.0),
    "positive_cubic_knot_C2": (lambda x: np.maximum(np.asarray(x) - 0.8, 0) ** 3, 0.12),
    "odd_cubic": (lambda x: np.asarray(x) ** 3, 3.0),
}
PHASES = (0.0, math.pi / 7, math.pi / 4, math.pi / 2, -math.pi / 5)
TAUS = (0.5, 2.0, 4.125, 8.0, 8.375)
ORDERS = (96, 144)


def jerusalem_stamp() -> str:
    return subprocess.check_output(
        ["date", "+%Y-%m-%d %H:%M:%S %Z"],
        env={**os.environ, "TZ": "Asia/Jerusalem"}, text=True,
    ).strip()


def close(name, actual, expected, tolerance, **context) -> float:
    actual, expected = np.asarray(actual), np.asarray(expected)
    error = float(np.max(np.abs(actual - expected)))
    if not np.all(np.isfinite(actual)) or not np.all(np.isfinite(expected)) or error > tolerance:
        raise AssertionError({"check": name, "error": error,
                              "tolerance": tolerance, **context})
    return error


def require(name, condition, **context) -> None:
    if not bool(condition):
        raise AssertionError({"check": name, **context})


def strip_bound(length: float) -> float:
    return math.sqrt(2 * length) * quad(
        lambda u: math.exp(math.pi * length * u * u), 0, 1,
        epsabs=2e-13, epsrel=2e-13,
    )[0]


def jump_coefficient(function) -> float:
    knot = math.sqrt(2 * 0.8 ** 2 - 1)
    points = [(1 - knot) / 2, 0.5, (1 + knot) / 2]

    def integrand(t):
        r = math.sqrt(1 - 2 * t * (1 - t))
        return float(function(r) + function(-r) - function(1.0) - function(-1.0)) / (
            2 * math.pi ** 2 * t * (1 - t)
        )

    return quad(integrand, 0, 1, points=points,
                epsabs=2e-11, epsrel=2e-11, limit=120)[0]


def mellin_coefficient(function) -> float:
    def integrand(s):
        r = math.sqrt(1 - 0.5 / math.cosh(math.pi * s) ** 2)
        return float(function(r) + function(-r) - function(1.0) - function(-1.0)) / math.pi

    knot = math.acosh(1 / math.sqrt(2 * (1 - 0.8 ** 2))) / math.pi
    return quad(integrand, -6, 6, points=[-knot, 0, knot],
                epsabs=2e-11, epsrel=2e-11, limit=120)[0]


def check_coefficients():
    rows, coefficients = [], {}
    for name, (function, lipschitz) in TESTS.items():
        direct, mellin = jump_coefficient(function), mellin_coefficient(function)
        # The omitted Mellin tails have this analytic bound. The quadrature
        # error itself remains an ordinary-precision diagnostic, not an enclosure.
        tail_bound = 4 * lipschitz / math.pi ** 2 * math.exp(-12 * math.pi)
        error = close("jump versus Mellin coefficient", direct, mellin,
                      2e-10 + tail_bound, function=name)
        coefficients[name] = direct
        rows.append({"function": name, "kappa": direct,
                     "mellin_difference": error, "mellin_tail_bound": tail_bound})
    moments = []
    for m in range(1, 8):
        exact = sum((Fraction(math.comb(m, j) * (-2) ** j
                              * math.factorial(j - 1) ** 2, math.factorial(2 * j - 1))
                     for j in range(1, m + 1)), Fraction())
        numeric = jump_coefficient(lambda x, m=m: np.asarray(x) ** (2 * m)) * math.pi ** 2
        error = close("beta-integral moment coefficient", numeric, float(exact),
                      2e-10, power=2 * m)
        moments.append({"power": 2 * m, "pi2_kappa_exact": str(exact),
                        "quadrature_difference": error})
    close("odd coefficient", coefficients["odd_cubic"], 0.0, 2e-12)
    return coefficients, rows, moments


def trig_columns(n: int):
    frequencies = np.arange(n)
    a = math.pi / 2

    def exponential_integral(beta):
        return np.expm1(1j * beta) / (1j * beta)

    plus = exponential_integral(a - 2 * math.pi * frequencies)
    minus = exponential_integral(-a - 2 * math.pi * frequencies)
    return (plus + minus) / 2, (plus - minus) / (2j)


def hilbert_square(n: int) -> np.ndarray:
    j, k = np.indices((n, n))
    result = np.divide(digamma(k + 1) - digamma(j + 1), k - j,
                       out=np.zeros((n, n)), where=k != j)
    np.fill_diagonal(result, polygamma(1, np.arange(1, n + 1)))
    return result


def check_finite_models():
    rows = []
    for n in (4, 8, 16, 32):
        c_column, s_column = trig_columns(n)
        fourier_errors = []
        for order in (128, 192):
            theta, weights = roots_legendre(order)
            theta, weights = (theta + 1) / 2, weights / 2
            exponentials = np.exp(-2j * math.pi * np.outer(np.arange(n), theta))
            c_quad = exponentials @ (weights * np.cos(math.pi * theta / 2))
            s_quad = exponentials @ (weights * np.sin(math.pi * theta / 2))
            fourier_errors.append(max(
                close("cosine Fourier coefficient", c_column, c_quad, 2e-12, n=n, order=order),
                close("sine Fourier coefficient", s_column, s_quad, 2e-12, n=n, order=order),
            ))
        c, s = toeplitz(c_column), toeplitz(s_column)
        core = c @ c + s @ s
        zero = np.zeros((n, n))
        star = np.block([[zero, c, s], [c, zero, zero], [s, zero, zero]])
        core_values, star_values = eigvalsh(core), eigvalsh(star)
        require("positive star core", core_values.min() >= -2e-12 and core_values.max() <= 1 + 2e-12, n=n)
        singular = np.sqrt(np.maximum(core_values, 0))
        star_error = 0.0
        for name, (function, _) in TESTS.items():
            star_error = max(star_error, close(
                "star trace versus scalar core", np.sum(function(star_values)),
                np.sum(function(singular) + function(-singular)), 3e-11,
                n=n, function=name,
            ))
        indices = np.arange(n)
        v_column, v2_column = np.zeros(n, complex), np.zeros(n)
        v_column[1:] = 1 / (math.sqrt(2) * math.pi * 1j * indices[1:])
        v2_column[0] = 1 / 6
        v2_column[1:] = 1 / (math.pi ** 2 * indices[1:] ** 2)
        v = toeplitz(v_column)
        defect_fourier = toeplitz(v2_column) - v @ v
        h2 = hilbert_square(n)
        defect_hilbert = (h2 + h2[::-1, ::-1]) / (2 * math.pi ** 2)
        defect_error = close("sawtooth defect versus two Hilbert corners",
                             defect_fourier, defect_hilbert, 3e-12, n=n)
        values = eigvalsh(defect_hilbert)
        require("half-unit defect bound", values.min() >= -2e-12 and values.max() <= 0.5 + 2e-12, n=n)
        comparison = float(np.sum(np.abs(eigvalsh(core - np.eye(n) + defect_hilbert))))
        rows.append({"n": n, "fourier_max_differences": fourier_errors,
                     "star_trace_max_difference": star_error,
                     "hilbert_defect_max_difference": defect_error,
                     "defect_min_max": [float(values.min()), float(values.max())],
                     "observed_core_trace_norm_difference": comparison})
    return rows


def continuum_case(tau: float, order: int, coefficients):
    nodes, weights = roots_legendre(order)
    difference = nodes[:, None] - nodes[None, :]
    total = nodes[:, None] + nodes[None, :]
    gram = np.sqrt(weights[:, None] * weights[None, :])
    c = tau * np.sinc(tau * difference) * gram
    modulation = np.exp(2j * math.pi * tau * nodes)
    cm = c * modulation[None, :]
    c2, cm2 = c @ c, cm @ cm
    c3 = c2 @ c
    cprime = modulation[:, None] * c * modulation.conj()[None, :]
    cprime2 = cprime @ cprime
    omega2_complex = np.trace((c3 * modulation[None, :]) @ cm)
    omega3_complex = np.trace(cm2 @ cm2)
    x_complex = np.trace(c2 @ cprime2)
    close("real cyclic traces", [omega2_complex.imag, omega3_complex.imag, x_complex.imag],
          np.zeros(3), 3e-11, tau=tau, order=order)
    omega2, omega3, xtrace = map(float, (omega2_complex.real, omega3_complex.real, x_complex.real))
    c4 = float(np.trace(c2 @ c2))
    v = np.exp(1j * math.pi * tau * total) * c
    fourth = []
    kernel_error = cyclic_error = 0.0
    base_values = None
    for phase in PHASES:
        direct = 2 * tau * np.cos(math.pi * tau * total + phase) * np.sinc(tau * difference) * gram
        reconstructed = np.exp(1j * phase) * v + np.exp(-1j * phase) * v.conj().T
        kernel_error = max(kernel_error, close("modulated factorization", direct, reconstructed,
                                               3e-11, tau=tau, order=order, phase=phase))
        shifted = tau * nodes + phase / (2 * math.pi)
        physical_difference = shifted[:, None] - shifted[None, :]
        sine = np.sin(2 * math.pi * shifted)
        kernel = np.divide(sine[:, None] - sine[None, :], math.pi * physical_difference,
                           out=np.zeros_like(physical_difference), where=physical_difference != 0)
        np.fill_diagonal(kernel, 2 * np.cos(2 * math.pi * shifted))
        kernel_error = max(kernel_error, close("shifted divided-difference kernel", direct,
                                               tau * gram * kernel, 3e-11,
                                               tau=tau, order=order, phase=phase))
        # Nonzero phase does not preserve the original reflection sectors.
        values = eigvalsh(direct)
        if phase == 0:
            base_values = values
        value = float(np.sum(values ** 4))
        prediction = 2 * omega3 * math.cos(4 * phase) + 8 * omega2 * math.cos(2 * phase) + 4 * xtrace + 2 * c4
        cyclic_error = max(cyclic_error, close("phase cyclic expansion", value, prediction,
                                               3e-10, tau=tau, order=order, phase=phase))
        fourth.append(value)
    extracted2 = (fourth[0] - fourth[3]) / 16
    extracted3 = (fourth[0] + fourth[3] - 2 * fourth[2]) / 8
    close("Omega2 phase extraction", omega2, extracted2, 5e-11, tau=tau, order=order)
    close("Omega3 phase extraction", omega3, extracted3, 5e-11, tau=tau, order=order)
    bound2, bound3 = strip_bound(0.25), 4 * strip_bound(0.125)
    require("Omega2 strip bound", abs(omega2) <= bound2 + 3e-10, tau=tau, order=order)
    require("Omega3 strip bound", abs(omega3) <= bound3 + 3e-10, tau=tau, order=order)
    phase_ratio = 0.0
    for j in range(len(PHASES)):
        for k in range(j):
            budget = 16 * strip_bound(abs(PHASES[j] - PHASES[k]) / (2 * math.pi))
            observed = abs(fourth[j] - fourth[k])
            require("fourth-trace endpoint-strip bound", observed <= budget + 3e-10,
                    tau=tau, order=order, phases=[PHASES[j], PHASES[k]])
            phase_ratio = max(phase_ratio, observed / budget)
    traces, residuals = {}, {}
    for name, (function, _) in TESTS.items():
        traces[name] = float(np.sum(function(base_values)))
        residuals[name] = float(traces[name] - 2 * tau * (function(1.0) + function(-1.0))
                               - coefficients[name] * math.log(tau))
    return {"tau": tau, "order": order, "kernel_max_difference": kernel_error,
            "cyclic_max_difference": cyclic_error, "omega2": omega2, "omega3": omega3,
            "xtrace": xtrace, "c4": c4, "phase_fourth_traces": fourth,
            "phase_bound_max_ratio": phase_ratio, "traces": traces, "residuals": residuals}


def run():
    coefficients, coefficient_rows, moments = check_coefficients()
    models = check_finite_models()
    continuum = []
    for tau in TAUS:
        print(f"General trace diagnostic: tau={tau:g}, orders={ORDERS}", file=sys.stderr, flush=True)
        coarse, fine = [continuum_case(tau, order, coefficients) for order in ORDERS]
        refinement = max(
            close("continuum refinement", coarse[key], fine[key], 1e-8, tau=tau, quantity=key)
            for key in ("omega2", "omega3", "xtrace", "c4", "phase_fourth_traces")
        )
        for name in TESTS:
            refinement = max(refinement, close("C2 trace refinement", coarse["traces"][name],
                                              fine["traces"][name], 1e-8, tau=tau, function=name))
        continuum.append({"tau": tau, "orders": list(ORDERS),
                          "refinement_max_difference": refinement,
                          "kernel_max_difference": max(coarse["kernel_max_difference"], fine["kernel_max_difference"]),
                          "cyclic_max_difference": max(coarse["cyclic_max_difference"], fine["cyclic_max_difference"]),
                          "omega2": fine["omega2"], "omega3": fine["omega3"],
                          "phase_fourth_traces": fine["phase_fourth_traces"],
                          "phase_bound_max_ratio": fine["phase_bound_max_ratio"],
                          "observed_trace_residuals": fine["residuals"]})
    return {"coefficient_checks": coefficient_rows, "moment_checks": moments,
            "finite_models": models, "phases": list(PHASES),
            "omega2_bound": strip_bound(0.25), "omega3_bound": 4 * strip_bound(0.125),
            "continuum_checks": continuum}


if __name__ == "__main__":
    timer = time.monotonic()
    metadata = {"started": jerusalem_stamp(), "numerical_threads": 1,
                "python": platform.python_version(), "numpy": np.__version__,
                "scipy": scipy.__version__,
                "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "precision": "ordinary floating point; no interval enclosures",
                "limitations": "Finite identity/refinement diagnostics, not proofs of uniform trace remainders or counting asymptotics. Reported residuals and trace-norm differences are observations, not inferred global bounds."}
    try:
        report = {"status": "PASS", **metadata, **run()}
        exit_code = 0
    except Exception as error:
        report = {"status": "FAIL", **metadata, "error_type": type(error).__name__, "error": str(error)}
        exit_code = 1
    report.update({"completed": jerusalem_stamp(), "elapsed_seconds": time.monotonic() - timer})
    print(json.dumps(report, indent=2, allow_nan=False))
    sys.exit(exit_code)
