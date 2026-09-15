#!/usr/bin/env python3
"""Independent floating-point diagnostics for sine Loewner parity, coefficients and covariance transfer.

Requires only NumPy, SciPy, and the standard library.  Run from any directory:
  python verify_parity_and_covariance.py --source-dir .. --output checks.json

The SHA-256 manifest identifies the manuscript snapshot inspected. Equations
are implemented independently here; no manuscript diagnostic is imported.
Two Gauss-Legendre resolutions check convergence. These are error-detection
samples, not interval certificates or proofs of universal assertions.

Source anchors in the compiled unified manuscript:
  sections/07_quadratic.tex: reflected continuum and finite second moments.
  sections/02_exchange.tex: strict signed/parity monotonicity.
  sections/04_finite_trace.tex: both finite covariance/compression comparisons.
  sections/05_density.tex: all Schatten powers.
  sections/06_coefficients.tex: squared Fredholm determinant coefficient.

The global gauge and trace constants now have proved numerical bounds. Here,
the finite tests check the sharper local compression bound using the actual
computed central mass, and check min-max domination and actual trace-norm
loss separately. The analytic manuscript supplies the universal C_0 bound.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess

# Avoid competing with other research work; set before importing BLAS users.
for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
import scipy
from scipy.linalg import eigvalsh
from scipy.special import roots_legendre, sici


PI = np.pi
GAP = 1.0 / np.sqrt(2.0)
EULER = np.euler_gamma
SOURCES = {
    "main.tex": (),
    "sections/07_quadratic.tex": ("eq:parity-G", "eq:parity-finite"),
    "sections/02_exchange.tex": ("eq:monotonicity",),
    "sections/04_finite_trace.tex": ("v7:eq:sample-compression", "v10:eq:finite-lip-error", "v10:eq:finite-lipschitz"),
    "sections/05_density.tex": ("cor:schatten",),
    "sections/06_coefficients.tex": ("eq:det-square",),
}


def check(condition, message):
    if not bool(condition):
        raise AssertionError(message)


def close(actual, expected, tolerance, message):
    error = float(np.max(np.abs(np.asarray(actual) - np.asarray(expected))))
    check(error <= tolerance, f"{message}: error={error:.6g}, limit={tolerance:.6g}")
    return error


def quad_unit(order):
    x, w = roots_legendre(order)
    return (x + 1.0) / 2.0, w / 2.0


def panel_integral(fun, upper, order):
    """Independent direct convolution integration; panels at most pi long."""
    panels = max(1, int(np.ceil(upper / PI)))
    edges = np.linspace(0.0, upper, panels + 1)
    x, w = quad_unit(order)
    widths = np.diff(edges)
    points = edges[:-1, None] + widths[:, None] * x[None, :]
    return np.sum(widths[:, None] * w[None, :] * fun(points))


def kernel(tau, x, y):
    return 2.0 * tau * np.cos(PI * tau * (x + y)) * np.sinc(tau * (x - y))


def continuum(tau, order):
    x, w = roots_legendre(order)
    matrix = kernel(tau, x[:, None], x[None, :]) * np.sqrt(w[:, None] * w[None, :])
    return x, w, matrix, eigvalsh(matrix)


def coefficient_checks(order):
    t, weight = quad_unit(order)
    u = t * (1.0 - t)
    ps = np.array([1.0, 1.5, 2.0, 2.5, 3.0])
    kappas = np.array([
        np.dot(weight, np.expm1((p / 2.0) * np.log1p(-2.0 * u)) / u) / PI**2
        for p in ps
    ])
    v = np.log1p(np.sqrt(2.0))
    closed = [
        -2.0 / PI**2 * (np.sqrt(2.0) * v - np.log(2.0)),
        -2.0 / PI**2 * (0.5 - np.log(2.0) + 5.0 * v / (2.0 * np.sqrt(2.0))),
    ]
    power_error = close(kappas[[0, -1]], closed, 3e-13, "closed p=1,3 coefficients")
    check(np.all(np.diff(kappas) < 0), "sampled kappa_p decrease")
    check(np.all(np.diff(kappas, n=2) > 0), "sampled kappa_p convexity")
    check(np.all(np.abs(kappas) <= abs(closed[0]) * ps + 3e-13), "sharp Lipschitz coefficient bound")
    det = []
    for z in [0.0, -0.5, 0.25, 1.0, 2.0 + 1.0j, -2.0 + 3.0j]:
        z = complex(z)
        integral = np.dot(weight, np.log1p(-2.0 * z * u / (1.0 + z)) / u) / PI**2
        formula = -4.0 / PI**2 * np.arcsin(0.5 * np.sqrt(2.0 * z / (1.0 + z)))**2
        err = close(integral, formula, 3e-13, f"squared determinant coefficient w={z}")
        if z == 1:
            close(integral, -1.0 / 9.0, 3e-13, "determinant coefficient at w=1")
        det.append({"w": [z.real, z.imag], "error": err})
    return {"p": ps.tolist(), "kappa": kappas.tolist(), "closed_max_error": power_error,
            "determinant_checks": det}


def parity_checks(order, panel_order):
    rows = []
    for tau in [1.0, 1.125, 2.375, 4.125, 8.375]:
        x, w, matrix, eigenvalues = continuum(tau, order)
        reflected = float(np.sum(matrix * matrix[::-1, :]))
        upper = 4.0 * PI * tau
        conv = float(panel_integral(lambda s: np.sinc(s / PI) * sici(upper - s)[0], upper, panel_order))
        conv_prime = float(panel_integral(
            lambda s: np.sinc(s / PI) * np.sinc((upper - s) / PI), upper, panel_order))
        si, ci = sici(2.0 * upper)
        prime_formula = (np.cos(upper) * (ci - EULER - np.log(2.0 * upper))
                         + np.sin(upper) * si) / upper
        prime_error = close(conv_prime, prime_formula, 2e-11, "G' convolution versus Ci/Si")
        parity_error = close(reflected, 2.0 * conv / PI**2, 3e-10, "reflected 2D trace versus G convolution")
        positive = x > 0
        xp, wp = x[positive], w[positive]
        parity_moments = []
        for parity in (1, -1):
            block = (kernel(tau, xp[:, None], xp[None, :])
                     + parity * kernel(tau, xp[:, None], -xp[None, :]))
            block *= np.sqrt(wp[:, None] * wp[None, :])
            moment = float(np.sum(block * block))
            close(moment, 0.5 * (np.dot(eigenvalues, eigenvalues) + parity * reflected),
                  2e-10, "parity-sector second moment")
            parity_moments.append(moment)
        old_bound = (4.0 * np.log(2.0 * PI * tau) + 2.0 + PI) / (PI**3 * tau)
        check(abs(reflected - 0.5) <= old_bound + 3e-10, "explicit parity bound")
        rows.append({"tau": tau, "reflected_trace": reflected, "convolution": conv,
                     "identity_error": parity_error, "derivative_identity_error": prime_error,
                     "sector_second_moments_even_odd": parity_moments})
    asymptotic = []
    for tau in [1.0, 1.125, 2.125, 4.125, 8.125, 16.125, 32.125]:
        upper = 4.0 * PI * tau
        conv = float(panel_integral(lambda s: np.sinc(s / PI) * sici(upper - s)[0], upper, panel_order))
        reflected = 2.0 * conv / PI**2
        leading = (0.5 - (np.log(8.0 * PI * tau) + EULER) * np.sin(upper) / (2.0 * PI**3 * tau)
                   - np.cos(upper) / (4.0 * PI**2 * tau))
        bound = (2.0 * np.log(8.0 * PI * tau) + 2.0 * EULER + PI - 2.5) / (8.0 * PI**4 * tau**2)
        remainder = abs(reflected - leading)
        check(remainder <= bound + 3e-11, "explicit parity oscillation remainder")
        scale = -(np.log(8.0 * PI * tau) + EULER) / (2.0 * PI**3 * tau)
        asymptotic.append({"tau": tau, "reflected_trace": reflected,
                           "error_after_leading_oscillation": remainder, "explicit_bound": bound,
                           "error_to_bound_ratio": remainder / bound,
                           "sharp_subsequence_ratio": (reflected - 0.5) / scale if tau != 1 else None})
    return {"exact_identities": rows, "asymptotic_samples": asymptotic}


def monotonicity_checks(order):
    x, w = roots_legendre(order)
    positive = x > 0
    x, w = x[positive], w[positive]
    scales = [2.0, 2.25, 2.5, 2.75]
    result = []
    for parity in (1, -1):
        signed = {1: [], -1: []}
        for tau in scales:
            block = (kernel(tau, x[:, None], x[None, :])
                     + parity * kernel(tau, x[:, None], -x[None, :]))
            spectrum = eigvalsh(block * np.sqrt(w[:, None] * w[None, :]))
            signed[1].append(spectrum[-4:][::-1])
            signed[-1].append(-spectrum[:4])
        for sign in (1, -1):
            values = np.array(signed[sign])
            check(np.min(values) > 1e-10, "sampled monotonicity modes resolved above roundoff")
            minimum_step = float(np.min(np.diff(values, axis=0)))
            check(minimum_step > 1e-10, "strict sampled all-rank monotonicity")
            check(np.max(values) < 1.0 - 1e-10, "sampled strict contractivity")
            result.append({"parity": parity, "sign": sign, "ranks": [1, 2, 3, 4],
                           "tau": scales, "eigenvalue_magnitudes": values.tolist(),
                           "minimum_strict_increment": minimum_step})
    return result


def finite_matrix(d, epsilon):
    m = np.arange(-(d // 2), d // 2 + 1, dtype=float)
    diff = m[:, None] - m[None, :]
    sums = m[:, None] + m[None, :]
    matrix = 2.0 * epsilon * np.cos(PI * epsilon * sums) * np.sinc(epsilon * diff)
    # An independent literal divided-difference construction checks the diagonal.
    sine = np.sin(2.0 * PI * epsilon * m)
    direct = np.empty_like(matrix)
    np.divide(sine[:, None] - sine[None, :], PI * diff, out=direct, where=diff != 0)
    np.fill_diagonal(direct, 2.0 * epsilon * np.cos(2.0 * PI * epsilon * m))
    close(matrix, direct, 4e-14, "finite kernel constructions")
    return m, matrix


def finite_checks(order):
    rows = []
    for tau, dimensions in [(2.125, [33, 65]), (4.25, [65, 129]), (6.125, [129, 257])]:
        _, _, _, spectrum_k = continuum(tau, order)
        absolute_k = np.abs(spectrum_k)
        for d in dimensions:
            h, epsilon = 2.0 / d, 2.0 * tau / d
            m, matrix_a = finite_matrix(d, epsilon)
            spectrum_a = eigvalsh(matrix_a)
            x = h * m
            node, weight = quad_unit(order)
            frequency, frequency_weight = tau * node, tau * weight
            e_hi = np.sqrt(h) * np.exp(2j * PI * frequency[:, None] * x[None, :])
            e_lo = np.sqrt(h) * np.exp(2j * PI * (frequency[:, None] - tau) * x[None, :])
            sampled = e_hi.conj().T @ (frequency_weight[:, None] * e_lo)
            sampled += sampled.conj().T
            sample_identity = close(sampled, matrix_a, 2e-12, "frequency exchange reconstructs sampled A")
            e_hi *= np.sinc(h * frequency)[:, None]
            e_lo *= np.sinc(h * (frequency - tau))[:, None]
            t_matrix = e_hi.conj().T @ (frequency_weight[:, None] * e_lo)
            t_matrix += t_matrix.conj().T
            close(t_matrix.imag, 0.0, 2e-12, "midpoint compression is real")
            t_matrix = t_matrix.real
            spectrum_t = eigvalsh(t_matrix)
            sample_error = float(np.sum(np.abs(eigvalsh(matrix_a - t_matrix))))
            eta = float(4.0 * np.dot(frequency_weight, 1.0 - np.sinc(h * frequency)**2))
            eta_bar = 2.0 * PI**2 / 9.0 * d * epsilon**3
            check(0 <= eta <= eta_bar + 3e-12, "sampling leakage eta bound")
            check(sample_error <= 1.5 * eta_bar + 3e-11, "sampling nuclear bound")
            norm_loss = float(np.sum(absolute_k) - np.sum(np.abs(spectrum_t)))
            check(norm_loss >= -3e-10, "compression trace-norm loss is nonnegative")
            a = GAP / 2.0
            ratio = np.minimum(absolute_k / a, 1.0)
            residual_mass = float(np.sum(np.where(absolute_k <= a,
                                      absolute_k * (1.0 - ratio)**3, 0.0)))
            local_compression_bound = residual_mass + 6.0 * eta / a
            check(norm_loss <= local_compression_bound + 3e-10, "local compression loss bound")
            for sign in (1, -1):
                k_list = np.sort(np.maximum(sign * spectrum_k, 0.0))[::-1]
                t_list = np.sort(np.maximum(sign * spectrum_t, 0.0))[::-1]
                size = max(len(k_list), len(t_list))
                difference = np.pad(t_list, (0, size - len(t_list))) - np.pad(k_list, (0, size - len(k_list)))
                check(np.max(difference) <= 3e-10, "signed min-max list domination")
            nuclear_bound = 16.0 * PI * tau**2 / (3.0 * d)
            local_min_bound = min(nuclear_bound, residual_mass + 6.0 * eta / a + 1.5 * eta_bar)
            actual_pairing_bound = sample_error + max(0.0, norm_loss)
            traces = []
            for p in [1.0, 1.5, 2.0, 3.0]:
                value_a = float(np.sum(np.abs(spectrum_a)**p))
                value_k = float(np.sum(absolute_k**p))
                value_t = float(np.sum(np.abs(spectrum_t)**p))
                difference = abs(value_a - value_k)
                check(abs(value_t - value_k) <= p * max(0.0, norm_loss) + 3e-10,
                      "Lipschitz trace versus actual norm loss")
                check(difference <= p * actual_pairing_bound + 3e-10,
                      "finite Schatten trace versus actual pairing bound")
                check(difference <= p * local_min_bound + 3e-10,
                      "finite Schatten trace versus minimum analytic local bound")
                traces.append({"p": p, "finite_trace": value_a, "continuum_trace": value_k,
                               "absolute_difference": difference,
                               "actual_pairing_bound_ratio": difference / (p * actual_pairing_bound)})
            exact_delta = 2.0 * np.sin(PI * d * epsilon) * (epsilon / np.sin(PI * epsilon) - 1.0 / PI)
            linear_error = close(np.trace(matrix_a) - 2.0 * np.sin(2.0 * PI * tau) / PI,
                                 exact_delta, 5e-13, "exact finite linear correction")
            prolate = 2.0 * epsilon * np.sinc(2.0 * epsilon * (m[:, None] - m[None, :]))
            reflected_a = float(np.sum(matrix_a * matrix_a[::-1, :]))
            reflected_c = float(np.sum(prolate * prolate[::-1, :]))
            reflected_error = close(reflected_a, reflected_c, 2e-12, "finite reflected second-moment identity")
            rows.append({"tau": tau, "d": d, "epsilon": epsilon, "traces": traces,
                         "sampling_identity_error": sample_identity, "linear_identity_error": linear_error,
                         "finite_reflection_identity_error": reflected_error,
                         "sampling_trace_norm_error": sample_error, "eta": eta, "eta_bar": eta_bar,
                         "compression_trace_norm_loss": norm_loss,
                         "actual_central_residual_mass": residual_mass,
                         "local_compression_bound": local_compression_bound,
                         "minimum_local_lipschitz_bound": local_min_bound})
    return rows


def compare_resolutions(low, high):
    parity_error = close([r["reflected_trace"] for r in low["parity"]["exact_identities"]],
                         [r["reflected_trace"] for r in high["parity"]["exact_identities"]],
                         4e-10, "parity resolution agreement")
    asymptotic_error = close([r["reflected_trace"] for r in low["parity"]["asymptotic_samples"]],
                             [r["reflected_trace"] for r in high["parity"]["asymptotic_samples"]],
                             5e-11, "scalar asymptotic resolution agreement")
    eigen_error = close([r["eigenvalue_magnitudes"] for r in low["monotonicity"]],
                        [r["eigenvalue_magnitudes"] for r in high["monotonicity"]],
                        5e-11, "signed rank resolution agreement")
    finite_error = close([[v["continuum_trace"] for v in r["traces"]] for r in low["finite"]],
                         [[v["continuum_trace"] for v in r["traces"]] for r in high["finite"]],
                         5e-10, "finite-to-continuum resolution agreement")
    convergence = []
    for index in range(0, len(high["finite"]), 2):
        small, large = high["finite"][index:index+2]
        ratios = []
        for first, second in zip(small["traces"], large["traces"]):
            ratios.append(second["absolute_difference"] / first["absolute_difference"])
            check(second["absolute_difference"] < first["absolute_difference"],
                  "fixed-tau trace errors decrease when finite dimension grows")
        convergence.append({"tau": small["tau"], "dimensions": [small["d"], large["d"]],
                            "error_ratios_p_1_1p5_2_3": ratios})
    return {"parity_max_difference": parity_error, "asymptotic_max_difference": asymptotic_error,
            "signed_eigenvalue_max_difference": eigen_error, "trace_max_difference": finite_error,
            "fixed_tau_dimension_convergence": convergence}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--resolutions", nargs=2, type=int, default=[96, 144])
    args = parser.parse_args()
    check(all(n >= 48 and n % 2 == 0 for n in args.resolutions), "resolutions must be even and at least 48")
    manifest = {}
    for name, labels in SOURCES.items():
        content = (args.source_dir / name).read_bytes()
        source = content.decode("utf-8")
        for label in labels:
            check(f"\\label{{{label}}}" in source, f"missing source anchor {name}:{label}")
        manifest[name] = hashlib.sha256(content).hexdigest()
    env = dict(os.environ, TZ="Asia/Jerusalem")
    timestamp = subprocess.check_output(["date", "+%Y-%m-%d %H:%M:%S %Z"], env=env, text=True).strip()
    result = {"status": "PASS", "scope": "floating-point diagnostic samples; not a proof or interval certificate",
              "timestamp_jerusalem": timestamp,
              "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
              "source_sha256": manifest,
              "global_constant_scope": "The manuscript proves C_V<65 and C_0<138; these local finite checks use actual residual mass, not those conservative bounds",
              "resolutions": []}
    for index, n in enumerate(args.resolutions):
        panel_order = [12, 20][index]
        row = {"matrix_order": n, "scalar_panel_order": panel_order,
               "coefficients": coefficient_checks(n), "parity": parity_checks(n, panel_order),
               "monotonicity": monotonicity_checks(n), "finite": finite_checks(n)}
        result["resolutions"].append(row)
    result["convergence"] = compare_resolutions(*result["resolutions"])
    # Detect concurrent manuscript edits rather than falsely naming one snapshot.
    for name, digest in manifest.items():
        check(hashlib.sha256((args.source_dir / name).read_bytes()).hexdigest() == digest,
              f"source changed during diagnostic: {name}; rerun after integration")
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "timestamp_jerusalem": timestamp,
                      "versions": result["versions"], "convergence": result["convergence"],
                      "output": str(args.output) if args.output else None,
                      "scope": result["scope"]}, indent=2))


if __name__ == "__main__":
    main()
