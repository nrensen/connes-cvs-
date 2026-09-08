#!/usr/bin/env python3
"""Deterministic falsifier for the signed spectral census of the critical operator.

Reproduces the retained census fixture (not printed in the paper) and checks the
explicit interior-window count bound (Proposition 'Explicit bounds for an interior
window') and the bounded compact gap counts (Corollary 'Gap counts and signed
transition counts').  The checks use ordinary double precision.  They
test counts, normalizations and asymptotic trends; they are not interval
certificates and do not replace the proofs.

The spectrum of the limit operator K_tau on L^2[-1,1] is computed by a Nystrom
discretisation in the Gauss-Legendre nodes.  The kernel

    K_tau(x,y) = [sin(2 pi tau x) - sin(2 pi tau y)] / [pi (x-y)]
               = 2 tau cos(pi tau (x+y)) sinc(tau (x-y)),   sinc(z)=sin(pi z)/(pi z),

is entire, the singularity at x=y being removable, so the discretisation
converges rapidly; the second form is the one evaluated, since it is free of
the cancellation that the first form suffers near x=y.  K_tau commutes with
f -> f(-.), which halves the matrices.
"""

from __future__ import annotations

import os
for _thread_variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
                         "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ[_thread_variable] = "1"

import math
import argparse
import platform
import resource
import time
import numpy as np
from scipy.linalg import eigvalsh
from scipy.special import roots_legendre

GAP_WINDOW = (0.05, 0.69)
BAND_WINDOW = (0.72, 0.995)
NEAR_ONE = 0.99
QUARTIC_SLOPE = 4.0 / (3.0 * math.pi**2)
QUARTIC_CONSTANT = 32.0
EXPECTED_COUNTS = {
    8: (2, 0, 2, 3, 26), 16: (2, 0, 2, 3, 58),
    32: (2, 0, 3, 4, 122), 64: (2, 0, 3, 4, 248),
    128: (2, 0, 4, 4, 504), 256: (2, 0, 4, 5, 1016),
    512: (2, 0, 4, 5, 2039), 1024: (2, 0, 5, 5, 4086),
}

EXPECTED_GAP_PAIRS = {
    8: (0.37354349543665194, 0.46625550703239022),
    16: (0.38450325682264763, 0.45271768266236517),
    32: (0.39262607403457417, 0.4429261025692805),
    64: (0.39867131151834945, 0.43581912079618462),
    128: (0.40318056822411086, 0.43064449098741481),
    256: (0.4065471553829485, 0.4268663498235965),
    512: (0.40906068441731569, 0.42410110829236281),
    1024: (0.41093638023818846, 0.42207294963425968),
}
# Acceptance tolerances for ordinary-precision diagnostics, not error enclosures.
GAP_FIXTURE_TOLERANCE = 5e-8
REFINEMENT_TOLERANCE = 5e-8
_SPECTRUM_CACHE: dict[tuple[float, float], np.ndarray] = {}


def peak_rss_gib() -> float:
    scale = 1 if platform.system() == "Darwin" else 1024
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale / 2**30


def kernel(tau: float, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """K_tau(x,y) in the cancellation-free form, valid on and off the diagonal."""
    return 2.0 * tau * np.cos(np.pi * tau * (x + y)) * np.sinc(tau * (x - y))


def spectrum(tau: float, nodes_per_tau: float) -> np.ndarray:
    """One construction/solve per parity; cache only the resulting eigenvalues."""
    cache_key = (tau, nodes_per_tau)
    if cache_key in _SPECTRUM_CACHE:
        return _SPECTRUM_CACHE[cache_key]
    n = int(round(nodes_per_tau * tau))
    if n % 2:
        n += 1
    x, w = roots_legendre(n)
    xp, wp = x[n // 2:], w[n // 2:]
    scale = np.sqrt(wp)
    values = []
    for sign, parity in ((1, "even"), (-1, "odd")):
        started = time.perf_counter()
        print(f"  tau={tau:g} density={nodes_per_tau:g} {parity}: "
              f"building {len(xp)} x {len(xp)} block", flush=True)
        same = kernel(tau, xp[:, None], xp[None, :])
        opposite = kernel(tau, xp[:, None], -xp[None, :])
        same += sign * opposite
        del opposite
        same *= scale[:, None]
        same *= scale[None, :]
        matrix = np.array(same, order="F")
        del same
        print(f"  tau={tau:g} density={nodes_per_tau:g} {parity}: "
              f"solving; peak RSS {peak_rss_gib():.3f} GiB", flush=True)
        values.append(eigvalsh(matrix, driver="evr", overwrite_a=True,
                               check_finite=False))
        del matrix
        print(f"  tau={tau:g} density={nodes_per_tau:g} {parity}: "
              f"completed in {time.perf_counter()-started:.2f}s", flush=True)
    result = np.sort(np.concatenate(values))
    _SPECTRUM_CACHE[cache_key] = result
    return result


def window_count(values: np.ndarray, low: float, high: float) -> int:
    return int(np.count_nonzero((values > low) & (values < high)))


def in_window(values: np.ndarray, low: float, high: float) -> np.ndarray:
    return values[(values > low) & (values < high)]


def report_stability(taus=(8,)) -> dict[int, np.ndarray]:
    """Exactly the 8-versus-13 node-density comparison and fixed gap-pair fixtures."""
    print("Gap-pair refinement: exactly 8 and 13 nodes per unit of tau")
    print("No density-26 or interval-certification claim is made.")
    quoted = {}
    for tau in taus:
        rows = []
        for density in (8, 13):
            values = spectrum(tau, density)
            inside = in_window(values, *GAP_WINDOW)
            assert len(inside) == 2, (tau, density, inside)
            assert np.max(np.abs(inside - np.asarray(EXPECTED_GAP_PAIRS[tau]))) < GAP_FIXTURE_TOLERANCE, \
                "gap-pair fixture differs beyond the stated diagnostic tolerance"
            rows.append(inside)
            print("  tau=%d density=%d gap pair %.12f %.12f" %
                  (tau, density, inside[0], inside[1]), flush=True)
        spread = float(np.max(np.abs(rows[0] - rows[1])))
        print(f"  tau={tau}: maximum 8-versus-13 pair spread {spread:.12e}", flush=True)
        assert spread < REFINEMENT_TOLERANCE, "gap pair changed under refinement"
        quoted[tau] = rows[1]
    print()
    return quoted


def report_census(taus, density: float) -> dict[int, np.ndarray]:
    print("Census of the spectrum of K_tau  (%g nodes per unit of tau)" % density)
    print("  tau   (%.2f,%.2f)  (%.2f,%.2f)   (%.3f,%.3f)  (%.3f,%.3f)   |lambda|>%.2f    4 tau"
          % (GAP_WINDOW[0], GAP_WINDOW[1], -GAP_WINDOW[1], -GAP_WINDOW[0],
             BAND_WINDOW[0], BAND_WINDOW[1], -BAND_WINDOW[1], -BAND_WINDOW[0], NEAR_ONE))
    spectra = {}
    for tau in taus:
        print(f"Census tau={tau}, density={density}: checking cached spectrum", flush=True)
        values = spectrum(tau, density)
        spectra[tau] = values
        gap_pos = window_count(values, *GAP_WINDOW)
        gap_neg = window_count(values, -GAP_WINDOW[1], -GAP_WINDOW[0])
        band_pos = window_count(values, *BAND_WINDOW)
        band_neg = window_count(values, -BAND_WINDOW[1], -BAND_WINDOW[0])
        near = int(np.count_nonzero(np.abs(values) > NEAR_ONE))
        print("%5d %11d %12d %13d %12d %14d %8d"
              % (tau, gap_pos, gap_neg, band_pos, band_neg, near, 4 * tau))
        assert (gap_pos, gap_neg, band_pos, band_neg, near) == EXPECTED_COUNTS[tau], \
            "published census fixture differs"
        assert (gap_pos, gap_neg) == (2, 0), "gap-window count moved"
        assert abs(values).max() <= 1.0 + 1e-9, "norm bound violated"
        assert near <= 4 * tau, "count above 0.99 exceeds 4 tau"
        # The guard bands 0.69 < |lambda| < 0.72 straddling +- 1/sqrt(2) separate
        # the two windows.  They are empty at these integer parameters, so
        # the two windows and their
        # reflections account for every eigenvalue with 0.05 < |lambda| < 0.995.
        guard = int(np.count_nonzero((np.abs(values) > GAP_WINDOW[1])
                                     & (np.abs(values) < BAND_WINDOW[0])))
        assert guard == 0, "an eigenvalue sits in the guard band around 1/sqrt(2)"
        inside = int(np.count_nonzero((np.abs(values) > GAP_WINDOW[0])
                                      & (np.abs(values) < BAND_WINDOW[1])))
        assert inside == gap_pos + gap_neg + band_pos + band_neg, \
            "the two windows do not exhaust 0.05 < |lambda| < 0.995"
    print()
    return spectra


def report_gap_pair(quoted: dict[int, np.ndarray],
                    spectra: dict[int, np.ndarray]) -> None:
    low, mid512, high = quoted[8], quoted[512], quoted[1024]
    separations = (high[1] - high[0], low[1] - low[0])
    factor = (separations[0] / separations[1]) ** (1.0 / 7.0)
    midpoints = (0.5 * (low[0] + low[1]), 0.5 * (high[0] + high[1]))
    print("The pair inside the gap window")
    print("  tau=8    : %.6f  %.6f   separation %.6f  midpoint %.6f"
          % (low[0], low[1], separations[1], midpoints[0]))
    print("  tau=512  : %.6f  %.6f   separation %.6f  midpoint %.6f"
          % (mid512[0], mid512[1], mid512[1] - mid512[0],
             0.5 * (mid512[0] + mid512[1])))
    print("  tau=1024 : %.6f  %.6f   separation %.6f  midpoint %.6f"
          % (high[0], high[1], separations[0], midpoints[1]))
    ladder = [t for t in sorted(spectra) if 8 <= t <= 1024]
    pairs = [in_window(spectra[t], *GAP_WINDOW) for t in ladder]
    seps = [float(p[1] - p[0]) for p in pairs]
    mids = [0.5 * float(p[0] + p[1]) for p in pairs]
    ratios = [seps[i + 1] / seps[i] for i in range(len(seps) - 1)]
    print("  per doubling of tau the separation falls by")
    print("    " + "  ".join("%d->%d %.4f" % (ladder[i], ladder[i + 1], r)
                             for i, r in enumerate(ratios)))
    print("  monotone increasing: %s;  range %.4f to %.4f;  geometric mean %.4f"
          % (all(ratios[i] < ratios[i + 1] for i in range(len(ratios) - 1)),
             min(ratios), max(ratios), factor))
    print("  the midpoint decreases monotonically (%s) and moves by %.6f in all"
          % (all(mids[i] > mids[i + 1] for i in range(len(mids) - 1)),
             abs(midpoints[1] - midpoints[0])))
    assert all(ratios[i] < ratios[i + 1] for i in range(len(ratios) - 1)), \
        "sampled separation ratios are not increasing"
    assert all(mids[i] > mids[i + 1] for i in range(len(mids) - 1)), \
        "sampled pair midpoints are not decreasing"
    print()


def report_count_stability(taus, densities) -> None:
    """The census counts must not move with the node density."""
    print("Stability of the table counts under a change of node density")
    print("  tau   nodes/tau   gap+ gap-  band+ band-  |lambda|>%.2f" % NEAR_ONE)
    for tau in taus:
        reference = None
        for density in densities:
            print(f"Count stability tau={tau}, density={density}: checking cached spectrum", flush=True)
            values = spectrum(tau, density)
            row = (window_count(values, *GAP_WINDOW),
                   window_count(values, -GAP_WINDOW[1], -GAP_WINDOW[0]),
                   window_count(values, *BAND_WINDOW),
                   window_count(values, -BAND_WINDOW[1], -BAND_WINDOW[0]),
                   int(np.count_nonzero(np.abs(values) > NEAR_ONE)))
            print("%5d %10d %6d %4d %6d %5d %9d" % ((tau, density) + row))
            assert row == EXPECTED_COUNTS[tau], "published census fixture differs"
            guard = int(np.count_nonzero((np.abs(values) > GAP_WINDOW[1])
                                         & (np.abs(values) < BAND_WINDOW[0])))
            assert guard == 0, "guard band is occupied at a sampled integer tau"
            inside = int(np.count_nonzero((np.abs(values) > GAP_WINDOW[0])
                                          & (np.abs(values) < BAND_WINDOW[1])))
            assert inside == sum(row[:4]), "windows do not exhaust the interior census"
            if reference is None:
                reference = row
            else:
                assert row == reference, "counts moved with the node density"
    print()


def report_kappa() -> float:
    """kappa = (2 pi^2)^{-1} int_{r_t in [a,b]} dt/(t(1-t)), evaluated in closed form.

    With u = t - 1/2 one has r_t^2 = 2u^2 + 1/2 and t(1-t) = 1/4 - u^2, so
    r_t in [a,b] is |u| in [u_a,u_b] with u_c = sqrt((c^2-1/2)/2), and
    int du/(1/4-u^2) = log((1/2+u)/(1/2-u)).
    """
    a, b = BAND_WINDOW
    ua, ub = math.sqrt((a * a - 0.5) / 2.0), math.sqrt((b * b - 0.5) / 2.0)
    ell = lambda u: math.log((0.5 + u) / (0.5 - u))
    kappa = 2.0 * (ell(ub) - ell(ua)) / (2.0 * math.pi**2)
    print("Second-order coefficient for the band window (%.3f,%.3f)" % (a, b))
    print("  kappa = %.6f,  kappa log 128 = %.4f" % (kappa, kappa * math.log(128.0)))
    print()
    return kappa


def report_counting_bound(spectra: dict[int, np.ndarray]) -> None:
    print("The unconditional counting bound, over a grid of windows")
    print("  tau      a       b      count    m(a,b) * count      tr(K^2-K^4)      bound")
    edges = [0.05, 0.2, 0.4, 0.6, 0.72, 0.85, 0.95, 0.99]
    worst = 0.0
    for tau, values in sorted(spectra.items()):
        trace = float(np.sum(values**2 - values**4))
        bound = QUARTIC_SLOPE * math.log(tau) + QUARTIC_CONSTANT
        assert abs(trace - QUARTIC_SLOPE * math.log(tau)) <= QUARTIC_CONSTANT, \
            "fourth-moment bound violated"
        for i, a in enumerate(edges):
            for b in edges[i + 1:]:
                m = min(a * a * (1 - a * a), b * b * (1 - b * b))
                count = int(np.count_nonzero((np.abs(values) >= a) & (np.abs(values) <= b)))
                assert m * count <= trace + 1e-9, "trace inequality violated"
                assert count <= bound / m + 1e-9, "counting bound violated"
                worst = max(worst, m * count / trace)
        a, b = 0.05, 0.99
        m = min(a * a * (1 - a * a), b * b * (1 - b * b))
        count = int(np.count_nonzero((np.abs(values) >= a) & (np.abs(values) <= b)))
        print("%5d  %6.2f  %6.2f  %7d  %16.6f  %15.6f  %9.1f"
              % (tau, a, b, count, m * count, trace, bound / m))
    print("  worst ratio of the left side to tr(K^2-K^4) over all windows: %.4f" % worst)
    print()


def report_quartic_slope(spectra: dict[int, np.ndarray]) -> None:
    taus = sorted(spectra)
    traces = [float(np.sum(spectra[t]**2 - spectra[t]**4)) for t in taus]
    logs = [math.log(t) for t in taus]
    slope = np.polyfit(logs, traces, 1)[0]
    print("Logarithmic coefficient of tr(K^2-K^4)")
    print("  least squares over tau = %s : %.6f" % (taus, slope))
    print("  4/(3 pi^2)                    : %.6f" % QUARTIC_SLOPE)
    print("  relative difference           : %.2e" % abs(slope / QUARTIC_SLOPE - 1.0))
    assert abs(slope / QUARTIC_SLOPE - 1.0) < 0.01, "measured slope off by over one percent"
    print()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full', action='store_true',
                        help='census through tau=1024, exactly densities 8 and 13')
    args = parser.parse_args()
    taus = [8, 16, 32, 64, 128, 256, 512, 1024] if args.full else [8, 16, 32, 64, 128]
    started = time.perf_counter()
    largest_side = 13 * taus[-1] // 2
    print(f"Work scope: tau={taus}; densities=(8,13); one numerical thread; "
          f"largest parity matrix {largest_side} x {largest_side}.", flush=True)
    print("Each spectrum is computed once and reused by all checks; only eigenvalues are cached.", flush=True)
    if args.full:
        print("Resource reference: the independent tau=1024 replay used about 2.8 GiB peak RSS; "
              "time and memory depend on the platform.", flush=True)
    quoted = report_stability(taus)
    spectra = report_census(taus, 13.0)
    report_count_stability(taus, (8.0, 13.0))
    if args.full:
        report_gap_pair(quoted, spectra)
    report_kappa()
    report_counting_bound(spectra)
    report_quartic_slope(spectra)
    if not args.full:
        print("Normal subset completed; use --full for tau=256,512,1024.")
    print(f"Elapsed {time.perf_counter()-started:.2f}s; peak RSS {peak_rss_gib():.3f} GiB.")
    print("SPECTRAL GAP DIAGNOSTIC PASS (ordinary precision; not a certificate)")


if __name__ == "__main__":
    main()
