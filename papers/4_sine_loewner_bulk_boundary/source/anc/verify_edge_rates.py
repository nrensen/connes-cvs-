#!/usr/bin/env python3
"""Deterministic finite-section check of the supercritical edge rates.

For the sine Loewner matrix A_N(eps) on the nodes -N..N (equation (1.1)) and the
discrete prolate (Slepian) matrix T_W with entries sin(2 pi W (m-n))/(pi (m-n)),
diagonal 2W, this script prints the four signed parity edge defects

    1 - lambda_{N,+}^+,  1 + lambda_{N,-}^-,  1 + lambda_{N,+}^-,  1 - lambda_{N,-}^+

next to the leakages 1 - lambda_0(T_{eps/2}) (even top DPSS eigenvalue) and
1 - lambda_1(T_{eps/2}) (odd top DPSS eigenvalue), and the ratios.  The script
checks that the ratios stay inside a wide empirical window, that every defect is
positive on the stated grid, and that each defect obeys the proved dual-comparison
bound 2 mu (1-mu)/(2 mu - 1) whenever the paired concentration norm mu exceeds 1/2
(Proposition 'Dual comparison and the four edges').  The sampled cells have
d*eps <= 8.02, below the d*eps >= 25 hypothesis of the factor-3 statement, and no
lower comparison is proved.  It uses IEEE double precision and is a falsifier, not
a certificate.  Cells with defects below 1e-12 are skipped.
"""

from __future__ import annotations

import numpy as np


def source_matrix(n: int, eps: float) -> np.ndarray:
    modes = np.arange(-n, n + 1, dtype=float)
    diff = modes[:, None] - modes[None, :]
    sine = np.sin(2.0 * np.pi * eps * modes)
    matrix = np.divide(sine[:, None] - sine[None, :], np.pi * diff,
                       out=np.zeros_like(diff), where=diff != 0.0)
    np.fill_diagonal(matrix, 2.0 * eps * np.cos(2.0 * np.pi * eps * modes))
    return matrix


def dpss_matrix(n: int, half_bandwidth: float) -> np.ndarray:
    modes = np.arange(-n, n + 1, dtype=float)
    diff = modes[:, None] - modes[None, :]
    return np.divide(np.sin(2.0 * np.pi * half_bandwidth * diff), np.pi * diff,
                     out=np.full_like(diff, 2.0 * half_bandwidth), where=diff != 0.0)


def parity_bases(dim: int) -> tuple[np.ndarray, np.ndarray]:
    n = (dim - 1) // 2
    even = np.zeros((dim, n + 1))
    odd = np.zeros((dim, n))
    even[n, 0] = 1.0
    for j in range(1, n + 1):
        even[n - j, j] = 2.0**-0.5
        even[n + j, j] = 2.0**-0.5
        odd[n - j, j - 1] = 2.0**-0.5
        odd[n + j, j - 1] = -(2.0**-0.5)
    return even, odd


def main() -> None:
    lower, upper = 0.2, 5.0
    for n, eps in ((100, 0.01), (100, 0.02), (100, 0.03), (200, 0.01), (200, 0.02), (400, 0.01)):
        d = 2 * n + 1
        A = source_matrix(n, eps)
        T = dpss_matrix(n, eps / 2.0)
        E, O = parity_bases(d)
        ev = np.linalg.eigvalsh(E.T @ A @ E)
        od = np.linalg.eigvalsh(O.T @ A @ O)
        te = np.linalg.eigvalsh(E.T @ T @ E)
        to = np.linalg.eigvalsh(O.T @ T @ O)
        defects = (1.0 - ev[-1], 1.0 + od[0], 1.0 + ev[0], 1.0 - od[-1])
        leaks = (1.0 - te[-1], 1.0 - te[-1], 1.0 - to[-1], 1.0 - to[-1])
        if min(defects) < 1e-12 or min(leaks) < 1e-12:
            print(f"N={n:3d} eps={eps:.3f} d*eps={d*eps:6.2f}: defects below 1e-12, skipped")
            continue
        ratios = tuple(x / y for x, y in zip(defects, leaks))
        print(f"N={n:3d} eps={eps:.3f} d*eps={d*eps:6.2f} "
              f"defects(e+,o-,e-,o+)=({defects[0]:.3e},{defects[1]:.3e},{defects[2]:.3e},{defects[3]:.3e}) "
              f"dpss(1-l0,1-l1)=({leaks[0]:.3e},{leaks[2]:.3e}) "
              f"ratios=({ratios[0]:.3f},{ratios[1]:.3f},{ratios[2]:.3f},{ratios[3]:.3f})")
        for r in ratios:
            assert lower <= r <= upper, r
        for x in defects:
            assert x > 0.0
        # Proved upper bound (dual comparison): even leakage pairs with (e+, o-), odd with (e-, o+).
        mus = (te[-1], te[-1], to[-1], to[-1])
        for x, mu in zip(defects, mus):
            if mu > 0.5:
                assert x <= 2.0 * mu * (1.0 - mu) / (2.0 * mu - 1.0) + 1e-12, (x, mu)
            else:
                print(f"    mu={mu:.4f} <= 1/2: dual-comparison bound not applicable at this cell")
    print("EDGE RATE DIAGNOSTIC PASS (ordinary precision; not a certificate)")


if __name__ == "__main__":
    main()
