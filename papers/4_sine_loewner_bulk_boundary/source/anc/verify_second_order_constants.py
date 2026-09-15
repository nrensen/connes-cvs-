#!/usr/bin/env python3
"""Deterministic check of the exact second-order trace asymptotics.

For the limit kernel K_tau(x,y) = [sin(2 pi tau x) - sin(2 pi tau y)]/(pi (x-y))
on L^2[-1,1], the concentration operator C_{2 tau} with kernel
sin(2 pi tau (x-y))/(pi (x-y)), and the sign split C_pm = (C_{2 tau} pm K_tau)/2,
this script compares Gauss--Legendre Nystrom traces with the closed forms

    D_R = tr(C_{2tau} - C_{2tau}^2)
        = (1/pi^2) [log(8 pi tau) + 1 + gamma - Ci(8 pi tau) - 2 J(2 pi tau)]      (exact),
    D_J = 4 tau - tr K_tau^2
        = (2/pi^2) [log(4 pi tau) + 1 + gamma + (log 2) cos(4 pi tau)] + O(1/tau),
    tr(C_+ - C_+^2) - tr(C_- - C_-^2) = (4 log 2/pi^2) cos(2 pi tau) + O(1/tau),
    tr(C_pm - C_pm^2) = (1/(4 pi^2)) [3 log tau + log(2 pi^3) + 3(1+gamma)
                          + 16 (log 2) cos^4(pi tau)  (resp. sin^4(pi tau))] + O(1/tau),
    tr(C_+ C_-) = (1/(4 pi^2)) [log(2 pi tau) + 1 + gamma + 2 (log 2) cos(4 pi tau)] + O(1/tau),

and, for the finite matrix A_N(eps) with d = 2N+1, the exact finite identity for
tr A_N(eps)^2 and the discrete law

    2 d eps - tr A_N(eps)^2 = (2/pi^2) [log(4 N sin(pi eps)) + 1 + gamma
                               + (log 2) cos(2 pi eps d)] + O_delta(eps^2 + 1/(N eps)),  eps <= 1/2 - delta.

It uses IEEE double precision, two quadrature resolutions, and checks that the
residuals are small and decrease along a ladder in tau.  It is a falsifier,
not a certificate; the proofs are in the paper.
"""

from __future__ import annotations

import math

import numpy as np
from numpy.polynomial.legendre import leggauss

EULER = 0.57721566490153286
LOG2 = math.log(2.0)


def cosine_integral(x: float, terms: int = 400) -> float:
    """Ci(x) = gamma + log x + sum_{k>=1} (-1)^k x^{2k} / (2k (2k)!)  for moderate x,
    and the asymptotic form for large x."""
    if x < 40.0:
        total = 0.0
        term = 1.0
        for k in range(1, terms):
            term *= -x * x / ((2 * k - 1) * (2 * k))
            total += term / (2 * k)
            if abs(term) < 1e-18 * max(1.0, abs(total)):
                break
        return EULER + math.log(x) + total
    # asymptotic expansion, error below 1e-12 for x >= 40
    f = sum((-1) ** k * math.factorial(2 * k) / x ** (2 * k) for k in range(0, 8)) / x
    g = sum((-1) ** k * math.factorial(2 * k + 1) / x ** (2 * k + 1) for k in range(0, 8)) / x
    return f * math.sin(x) - g * math.cos(x)


def tail_integral_j(a: float, n: int = 200000) -> float:
    """J(a) = int_2^infty cos(2 a u)/u^2 du by the integrated-by-parts form
    -sin(4a)/(8a) + (1/a) int_2^infty sin(2 a u)/u^3 du, the last integral on a fine grid."""
    u = np.linspace(2.0, 2.0 + 400.0 / a, n)
    integrand = np.sin(2.0 * a * u) / u ** 3
    inner = np.trapezoid(integrand, u)
    return -math.sin(4.0 * a) / (8.0 * a) + inner / a


def nystrom(tau: float, n: int):
    x, w = leggauss(n)
    sw = np.sqrt(w)
    X, Y = x[:, None], x[None, :]
    D = X - Y
    with np.errstate(divide="ignore", invalid="ignore"):
        K = np.where(np.abs(D) < 1e-14, 2.0 * tau * np.cos(2.0 * np.pi * tau * X),
                     (np.sin(2.0 * np.pi * tau * X) - np.sin(2.0 * np.pi * tau * Y)) / (np.pi * D))
        C = np.where(np.abs(D) < 1e-14, 2.0 * tau, np.sin(2.0 * np.pi * tau * D) / (np.pi * D))
    K = sw[:, None] * K * sw[None, :]
    C = sw[:, None] * C * sw[None, :]
    return K, C


def continuum_checks() -> None:
    print("continuum: tau | D_R exact-formula residual | D_J residual | cross residual | C+ residual | C- residual | tr(C+C-) residual")
    worst_by_tau = []
    for tau in (2.0, 4.0, 4.25, 8.0, 8.125, 8.5, 16.0, 16.25, 32.0, 32.75):
        n = int(80 + 40 * tau)
        K, C = nystrom(tau, n)
        K2, C2 = nystrom(tau, int(1.5 * n))
        # resolution check on the two basic traces
        assert abs(np.trace(K @ K) - np.trace(K2 @ K2)) < 1e-8
        assert abs(np.trace(C @ C) - np.trace(C2 @ C2)) < 1e-8
        Cp, Cm = (C + K) / 2.0, (C - K) / 2.0
        d_r = np.trace(C) - np.trace(C @ C)
        d_j = 4.0 * tau - np.trace(K @ K)
        dp = np.trace(Cp) - np.trace(Cp @ Cp)
        dm = np.trace(Cm) - np.trace(Cm @ Cm)
        cross = np.trace(Cp @ Cm)
        a = 2.0 * np.pi * tau
        d_r_exact = (math.log(4 * a) + 1 + EULER - cosine_integral(4 * a) - 2 * tail_integral_j(a)) / np.pi ** 2
        c4, c2 = math.cos(4 * np.pi * tau), math.cos(2 * np.pi * tau)
        d_j_pred = 2 / np.pi ** 2 * (math.log(4 * np.pi * tau) + 1 + EULER + LOG2 * c4)
        diff_pred = 4 * LOG2 / np.pi ** 2 * c2
        base = 3 * math.log(tau) + math.log(2 * np.pi ** 3) + 3 * (1 + EULER)
        dp_pred = (base + 16 * LOG2 * math.cos(np.pi * tau) ** 4) / (4 * np.pi ** 2)
        dm_pred = (base + 16 * LOG2 * math.sin(np.pi * tau) ** 4) / (4 * np.pi ** 2)
        cross_pred = (math.log(2 * np.pi * tau) + 1 + EULER + 2 * LOG2 * c4) / (4 * np.pi ** 2)
        res = (abs(d_r - d_r_exact), abs(d_j - d_j_pred), abs((dp - dm) - diff_pred),
               abs(dp - dp_pred), abs(dm - dm_pred), abs(cross - cross_pred))
        print(f"tau={tau:6.3f}: {res[0]:.1e} | {res[1]:.1e} | {res[2]:.1e} | {res[3]:.1e} | {res[4]:.1e} | {res[5]:.1e}")
        assert res[0] < 1e-6, "exact defect formula failed"
        # Proved remainder constants, in units of 1/(pi^3 tau): D_J 4 (Theorem 'Hilbert--Schmidt norm'),
        # signed cross term S 3 (Theorem 'Signed cross term'), tr(C_pm - C_pm^2) 21/8 and
        # tr(C_+ C_-) 9/8 (Corollary 'Signed second-order constants'); 1e-7 covers quadrature noise.
        proved = (4.0, 3.0, 21.0 / 8.0, 21.0 / 8.0, 9.0 / 8.0)
        for r, c in zip(res[1:], proved):
            assert r < c / (np.pi ** 3 * tau) + 1e-7, ("asymptotic law failed beyond its proved remainder", r, c)
        worst_by_tau.append((tau, max(res[1:])))
    # residuals at the integer ladder must decrease
    ladder = [w for t, w in worst_by_tau if float(t).is_integer()]
    assert all(ladder[i + 1] < ladder[i] for i in range(len(ladder) - 1)), ladder


def source_matrix(n: int, eps: float) -> np.ndarray:
    modes = np.arange(-n, n + 1, dtype=float)
    diff = modes[:, None] - modes[None, :]
    sine = np.sin(2.0 * np.pi * eps * modes)
    matrix = np.divide(sine[:, None] - sine[None, :], np.pi * diff, out=np.zeros_like(diff), where=diff != 0.0)
    np.fill_diagonal(matrix, 2.0 * eps * np.cos(2.0 * np.pi * eps * modes))
    return matrix


def finite_identity(n: int, eps: float) -> float:
    d = 2 * n + 1
    k = np.arange(1, 2 * n + 1, dtype=float)
    s = np.sin(np.pi * eps * k) ** 2 / (np.pi ** 2 * k ** 2)
    total = np.sum(2.0 * s * (2.0 * (d - k) + 2.0 * np.sin(2.0 * np.pi * eps * (d - k)) / np.sin(2.0 * np.pi * eps)))
    return total + 2.0 * eps ** 2 * d + 2.0 * eps ** 2 * np.sin(2.0 * np.pi * eps * d) / np.sin(2.0 * np.pi * eps)


def discrete_checks() -> None:
    print("discrete: N, eps | |tr A^2 - finite identity| | deficit residual vs law | allowed")
    for n, eps in ((100, 0.1), (400, 0.1), (1000, 0.01), (4000, 0.0025), (300, 0.3), (500, 0.07), (50, 0.2)):
        d = 2 * n + 1
        A = source_matrix(n, eps)
        tr2 = float(np.sum(A * A))
        ident = finite_identity(n, eps)
        law = 2.0 / np.pi ** 2 * (math.log(4 * n * math.sin(np.pi * eps)) + 1 + EULER + LOG2 * math.cos(2 * np.pi * eps * d))
        resid = abs(2 * d * eps - tr2 - law)
        # Theorem 'Discrete versions' proves O_delta(eps^2 + 1/(N eps)) without an explicit constant;
        # 0.5 is an empirical acceptance envelope, not a proved value.
        allowed = 0.5 * eps ** 2 + 0.5 / (n * eps)
        print(f"N={n:5d} eps={eps:.4f}: {abs(tr2 - ident):.1e} | {resid:.2e} | {allowed:.2e}")
        assert abs(tr2 - ident) < 1e-9 * max(1.0, tr2)
        assert resid < allowed


def main() -> None:
    continuum_checks()
    discrete_checks()
    print("SECOND-ORDER CONSTANTS DIAGNOSTIC PASS (ordinary precision; not a certificate)")


if __name__ == "__main__":
    main()
