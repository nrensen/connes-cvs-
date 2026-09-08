#!/usr/bin/env python3
"""Replay logarithm-centered bounds and refine the quartic trace itself.

The algebraic identity with W=C-K^2 checks implementation consistency. A separate
external-integral construction of W supplies a non-tautological boundary check.
Ordinary precision; diagnostics, not certificates. Both integer and off-integer
phases are sampled.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import roots_legendre, sici

EULER = 0.57721566490153286
CROSS_BOUND, W2_BOUND, QUARTIC_BOUND = 7.33, 15.676, 32.0


def operators(tau: float, n: int):
    x, w = roots_legendre(n)
    d = x[:, None] - x[None, :]
    gram = np.sqrt(w[:, None] * w[None, :])
    k = 2*tau*np.cos(np.pi*tau*(x[:, None]+x[None, :]))*np.sinc(tau*d)
    c = 2*tau*np.sinc(2*tau*d)
    return gram*k, gram*c


def traces(tau: float, n: int):
    k, c = operators(tau, n)
    k2 = k@k
    w = c-k2
    dr = float(np.trace(c)-np.sum(c*c))
    dj = float(np.trace(w))
    cw, w2 = float(np.sum(c*w)), float(np.sum(w*w))
    quartic = float(np.trace(k2)-np.sum(k2*k2))
    return dr, dj, cw, w2, quartic, float(np.trace(k2))


def exact_prolate_defect(tau: float) -> float:
    a = 2*math.pi*tau
    si, ci = sici(4*a)
    return (2*a*math.pi-4*a*si+1-math.cos(4*a)
            +math.log(4*a)+EULER-ci)/math.pi**2


def boundary_check() -> None:
    tau, n, zmax = 2.0, 96, 80.0
    x, w = roots_legendre(n)
    z0, zw0 = roots_legendre(10)
    edges = np.linspace(1.0, zmax, int(8*tau*(zmax-1))+1)
    mid, half = (edges[1:]+edges[:-1])/2, np.diff(edges)/2
    z = (mid[:, None]+half[:, None]*z0).ravel()
    zw = (half[:, None]*np.broadcast_to(zw0, (len(half), len(zw0)))).ravel()
    z, zw = np.concatenate((z, -z)), np.concatenate((zw, zw))
    kernel = 2*tau*np.cos(np.pi*tau*(x[:, None]+z))*np.sinc(tau*(x[:, None]-z))
    u = np.sqrt(w)[:, None]*kernel*np.sqrt(zw)[None, :]
    external_w = u@u.T
    k, c = operators(tau, n)
    difference = c-k@k-external_w
    # The omitted exterior tail is positive with trace <=16/[pi^2(zmax-1)].
    budget = 16/(math.pi**2*(zmax-1))
    error = float(np.linalg.norm(difference, 2))
    assert error <= budget+1e-9, (error, budget)
    assert float(np.linalg.eigvalsh(difference)[0]) >= -1e-9
    print(f"Independent exterior W check: norm error={error:.3e}, tail budget={budget:.3e}", flush=True)


def main() -> None:
    boundary_check()
    print("tau | identity | exact D_R | log-centered CW | log-centered W2 | quartic remainder | refinement", flush=True)
    for tau in (2., 4., 8., 8.125, 8.25, 8.375, 8.5, 16., 32.):
        n = max(128, int(math.ceil(12*tau)))
        dr, dj, cw, w2, quartic, quadratic = traces(tau, n)
        fine = traces(tau, max(n+64, int(math.ceil(16*tau))))
        identity = quartic-(dr-dj+2*cw-w2)
        rcw = cw-3/(2*math.pi**2)*math.log(tau)
        rw2 = w2-2/(3*math.pi**2)*math.log(tau)
        rq = quartic-4/(3*math.pi**2)*math.log(tau)
        refinement = max(abs(quartic-fine[4]), abs(quadratic-fine[5]),
                         abs(cw-fine[2]), abs(w2-fine[3]))
        assert abs(identity)<1e-9
        assert abs(dr-exact_prolate_defect(tau))<1e-9
        assert refinement<1e-8
        assert abs(rcw)<=CROSS_BOUND and abs(rw2)<=W2_BOUND
        assert abs(rq)<=QUARTIC_BOUND
        print(f"{tau:7.3f} | {identity:+.2e} | {dr-exact_prolate_defect(tau):+.2e} | {rcw:+.6f} | {rw2:+.6f} | {rq:+.6f} | {refinement:.2e}", flush=True)
        # The identity is an algebraic tautology (implementation consistency); the bounds are proved
        # constants, not estimates of the remainders. Ratios observed/bound are printed as diagnostics.
        print(f"          ratios observed/bound: CW {abs(rcw)/CROSS_BOUND:.3f}, W2 {abs(rw2)/W2_BOUND:.3f}, quartic {abs(rq)/QUARTIC_BOUND:.3f}", flush=True)
    print("QUARTIC LAW DIAGNOSTIC PASS (ordinary precision; not a certificate)")


if __name__ == '__main__': main()
