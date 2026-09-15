#!/usr/bin/env python3
"""Reproduce higher-moment and nonpolynomial numerical observations.

These fixtures and sampled octave increments do not prove a general trace law.
"""
from __future__ import annotations
import json
import math
import numpy as np
from scipy.special import roots_legendre
from verify_spectral_gap import spectrum

OCTAVES = np.array([
    [-2.000094, -2.000023, -2.000006], [-3.335994, -3.335107, -3.334370],
    [-4.277271, -4.272305, -4.269480], [-4.971933, -4.961245, -4.956267],
    [-5.502697, -5.488457, -5.482907], [-5.922014, -5.907692, -5.903234],
    [-6.264565, -6.253194, -6.250783]])


def concentration_spectrum(tau: float, density: float) -> np.ndarray:
    n = int(math.ceil(density*tau)); n += n % 2
    x, w = roots_legendre(n)
    p = x>0; x, w = x[p], w[p]
    gram = np.sqrt(w[:, None]*w[None, :])
    same = 2*tau*np.sinc(2*tau*(x[:, None]-x[None, :]))
    flip = 2*tau*np.sinc(2*tau*(x[:, None]+x[None, :]))
    return np.concatenate([np.linalg.eigvalsh(gram*(same+flip)),
                           np.linalg.eigvalsh(gram*(same-flip))])


def entropy_square(v):
    x=v*v; out=np.zeros_like(x); p=x>0
    out[p]=-x[p]*np.log(x[p]); return out


def entropy_positive(v):
    assert v.min()>-1e-10
    x=np.maximum(v,0); out=np.zeros_like(x); p=x>0
    out[p]=-x[p]*np.log(x[p]); return out


def quantities(v, c, tau):
    return np.array([*[float(np.sum(v**(2*m))-4*tau) for m in range(1,8)],
                     float(np.sum(abs(v))-4*tau),
                     float(np.sum(2*np.sin(np.pi*v/4)**2)-4*tau),
                     float(np.sum(entropy_square(v))),
                     float(np.sum(c*c)-4*tau),
                     float(np.sum(entropy_positive(c)))])



def kappa_power(p: float) -> float:
    """kappa(|x|^p) by its integral form."""
    from scipy.integrate import quad
    r = lambda s: math.sqrt(s*s+(1-s)*(1-s))
    return quad(lambda s: (r(s)**p-1)/(s*(1-s)), 0, 1,
                epsabs=1e-13, epsrel=1e-13)[0]/math.pi**2


def schatten_and_parity() -> None:
    """Coefficients of the Schatten norms, and the parity-resolved second moment."""
    from scipy.special import sici
    k1 = -2/math.pi**2*(math.sqrt(2)*math.log(1+math.sqrt(2))-math.log(2))
    k3 = -2/math.pi**2*(0.5-math.log(2)+5/(2*math.sqrt(2))*math.log(1+math.sqrt(2)))
    assert abs(kappa_power(1)-k1) < 1e-12 and abs(kappa_power(3)-k3) < 1e-12
    assert abs(kappa_power(2)+2/math.pi**2) < 1e-12
    assert abs(kappa_power(4)+10/(3*math.pi**2)) < 1e-12
    ps = [1+0.25*i for i in range(17)]
    ks = [kappa_power(p) for p in ps]
    assert all(ks[i] > ks[i+1] for i in range(len(ks)-1))
    assert all(ks[i]-2*ks[i+1]+ks[i+2] > 0 for i in range(len(ks)-2))
    print(f"kappa_1 {k1:+.12f} | kappa_2 {kappa_power(2):+.12f} | "
          f"kappa_3 {k3:+.12f} | kappa_4 {kappa_power(4):+.12f}", flush=True)
    # the sharp Lipschitz bound |kappa(f)| <= |kappa_1| Lip(f), attained at f(x)=|x|
    r = lambda s: math.sqrt(s*s+(1-s)*(1-s))
    for f, lip in ((lambda x: abs(x), 1.0), (lambda x: x*x, 2.0),
                   (lambda x: math.sin(math.pi*x/2), math.pi/2),
                   (lambda x: x*abs(x), 2.0)):
        from scipy.integrate import quad
        v = quad(lambda s: (f(r(s))+f(-r(s))-f(1.0)-f(-1.0))/(s*(1-s)), 0, 1,
                 epsabs=1e-13, epsrel=1e-13)[0]/(2*math.pi**2)
        assert abs(v) <= abs(k1)*lip+1e-12, (v, lip)
    # tr(R K_tau^2) -> 1/2 with the bound of the parity proposition
    for tau in (1.0, 2.0, 4.0, 8.0, 16.0):
        n = max(400, int(80*tau))
        x, w = roots_legendre(n)
        k = 2*tau*np.cos(np.pi*tau*(x[:, None]+x[None, :]))*np.sinc(tau*(x[:, None]-x[None, :]))
        a = np.sqrt(w[:, None]*w[None, :])*k
        flip = np.eye(n)[::-1]
        val = float(np.trace(flip@a@a))
        bound = (4*math.log(2*math.pi*tau)+2+math.pi)/(math.pi**3*tau)
        assert abs(val-0.5) <= bound, (tau, val, bound)
        print(f"tau {tau:6.2f} | tr(R K^2) {val:+.10f} | deviation {abs(val-0.5):.3e} "
              f"| bound {bound:.3e}", flush=True)
    # the exact finite identity tr(R_N A_N^2) = tr(R_N C_N(2 eps)^2)
    for N, eps in ((6, 0.17), (9, 0.23), (12, 0.4), (20, 0.05)):
        m = np.arange(-N, N+1); dd = m[:, None]-m[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            amat = (np.sin(2*np.pi*m[:, None]*eps)-np.sin(2*np.pi*m[None, :]*eps))/(np.pi*dd)
        np.fill_diagonal(amat, 2*eps*np.cos(2*np.pi*m*eps))
        c = np.where(dd == 0, 2*eps, np.sin(2*np.pi*eps*np.where(dd == 0, 1, dd))/(np.pi*np.where(dd == 0, 1, dd)))
        flip = np.eye(2*N+1)[::-1]
        lhs = float(np.trace(flip@amat@amat)); rhs = float(np.trace(flip@c@c))
        assert abs(lhs-rhs) < 1e-11, (N, eps, lhs, rhs)
    print("parity-resolved second moment and finite identity checked", flush=True)


def main():
    taus=[16,32,64,128]; rows=[]; refinements=[]
    for tau in taus:
        print(f"Extended traces tau={tau}: primary rule",flush=True)
        row=quantities(spectrum(tau,12),concentration_spectrum(tau,12),tau)
        print(f"Extended traces tau={tau}: refinement",flush=True)
        fine=quantities(spectrum(tau,16),concentration_spectrum(tau,16),tau)
        discrepancy=float(max(abs(row-fine)))
        assert discrepancy<1e-8, (tau,discrepancy)
        rows.append(row); refinements.append(discrepancy)
    rows=np.array(rows); increments=np.diff(rows,axis=0).T*math.pi**2/math.log(2)
    assert np.max(abs(increments[:7]-OCTAVES))<5.1e-7
    observations=np.array([[-1.09807,-1.10227,-1.10441], [-1.70608,-1.70613,-1.70624],
                           [1.61400,1.61173,1.61046], [-1.000007,-1.000002,-1.000000],
                           [1.642886,1.643907,1.644420]])
    tolerances=np.array([5.1e-6,5.1e-6,5.1e-6,5.1e-7,5.1e-7])[:,None]
    assert np.all(abs(increments[7:]-observations)<tolerances), increments[7:]-observations
    assert max(abs(rows[:,9]-np.array([1.362492,1.475844,1.589036,1.702139])))<5.1e-7
    print(json.dumps({'taus':taus,'refinement_max':refinements,
                      'normalized_octaves':increments.tolist(),
                      'entropy_square_traces':rows[:,9].tolist()},indent=2))
    schatten_and_parity()
    print("EXTENDED TRACE DIAGNOSTIC PASS (ordinary precision; not a certificate)")


if __name__=='__main__':main()
