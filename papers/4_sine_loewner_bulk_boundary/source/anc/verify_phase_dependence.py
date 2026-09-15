#!/usr/bin/env python3
"""Independent quadrature replay of the phase-resolved spectral illustration.

Gauss-Legendre uses divided differences; Clenshaw-Curtis uses the stable
trigonometric form. Agreement is diagnostic, not a certified operator enclosure.
"""
import json
import numpy as np
from scipy.special import roots_legendre

FIXTURES = {
    8.: [.373543495437,.466255507032],
    8.125: [.605496937515,.716074897848],
    8.25: [-.098189654522,-.064436120477],
    8.375: [-.236090353158,-.172943608256],
    8.5: [-.464887389540,-.374640062374],
    16.125: [.612028524251,.701282374557],
}


def clenshaw_curtis(n):
    assert n%2==0
    theta=np.pi*np.arange(n+1)/n
    w=np.zeros(n+1); v=np.ones(n-1)
    w[0]=w[-1]=1/(n*n-1)
    for k in range(1,n//2):
        v-=2*np.cos(2*k*theta[1:-1])/(4*k*k-1)
    v-=np.cos(n*theta[1:-1])/(n*n-1)
    w[1:-1]=2*v/n
    return np.cos(theta),w


def eigenvalues(tau,n,method):
    x,w=roots_legendre(n) if method=='GL' else clenshaw_curtis(n)
    X=x[:,None];Y=x[None,:]
    if method=='GL':
        d=X-Y
        with np.errstate(divide='ignore',invalid='ignore'):
            k=(np.sin(2*np.pi*tau*X)-np.sin(2*np.pi*tau*Y))/(np.pi*d)
        np.fill_diagonal(k,2*tau*np.cos(2*np.pi*tau*x))
    else:
        k=2*tau*np.cos(np.pi*tau*(X+Y))*np.sinc(tau*(X-Y))
    return np.linalg.eigvalsh(k*np.sqrt(w[:,None]*w[None,:]))


def main():
    rows=[]
    for tau,expected in FIXTURES.items():
        print(f'Phase comparison tau={tau}',flush=True)
        values=[eigenvalues(tau,n,m) for m,n in [('GL',256),('GL',384),('CC',384)]]
        chosen=[v[(abs(v)>.05)&(abs(v)<.72)] for v in values]
        assert all(len(v)==len(expected) for v in chosen)
        spread=max(float(max(abs(v-chosen[-1]))) for v in chosen)
        assert spread<1e-10
        assert max(abs(chosen[-1]-expected))<1e-10
        rows.append({'tau':tau,'selected_eigenvalues':chosen[-1].tolist(),'agreement':spread})
    print(json.dumps(rows,indent=2))
    print('PHASE DEPENDENCE DIAGNOSTIC PASS (ordinary precision; not a certificate)')


if __name__=='__main__':main()
