#!/usr/bin/env python3
"""Portable rigorous cosine-defect midpoint bound: stdlib + python-flint.

See "Proof of the certified root enclosure" in the manuscript for every
scalar inequality (sections/24_root_certificate.tex). This
validates the analytic midpoint interpolation error and the complete output
tail. It does NOT include rounding of stored midpoint kernel coefficients.

Public function:
    cosine_defect_midpoint_bound(source_edges, target_edges) -> dict of arb
Inputs must be exact rationals (Fraction, int, rational string or [n,d]).
Binary floats are rejected. Source edges must be exactly reflection-symmetric;
target edges must start at 0, include 1, and finish above 1. No source/root solve.
"""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import subprocess
import time
from flint import arb, fmpq, ctx


def rational(value):
    if isinstance(value,Fraction):
        return value
    if isinstance(value,float):
        raise TypeError('Use exact rational endpoints; binary floats are rejected')
    if isinstance(value,(list,tuple)) and len(value)==2:
        return Fraction(int(value[0]),int(value[1]))
    if isinstance(value,(int,str)):
        return Fraction(value)
    if isinstance(value,fmpq):
        return Fraction(str(value))
    raise TypeError('Unsupported exact rational endpoint')


def ball(value):
    value=rational(value)
    return arb(fmpq(value.numerator,value.denominator))


def exact_mesh(nsource=48,outer=6,octaves=40):
    if nsource%2 or nsource<2 or outer<1 or octaves<2:
        raise ValueError('Invalid mesh parameters')
    half=nsource//2
    low=[Fraction(j**4,2*half**4) for j in range(half+1)]
    source=low+[1-x for x in low[-2::-1]]
    target=set(source)
    target.update(1+x for x in source[1:])
    for k in range(1,octaves):
        target.update(Fraction(2**k)*(1+Fraction(j,outer)) for j in range(1,outer+1))
    return source,sorted(target)


def _midpoint_cell(ht,hs,bt,bs):
    return ht*hs*(bt*bt*ht*ht/12+bs*bs*hs*hs/12+bt*bs*ht*hs/8)


def _corner_derivatives(a0,a1,b0,b1,c):
    d=a0+b0
    if d==0:
        return None
    denominator=ball(d)**2
    aa,bb=ball(a1),ball(b1)
    return (c*bb+c**3*aa**3/3)/denominator,c*aa/denominator


def cosine_defect_midpoint_bound(source_edges,target_edges):
    """Return rigorous ball upper-bound expressions for midpoint error norms.

    Each returned expression encloses an explicit analytical upper bound.
    `hs_bound` is the bound to add in norm to any midpoint-value rounding
    error, and then to the independently certified model-trial residual.
    The caller owns ctx.prec; 128 bits or more is recommended.
    """
    source=[rational(x) for x in source_edges]
    target=[rational(x) for x in target_edges]
    if source[0]!=0 or source[-1]!=1 or target[0]!=0 or target[-1]<=1 or 1 not in target:
        raise ValueError('Require source[0,1] and target[0,T] with 1 a boundary, T>1')
    if any(a>=b for a,b in zip(source,source[1:])) or any(a>=b for a,b in zip(target,target[1:])):
        raise ValueError('Edges must increase strictly')
    if any(x+y!=1 for x,y in zip(source,source[::-1])):
        raise ValueError('Source grid must be exactly invariant under s -> 1-s')
    c=arb(11)/7             # pi/2 < 11/7, a rational upper bound.
    c2=c*c
    e_d,e_g,e_f=arb(0),arb(0),arb(0)
    for t0,t1 in zip(target,target[1:]):
        ht=ball(t1-t0)
        for s0,s1 in zip(source,source[1:]):
            hs=ball(s1-s0)
            if t1<=1:
                e_d+=_midpoint_cell(ht,hs,c2/2,c2/2)
                e_g+=_midpoint_cell(ht,hs,3*c2/2,3*c2/2)
                derivatives=_corner_derivatives(1-t1,1-t0,s0,s1,c)
                if derivatives is None:
                    e_f+=c2*ht*hs
                else:
                    ba,bb=derivatives
                    e_f+=_midpoint_cell(ht,hs,ba,bb)
            else:
                if t0<1:
                    raise ValueError('Target cell straddles 1')
                # D=-sin(c alpha)/(alpha+beta), alpha=1-s,beta=t-1.
                derivatives=_corner_derivatives(1-s1,1-s0,t0-1,t1-1,c)
                if derivatives is None:
                    e_d+=c2*ht*hs
                else:
                    ba,bb=derivatives
                    e_d+=_midpoint_cell(ht,hs,bb,ba)
                denominator=ball(t0+s0)
                e_g+=_midpoint_cell(ht,hs,1/denominator**2,c/denominator+1/denominator**2)
                # F is identically 0 above 1.
    # 1/(4 pi^2) < 1/36, 1/(8 pi^2) < 1/72, using pi>3.
    finite=(4*e_d+2*e_g+2*e_f)/36
    top=ball(target[-1])
    tail=(4/(top-1)+2/top)/72
    return {'scalar_D_squared_bound':e_d,'scalar_G_squared_bound':e_g,
            'scalar_F_squared_bound':e_f,'finite_HS_squared_bound':finite,
            'tail_HS_squared_bound':tail,'total_HS_squared_bound':finite+tail,
            'hs_bound':(finite+tail).sqrt()}

