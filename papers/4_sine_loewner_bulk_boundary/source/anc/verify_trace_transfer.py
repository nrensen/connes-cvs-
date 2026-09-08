#!/usr/bin/env python3
"""Bounded diagnostics for sampling, local regularity and transition profiles.

Independent finite divided differences, two continuum quadrature orders,
frequency-factor reconstruction, and exact scalar coefficient comparisons.
Only NumPy/SciPy are required. These are error-detection checks, not interval
enclosures or proofs of uniform asymptotic laws. Progress goes to stderr;
the complete reproducible record goes to stdout as JSON.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
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
from scipy.special import roots_legendre


def stamp():
    return subprocess.check_output(['date','+%Y-%m-%d %H:%M:%S %Z'],
        env={**os.environ,'TZ':'Asia/Jerusalem'},text=True).strip()


def require(ok, name, **data):
    if not ok:
        raise AssertionError({'check':name,**data})


def finite_matrix(n, epsilon):
    m=np.arange(-n,n+1,dtype=float)
    diff=m[:,None]-m[None,:]
    numerator=np.sin(2*np.pi*epsilon*m[:,None])-np.sin(2*np.pi*epsilon*m[None,:])
    a=np.divide(numerator,np.pi*diff,out=np.zeros_like(diff),where=diff!=0)
    np.fill_diagonal(a,2*epsilon*np.cos(2*np.pi*epsilon*m))
    return a


def continuum_spectrum(t, order):
    x,w=roots_legendre(order)
    k=2*t*np.cos(np.pi*t*(x[:,None]+x[None,:]))*np.sinc(t*(x[:,None]-x[None,:]))
    return eigvalsh(np.sqrt(w[:,None]*w[None,:])*k,check_finite=False)


FUNCTIONS={
    'linear':(lambda x:x,1.,0.),
    'square':(lambda x:x*x,0.,2.),
    'cube':(lambda x:x**3,0.,6.),
    'fourth':(lambda x:x**4,0.,12.),
    'sine':(np.sin,1.,math.sin(1.)),
    'cosine_minus_one':(lambda x:-2*np.sin(x/2)**2,0.,1.),
    'exponential_minus_one':(np.expm1,1.,math.e),
}


def finite_checks():
    cases=[(1,.001),(2,.03),(4,.2),(5,.49),(9,.07),(12,.35),
           (20,.1),(25,.25),(40,.02),(50,.49)]
    rows=[]
    for k,(n,eps) in enumerate(cases,1):
        d=2*n+1;h=2/d;t=d*eps/2
        a=finite_matrix(n,eps);ev=eigvalsh(a)
        order=max(96,int(math.ceil(10*t)))
        low=continuum_spectrum(t,order);high=continuum_spectrum(t,order+48)
        delta=2*math.sin(math.pi*d*eps)*(eps/math.sin(math.pi*eps)-1/math.pi)
        # Independent feature integral: pair xi in (0,t) with xi-t.
        z,w=roots_legendre(max(96,int(8*t)));xi=t*(z+1)/2;w=t*w/2
        centres=np.arange(-n,n+1)*h
        f0=np.exp(2j*np.pi*centres[:,None]*xi[None,:])
        f1=np.exp(2j*np.pi*centres[:,None]*(t-xi)[None,:])
        sampled=2*h*np.real((f1*w)@f0.T)
        s0=np.sinc(h*xi);s1=np.sinc(h*(t-xi))
        averaged=2*h*np.real((f1*(w*s0*s1))@f0.T)
        feature_error=float(np.max(np.abs(a-sampled)))
        eta=4*float(np.sum(w*(1-s0*s0)))
        eta_bound=2*np.pi**2*d*eps**3/9
        trace_distance=float(np.sum(np.abs(eigvalsh(a-averaged))))
        require(feature_error<2e-11,'finite feature identity',n=n,error=feature_error)
        require(-1e-12<=eta<=eta_bound+1e-11,'cell leakage',n=n,eta=eta,bound=eta_bound)
        require(trace_distance<=1.5*eta_bound+2e-11,'sample versus compression',
                n=n,error=trace_distance,bound=1.5*eta_bound)
        checks=[]
        for name,(fun,linear,m2) in FUNCTIONS.items():
            observed=float(np.sum(fun(ev))-np.sum(fun(high)))
            refinement=float(abs(np.sum(fun(low))-np.sum(fun(high))))
            residual=abs(observed-linear*delta)
            bound=5*np.pi**2*d*eps**3*m2/9
            require(refinement<3e-10,'continuum refinement',n=n,function=name,error=refinement)
            require(residual<=bound+3e-10,'affine-centered trace bound',n=n,function=name,
                    residual=residual,bound=bound)
            checks.append({'function':name,'trace_difference':observed,
                           'linear_correction':linear*delta,'residual':residual,
                           'proved_upper_bound':bound,'quadrature_refinement':refinement})
        rows.append({'N':n,'epsilon':eps,'t':t,'orders':[order,order+48],
                     'feature_identity_error':feature_error,'cell_leakage':eta,
                     'cell_leakage_bound':eta_bound,'sampling_trace_distance':trace_distance,
                     'functional_checks':checks})
        print(f'finite checks {k}/{len(cases)}',file=sys.stderr,flush=True)
    return rows


def star_checks():
    theta,w=roots_legendre(384);theta=(theta+1)/2;w=w/2
    rows=[]
    for n in (1,2,5,12,30,60,120,240):
        modes=np.arange(n)
        c=np.cos(np.pi*theta/2);s=np.sin(np.pi*theta/2)
        expo=np.exp(-2j*np.pi*modes[:,None]*theta[None,:])
        tc=toeplitz(expo@(w*c));ts=toeplitz(expo@(w*s))
        hn=tc@tc+ts@ts
        edge=float(eigvalsh(hn)[0]);identity=float(np.linalg.norm(
            hn-.5*((tc+ts)@(tc+ts)+(tc-ts)@(tc-ts)),ord=2))
        require(edge>=.5-2e-12,'actual star squared gap',n=n,edge=edge)
        require(identity<2e-12,'polarization identity (implementation consistency)',n=n,error=identity)
        rows.append({'n':n,'minimum_squared_singular_value':edge,'polarization_consistency_error':identity})
    return rows


def coefficient_and_profiles():
    g=1/math.sqrt(2)
    j=quad(lambda s:1/(1+math.sqrt((1+s*s)/2)),0,1,epsabs=1e-13)[0]
    closed=math.sqrt(2)*math.log(1+math.sqrt(2))-math.log(2)
    require(abs(j-closed)<2e-14,'trace norm coefficient',integral=j,closed=closed)
    gamma=2*closed/math.pi**2
    rows=[]
    # Observations, not pass/fail assertions about an unknown asymptotic constant.
    for t in (4.,8.125,16.,32.,64.):
        order=max(120,int(10*t));a=continuum_spectrum(t,order);b=continuum_spectrum(t,order+48)
        norms=[float(np.sum(np.abs(v))) for v in (a,b)]
        require(abs(norms[1]-norms[0])<5e-9,'trace norm refinement',t=t,values=norms)
        signs={}
        for name,v in [('positive',b[b>0][::-1]),('negative',-b[b<0])]:
            ranks=[]
            for nu in (.1,.2,.4):
                index=int(math.floor(2*t-nu*math.log(t)))
                if 1<=index<=len(v):
                    prediction=math.sqrt((1+math.tanh(math.pi**2*nu/2)**2)/2)
                    ranks.append({'nu':nu,'one_based_rank':index,'eigenvalue':float(v[index-1]),
                                  'asymptotic_prediction':prediction})
            signs[name]={'count_above_0_5':int(np.sum(v>.5)),
                         'count_above_0_8':int(np.sum(v>.8)),'rank_observations':ranks}
        rows.append({'t':t,'orders':[order,order+48],'trace_norm':norms[1],
                     'trace_norm_refinement':abs(norms[1]-norms[0]),
                     'trace_norm_remainder':norms[1]-4*t+gamma*math.log(t),
                     'gap_mass_b_0_5':float(np.sum(np.abs(b[np.abs(b)<=.5]))),'signed':signs})
        print(f'profile t={t:g} complete',file=sys.stderr,flush=True)
    return {'trace_norm_integral':j,'closed_integral':closed,'logarithmic_deficit':gamma,
            'gap_edge':g,'observations':rows}


def main():
    start=time.monotonic()
    result={'started_at':stamp(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
            'scope':'Ordinary-precision diagnostics; profile observations are not asymptotic certificates.'}
    result['finite']=finite_checks();result['star']=star_checks()
    result['profiles']=coefficient_and_profiles()
    result['status']='PASS';result['elapsed_seconds']=time.monotonic()-start
    result['completed_at']=stamp()
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
