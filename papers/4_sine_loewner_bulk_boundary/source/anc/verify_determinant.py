#!/usr/bin/env python3
"""Independent determinant/coefficient diagnostics; not interval certificates."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.linalg import eigvalsh
from scipy.special import roots_legendre


def stamp():
    return subprocess.check_output(['date','+%Y-%m-%d %H:%M:%S %Z'],
        env={**os.environ,'TZ':'Asia/Jerusalem'},text=True).strip()


def coefficient(z):
    return 2/math.pi**2*math.atanh(z/math.sqrt(2-z*z))**2


def logdet(t,phi,z,order):
    x,w=roots_legendre(order)
    k=2*t*np.cos(np.pi*t*(x[:,None]+x[None,:])+phi)*np.sinc(t*(x[:,None]-x[None,:]))
    a=np.sqrt(w[:,None]*w[None,:])*k
    ev=eigvalsh(a,check_finite=False)
    spectral=float(np.log1p(-z*ev).sum())
    sign,direct=np.linalg.slogdet(np.eye(order)-z*a)
    assert sign==1 and abs(spectral-direct)<5e-11,(t,phi,z,order,spectral,direct)
    return spectral,abs(spectral-direct)


def main():
    record={'started_at':stamp(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
            'scope':'Ordinary-precision identities and refinement; bounded residual observations do not prove an asymptotic theorem.'}
    rows=[]
    for z in (-.95,-.8,-.3,0.,.001,.3,.8,.95):
        a=2*z*z/(1-z*z)
        val=quad(lambda t:math.log1p(a*t*(1-t))/(t*(1-t)),0,1,epsabs=1e-13,epsrel=1e-13)[0]/(2*math.pi**2)
        closed=coefficient(z)
        assert abs(val-closed)<3e-13,(z,val,closed)
        rows.append({'z':z,'integral':val,'closed':closed,'discrepancy':abs(val-closed)})
    record['coefficients']=rows
    small=coefficient(1e-4)/1e-8
    assert abs(small-1/math.pi**2)<2e-9,small
    record['quadratic_coefficient_limit']={'at_z':1e-4,'ratio':small,'limit':1/math.pi**2}
    # square form: -(4/pi^2) arcsin^2( sqrt(2w/(1+w))/2 ) equals twice the coefficient at z=i sqrt(w)
    rows=[]
    for w in (.05,.25,.5,1.,2.,5.):
        square=-4/math.pi**2*math.asin(math.sqrt(2*w/(1+w))/2)**2
        viaartanh=-4/math.pi**2*math.atan(math.sqrt(w/(2+w)))**2
        direct=quad(lambda t:math.log((1+w*(1-2*t*(1-t)))/(1+w))/(t*(1-t)),0,1,
                    epsabs=1e-13,epsrel=1e-13)[0]/math.pi**2
        assert abs(square-viaartanh)<1e-14 and abs(square-direct)<3e-12,(w,square,viaartanh,direct)
        rows.append({'w':w,'arcsin_form':square,'arctan_form':viaartanh,'direct_integral':direct})
    assert abs(rows[3]['arcsin_form']+1/9)<1e-14,rows[3]
    record['square_coefficients']=rows
    rows=[]
    for t in (2.,4.125,8.,16.):
        for phi in (0.,.7,math.pi/2):
            for z in (-.8,.3,.8):
                order=max(96,int(10*t))
                low,d0=logdet(t,phi,z,order);high,d1=logdet(t,phi,z,order+48)
                error=abs(high-low)
                assert error<2e-10,(t,phi,z,error)
                rows.append({'tau':t,'phase':phi,'z':z,'orders':[order,order+48],
                    'logdet':high,'refinement':error,'matrix_identity_error':max(d0,d1),
                    'centered_remainder':high-2*t*math.log1p(-z*z)-coefficient(z)*math.log(t)})
        print(f'determinants tau={t:g} complete',file=sys.stderr,flush=True)
    record['determinants']=rows;record['status']='PASS';record['completed_at']=stamp()
    print(json.dumps(record,indent=2))


if __name__=='__main__':
    main()
