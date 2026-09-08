#!/usr/bin/env python3
"""Bounded independent physical-kernel and scalar-transform checks; not certification."""
from pathlib import Path
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import numpy as np
from scipy.special import roots_legendre,hyp2f1,gamma,digamma,expit
from scipy.integrate import quad
import scipy,subprocess,json,hashlib,time,platform
OUT=Path(os.environ.get('BOUNDARY_DIAGNOSTIC_OUTPUT','harmonic_homogenization.json')).resolve()
def stamp():return subprocess.check_output(['env','TZ=Asia/Jerusalem','date','+%Y-%m-%d %H:%M:%S %Z'],text=True).strip()
def poly(x,c):return np.polynomial.polynomial.polyval(x,c)
def hpoly(x,c):
 f=poly(x,c);g=f*np.log(x/(1-x))
 for k,ck in enumerate(c):
  for j in range(k):g-=ck*x**j/(k-j)
 return g/np.pi
# Endpoint-vanishing polynomials: independent analytic principal-value integral.
coeff=[np.array([0,0,1,-2,1.]),np.array([0,1,-3,2.]),np.array([0,0,2,-5,3.])]
def compute(length,order):
 x,w=roots_legendre(order);x=(x+1)/2;w=w/2
 f=np.array([poly(x,c) for c in coeff]);hf=np.array([hpoly(x,c) for c in coeff]);q=(f+1j*hf)/2;nq=(f-1j*hf)/2
 cf=np.array([q[1],q[0]+nq[2],nq[1]])
 phases=np.exp(2j*np.pi*length*np.arange(-1,2)[:,None]*x)
 u=np.sum(phases*f,axis=0);expected=np.sum(phases*cf,axis=0)
 # Exact physical kernel after dilation, including the diagonal continuously.
 xx=x[:,None];yy=x[None,:]
 kernel=2*length*np.cos(np.pi*length*(xx+yy))*np.sinc(length*(xx-yy))
 actual=kernel@(w*u)
 return {'length':length,'order':order,'norm_squared':float(np.sum(w*abs(u)**2)),'envelope_norm_squared':float(np.sum(w*np.sum(abs(f)**2,axis=0))),'intertwining_error':float(np.sqrt(np.sum(w*abs(actual-expected)**2)))}
def ell_logit(x,nu):
 if x>25: return (x+digamma(1)-digamma(nu))/np.pi
 if x < -25:return (1/nu+np.exp(x)/(nu+1))/np.pi
 return hyp2f1(1,nu,nu+1,expit(x))/(np.pi*nu)
def transform(nu,omega):
 def integrand(x):return ell_logit(x,nu)/(2*np.cosh(x/2))*np.exp(-1j*omega*x)/np.sqrt(2*np.pi)
 a,ea=quad(lambda x:integrand(x).real,-65,65,epsabs=2e-11,epsrel=2e-11,limit=500,points=[-25,0,25])
 b,eb=quad(lambda x:integrand(x).imag,-65,65,epsabs=2e-11,epsrel=2e-11,limit=500,points=[-25,0,25])
 direct=a+1j*b;z=.5+1j*omega;closed=gamma(nu)*gamma(z)/gamma(nu+z)/(np.sqrt(2*np.pi)*np.cosh(np.pi*omega))
 return {'nu':nu,'omega':omega,'direct':[a,b],'beta_formula':[closed.real,closed.imag],'difference':float(abs(direct-closed)),'quadrature_reported_error':ea+eb,'scope':'ordinary quadrature plus asymptotic tail truncation, not an interval enclosure'}
def main():
 assert not OUT.exists();record={'started':stamp(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'threads':1,'scope':'Independent finite physical kernel versus analytic polynomial Hilbert model, and logit Fourier integral versus beta expression; diagnostic, not a proof or certificate.','homogenization':[],'transforms':[]};start=time.monotonic()
 def save():
  tmp=OUT.with_suffix('.partial.tmp');tmp.write_text(json.dumps(record,indent=2)+'\n');os.replace(tmp,OUT.with_suffix('.partial.json'))
 for length in (8,16,32,64):
  cell=[compute(length,int(8*length+80)),compute(length,int(12*length+120))]
  assert abs(cell[0]['intertwining_error']-cell[1]['intertwining_error'])<2e-8
  record['homogenization'].append(cell);save();print(stamp(),'physical case',length,'error',cell[-1]['intertwining_error'],flush=True)
 for nu in (.6,.72):
  for omega in (0.,.5,2.):
   row=transform(nu,omega);record['transforms'].append(row);save();assert row['difference']<2e-9;print(stamp(),'transform',nu,omega,'difference',row['difference'],flush=True)
 errors=[x[-1]['intertwining_error'] for x in record['homogenization']];assert all(a>b for a,b in zip(errors,errors[1:]))
 record.update(status='PASS',finished=stamp(),elapsed_seconds=time.monotonic()-start);save();os.rename(OUT.with_suffix('.partial.json'),OUT);print(stamp(),'PASS',flush=True)
if __name__=='__main__':main()
