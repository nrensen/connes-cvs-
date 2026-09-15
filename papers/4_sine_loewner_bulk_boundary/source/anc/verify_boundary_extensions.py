"""Independent finite/covariance/multifrequency falsification; floating diagnostics only."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name]='1'
import json, hashlib, sys, time, argparse
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import numpy as np
import scipy
from scipy.linalg import eigvalsh, eigh
from scipy.special import roots_legendre
import mpmath as mp

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path)
args=parser.parse_args()
DEST=args.output
def stamp(): return datetime.now(ZoneInfo('Asia/Jerusalem')).strftime('%Y-%m-%d %H:%M:%S %Z')
def norm1(a): return float(np.linalg.svd(a,compute_uv=False).sum())
def sqrtpos(a):
    e,v=eigh(a);return (v*np.sqrt(np.maximum(e,0)))@v.T
def kernel(tau,order,coeff=(1.,)):
    x,w=roots_legendre(order);s=x[:,None]+x[None,:];d=x[:,None]-x[None,:]
    a=sum(c*2*m*tau*np.cos(np.pi*m*tau*s)*np.sinc(m*tau*d) for m,c in enumerate(coeff,1))
    return a*np.sqrt(w[:,None]*w[None,:])
def add(row):
    row['timestamp']=stamp();data['checks'].append(row)
    if DEST is not None:
        tmp=DEST.with_suffix('.json.tmp');tmp.write_text(json.dumps(data,indent=2)+'\n');os.replace(tmp,DEST)
    print(stamp(),row['kind'],row.get('N',row.get('coeff',row.get('dimension',''))),flush=True)

if DEST is not None and DEST.exists(): raise SystemExit('Finalized result already exists; refusing overwrite')
data={'started':stamp(),'status':'ordinary floating-point and arbitrary-precision falsification, not proof',
      'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mp.__version__,'threads':1,'checks':[]}
t0=time.monotonic(); rng=np.random.default_rng(729301)
for n in (3,7,19,41):
    z=rng.normal(size=(n,n));c0=z@z.T;c0/=np.linalg.norm(c0,2)*1.2
    z=rng.normal(size=(n,n));c1=z@z.T;c1/=np.linalg.norm(c1,2)*1.2
    j=rng.normal(size=(n,n));j=(j+j.T)/2;j/=np.linalg.norm(j,2)
    s0,s1=sqrtpos(c0),sqrtpos(c1);a0=s0@j@s0;a1=s1@j@s1
    diff=c1-c0;bound=norm1(diff);ev0,ev1=eigvalsh(a0),eigvalsh(a1)
    # Noncommuting covariance path; exponential checks all powers at once.
    value=float(np.expm1(ev1).sum()-np.expm1(ev0).sum()-np.trace(j@diff))
    sqnorm=norm1(s1-s0);sqbound=float(np.sqrt(np.abs(eigvalsh(diff))).sum())
    assert abs(value)<=np.e*bound+1e-12 and sqnorm<=sqbound+1e-12
    add({'kind':'noncommuting covariance and square-root inequality','dimension':n,
         'affine_removed_exp_trace':value,'C2_bound':np.e*bound,'sqrt_nuclear':sqnorm,'sqrt_bound':sqbound})
for N in (8,32,128):
  for eps in (.015,.1,.3,.45):
    d=2*N+1;tau=d*eps/2;x=np.arange(-N,N+1,dtype=float)
    a=2*eps*np.cos(np.pi*eps*(x[:,None]+x[None,:]))*np.sinc(eps*(x[:,None]-x[None,:]))
    ea=eigvalsh(a);order=2*int(np.ceil(4*tau+64));e0=eigvalsh(kernel(tau,order));e1=eigvalsh(kernel(tau,order+64))
    delta=2*np.sin(np.pi*d*eps)*(eps/np.sin(np.pi*eps)-1/np.pi)
    L=2*eps/np.sin(2*np.pi*eps)-1/np.pi
    gamma=16*np.sqrt(2/np.pi)*eps/((1-2*eps)*np.sqrt(1-4*eps**2))
    tests={}
    for name,f,fp0,m2 in [('x2',lambda x:x*x,0.,2.),('x3',lambda x:x**3,0.,6.),
                         ('x4',lambda x:x**4,0.,12.),('exp',np.expm1,1.,np.e)]:
      err=float(f(ea).sum()-f(e1).sum()-fp0*delta);quad=float(abs(f(e0).sum()-f(e1).sum()))
      assert abs(err)<=L*m2+1e-9
      tests[name]={'affine_removed_error':err,'bound':L*m2,'quadrature_difference':quad}
    liperr=float(abs(ea).sum()-abs(e1).sum());assert abs(liperr)<=gamma+1e-9
    add({'kind':'direct finite-continuum comparison','N':N,'epsilon':eps,'tau':tau,'d_epsilon_cubed':d*eps**3,
         'orders':[order,order+64],'tests':tests,'abs_trace_error':liperr,'Lipschitz_bound':gamma})
# Exact algebra at independent high precision, without using the Toeplitz asymptotic.
for digits in (50,90):
  mp.mp.dps=digits
  c=list(map(mp.mpf,['0.7','-1.2','0.4','0.9']));M=len(c);t=mp.mpf('0.371')
  h0=mp.matrix(2*M+1);h1=mp.matrix(2*M+1)
  for m,cm in enumerate(c,1):
    for j in range(1,m+1):
      u,v=j+M,j-m+M;h0[u,v]+=cm;h0[v,u]+=cm;h1[u-1,v-1]+=cm;h1[v-1,u-1]+=cm
  ht=(1-t)*h0+t*h1
  cubic=sum((ht**3)[i,i] for i in range(2*M+1))
  resonance=sum(c[j-1]*c[k-1]*c[j+k-1] for j in range(1,M+1) for k in range(1,M+1) if j+k<=M)
  quadratic=sum((ht**2)[i,i] for i in range(2*M+1))
  predicted2=2*sum(m*cm**2 for m,cm in enumerate(c,1))-4*t*(1-t)*sum(cm**2 for cm in c)
  assert abs(cubic-6*t*(1-t)*resonance)<mp.mpf(10)**(-digits+4)
  add({'kind':'arbitrary-precision seam triangle identity','dps':digits,'cubic':mp.nstr(cubic,digits),
       'cubic_error':mp.nstr(cubic-6*t*(1-t)*resonance,digits),'quadratic_error':mp.nstr(quadratic-predicted2,digits)})
for coeff in ((1.,1.,1.),(0.,1.,1.),(.7,-1.2,.4)):
  M=len(coeff);res=sum(coeff[j-1]*coeff[k-1]*coeff[j+k-1] for j in range(1,M+1) for k in range(1,M+1) if j+k<=M)
  rows=[]
  for tau in (8.125,16.125,32.125):
    order=2*int(np.ceil(4*M*tau+64));vals=[]
    for q in (order,order+64):
      a=kernel(tau,q,coeff);a2=a@a
      vals.append((float(np.sum(a*a)),float(np.sum(a2*a.T))))
    rows.append({'tau':tau,'orders':[order,order+64],'traces':vals,
      'quadratic_residual':vals[1][0]-4*tau*sum(m*c*c for m,c in enumerate(coeff,1))+2/np.pi**2*sum(c*c for c in coeff)*np.log(tau),
      'cubic_residual':vals[1][1]-3*res/np.pi**2*np.log(tau)})
  add({'kind':'general sine-family physical moments','coeff':coeff,'cubic_log_coefficient':3*res/np.pi**2,'scales':rows})

for tau in (8.125,16.125,32.125):
  order=2*int(np.ceil(4*tau+80));rows=[]
  for q in (order,order+64):
    a=kernel(tau,q);h=q//2;aa=a[:h,:h];ab=a[:h,h:][:,::-1]
    ev,od=eigvalsh(aa+ab),eigvalsh(aa-ab)
    squares={str(w):float(np.sum(np.log1p(w*ev**2))-np.sum(np.log1p(w*od**2))) for w in (.25,1.,4.)}
    z=.5
    linear=float(np.log1p(-z*ev).sum()-np.log1p(-z*od).sum())
    rows.append({'order':q,'squared_log_ratios':squares,'linear_log_ratio':linear,
      'reflected_powers':{str(m):float((ev**m).sum()-(od**m).sum()) for m in range(1,9)}})
  discrepancy=max(abs(rows[0]['squared_log_ratios'][str(w)]-rows[1]['squared_log_ratios'][str(w)]) for w in (.25,1.,4.))
  assert discrepancy<1e-7
  add({'kind':'parity determinant and moment limits','tau':tau,'resolutions':rows,
       'quadrature_discrepancy':discrepancy,'predicted_squared_log_ratios':{str(w):float(np.log1p(w/2)) for w in (.25,1.,4.)},
       'predicted_linear_log_ratio':float(.5*np.log((1-z*z/2)*(1-z)/(1+z)))})

data['finished']=stamp();data['elapsed_seconds']=time.monotonic()-t0
data['status']='PASS'
data['interpretation']='Numerical diagnostics only; no asymptotic limit or interval enclosure is certified.'
if DEST is not None:
    tmp=DEST.with_suffix('.json.tmp');tmp.write_text(json.dumps(data,indent=2)+'\n');os.replace(tmp,DEST)
print(json.dumps(data,indent=2),flush=True)
