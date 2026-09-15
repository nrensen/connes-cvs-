"""Parity determinant and moment diagnostics; floating-point checks only."""
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
