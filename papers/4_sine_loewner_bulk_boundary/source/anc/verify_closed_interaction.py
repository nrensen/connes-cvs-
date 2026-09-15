#!/usr/bin/env python3
"""Independent scalar integrals and matrix algebra; floating diagnostics only."""
from pathlib import Path
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import mpmath as mp
import subprocess,json,hashlib,time,platform
OUT=Path(os.environ.get('BOUNDARY_DIAGNOSTIC_OUTPUT','closed_interaction.json')).resolve()
def stamp():
    return subprocess.check_output(['env','TZ=Asia/Jerusalem','date','+%Y-%m-%d %H:%M:%S %Z'],text=True).strip()
def ell_reflected(s,nu):
    if s < mp.mpf('.01'):
        # Convergent incomplete-beta series, with a precision-sized truncation.
        term=mp.mpf(1);series=mp.mpf(0)
        for k in range(1,mp.mp.dps+20):
            term *= (k-nu)*s/k
            series += term/k
        return (1-s)**(-nu)*(-mp.log(s)-mp.euler-mp.digamma(nu)-series)/mp.pi
    return mp.hyp2f1(1,nu,nu+1,1-s)/(mp.pi*nu)
def row(alpha_text,dps):
    mp.mp.dps=dps
    alpha=mp.mpf(alpha_text);a=mp.sqrt(1-2*alpha**2)
    beta=mp.atan(a)/mp.pi;nu=mp.mpf('.5')+beta
    def lminus(w):
        return mp.beta(nu,mp.mpf('.5')-1j*w)/(mp.sqrt(2*mp.pi)*mp.cosh(mp.pi*w))
    # Independent real-space direct term and frequency-space inverse term.
    j=mp.quad(lambda t:4*t**(3-4*nu)*ell_reflected(t**4,nu),[0,mp.mpf('.1'),mp.mpf('.5'),1])
    integral=2*mp.quad(lambda w:mp.re(lminus(w)**2/(a-1j*mp.tanh(mp.pi*w))),[0,mp.mpf('.25'),1,3,8,mp.inf])
    closed=mp.pi/mp.sin(mp.pi*nu)**2
    errors=[]
    S=mp.matrix([[0,1,0],[1,0,-1],[0,-1,0]])
    A=mp.matrix([[0,1,0],[1,0,1],[0,1,0]])
    q=(1-1j*a)/2;r=mp.conj(q)
    v=mp.matrix([r/alpha,1,q/alpha]);u=S*v
    errors.append(abs((v.H*(A/2-alpha*mp.eye(3))*v)[0]))
    for w in (mp.mpf('-1.7'),mp.mpf('-.2'),mp.mpf(0),mp.mpf('.9')):
        x=mp.tanh(mp.pi*w);qw=(1+x)/2
        C=mp.matrix([[0,qw,0],[qw,0,1-qw],[0,1-qw,0]])
        # Direct arbitrary-precision matrix inversion, not the proposed polynomial.
        R=(C-alpha*mp.eye(3))**-1
        errors.append(abs((u.H*R*u)[0]))
        errors.append(abs((u.H*R*u.conjugate())[0]+4*a/(alpha*(a-1j*x))))
    discrepancy=abs(j+integral-closed)
    assert discrepancy < mp.mpf(10)**(-min(28,dps-12)), mp.nstr(discrepancy,12)
    assert max(errors)<mp.mpf(10)**(-dps+8)
    return {k:mp.nstr(v,dps) if isinstance(v,(mp.mpf,mp.mpc)) else v for k,v in dict(alpha=alpha_text,dps=dps,beta=beta,nu=nu,direct_integral=j,complement_integral=integral,sum=j+integral,closed=closed,absolute_discrepancy=discrepancy,max_matrix_identity_error=max(errors),universal_offdiagonal=-a*closed/alpha).items()}
def main():
    assert not OUT.exists()
    record=dict(started=stamp(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),python=platform.python_version(),mpmath=mp.__version__,threads=1,scope='Floating mpmath quadratures at two precisions plus direct 3x3 matrix inversions; neither interval enclosure nor numerical proof of actual boundary residues.',rows=[])
    started=time.monotonic()
    def checkpoint():
        tmp=OUT.with_suffix('.tmp');tmp.write_text(json.dumps(record,indent=2)+'\n');os.replace(tmp,OUT.with_suffix('.partial.json'))
    for dps in (40,65):
        for alpha in ('0.2','0.416','0.6','-0.416'):
            item=row(alpha,dps);record['rows'].append(item);checkpoint()
            print(stamp(),alpha,dps,'error',item['absolute_discrepancy'],flush=True)
    for i in range(4):
        mp.mp.dps=65
        assert abs(mp.mpf(record['rows'][i]['sum'])-mp.mpf(record['rows'][i+4]['sum']))<mp.mpf('1e-28')
    record.update(status='PASS',finished=stamp(),elapsed_seconds=time.monotonic()-started);checkpoint()
    os.rename(OUT.with_suffix('.partial.json'),OUT)
    print(stamp(),'PASS',flush=True)
if __name__=='__main__':main()
