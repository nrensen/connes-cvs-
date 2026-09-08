#!/usr/bin/env python3
"""Exact-input interval gates for the physical residue functional.

Checks elementary global row/derivative bounds and the exact step-witness
functional by special-function primitives.  No quadrature or eigensolver.
The negative-window inverse certificate and analytic proof remain separate
inputs to the nonvanishing theorem; this file certifies their residue gates.
"""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import time
import flint
from flint import arb, acb, fmpq, ctx

WITNESS_SHA = '61946f85902ba95c0b695fb7c43e8e814b84c081918e7602c145b739dbfea57b'


def stamp():
    return subprocess.check_output(['env','TZ=Asia/Jerusalem','date','+%Y-%m-%d %H:%M:%S %Z'],text=True).strip()


def rat(a,b=1):
    return arb(fmpq(a,b))


def row_bounds():
    pi=arb.pi()
    lo,hi=rat(2,5),rat(11,25)
    amin=(1-2*hi*hi).sqrt();amax=(1-2*lo*lo).sqrt()
    cmin=((1-hi*hi)/2).sqrt()
    nulo=rat(1,2)+amin.atan()/pi;nuhi=rat(1,2)+amax.atan()/pi
    # Fixed rational bounds used in the analytic logarithmic envelopes.
    aa,cc=rat(39,50),rat(63,100)
    tp=2*hi/(aa*(1+aa*aa))
    cp=hi/(2*cc)
    kp=1/cc+hi*hi/(2*cc**3)
    pref=rat(7,10).gamma()/(2*pi)**rat(7,10)/(8*pi*pi*aa*cc)
    logder=rat(71,100)/pi*(rat(13,10)+(2*pi).log()) + rat(88,100)/aa**2 + rat(35,100)/cc + rat(3,2)*rat(71,100)
    jnormsq=rat(186656,1000)
    jprimenormsq=rat(2055792,1000)
    fnorm=pref*jnormsq.sqrt()
    fprime=4*rat(13,100)+pref*jprimenormsq.sqrt()
    c1=rat(5,2)*rat(47,10)+(3+rat(214,100))*rat(10,7) + rat(23,100)*(rat(7,10)*18+rat(17,10)*rat(15,4)+rat(21,10))
    c2=3*rat(47,10)+(rat(5,2)+rat(214,100))*rat(10,7) + rat(23,100)*(rat(7,10)*18+rat(17,10)*rat(15,4)+rat(21,10))
    checks={
        'a_lower':bool(amin>aa),'c_lower':bool(cmin>cc),
        'nu_range':bool(nulo>rat(7,10) and nuhi<rat(3,4)),
        'kappa_bound':bool(hi/cc<rat(7,10)),
        'tprime_bound':bool(tp<rat(71,100)),
        'nuprime_bound':bool(rat(71,100)/pi<rat(23,100)),
        'cprime_bound':bool(cp<rat(35,100)),
        'kappaprime_bound':bool(kp<2),
        'k2prime_bound':bool(2+rat(7,10)*rat(71,100)<rat(5,2)),
        'k3prime_bound':bool(2+rat(14,10)*rat(71,100)<3),
        'k4prime_bound':bool(3*rat(71,100)<rat(214,100)),
        'digamma_bounds':bool(rat(7,10).digamma()>-rat(13,10) and rat(3,4).digamma()<0),
        'log_two_bound':bool(rat(2).log()<rat(7,10)),
        'positive_ell_derivative':bool(rat(100,49)+pi*pi/6<rat(15,4)),
        'negative_ell_derivative':bool(rat(100,49)<rat(21,10)),
        'A_derivative':bool(16+pi*pi/6<18),
        'log_cross_moment':bool(2-pi*pi/6<rat(36,100)),
        'row_constant':bool(rat(7,10)*rat(47,10)+rat(27,10)*rat(10,7)<rat(36,5)),
        'row_derivative_constants':bool(c1<24 and c2<26),
        'prefactor_bound':bool(pref<rat(93,10000)),
        'prefactor_log_derivative':bool(logder<4),
        'global_row_norm':bool(fnorm<rat(13,100)),
        'global_row_derivative_norm':bool(fprime<1),
    }
    assert all(checks.values()),checks
    return {'checks':checks,'row_norm_upper':str(fnorm),'row_derivative_upper':str(fprime),
            'prefactor_upper':str(pref),'prefactor_log_derivative_upper':str(logder)}


def anchor(witness):
    alpha=arb(fmpq(witness['lambda_exact_dyadic']))
    a=(1-2*alpha*alpha).sqrt(); t=a.atan();nu=rat(1,2)+t/arb.pi()
    c=1/(2*t.cos());pi=arb.pi();eg=arb.const_euler()
    phase=lambda x:acb(0,x).exp()
    k=[acb(1),alpha/c*phase(t),acb(0,1)*alpha/c*phase(2*t),acb(0,1)*phase(3*t)]
    pref=acb(0,-1)*nu.gamma()*phase(-(pi/4+3*t/2))/(8*pi*pi*a*c*(2*pi)**nu)
    ecache={};dcache={}

    def ell(p,z):
        # The two flags certify the exact relations b-c=-1, a+b-c=0,
        # even though p is represented by an outward-rounded ball.
        return z.hypgeom_2f1(arb(1),p,p+1,bc=True,abc=True)/(pi*p)

    def e(zq):
        key=str(zq)
        if key not in ecache:
            if zq==0: value=arb(0)
            elif zq==1:value=-(nu.digamma()+eg)/(pi*(1-nu))
            else:
                z=arb(zq);value=(z*ell(nu,z)+(1-z).log()/pi)/(1-nu)
            ecache[key]=value
        return ecache[key]

    def d(sq):
        key=str(sq)
        if key not in dcache:
            if sq==0:value=((1-nu).digamma()+eg)/(1-nu)
            else:
                s=arb(sq)
                av=-pi/(1+s)*ell(1-nu,1/(1+s))
                value=((1+s)*av-s.log())/(1-nu)
            dcache[key]=value
        return dcache[key]

    def primitive(s):
        j1=k[1]*d(s)+pi*(-e(s)-k[2]*e(1-s)-k[3]*e(-s))
        j2=k[2]*d(1-s)+pi*(-e(s-1)-k[1]*e(s)-k[3]*e(1-s))
        return pref*j1,pref*j2

    mesh=[fmpq(x) for x in witness['positive_mesh_exact_dyadics']]
    assert mesh[0]==0 and mesh[-1]==1 and all(x<y for x,y in zip(mesh,mesh[1:]))
    n=len(mesh)-1
    coeff=[acb(arb(fmpq(x)),arb(fmpq(y))) for x,y in zip(
        witness['coeff_real_exact_dyadics'],witness['coeff_imag_exact_dyadics'])]
    assert len(coeff)==2*n
    norm=sum((z.real*z.real+z.imag*z.imag for z in coeff),arb(0))
    result=acb(0);previous=primitive(mesh[0])
    for i,(left,right) in enumerate(zip(mesh,mesh[1:])):
        current=primitive(right);width=arb(right-left)
        result+=(current[0]-previous[0])*coeff[n+i]/width.sqrt()
        result+=(current[1]-previous[1])*coeff[i]/width.sqrt()
        previous=current
        if (i+1)%64==0:print(stamp(),ctx.prec,'bits',i+1,'/',n,'exact residue cells',flush=True)
    checks={'frequency_norm_bound':bool(norm<rat(1001,1000)**2),
            'trial_center_lower':bool(alpha>rat(41635,100000)),
            'trial_center_upper':bool(alpha<rat(41637,100000)),
            'anchor_lower':bool(abs(result)>rat(69,1000))}
    assert all(checks.values()),(checks,str(result))
    return {'checks':checks,'frequency_norm_squared':str(norm),'functional_value':str(result),
            'functional_absolute_value':str(abs(result)),'exact_cells_per_half':n}


def transfer_margin():
    error=rat(311,50000);unorm=rat(1001,1000);step=rat(3,250)
    lower=rat(69,1000)-step*unorm-rat(13,100)/rat(81,200)*(error+step*unorm+error/rat(49,1000))
    checks={'projected_functional_lower':bool(lower>rat(103,10000)),
            'normalized_amplitude_lower':bool(lower/unorm>rat(1,100)),
            'frequency_residual_from_existing_certificate':bool(rat(38641,1000000000)<error*error)}
    assert all(checks.values()),checks
    return {'checks':checks,'projected_functional_lower':str(lower),
            'normalized_amplitude_lower':str(lower/unorm),
            'conditional_inputs':'Existing exact witness residual and a unique eigenvalue with outside spectrum distance at least49/1000 from a0.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--bits',type=int,nargs='+',default=[192,256])
    args=parser.parse_args()
    assert args.output is None or not args.output.exists(),'Output exists'
    raw=Path(__file__).with_name('boundary_gap_witness.json').read_bytes();assert hashlib.sha256(raw).hexdigest()==WITNESS_SHA
    witness=json.loads(raw)
    record={'started':stamp(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'witness_sha256':WITNESS_SHA,'python_flint':flint.__version__,
            'scope':'Rigorous residue row, anchor and conditional transfer gates; negative inverse certificate is a separate required input.',
            'runs':[]}
    start=time.monotonic()
    for bits in args.bits:
        ctx.prec=bits;ctx.threads=1
        run={'precision_bits':bits,'row_bounds':row_bounds(),'anchor':anchor(witness),'transfer':transfer_margin()}
        record['runs'].append(run)
        print(stamp(),bits,'bits: residue interval gates PASS',json.dumps(run),flush=True)
    record.update(finished=stamp(),elapsed_seconds=time.monotonic()-start,status='PASS_RESIDUE_GATES')
    if args.output:
        with args.output.open('x') as stream:json.dump(record,stream,indent=2);stream.write('\n')
    print(json.dumps(record,indent=2),flush=True)


if __name__=='__main__':main()
