"""Rigorous full-operator negative-window certificate from an exact step witness.

Only stdlib and python-flint are used. No eigenvalue solve, FFT, numerical
quadrature or floating-point bound enters this checker.
"""
from pathlib import Path
from fractions import Fraction
import json,hashlib,subprocess,time,argparse,importlib.util
from flint import arb,acb,acb_mat,fmpq,ctx

def stamp():return subprocess.check_output(['env','TZ=Asia/Jerusalem','date','+%Y-%m-%d %H:%M:%S %Z'],text=True).strip()
def ball(x):
    x=Fraction(x);return arb(fmpq(x.numerator,x.denominator))
def primitive(x):
    x=Fraction(x)
    if x==0:return arb(0)
    y=ball(x);return y*abs(y).log()
def square_norm(a):
    ans=arb(0)
    for i in range(a.nrows()):
      for j in range(a.ncols()):
        z=a[i,j];ans+=z.real*z.real+z.imag*z.imag
    return ans
def real_inner(a,b):
    ans=arb(0)
    for i in range(a.nrows()):
      for j in range(a.ncols()):
        z,w=a[i,j],b[i,j];ans+=z.real*w.real+z.imag*w.imag
    return ans
def model_matrices(tr,emit):
    edge=[Fraction(x) for x in tr];width=[ball(b-a) for a,b in zip(edge,edge[1:])];sq=[x.sqrt() for x in width]
    n=len(width);h=acb_mat(n,n);c=acb_mat(n,n);pi=arb.pi()
    for i in range(n):
      for j in range(i+1):
        a,b,u,v=edge[i],edge[i+1],edge[j],edge[j+1]
        cc=(primitive(b+v)-primitive(a+v)-primitive(b+u)+primitive(a+u))/(pi*sq[i]*sq[j])
        c[i,j]=cc;c[j,i]=cc
        if i!=j:
          hh=(primitive(b-u)-primitive(a-u)-primitive(b-v)+primitive(a-v))/(pi*sq[i]*sq[j])
          h[i,j]=hh;h[j,i]=-hh
      if (i+1)%64==0:emit('exact model rows',i+1,'/',n)
    return h,c,width
def midpoint_defect(sr,tr):
    nt=len(tr)-1;ns=len(sr)-1;d=[acb_mat(nt,2*ns) for _ in range(4)];c=arb.pi()/2
    def chi(x):return (c*ball(x)).cos() if x<1 else arb(0)
    def D(t,s):return -c*(c*ball(t)).sin() if t==s and t<1 else (chi(t)-chi(s))/ball(t-s)
    def G(t,s):return (chi(t)-chi(s))/ball(t+s)
    def F(t,s):return chi(t)/ball(1-t+s) if t<1 else arb(0)
    for i,(t0,t1) in enumerate(zip(tr,tr[1:])):
      t=(t0+t1)/2
      for j,(s0,s1) in enumerate(zip(sr,sr[1:])):
        s=(s0+s1)/2;weight=(ball(t1-t0)*ball(s1-s0)).sqrt()/(2*arb.pi())
        rows=((D(t,s),G(t,1-s)),(-F(t,s),D(t,s)),(-D(t,1-s),F(t,1-s)),(-G(t,s),-D(t,1-s)))
        for ch in range(4):
          for side in range(2):d[ch][i,j+side*ns]=acb(0,weight*rows[ch][side])
    return d
def projected_norm_squared(sr,tr,w):
    ns=len(sr)-1;nt=len(tr)-1;cols=2*ns
    assert tr[:ns+1]==sr
    answer=arb(0);pi=arb.pi()
    for i,(a,b) in enumerate(zip(sr,sr[1:])):
      ht=ball(b-a);ca,cb=pi*ball(a),pi*ball(b)
      wc=ht/2+(cb.sin()-ca.sin())/(2*pi);wi=ht-wc;cross=(ca.cos()-cb.cos())/(2*pi)
      assert sr[ns-i]==1-a and sr[ns-1-i]==1-b
      for ch in (0,1):
        for j in range(cols):
          x=w[ch][i,j];y=w[ch+2][ns-1-i,j]
          xx=x.real*x.real+x.imag*x.imag;yy=y.real*y.real+y.imag*y.imag
          xy=x.real*y.real+x.imag*y.imag
          answer+=(wc*xx+wi*yy+2*cross*xy)/ht
    return answer
def certify(data,bits,bound_function,emit):
    ctx.prec=bits;ctx.threads=1;start=time.monotonic()
    sr=[Fraction(x) for x in data['source_edges']];tr=[Fraction(x) for x in data['model_edges']]
    assert all(x+y==1 for x,y in zip(sr,sr[::-1]))
    assert all(a<b for a,b in zip(sr,sr[1:])) and all(a<b for a,b in zip(tr,tr[1:]))
    ns=len(sr)-1;nt=len(tr)-1;nc=2*ns;den=1<<data['coefficient_denominator_power']
    assert data['shape']==[4*nt,nc]
    re,im=data['W_real_integers'],data['W_imag_integers'];assert len(re)==len(im)==4*nt*nc
    w=[acb_mat(nt,nc) for _ in range(4)]
    for ch in range(4):
      for i in range(nt):
        for j in range(nc):
          k=(ch*nt+i)*nc+j;w[ch][i,j]=acb(arb(fmpq(re[k],den)),arb(fmpq(im[k],den)))
    emit('loaded exact witness',bits,'bits')
    h,c,_=model_matrices(tr,emit);d=midpoint_defect(sr,tr)
    emit('exact model action')
    hw=[h*x for x in w];cw=[c*x for x in w]
    iu=acb(0,arb(1)/2);half=arb(1)/2
    qw=[w[0]*half+hw[0]*iu+cw[3]*iu,w[1]*half+hw[1]*iu,
        w[2]*half-hw[2]*iu,w[3]*half-hw[3]*iu-cw[0]*iu]
    assert Fraction(data['a'])==Fraction(-833,2000)
    a=ball(data['a']);exchange=(1,0,3,2)
    normw=sum((square_norm(x) for x in w),arb(0));normd=sum((square_norm(x) for x in d),arb(0))
    qenergy=sum((real_inner(w[j],qw[j]) for j in range(4)),arb(0))
    qcross=sum((real_inner(qw[j],w[exchange[j]]) for j in range(4)),arb(0))
    carleman=square_norm(cw[1])+square_norm(cw[2])
    mw_upper=qenergy-carleman/4-2*a*qcross+a*a*normw
    cross=sum((real_inner(qw[j]-w[exchange[j]]*a,d[j]) for j in range(4)),arb(0))
    residual_squared=mw_upper-2*cross+normd
    projected_squared=projected_norm_squared(sr,tr,w)
    assert residual_squared>0 and projected_squared>0
    err=bound_function(sr,tr)['hs_bound']
    k_bound=projected_squared.sqrt()+5*(residual_squared.sqrt()+err)
    checks={'full_correction_norm_below_three_fifths':bool(k_bound<arb(3)/5),
            'model_residual_below_9_over_500':bool(residual_squared.sqrt()<arb(9)/500),
            'projected_trial_HS_below_27_over_100':bool(projected_squared.sqrt()<arb(27)/100),
            'global_defect_error_below_37_over_1000':bool(err<arb(37)/1000)}
    assert all(checks.values()),checks
    # From k<3/5, ||P_a^-1||<25/2. Thus |lambda-a|<=3/50
    # has Neumann product strictly below 3/4; the entire closed interval is excluded.
    inverse_bound=5/(1-k_bound)
    return {'precision_bits':bits,'model_cells':nt,'source_cells_per_half':ns,'model_trial_HS_squared':str(normw),
      'midpoint_defect_HS_squared':str(normd),'model_residual_squared_upper':str(residual_squared),
      'model_residual_HS_upper':str(residual_squared.sqrt()),'projected_trial_HS_squared':str(projected_squared),
      'projected_trial_HS_upper':str(projected_squared.sqrt()),'global_midpoint_defect_HS_upper':str(err),
      'full_correction_norm_upper':str(k_bound),'full_pencil_inverse_norm_upper':str(inverse_bound),
      'certified_closed_negative_interval':['-953/2000','-713/2000'],'checks':checks,'elapsed_seconds':time.monotonic()-start}
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--witness',type=Path);parser.add_argument('--output',type=Path);parser.add_argument('--bits',type=int,nargs='+',default=[128,192]);args=parser.parse_args()
    assert args.output is None or not args.output.exists();here=Path(__file__).resolve();dep=here.with_name('parametrix_defect_step_bound.py')
    expected='f6f952f45b4eeb89cd73f2dbabf11044295bf28d5a69c9ec9d0641b871b20200';assert hashlib.sha256(dep.read_bytes()).hexdigest()==expected
    spec=importlib.util.spec_from_file_location('exact_defect_bound',dep);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    fixture=args.witness or here.with_name('step_parametrix_witness.json')
    data=json.loads(fixture.read_text());record={'started':stamp(),'source_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),'bound_source_sha256':expected,'witness_sha256':hashlib.sha256(fixture.read_bytes()).hexdigest(),'scope':'Full infinite-dimensional negative-pencil exclusion and consequent positive-root uniqueness; not a finite-matrix spectral count.','runs':[]}
    def emit(*args):print(stamp(),*args,flush=True)
    for bits in args.bits:
      result=certify(data,bits,module.cosine_defect_midpoint_bound,emit);record['runs'].append(result);emit('certificate PASS',result)
    record.update(finished=stamp(),status='PASS')
    if args.output:
      with args.output.open('x') as f:json.dump(record,f,indent=2);f.write('\n')
    print(json.dumps(record,indent=2),flush=True)
if __name__=='__main__':main()
