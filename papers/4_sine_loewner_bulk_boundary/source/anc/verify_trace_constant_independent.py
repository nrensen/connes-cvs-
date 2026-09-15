#!/usr/bin/env python3
"""Independent exact scalar verification for the new admissible trace bound.

The operator inequalities are proved separately. This checks the product-rule
combinatorics and arithmetic without importing the proposing implementation.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import argparse
import hashlib
import json


def atan_bounds(x, n):
    s=sum(((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(n)),F(0))
    next_term=(-1)**n*x**(2*n+1)/F(2*n+1)
    return min(s,s+next_term),max(s,s+next_term)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    a=atan_bounds(F(1,5),24);b=atan_bounds(F(1,239),8)
    pi_lo=16*a[0]-4*b[1];pi_hi=16*a[1]-4*b[0]
    checks={'machin_pi_enclosure':F(31,10)<pi_lo<pi_hi<F(22,7)}
    derived=[]
    for n in (1,2,3):
        coefficients=[F(0),F(0),F(0)]
        for i in range(n+1):
            for j in range(n-i+1):
                k=n-i-j
                m=F(factorial(n),factorial(i)*factorial(j)*factorial(k))
                if i:
                    coefficients[0]+=m*2**(i-1)*F(1,2)**k
                elif j:
                    coefficients[1]+=m*F(1,2)**k
                else:
                    coefficients[2]+=m*F(1,2)**(k-1)
        derived.append(coefficients)
    expected=[[F(1),F(1),F(1)],[F(5),F(2),F(1,2)],[F(79,4),F(13,4),F(1,4)]]
    checks['independent_multinomial_derivatives']=derived==expected
    p=F(22,7);k=F(67,10);ell=F(15)
    checks['rank_four_hs_trace_bound']=4*(4*p*p/3-2)<k*k
    checks['infinite_rank_log_bound']=p*(k+2*F(99,70))/2<ell
    v1=k+ell+p
    v2=5*p*k+2*p*ell+p*p/2
    v3=F(79,4)*p*p*k+F(13,4)*p*p*ell+p**3/4
    cv=v1/4+v2/(6*p)+v3/(12*p)
    checks['seam_below_65']=cv<65
    lower=F(31,10)
    scalar=4*F(191,100)+2+2/lower**2+F(3,5)+F(19,10)+F(1,2)
    checks['trace_143']=2*cv+scalar<143
    band=130+4*F(191,100)+2+2/lower**2+F(99,70)*(F(3,5)+F(19,10)+F(1,4))
    checks['band_144']=band<144
    checks['plateau_138']=130+4*F(191,100)<138
    checks['phase_increment_below_8']=4*F(191,100)<8
    assert all(checks.values()),checks
    stamp=datetime.now(ZoneInfo('Asia/Jerusalem')).strftime('%Y-%m-%d %H:%M:%S %Z')
    record={'recorded':stamp,'status':'PASS_INDEPENDENT_EXACT_ARITHMETIC',
            'scope':'Operator proof is separate; this checks scalar arithmetic and differentiated-product coefficients.',
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'checks':checks,'derived_coefficients':[[str(x) for x in row] for row in derived],
            'gauge_upper':str(cv),'trace_upper':str(2*cv+scalar),'band_upper':str(band)}
    with args.output.open('x') as f: json.dump(record,f,indent=2);f.write('\n')
    print(json.dumps(record))

if __name__=='__main__':main()
