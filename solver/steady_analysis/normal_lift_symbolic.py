"""Equation-generated normal jets, exact cubic benchmark and resummation gates."""
import argparse
import json
from pathlib import Path
import sympy as sp
from normal_lift_algebra import curl_coefficients


def cubic_jets(count,k):
    s=sp.symbols('s',real=True)
    J=1+s*s;mu=1/(4*k);lam=3*mu
    coeff=[];results=[]
    for n in range(1,count+1):
        j=[J]+[sp.diff(a,s) for a in coeff]
        js=[sp.diff(v,s) for v in j]
        jss=[sp.diff(v,s) for v in js]
        kk=[(i+1)*a for i,a in enumerate(coeff)]
        ks=[sp.diff(v,s) for v in kk]
        ky=[(i+1)*kk[i+1] for i in range(len(kk)-1)]
        E=curl_coefficients(j,js,jss,kk,ks,ky,s,lam,mu,n-1,sp.cancel)[-1]
        source=sp.cancel(E/(n*J**3))
        rhs=sp.cancel(source*J**(2*n-1)/mu)
        assert sp.denom(rhs)==1 or not sp.denom(rhs).has(s)
        rhs=sp.Poly(rhs,s)
        assert rhs.degree()<=n-1
        P=sp.Integer(0);previous={}
        for i in range(n):
            previous[i]=sp.cancel((rhs.nth(i)-(i-n-1)*previous.get(i-2,0))/(i+5*n-1))
            P+=previous[i]*s**i
        gate=sp.expand(s*J*sp.diff(P,s)+(5*n-1+(1-n)*s*s)*P-rhs.as_expr())
        assert gate==0
        aa=sp.factor(P/J**(2*n-2))
        coeff.append(aa)
        results.append({'n':n,'a_n':str(aa),'numerator_degree':int(sp.degree(P,s)),
                        'ODE_residual':'0','source_degree':int(rhs.degree())})
        print(f'exact normal coefficient n={n} verified',flush=True)
    if count>=4:
        pade_gate=sp.factor(coeff[1]*coeff[3]-coeff[2]**2)
        binomial_gate=sp.factor(16*coeff[1]*coeff[3]-21*coeff[2]**2)
    else:pade_gate=binomial_gate=None
    return {'k':str(k),'rates':{'lambda':str(lam),'mu':str(mu)},'coefficients':results,
            'one_pole_curvature_gate_a2a4_minus_a3_squared':str(pade_gate),
            'height_power_2over3_gate_16a2a4_minus_21a3_squared':str(binomial_gate),
            'scope':'Cubic symmetric benchmark; finite-order exact generation, not an all-order convergence or summation theorem.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--order',type=int,default=7)
    parser.add_argument('--symbolic-k',action='store_true');args=parser.parse_args()
    out=Path(__file__).resolve().parents[2]/'runs/normal_lift_20260928';out.mkdir(exist_ok=True)
    k=sp.symbols('k',positive=True) if args.symbolic_k else sp.Integer(1)
    result=cubic_jets(args.order,k)
    name='cubic_symbolic_k.json' if args.symbolic_k else 'cubic_exact_jets.json'
    (out/name).write_text(json.dumps(result,indent=2)+'\n')
