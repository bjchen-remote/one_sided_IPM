"""A bounded exact fourth-order gate, reusing the saved first three orders.

Uses polynomial-domain arithmetic to avoid repeated general expression
simplification. Stops after a short local wall-clock budget, saving the stage.
This is not a long-running numerical job.
"""
import json
from pathlib import Path
import signal
import time
import argparse

import sympy as sp
from sympy.polys.fields import field
from normal_lift_algebra import curl_coefficients

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/ridge_density_closure_20260928'


class BudgetExpired(Exception):
    pass


def timeout(signum,frame):
    raise BudgetExpired('45-second bounded local symbolic check expired')


def main(order=4):
    signal.signal(signal.SIGALRM,timeout)
    signal.alarm(45)
    start=time.monotonic()
    result=dict(stage='load',status='running',order=order,scope=f'lambda=3, mu=1, v=1, exact cubic wall; normal coefficient {order} only')
    try:
        r,L=sp.symbols('r L',real=True)
        F,rr,ll=field((r,L),sp.QQ)
        J=1+rr*rr;s=rr+1
        def D(f):
            f=F(f)
            return f.diff(rr)+f.diff(ll)/J
        old=json.loads((ROOT/'runs/ridge_symbolic_20260928/elementary_lift.json').read_text())
        aa=[F.from_expr(sp.sympify(row['coefficient'],locals={'r':r,'L':L})) for row in old['rows']]
        for previous in range(4,order):
            filename='fourth_elementary_gate.json' if previous==4 else f'elementary_gate_{previous}.json'
            stored=json.loads((OUT/filename).read_text())
            assert stored['status']=='closed_in_same_ring'
            aa.append(F.from_expr(sp.sympify(stored['coefficient'],locals={'r':r,'L':L})))
        assert len(aa)==order-1
        result['stage']='source'
        j=[J]+[D(a) for a in aa];k=[(i+1)*a for i,a in enumerate(aa)]
        E=curl_coefficients(j,[D(a) for a in j],[D(D(a)) for a in j],
            k,[D(a) for a in k],[(i+1)*k[i+1] for i in range(len(k)-1)],s,F(3),F(1),order-1,F)[-1]
        source=E/(order*J**3)
        gamma=5*order-1
        integrand=s**(gamma-1)*source/J**(order+1)
        result['sourceSeconds']=time.monotonic()-start
        print(f'order {order} source generated',result['sourceSeconds'],flush=True)
        # Preserve the exact source even if the following integration needs more time.
        expr=integrand.as_expr().subs({sp.Symbol('r'):r,sp.Symbol('L'):L})
        result['integrand']=str(expr)
        pending=sp.Poly(expr,L)
        primitive=sp.Integer(0)
        stages=[]
        for degree in range(pending.degree(),-1,-1):
            result['stage']=f'integrate_L_power_{degree}'
            coef=sp.cancel(pending.nth(degree))
            anti=sp.expand(sp.integrate(coef,r)).subs(sp.atan(r),L)
            if anti.has(sp.log,sp.Integral,sp.atan):
                result.update(status='new_function_required_by_this_integration_test',
                              obstructedPower=int(degree),coefficient=str(coef),integral=str(anti))
                break
            pp=sp.Poly(anti,L)
            rat=sp.cancel(pp.nth(0));angle=sp.cancel(pp.nth(1))
            assert pp.degree()<=1 and not angle.has(r)
            term=rat*L**degree+angle*L**(degree+1)/sp.Integer(degree+1)
            primitive+=term
            derivative=sp.diff(term,r)+sp.diff(term,L)/(1+r*r)
            pending=sp.Poly(sp.cancel(pending.as_expr()-derivative),L)
            assert pending.nth(degree)==0
            stages.append(dict(power=int(degree),angleCoefficient=str(angle)))
            result['integrationStages']=stages
            print('integrated angle power',degree,flush=True)
        else:
            assert pending.as_expr()==0
            result['stage']='verify'
            primitive=sp.cancel(primitive)
            base=sp.cancel(primitive.subs({r:-1,L:0}))
            coeff=sp.cancel((1+r*r)**(order+1)/(r+1)**gamma*(primitive-base))
            a4=F.from_expr(coeff)
            defect=s*D(a4)+(F(gamma)-(order+1)*s*2*rr/J)*a4-source
            assert defect==0
            result.update(status='closed_in_same_ring',coefficient=str(coeff),primitive=str(primitive),
                          angleDegree=int(sp.Poly(coeff,L).degree()),odeResidual=0)
    except BudgetExpired as exc:
        result.update(status='bounded_check_timed_out',message=str(exc))
    finally:
        signal.alarm(0)
        result['seconds']=time.monotonic()-start
        OUT.mkdir(parents=True,exist_ok=True)
        filename='fourth_elementary_gate.json' if order==4 else f'elementary_gate_{order}.json'
        (OUT/filename).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:value for key,value in result.items() if key not in ['integrand','coefficient','primitive','integral']},indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--order',type=int,default=4)
    main(parser.parse_args().order)
