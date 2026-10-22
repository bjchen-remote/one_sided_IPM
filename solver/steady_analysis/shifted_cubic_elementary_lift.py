"""Bounded exact test of a rational-plus-arctangent normal-lift algebra.

The benchmark lambda=3, mu=1 has d=1 by default; --symbolic-d retains d as
a parameter. Finite-order success is not all-order closure or a closed sum.
"""
import argparse
import json
from pathlib import Path
import time

import sympy as sp
from normal_lift_algebra import curl_coefficients

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/ridge_symbolic_20260928'


def run(order, symbolic_d=False):
    r,L=sp.symbols('r L',real=True)
    J=1+r*r
    d=sp.symbols('d',real=True) if symbolic_d else sp.Integer(1)
    shifted=r+d
    def deriv(f):
        return sp.cancel(sp.diff(f,r)+sp.diff(f,L)/J)
    def integrate_ring(f):
        """Triangular Hermite integration of a polynomial in L, L'=1/J."""
        pending=sp.Poly(sp.cancel(f),L)
        primitive=sp.Integer(0)
        stages=[]
        for degree in range(pending.degree(),-1,-1):
            coeff=sp.cancel(pending.nth(degree))
            anti=sp.integrate(coeff,r)
            anti=sp.expand(anti).subs(sp.atan(r),L)
            if anti.has(sp.log,sp.Integral,sp.atan):
                return None,dict(failedDegree=degree,coefficient=str(coeff),
                                 rationalIntegral=str(anti),stages=stages)
            poly=sp.Poly(anti,L)
            if poly.degree()>1:
                raise AssertionError('Unexpected rational-integration output')
            rational=sp.cancel(poly.nth(0))
            angle=sp.cancel(poly.nth(1))
            assert not angle.has(r)
            assert sp.cancel(sp.diff(rational,r)+angle/J-coeff)==0
            term=rational*L**degree+angle*L**(degree+1)/sp.Integer(degree+1)
            primitive+=term
            pending=sp.Poly(sp.cancel(pending.as_expr()-deriv(term)),L)
            assert pending.nth(degree)==0
            stages.append(dict(degree=int(degree),addedAngleConstant=str(angle)))
        assert pending.as_expr()==0
        return sp.cancel(primitive),dict(stages=stages)

    coefficients=[]
    rows=[]
    for n in range(1,order+1):
        start=time.monotonic()
        js0=[J]+[deriv(x) for x in coefficients]
        js1=[deriv(x) for x in js0]
        js2=[deriv(x) for x in js1]
        kk=[(i+1)*x for i,x in enumerate(coefficients)]
        ks=[deriv(x) for x in kk]
        ky=[(i+1)*kk[i+1] for i in range(len(kk)-1)]
        E=curl_coefficients(js0,js1,js2,kk,ks,ky,shifted,sp.Integer(3),
                            sp.Integer(1),n-1,sp.cancel)[-1]
        source=sp.cancel(E/(n*J**3))
        gamma=5*n-1
        integrand=sp.cancel(shifted**(gamma-1)*source/J**(n+1))
        primitive,details=integrate_ring(integrand)
        if primitive is None:
            rows.append(dict(n=n,closedInRing=False,**details))
            break
        base=sp.cancel(primitive.subs({r:-d,L:0}))
        if base.has(sp.zoo,sp.nan):
            rows.append(dict(n=n,closedInRing=False,reason='Primitive has an apparent stationary singularity; needs a separate limit calculation.'))
            break
        coeff=sp.cancel(J**(n+1)/shifted**gamma*(primitive-base))
        defect=sp.cancel(shifted*deriv(coeff)+(gamma-(n+1)*shifted*2*r/J)*coeff-source)
        assert defect==0
        coefficients.append(coeff)
        row=dict(n=n,closedInRing=True,angleDegree=int(sp.Poly(coeff,L).degree()),
                 coefficient=str(coeff),primitive=str(primitive),
                 odeResidual=0,seconds=time.monotonic()-start,**details)
        rows.append(row)
        print(f'n={n}: exact ODE verified, degree in angle={row["angleDegree"]}, seconds={row["seconds"]:.3f}',flush=True)
        # Save completed orders immediately; the experiment is safely interruptible.
        OUT.mkdir(parents=True,exist_ok=True)
        filename='elementary_lift_symbolic_d.json' if symbolic_d else 'elementary_lift.json'
        (OUT/filename).write_text(json.dumps(dict(
            benchmark=f'lambda=3, mu=1, v={d}; x(r,0)=r+r^3/3',
            angle=f'L=atan(r)+atan({d}); L(-{d})=0',rows=rows,
            scope='Exact finite-order coefficient checks only. No all-order closure, convergence, or closed sum is asserted.'),indent=2)+'\n')
    if rows and not rows[-1]['closedInRing']:
        filename='elementary_lift_symbolic_d_failure.json' if symbolic_d else 'elementary_lift_failure.json'
        (OUT/filename).write_text(json.dumps(rows[-1],indent=2)+'\n')
    if coefficients and len(coefficients)==order:
        # An independent substitution into the complete curl polynomial, after
        # all newly generated coefficients have been included.
        j=[J]+[deriv(x) for x in coefficients]
        kk=[(i+1)*x for i,x in enumerate(coefficients)]
        full=curl_coefficients(j,[deriv(x) for x in j],
            [deriv(deriv(x)) for x in j],kk,[deriv(x) for x in kk],
            [(i+1)*kk[i+1] for i in range(len(kk)-1)],shifted,
            sp.Integer(3),sp.Integer(1),order-1,sp.cancel)
        assert all(x==0 for x in full)
        filename='full_substitution_symbolic_d.json' if symbolic_d else 'full_substitution.json'
        (OUT/filename).write_text(json.dumps(dict(order=order,
            fullCurlPolynomialCoefficients=list(map(str,full)),
            scope='Coefficients of y^0 through y^(order-1) vanish exactly; no claim about omitted powers.'),indent=2)+'\n')
        print(f'complete curl substitution through y^{order-1}: all zero',flush=True)
    return rows


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--order',type=int,default=3)
    parser.add_argument('--symbolic-d',action='store_true')
    args=parser.parse_args()
    run(args.order,args.symbolic_d)
