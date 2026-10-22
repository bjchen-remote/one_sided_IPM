"""Exact certificates supporting the formal all-order ring argument.

The induction proof is in REPORT_ZH.md. These checks verify its two scaling
identities, the pole/growth hypotheses through order five, full substitution,
and strict exclusions of two simplest finite resummation classes.
"""
import json
from pathlib import Path
import sympy as sp
from sympy.polys.fields import field
from normal_lift_algebra import curl_coefficients

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/ridge_density_closure_20260928'


def main():
    s,y,mu=sp.symbols('s y mu',nonzero=True)
    cs=sp.symbols('c0:4')
    homogeneous=[]
    for p,b in [(1,5),(3,3)]:
        X=sum(cs[i]*s**(p-b*i)*y**i for i in range(4))
        J=sp.diff(X,s);K=sp.diff(X,y)
        H=sp.integrate(J,(y,0,y))
        flux=5*mu*H-mu*s*sp.diff(H,s)
        U=sp.cancel(flux/J)
        V1=sp.cancel(mu*s*J+K*U)
        assert sp.cancel(U-b*mu*y)==0
        assert sp.cancel(V1-p*mu*X)==0
        curl=sp.cancel(sp.diff(U,s)/J-sp.diff(V1,y)+K*sp.diff(V1,s)/J)
        assert curl==0
        homogeneous.append(dict(horizontalWeight=p,verticalWeight=b,
                                 V1=f'{p}*mu*X',V2=f'{b}*mu*y',curlResidual=0))

    r,L=sp.symbols('r L',real=True);J=1+r*r
    old=json.loads((ROOT/'runs/ridge_symbolic_20260928/elementary_lift.json').read_text())
    a=[];sources=[]
    for row in old['rows']:
        a.append(sp.sympify(row['coefficient'],locals={'r':r,'L':L}))
        prim=sp.sympify(row['primitive'],locals={'r':r,'L':L})
        sources.append(sp.cancel(sp.diff(prim,r)+sp.diff(prim,L)/J))
    for name in ['fourth_elementary_gate.json','elementary_gate_5.json']:
        row=json.loads((OUT/name).read_text())
        assert row['status']=='closed_in_same_ring'
        a.append(sp.sympify(row['coefficient'],locals={'r':r,'L':L}))
        sources.append(sp.sympify(row['integrand'],locals={'r':r,'L':L}))
    poles=[]
    for n,source in enumerate(sources,1):
        pp=sp.Poly(source,L)
        rows=[]
        for k in range(pp.degree()+1):
            f=sp.cancel(pp.nth(k))
            if f==0:
                continue
            num,den=sp.fraction(f)
            N=sp.Poly(num,r);D=sp.Poly(den,r);jpower=0
            while D.degree()>0:
                Q,rem=sp.div(D,sp.Poly(J,r))
                assert rem.is_zero
                D=Q;jpower+=1
            decay=2*jpower-N.degree()
            assert decay>=2
            rows.append(dict(anglePower=k,denominatorJpower=jpower,decayPower=int(decay)))
        poles.append(dict(n=n,coefficientChecks=rows))

    F,rr,ll=field((r,L),sp.QQ);jj=1+rr*rr
    D=lambda f:F(f).diff(rr)+F(f).diff(ll)/jj
    aa=[F.from_expr(f) for f in a]
    j=[jj]+[D(f) for f in aa];k=[(i+1)*f for i,f in enumerate(aa)]
    E=curl_coefficients(j,[D(f) for f in j],[D(D(f)) for f in j],
        k,[D(f) for f in k],[(i+1)*k[i+1] for i in range(len(k)-1)],
        rr+1,F(3),F(1),4,F)
    assert all(f==0 for f in E)

    # L is transcendental over C(r). A nonzero determinant at this algebraic
    # specialization proves the corresponding rational determinant is not zero.
    # This is NOT a physical evaluation of atan(0)+pi/4.
    vals=[sp.factor(f.subs({r:0,L:0})) for f in a]
    w=sum(f*y**(i+1) for i,f in enumerate(vals))
    ranks=[]
    for name,basis in [('total_degree_two_implicit',[1,y,w,y*y,y*w,w*w]),
                       ('rational_pade_2_2',[1,y,y*y,w,y*w,y*y*w])]:
        matrix=sp.Matrix([[sp.expand(f).coeff(y,n) for f in basis] for n in range(6)])
        det=sp.factor(matrix.det())
        assert det!=0
        ranks.append(dict(family=name,
                          monomials=(['1','y','w','y^2','y*w','w^2'] if name.startswith('total') else ['1','y','y^2','w','y*w','y^2*w']),
                          specialization=dict(r=0,L=0),determinant=str(det),rank=6,
                          conclusion='No nonzero relation in this monomial space over C(r,L).'))
    result=dict(homogeneousCancellations=homogeneous,polesThroughOrderFive=poles,
                fullCurlCoefficients=list(map(str,E)),finiteResummationExclusions=ranks,
                scope='Exact algebraic checks. The general proof concerns formal coefficients, not convergence or a finite closed sum.')
    (OUT/'ring_certificates.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(homogeneousCancellations=homogeneous,
                         verifiedOrders=len(poles),fullCurlCoefficients=result['fullCurlCoefficients'],
                         exclusions=ranks),indent=2))


if __name__=='__main__':
    main()
