"""Globally positive analytic wall geometries with free alpha=1/(r+1).

These solve the FIRST wall-normal compatibility equation exactly. They are not
claimed to solve the full two-dimensional IPM or its Biot-Savart closure.
"""
from pathlib import Path
import json
import numpy as np
import sympy as sp


def counterterm_wall(s, r, b=1., c=1., eps=0., L=1., A=1., mu=1., delta=0.):
    """Quadratic-radical family; positivity bound is documented below.

    Uniform real analyticity across s=0 is built into the formula even when r
    is noninteger. For 0<r<4 use counterterm_positivity_bound for b/c.
    """
    x=np.asarray(s)/L
    w=1+x*x
    chi=x/np.sqrt(w)
    p=r/2
    B=b*(1+delta*chi)+c*(w**p-p*w**(p-1)+eps*x*w**((r-1)/2)
           -eps*(r-1)/2*x*w**((r-3)/2))
    H=r*b+b*delta*chi*(r-1/w)-c/2*w**(p-2)*(r*(r-2)+eps*(r-1)*(r-3)*chi)
    K=A*H**(-2/3)
    J=B/(2*mu*A*H**(1/3))
    return {'B':B,'H':H,'K':K,'J':J}


def counterterm_positivity_bound(r,c,eps,delta=0.):
    """Sufficient b>return value; 0<r<4, |eps|,|delta|<1.

    For r>1 the denominator is exactly 2*r*(1-abs(delta)).
    """
    base=r-abs(delta)*max(r,abs(r-1))
    if base<=0:return float('inf')
    return c*(abs(r*(r-2))+abs(eps*(r-1)*(r-3)))/(2*base)


def soft_power_wall(s,r,b=1.,c=1.,eps=0.,L=1.,A=1.,mu=1.,m=2):
    """Parent agent's simpler soft-power family; 2*m>r required for H->r*b.

    Positivity: b,c,A,L,mu>0, |eps|<1, r>|eps|*|r-1|.
    """
    x=np.asarray(s)/L
    R=(1+x**(2*m))**(1/(2*m))
    chi=x/R
    B=b+c*R**r*(1+eps*chi)
    H=r*b+c*R**(r-2*m)*(r+eps*(r-1)*chi)
    K=A*H**(-2/3)
    J=B/(2*mu*A*H**(1/3))
    return {'B':B,'H':H,'K':K,'J':J}


def symbolic_checks():
    x,r,b,c,e=sp.symbols('x r b c e',real=True)
    w=1+x*x
    p=r/2
    B=b+c*(w**p-p*w**(p-1)+e*x*w**((r-1)/2)
            -e*(r-1)/2*x*w**((r-3)/2))
    H=r*b-c*w**(p-2)*(r*(r-2)+e*(r-1)*(r-3)*x/sp.sqrt(w))/2
    # Positive w permits collection of real powers without branch ambiguity.
    residual=sp.powsimp(sp.expand_power_base(r*B-x*sp.diff(B,x)-H,force=True),force=True)
    assert sp.simplify(residual)==0
    delta=sp.symbols('delta',real=True)
    chi=x/sp.sqrt(w)
    Bd=B+b*delta*chi
    Hd=H+b*delta*chi*(r-1/w)
    delta_residual=sp.powsimp(sp.expand_power_base(r*Bd-x*sp.diff(Bd,x)-Hd,force=True),force=True)
    assert sp.simplify(delta_residual)==0
    B0=sp.simplify(B.subs({e:0,x:0}))
    b2=sp.simplify(sp.diff(B,x,2).subs({e:0,x:0})/2)
    H0=sp.simplify(H.subs({e:0,x:0}))
    h2=sp.simplify(sp.diff(H,x,2).subs({e:0,x:0})/2)
    j_rel2=sp.factor(b2/B0-h2/(3*H0))
    expected=c*(4-r)*(r+1)/(6*B0)
    assert sp.factor(j_rel2-expected)==0

    soft_results=[]
    for m in [1,2,3]:
        R=(1+x**(2*m))**sp.Rational(1,2*m)
        chi=x/R
        Bs=b+c*R**r*(1+e*chi)
        Hs=r*b+c*R**(r-2*m)*(r+e*(r-1)*chi)
        soft_residual=sp.powsimp(sp.expand_power_base(r*Bs-x*sp.diff(Bs,x)-Hs,force=True),force=True)
        assert sp.simplify(soft_residual)==0
        soft_results.append(m)
    return {'countertermIdentity':'r*B-s*B_prime-H = exact zero',
            'deltaExtensionIdentity':'exact symbolic zero',
            'softPowerIdentityVerifiedM':soft_results,
            'countertermSymmetricRelativeQuadraticCoefficient':str(j_rel2),
            'positiveRange':'0<r<4, |eps|<1, b>counterterm_positivity_bound(r,c,eps,delta)',
            'deltaExtensionForRGreaterOne':'|delta|<1 and r*b*(1-|delta|)>c*(abs(r*(r-2))+abs(eps*(r-1)*(r-3)))/2',
            'softPowerRange':'r>0, 2*m>r, b,c>0, |eps|<1, r>|eps|*abs(r-1)',
            'warning':'soft-power m>=2 and eps=0 has a degenerate wall-gradient maximum at s=0'}


def numerical_examples():
    # Numerical values supplement the analytic positivity proof; no grid test
    # is used as a substitute for global positivity.
    s=np.r_[0.,np.geomspace(1e-8,1e6,4000),-np.geomspace(1e-8,1e6,4000)]
    result=[]
    for r in [1.8,1.9,2.,2.1,2.2]:
        eps=.6;b=c=1.;A=mu=L=1.
        bound=counterterm_positivity_bound(r,c,eps)
        assert b>bound
        f=counterterm_wall(s,r,b,c,eps,L,A,mu)
        for name in ['B','H','K','J']:assert np.all(f[name]>0)
        fplus=counterterm_wall(1e6,r,b,c,eps)
        fminus=counterterm_wall(-1e6,r,b,c,eps)
        expected_plus=c*(1+eps)/(2*mu*A*(r*b)**(1/3)*L**r)
        expected_minus=c*(1-eps)/(2*mu*A*(r*b)**(1/3)*L**r)
        result.append({'r':r,'alpha':1/(r+1),'b':b,'c':c,'eps':eps,
            'sufficient_b_lower_bound':bound,
            'analytic_H_lower_bound':r*(b-bound),
            'sampled_min_B':float(f['B'].min()),'sampled_min_H':float(f['H'].min()),
            'sampled_min_J':float(f['J'].min()),
            'J_plus_tail_coefficient':expected_plus,
            'J_minus_tail_coefficient':expected_minus,
            'J_plus_tail_relative_check':float(fplus['J']/(1e6**r*expected_plus)-1),
            'J_minus_tail_relative_check':float(fminus['J']/(1e6**r*expected_minus)-1)})
    return result


def main():
    out={'symbolic':symbolic_checks(),'examples':numerical_examples(),
         'scope':'Exact analytic positive wall-geometry families, not yet full 2D steady states.'}
    path=Path(__file__).resolve().parents[2]/'runs/structure_free_alpha_20260928/free_wall_checks.json'
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(out,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
