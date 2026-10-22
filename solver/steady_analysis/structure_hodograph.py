"""Exact hodograph structure probes, deliberately independent of Newton solves.

Requires sympy. All checks are symbolic; no claim of full half-plane closure is made.
"""
from pathlib import Path
import json
import sympy as sp

z, y = sp.symbols('z y', real=True)
mu, k, L, A = sp.symbols('mu k L A', positive=True)
eps, C = sp.symbols('eps C', real=True)


def wall_operator(J, K, lam=3*mu):
    return sp.factor(mu*z*sp.diff(K,z)*J +
        2*(lam-mu)*K*J - 2*mu*z*K*sp.diff(J,z) - 1)


def exact_wall_checks():
    R=sp.sqrt(L**2+z**2)
    H=1+eps*z/R
    K=A*H**sp.Rational(-2,3)
    J=(1+2*eps*z*R/L**2+C*z**2/L**2)/(4*mu*A*H**sp.Rational(1,3))
    # Avoid fractional-power branch rewriting: eliminate J'=J(P'/P-H'/3H)
    # and K'=-(2/3)K H'/H. Cross terms cancel identically.
    P=1+2*eps*z*R/L**2+C*z**2/L**2
    reduced=sp.factor((4*P-2*z*sp.diff(P,z))/(4*H)-1)
    assert sp.simplify(reduced)==0
    assert sp.simplify(sp.diff(H,z,2).subs(z,0))==0
    # Hyperbolic positivity of the C=2 numerator.
    t=sp.symbols('t',real=True)
    positivity_identity=sp.trigsimp(P.subs(C,2).subs({z:L*sp.sinh(t),R:L*sp.cosh(t)},simultaneous=True).rewrite(sp.exp)
        -(sp.cosh(2*t)+eps*sp.sinh(2*t)).rewrite(sp.exp))
    assert sp.simplify(positivity_identity)==0
    return {'family': 'K=A*H^(-2/3), J=(1+2*eps*z*sqrt(L^2+z^2)/L^2+C*z^2/L^2)/(4*mu*A*H^(1/3))',
            'H': '1+eps*z/sqrt(L^2+z^2)',
            'wallCompatibility': 'exact symbolic zero',
            'resonance': "H''(0)=0 exactly",
            'positiveBranch': 'A,L,mu>0, |eps|<1, C=2 imply J>0 and K>0 for every real z'}


def polynomial_and_jet_checks():
    r,d=sp.symbols('r d',real=True)
    b=sp.symbols('b0:4')
    q=sp.symbols('q',real=True)
    B=sum(b[i]*q**i for i in range(4))
    polynomial=sp.Poly(sp.expand(mu*(q+d)*(1+q*q)*sp.diff(B,q)
        +(2*mu*r*(1+q*q)-4*mu*q*(q+d))*B-1),q)
    # r=lambda/mu-1. These are successive leading coefficients once higher b vanish.
    leading=[]
    for n in range(3,-1,-1):
        coeff=sp.factor(polynomial.nth(n+2).subs({b[j]:0 for j in range(n+1,4)}))
        assert sp.simplify(coeff-mu*(n+2*r-4)*b[n])==0
        leading.append(str(coeff))
    assert sp.factor(polynomial.as_expr().subs({r:2,b[3]:0,b[2]:0,b[1]:0,b[0]:1/(4*mu)}))==-d*q

    J0=1+z*z
    a2=-(1+k*k)*z/(5*J0**2)
    a3=2*k*(1+k*k)*(17*z*z-4)/(105*J0**4)
    a4=-z*(1+k*k)*((3767*k*k-910)*z*z-1903*k*k+770)/(5775*J0**6)
    S2=-(1+k*k)*z/(2*k*J0**3)
    S3=4*(1+k*k)*(5*z*z-1)/(15*J0**5)
    S4=-z*(1+k*k)*((394*k*k-98)*z*z-173*k*k+70)/(105*k*J0**7)
    for m,am,Sm in [(2,a2,S2),(3,a3,S3),(4,a4,S4)]:
        op=(z*sp.diff(am,z)+(5*m-1-(m+1)*z*sp.diff(J0,z)/J0)*am)/(4*k)
        assert sp.factor(op-Sm)==0
    return {'cubicFirstNormalLeadingCoefficients':leading,
            'a2':str(a2),'a3':str(a3),'a4':str(a4),
            'recurrenceOrdersVerified':[2,3,4],
            'note':'r here denotes lambda/mu-1; source S2-S4 derived from the cleared Poisson numerator.'}


def affine_tilt_checks():
    lam,v=sp.symbols('lam v',positive=True)
    kval=1/(2*(lam-mu)+4*v*v/(2*lam-mu))
    cval=4*v*kval/(2*lam-mu)
    d=v/mu
    assert sp.factor(mu*d*cval+2*(lam-mu)*kval-1)==0
    assert sp.factor(4*mu*d*kval-(2*lam-mu)*cval)==0
    numeric={lam:0.54354325,mu:0.18044829,v:0.16458512}
    return {'normalizedPeakTilt':str(kval),'normalizedPeakTiltDerivative':str(cval),
            'L4FineIllustration':{'K_at_peak':float(kval.subs(numeric)),
                                 'K_derivative_at_peak':float(cval.subs(numeric))},
            'hypergeometricAnalyticBranch':
            'J=(2*mu*r*a)^(-1)*2F1(1,-r-1/2;1-r;-c*z/a), r=lambda/mu-1 not a positive integer'}


def main():
    out={'wallFamily':exact_wall_checks(),
         'polynomialAndJets':polynomial_and_jet_checks(),
         'affineTilt':affine_tilt_checks(),
         'scope':'Exact wall/normal compatibility only, not a full 2D or standard Biot-Savart solution.'}
    dest=Path(__file__).resolve().parents[2]/'runs/structure_20260928/hodograph_checks.json'
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(out,indent=2,ensure_ascii=False))

if __name__=='__main__':main()
