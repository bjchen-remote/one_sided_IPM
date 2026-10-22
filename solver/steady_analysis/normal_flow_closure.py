"""Symbolic screening of a curved, moving-density geometry.

Read-only wrt existing computations.  No sampled PDE residual minimization.
"""
from pathlib import Path
import json
import sympy as sp


def screen():
    t, k, c, lam, C = sp.symbols('t k c lambda C', real=True)
    results = []
    for alpha in [sp.Rational(1,3), sp.Rational(2,5), sp.Rational(3,8),
                  sp.Rational(1,4), sp.Rational(1,5)]:
        r = (1-alpha)/alpha
        K = k+(1-alpha)*c*t
        L = k+c*t
        M = K+L
        q = 2*lam*(1-alpha)*k + lam*(1-alpha)**2*(2-alpha)*c*t/(1-2*alpha)
        if r.is_integer:
            q += C*t**r
        defect = sp.expand((1+K*K)*(M*sp.diff(q,t)-sp.diff(M,t)*q)
                           +(2-alpha)*c*q*M*M-lam*(1-alpha)**2*c*M**3)
        coeffs = [sp.factor(x) for x in sp.Poly(defect,t).all_coeffs()]
        results.append(dict(alpha=str(alpha),regular_exponent=str(r),
                            q_wall=str(q),second_gate_coefficients=list(map(str,coeffs))))
    return results


def variable_slope_gate():
    alpha,k,k1,lam=sp.symbols('alpha k k1 lambda',real=True)
    q0=2*lam*(1-alpha)*k
    q1=lam*(1-alpha)*(2-alpha)*k1/(1-2*alpha)
    M0=2*k
    M1=(2-alpha)*k1/(1-alpha)
    gate=sp.factor((1+k*k)*(M0*q1-M1*q0)+M1*q0*M0*M0
                   -lam*(1-alpha)*k1*M0**3)
    expected=2*lam*k*k1*((alpha*(2-alpha)/(1-2*alpha))*(1+k*k)+4*k*k)
    assert sp.factor(gate-expected)==0
    return {'stationary_second_gate':str(gate),
            'positive_factor_form':str(expected),
            'consequence':'For 0<alpha<1/2, lambda>0, k>0, smoothness forces k1=0 and q1=0.',
            'scope':'K=Kappa(s*q**(alpha/(1-alpha))), without a shifted invariant.'}


def required_stationary_curvature():
    lam,mu,q,qs=sp.symbols('lambda mu q q_s',positive=True)
    J=1/q;Jp=-qs/q**2
    K=q/(2*(lam-mu))
    Kp=-2*(lam-2*mu)/(2*lam-mu)*K*Jp/J
    S2=-mu*(1+K*K)*Jp/(2*J*J)-(lam-2*mu)*K*Kp/J+Kp/J**2
    curvature=sp.factor(2*S2/(4*lam-3*mu))
    compact=qs*(mu*(2*lam-mu)+(4*lam*lam-6*lam*mu-mu*mu)*K*K)/((2*lam-mu)*(4*lam-3*mu))
    assert sp.factor(curvature-compact)==0
    kval=sp.symbols('k')
    third=sp.factor(curvature.subs(lam,3*mu).subs(q,4*mu*kval))
    assert sp.factor(third-qs*(5+17*kval*kval)/45)==0
    # Independent coefficient extraction from the complete hodograph curl identity.
    s,y=sp.symbols('s y')
    j0,j1,k0,k1,k2,a0,a1=sp.symbols('j0 j1 k0 k1 k2 a0 a1')
    X=j0*s+j1*s*s/2+(k0+k1*s+k2*s*s/2)*y+(a0+a1*s)*y*y
    JJ=sp.diff(X,s);KK=sp.diff(X,y)
    area=sp.integrate(JJ,(y,0,y))
    FF=(2*lam-mu)*area-mu*s*sp.diff(area,s)
    E=JJ*(1+KK*KK)*sp.diff(FF,s)
    E+=FF*(2*JJ*KK*sp.diff(KK,s)-(1+KK*KK)*sp.diff(JJ,s)-JJ*JJ*sp.diff(KK,y))
    E-=2*(lam-mu)*JJ**3*KK
    E+=2*mu*s*JJ*JJ*KK*sp.diff(JJ,s)-mu*s*JJ**3*sp.diff(KK,s)+JJ*JJ
    first=sp.diff(E,y).subs({s:0,y:0})
    solved=sp.solve(first,a0)[0]
    reference=(-mu*(1+k0*k0)*j1/(2*j0*j0)-(lam-2*mu)*k0*k1/j0+k1/j0**2)/(4*lam-3*mu)
    assert sp.factor(solved-reference)==0
    return {'X_yy_at_s0_wall':str(curvature),
            'alpha_one_third':str(third),
            'scope':'Exact necessary second normal jet, regular steady solution; q_s is the wall derivative with respect to shifted density s.'}


if __name__ == '__main__':
    result = {'constant_C_trials':screen(),'variable_slope_gate':variable_slope_gate(),
              'required_stationary_curvature':required_stationary_curvature()}
    out = Path(__file__).resolve().parents[2]/'runs/normal_flow_closure_20260928'
    out.mkdir(exist_ok=True)
    (out/'symbolic_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
