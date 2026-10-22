"""Exact algebra and numerical audits for the nonzero-amplitude hodograph.

These checks do not establish existence, series convergence or a global root.
Run with numpy, scipy and sympy installed; output is tied to the saved seeds.
"""
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from fixed_point import fd_matrix


def symbolic_checks():
    z,y,e=sp.symbols('z y e',real=True)
    lam,mu,v=sp.symbols('lambda mu v',nonzero=True)
    F=sp.Function('F')(z);h=sp.Function('h')(z)

    def laplacian(X,phi):
        J=sp.diff(X,z);K=sp.diff(X,y)
        dx=lambda expr:sp.diff(expr,z)/J
        dy=lambda expr:sp.diff(expr,y)-K*sp.diff(expr,z)/J
        return dx(dx(phi))+dy(dy(phi))

    checked=[]
    for m in range(1,5):
        X=F+e*h*y**m
        phi=(lam*F+v-mu*z*sp.diff(F,z))*y
        phi+=e*(lam*(1-m)*h-mu*z*sp.diff(h,z))*y**(m+1)/(m+1)
        J=sp.diff(X,z);K=sp.diff(X,y)
        transport=-sp.diff(phi,y)+lam*(X-y*K)+v-mu*z*J
        assert sp.simplify(transport)==0
        residual=-laplacian(X,phi)-1/J
        linearized=sp.diff(residual,e).subs(e,0).expand().coeff(y,m-1)
        expected=m*(mu*z*sp.diff(h,z)+(2*m*lam-(m+1)*mu
                      -(m+1)*mu*z*sp.diff(F,z,2)/sp.diff(F,z))*h)
        assert sp.simplify(linearized-expected)==0
        checked.append(m)

    # A cubic inverse wall profile gives an elementary wall trace, but its
    # simplest tilted extension has a nonzero Poisson residual.
    k=sp.symbols('k',positive=True)
    F3=z+z**3/3;J3=sp.diff(F3,z)
    X1=F3+k*y;phi1=2*mu*z*y
    residual1=sp.factor(-laplacian(X1,phi1)-1/J3)
    expected1=(4*mu*k-1)/J3+4*mu*(1+k*k)*z*y/J3**3
    assert sp.simplify(residual1-expected1)==0
    # The bounded second normal coefficient for v=0 and lambda=3*mu.
    a2=-(1+k*k)*z/(5*J3**2)
    X2=X1+a2*y*y
    phi2=z*y/(2*k)+(-3*a2-z*sp.diff(a2,z))*y**3/(12*k)
    residual2=-laplacian(X2,phi2)-1/sp.diff(X2,z)
    assert sp.simplify(residual2.subs(y,0))==0
    assert sp.simplify(sp.diff(residual2,y).subs(y,0))==0
    return dict(transportIdentities='exact symbolic zero',
        highestNewCoefficientVerified=checked,
        cubicTiltPoissonResidual=str(sp.factor(expected1)),
        cubicTiltResidualAfterKCondition=str(sp.factor(expected1.subs(k,1/(4*mu)))),
        cubicNormalCoefficient2=str(a2),
        secondOrderJetPoissonOrders0and1='exact symbolic zero',
        claimsExcluded=['global existence','convergence of the series','identification with the time-evolution limit'])


def numerical_checks(folder):
    results=[]
    s=np.linspace(0,2,1001)
    cubic_f=2*np.sinh(np.arcsinh(1.5*s)/3)
    cubic_q=1/(1+cubic_f**2)
    for name in ('L4_129x65','L4_257x129','L6_193x97','L8_257x129_continued'):
        d=np.load(folder/name/'state.npz')
        x,y,Q=d['x'],d['y'],d['Q'];lam,mu,v,b=d['rates']
        f=CubicSpline(x,Q[0]);q=CubicSpline(x,d['qx'][0])
        g=CubicSpline(x,(fd_matrix(y,1,7)@Q)[0])
        xs=brentq(lambda xx:f(xx)+v/mu,x[0],0)
        inverse=lambda z:brentq(lambda xx:f(xx)+v/mu-z,x[0],x[-1])
        gamma=2*(lam-mu)/mu
        predicted_k=quad(lambda t:t**(gamma-1)*q(inverse(t*v/mu))**3,
                         0,1,epsabs=1e-10)[0]/mu
        gs=-q(xs)**2/(2*(lam-mu))
        even=(q(s)+q(-s))/2;odd=(q(s)-q(-s))/2
        results.append(dict(name=name,stagnationX=xs,
            QyAtStagnation=float(g(xs)),QyAtStagnationFromIdentity=float(gs),
            QyAtPeak=float(g(0)),QyAtPeakFromHodographIntegral=-predicted_k,
            wallEvenCubicTrialLinf=float(np.max(abs(even-cubic_q))),
            wallOddLinf=float(np.max(abs(odd))),
            alpha=float(mu/lam),v=float(v)))
    return results


if __name__=='__main__':
    folder=Path(sys.argv[1])
    result=dict(symbolic=symbolic_checks(),numerical=numerical_checks(folder))
    (folder/'analytical_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
