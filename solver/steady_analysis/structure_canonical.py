"""Exact area-coordinate identities and a read-only probe of saved profiles.

This is structure discovery, not a steady-state solver or error-correction step.
All fields use (y,x) layout. The numerical checks test identities locally and
do not identify a continuum limit or impose a far-field condition.
"""
from pathlib import Path
import json
import numpy as np
import sympy as sp
from scipy.interpolate import RectBivariateSpline
from scipy.integrate import cumulative_simpson
from scipy.optimize import brentq


def exact_checks():
    s, y, lam, mu, v = sp.symbols('s y lambda mu v', real=True)
    I = sp.Function('I')(s, y)
    X = sp.diff(I, y)
    J, K = sp.diff(X, s), sp.diff(X, y)
    T = sp.diff(I, s)
    phi = lam*(2*I-y*X)+v*y-mu*s*sp.diff(I, s)
    Vy = sp.diff(phi, s)/J + lam*y
    Vx = -sp.diff(phi, y)+K*sp.diff(phi, s)/J+lam*X+v
    sdot = sp.simplify((Vx-K*Vy)/J)
    tdot = sp.simplify(sp.diff(T, s)*sdot+sp.diff(T, y)*Vy)
    assert sp.simplify(sdot-mu*s) == 0
    assert sp.simplify(tdot-(2*lam-mu)*T) == 0
    assert sp.simplify(sp.diff(T, y)-J) == 0

    # Pullback of curl(V)=-s_x for an arbitrary area-preserving map Z.
    t, nu = sp.symbols('t nu', real=True)
    A, B = sp.Function('A')(s,t), sp.Function('B')(s,t)
    W1=mu*s*sp.diff(A,s)+nu*t*sp.diff(A,t)
    W2=mu*s*sp.diff(B,s)+nu*t*sp.diff(B,t)
    curl_pullback = (sp.diff(W1*sp.diff(A,t)+(W2+s)*sp.diff(B,t),s)
                     -sp.diff(W1*sp.diff(A,s)+(W2+s)*sp.diff(B,s),t))
    Zs=sp.Matrix([sp.diff(A,s),sp.diff(B,s)])
    Zt=sp.Matrix([sp.diff(A,t),sp.diff(B,t)])
    metric_form=((mu-nu)*Zs.dot(Zt)
        +mu*s*(Zs.diff(s).dot(Zt)-Zt.diff(s).dot(Zs))
        +nu*t*(Zs.diff(t).dot(Zt)-Zt.diff(t).dot(Zs))
        +sp.diff(B,t))
    assert sp.simplify(curl_pullback-metric_form)==0
    return dict(area_jacobian='exact zero defect',
                density_transport='ds/dtau = mu*s',
                area_transport='dt/dtau = (2*lambda-mu)*t',
                curl_pullback='exact zero defect',
                caveat='Identities; no global existence or exact-profile claim.')


def numerical_checks(root):
    results=[]
    for name in ['L4_129x65','L4_257x129','L6_193x97','L8_257x129_continued']:
        d=np.load(root/name/'state.npz')
        x,y,Q=d['x'],d['y'],d['Q']
        lam,mu,v,b=d['rates']; nu=2*lam-mu
        qs=RectBivariateSpline(y,x,Q,kx=3,ky=3)
        ps=RectBivariateSpline(y,x,d['psi'],kx=3,ky=3)
        levels=np.linspace(-.55,.4,41)
        hh=np.linspace(0,1.5,193)
        xx=np.empty((len(hh),len(levels)))
        for j,h in enumerate(hh):
            for i,q in enumerate(levels):
                xx[j,i]=brentq(lambda u:float(qs(h,u)[0,0])-q,x[0],x[-1])
        yy=np.broadcast_to(hh[:,None],xx.shape)
        qx=qs.ev(yy,xx,dx=0,dy=1)
        qxx=qs.ev(yy,xx,dx=0,dy=2)
        J=1/qx
        T=cumulative_simpson(J,x=hh,axis=0,initial=0)
        Ts=cumulative_simpson(-qxx/qx**3,x=hh,axis=0,initial=0)
        Vy=ps.ev(yy,xx,dx=0,dy=1)+lam*yy
        ss=np.broadcast_to(levels[None,:]+v/mu,xx.shape)
        defect=mu*ss*Ts+Vy*J-nu*T
        xs=brentq(lambda u:float(qs(0,u)[0,0])+v/mu,x[0],0)
        q0=float(qs(0,xs,dx=0,dy=1)[0,0])
        g0=float(qs(0,xs,dx=1,dy=0)[0,0])
        results.append(dict(name=name,alpha=float(mu/lam),
            transverse_to_density_rate=float(nu/mu),
            predicted_physical_space_exponent=float(lam/(lam-mu)),
            predicted_physical_density_exponent=float(mu/(lam-mu)),
            stagnation_x=xs,stagnation_q=q0,stagnation_Qy=g0,
            stagnation_Qy_identity=-q0*q0/(nu-mu),
            canonical_transport_abs_linf=float(np.max(abs(defect))),
            canonical_transport_relative_linf=float(np.max(abs(defect))/np.max(abs(nu*T))),
            density_levels=levels[[0,-1]].tolist(),height=float(hh[-1]),
            interpretation='Local spline/quadrature identity check on finite-box data.'))
    return results


if __name__=='__main__':
    root=Path(__file__).resolve().parents[2]
    output=root/'runs/structure_20260928'
    output.mkdir(exist_ok=True)
    result=dict(exact=exact_checks(),data=numerical_checks(root/'runs/steady_20260928_tau20'))
    (output/'canonical_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
