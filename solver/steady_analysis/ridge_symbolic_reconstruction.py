"""Local ridge Cauchy reconstruction and a full-profile symbolic closure gate.

The analytic Cauchy statement is local and is not a Sobolev well-posedness
claim. The elementary first normal lift is not a full two-dimensional solution.
Existing fields and earlier agents' outputs are only read.
"""
from pathlib import Path
import json

import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.interpolate import RectBivariateSpline, make_interp_spline
from scipy.optimize import brentq

from ridge_trace_equation import load, peak

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'runs/ridge_symbolic_20260928'


def cauchy_checks():
    A, g, b, bp, m1, m2, U, Up, S, lam, mu = sp.symbols(
        'A g beta beta_prime m1 m2 U Up S lambda mu')
    k = m1-b
    kp = m2-bp
    D = 1+m1*m1
    T = S-b*U
    Tp = mu*b-S*g-bp*U-b*Up
    r = (Tp+2*m1*Up+m2*U-2*lam*m1-A)/D
    W1 = 2*lam-Up
    R2 = ((mu-W1)*A-r*A*b-U*A*g)/T
    r1 = ((m1+k)*g+kp)/(1+k*k)
    r0 = (2*(lam-mu)*k-A-S*g)/(1+k*k)
    a = g+b*r1
    f = 2*lam-mu+b*r0
    factor = A*(1+k*k)/(T*D)
    assert sp.factor(R2-factor*(Up-a*U-f)) == 0

    pp, rp, Aprime = sp.symbols('p_prime r_prime A_prime')
    U2 = (pp+m1*rp)/D
    h_from_recurrence = (A*rp-2*r*Aprime-U2*A*b)/T
    h_from_old_gate = -(2*r*Aprime+A*(b*pp-(1+k*m1)*rp)/D)/T
    assert sp.factor(h_from_recurrence-h_from_old_gate) == 0

    # Directly verify the change of variables in divergence and curl.
    Wn, Wy, Un, Uy, Rn, W, Wyy = sp.symbols('Wn Wy Un Uy Rn W Wyy')
    V1x = Wn+m1*Un
    V2y = Uy-m1*Un
    V1y = Wy+m2*U+m1*Uy-m1*Wn-m1*m1*Un
    assert sp.expand(V1x+V2y-(Wn+Uy)) == 0
    curl = sp.expand(Un-V1y)
    expected = D*Un-Wy+m1*Wn-m1*Uy-m2*U
    assert sp.expand(curl-expected) == 0
    solved = curl.subs({Wn:2*lam-Uy,
                        Un:(Wy+2*m1*Uy+m2*U-2*lam*m1-Rn)/D})
    assert sp.factor(solved+Rn) == 0
    return dict(
        ridgeConditionFactor=str(factor),
        ridgeConditionIdentity='R2 = factor * (U_prime-a*U-f)',
        residuals=dict(divergence=0, curl=0, ridgeODE=0, curvature=0),
        scope='Analytic local Cauchy reconstruction when T != 0. Neither global boundary conditions nor stability in smooth norms follows.')


def cubic_gate():
    r, d, mu = sp.symbols('r d mu', real=True, nonzero=True)
    J = 1+r*r
    atan_coefficient = 3*d*(1+d*d)/8
    rational_part = (-6*d*d+3*d*r**3*(1+d*d)+d*r*(5*d*d-3)
                     -4*r*r-2)/(8*J*J)
    primitive = rational_part+atan_coefficient*sp.atan(r)
    assert sp.factor(sp.diff(primitive, r)-(r+d)**3/J**3) == 0
    base = sp.simplify(primitive.subs(r, -d))
    K = J*J*(primitive-base)/(mu*(r+d)**4)
    defect = mu*(r+d)*J*sp.diff(K,r)+(4*mu*J-2*mu*(r+d)*sp.diff(J,r))*K-1
    assert sp.factor(defect) == 0
    # The logarithmic residues are invariant under adding the homogeneous solution.
    residues = [sp.simplify(sp.residue((r+d)**3/J**3,r,root))
                for root in [sp.I, -sp.I]]
    assert sp.simplify(residues[0]-atan_coefficient/(2*sp.I)) == 0
    assert sp.simplify(residues[1]+atan_coefficient/(2*sp.I)) == 0

    mp.mp.dps = 65
    Kfun = sp.lambdify((r,d,mu),K,'mpmath')
    rows=[]
    dd=mp.mpf('0.9'); mm=mp.mpf('0.18')
    for rr in [mp.mpf('-2'),mp.mpf('-.90001'),mp.mpf('-.5'),mp.mpf('0'),mp.mpf('1')]:
        direct=(1+rr*rr)**2/mm*mp.quad(
            lambda t:t**3/(1+((rr+dd)*t-dd)**2)**3,[0,1])
        value=Kfun(rr,dd,mm)
        rows.append(dict(r=float(rr),closedK=float(value),
                         differenceFromRegularIntegral=float(abs(value-direct))))
    # At the shifted stationary density r=-d the apparently singular formula is regular.
    stationary = 1/(4*mm*(1+dd*dd))
    direct_stationary=(1+dd*dd)**2/mm*mp.quad(lambda t:t**3/(1+dd*dd)**3,[0,1])
    assert abs(stationary-direct_stationary) < mp.mpf('1e-60')
    slopes = dict(K0=sp.factor(K.subs(r,0)),
                  A_y=-sp.factor(sp.diff(K,r).subs(r,0)),
                  m_y=sp.factor((K-sp.diff(K,r,2)/2).subs(r,0)))
    return dict(
        wall='x=r+r^3/3, r=Q(x,0), s=r+d, d=v/mu, lambda=3*mu',
        primitive=str(primitive),atanCoefficient=str(atan_coefficient),
        complexResidues=list(map(str,residues)),K=str(K),
        stationaryValue='1/(4*mu*(1+d^2))',
        wallSlopes={key:str(value) for key,value in slopes.items()},
        exactResidual=0, checks=rows,
        conclusion='For real d != 0, the unique regular first lift is elementary but not algebraic over C(r). A finite algebraic hodograph family through this exact cubic wall trace cannot solve the full PDE.',
        exclusions='Does not exclude elementary/logarithmic closure, other wall traces, or the empirical cubic relation m(A).')


def wall_seed_check():
    Ay, my, cy, K, lam, mu, v, Ms = sp.symbols('Ay my cy K lambda mu v Mstar')
    nu=2*lam-mu
    relations=[v*(-Ay)+2*(lam-mu)*K-1,
               my+K+nu*Ay/(2*v), my+sp.Rational(5,2)*Ms*Ay]
    sol=sp.solve(relations,[Ay,my,K],dict=True)[0]
    expected=-v/(v*v+(lam-mu)*(nu-5*v*Ms))
    assert sp.factor(sol[Ay]-expected) == 0
    return dict(A_y=str(sp.factor(sol[Ay])),m_y=str(sp.factor(sol[my])),
                C_y=str(sp.factor(sol[my]-sol[K])),
                Y_t_at_1=str(sp.factor(-2/expected)),
                assumptions='Q(0,0)=0, q=1, q_x=0, q_xx=-2; lambda>mu>0; drift and density offset both v; m=Mstar*(7t^3-8t^2+1), A=t^-2.')


def data_audit(name):
    p,_=load(name)
    q=RectBivariateSpline(p['y'],p['x'],p['q'])
    Q=RectBivariateSpline(p['y'],p['x'],p['Q'])
    lam,mu,drift,offset=p['rates']
    A0=peak(q,p['x'],0)[1]
    heights=p['y'][(p['y']>=0)&(p['y']<=min(3,float(p['y'][-1])))]
    pa=np.array([peak(q,p['x'],float(y)) for y in heights])
    crossings=np.flatnonzero((pa[:-1,1]>A0/2)&(pa[1:,1]<=A0/2))
    half=None
    if len(crossings):
        i=crossings[0]
        yh=brentq(lambda y:peak(q,p['x'],y)[1]-A0/2,heights[i],heights[i+1])
        mh=peak(q,p['x'],yh)[0]
        Ms=(mh-pa[0,0])/(14*np.sqrt(2)-15)
        denom=offset**2+(lam-mu)*(2*lam-mu-5*offset*Ms)
        half=dict(yHalf=float(yh),mHalf=float(mh),Mstar=float(Ms),
                  predictedA_y=float(-offset/denom),
                  predictedm_y=float(2.5*Ms*offset/denom),
                  predictedY_t=float(2*denom/offset))
    slopes=[]
    for stride,degree in [(1,3),(1,5),(2,5)]:
        ys=p['y'][(p['y']>=0)&(p['y']<=.75)][::stride]
        if len(ys)<=degree:
            continue
        mm,aa=np.array([peak(q,p['x'],float(y)) for y in ys]).T
        M=make_interp_spline(ys,mm,k=degree)
        A=make_interp_spline(ys,aa,k=degree)
        C=make_interp_spline(ys,Q.ev(ys,mm),k=degree)
        slopes.append(dict(stride=stride,degree=degree,
                           A_y=float(A(0,1)),m_y=float(M(0,1)),C_y=float(C(0,1))))
    # Test the EXACT cubic wall trace at the saved rates, without assuming lambda=3*mu.
    d=offset/mu;gamma=2*lam/mu-2
    K=quad(lambda t:t**(gamma-1)/(1+d*d*(t-1)**2)**3,0,1,
           epsabs=1e-13,epsrel=1e-12)[0]/mu
    K1=(1-2*(lam-mu)*K)/offset
    K2=4*K-(2*lam-mu)*K1/offset
    cubic=dict(A_y=float(-K1),m_y=float(K-K2/2),C_y=float(-K2/2))
    return dict(name=name,rates=p['rates'].tolist(),measuredWallSlopes=slopes,
                cubicWallPrediction=cubic,empiricalRidgeWallSeed=half,
                sourceKind='finite-time trajectory' if name=='tau20' else 'conditional finite-box root',
                caveat='Near-wall derivatives are interpolation-sensitive. The ridge-seed calculation uses ideal normalizations q=1,q_xx=-2; the exact-cubic wall is a trial, not a fit to the full measured wall.')


if __name__=='__main__':
    OUT.mkdir(parents=True,exist_ok=True)
    result=dict(cauchy=cauchy_checks(),cubicGate=cubic_gate(),
                wallSeed=wall_seed_check(),
                data=[data_audit(name) for name in ['L4_257x129','L8_257x129_continued','tau20']])
    (OUT/'checks.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(dict(cauchy=result['cauchy'],wallSeed=result['wallSeed'],
        cubicQuadratureMaxError=max(row['differenceFromRegularIntegral'] for row in result['cubicGate']['checks']),
        data=result['data']),indent=2))
