"""Positive one-parameter asymmetric wall and exact normal-lift ring test.

Data fit does not select lambda/mu=3. The separate symbolic experiment uses
that resonance explicitly. The wall primitive is generally hypergeometric;
the normal corrections use a fixed elementary alphabet.
"""
import argparse
import json
from pathlib import Path
import signal
import time
import numpy as np
import sympy as sp
from sympy.polys.fields import field
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import RectBivariateSpline, PchipInterpolator
from scipy.optimize import minimize_scalar
from scipy.integrate import quad
from ridge_trace_equation import load
from normal_lift_algebra import curl_coefficients

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/ridge_density_closure_20260928'
GN,GW=leggauss(64);GN=(GN+1)/2;GW=GW/2


def jacobian(r,b):
    u0=-b/2;h=np.sqrt(1+u0*u0);u=u0+h*np.asarray(r)
    return (1+u*u)/(1+u0*u0)*np.exp(b*(np.arctan(u)-np.arctan(u0)))


def wall_certificate():
    r,b,d=sp.symbols('r b d',real=True)
    u0=-b/2;h=sp.sqrt(1+u0*u0);u=u0+h*r
    J=(1+u*u)/(1+u0*u0)*sp.exp(b*(sp.atan(u)-sp.atan(u0)))
    normalized=[sp.simplify(sp.diff(J,r,n).subs(r,0)) for n in range(3)]
    assert normalized==[1,0,2]
    coeffs=sp.symbols('p0:5');P=sum(c*r**i for i,c in enumerate(coeffs))
    residual=(1+r*r)*sp.diff(P,r)-(4*r+3*b)*P-(r+d)**3
    matrix,_=sp.linear_eq_to_matrix(sp.Poly(sp.expand(residual),r).all_coeffs(),coeffs)
    determinant=sp.factor(matrix.det())
    assert sp.simplify(determinant+3*b*(9*b*b+4)*(9*b*b+16))==0
    return dict(normalizedJacobianDerivatives=list(map(str,normalized)),
                firstLiftPolynomialDeterminant=str(determinant),
                positivity='J>0 on the real wall. This does not assert positivity of a summed two-dimensional solution.')


def data_audit(name):
    p,_=load(name)
    inv=PchipInterpolator(p['Q'][0],p['x'],extrapolate=False)
    q=RectBivariateSpline(p['y'],p['x'],p['q'])
    r=np.linspace(-.95,.7,171);x=inv(r)
    assert np.all(np.isfinite(x))
    measuredJ=1/q.ev(0*r,x)
    train=(r>=-.5)&(r<=.35);hold=~train
    fit=minimize_scalar(lambda b:np.linalg.norm((np.log(jacobian(r,b))-np.log(measuredJ))[train])**2,
        bounds=(-5,5),method='bounded',options={'xatol':1e-12})
    def metrics(b):
        j=jacobian(r,b);X=r*np.sum(GW*jacobian(r[:,None]*GN,b),axis=1)
        return dict(b=float(b),jacobianTrainRelativeL2=float(np.linalg.norm((j-measuredJ)[train])/np.linalg.norm(measuredJ[train])),
            jacobianWingRelativeL2=float(np.linalg.norm((j-measuredJ)[hold])/np.linalg.norm(measuredJ[hold])),
            inversePositionWingRelativeL2=float(np.linalg.norm((X-x)[hold])/np.linalg.norm(x[hold])))
    lam,mu,drift,off=p['rates'];d=off/mu;gamma=2*lam/mu-2
    K=quad(lambda t:t**(gamma-1)/jacobian(d*(t-1),fit.x)**3,0,1,epsabs=1e-12)[0]/mu
    kp=(1-2*(lam-mu)*K)/off;kpp=4*K-(2*lam-mu)*kp/off
    return dict(name=name,rates=p['rates'].tolist(),fit=metrics(fit.x),symmetricBaseline=metrics(0),
                wallRidgePrediction=dict(A_y=float(-kp),m_y=float(K-kpp/2),C_y=float(-kpp/2)),
                protocol='Fit only b to log(J) on Q in [-.5,.35]; wings to [-.95,.7] held out. Predict wall ridge slopes using actual saved rates, not by imposing lambda=3*mu.',
                sourceKind='finite-time trajectory' if name=='tau20' else 'conditional finite-box root')


def symbolic(order=2):
    start=time.monotonic()
    r,e,L=sp.symbols('r e L')
    F,rr,ee,ll=field((r,e,L),sp.QQ)
    b=sp.Rational(-4,5);jbase=1+rr*rr;J=jbase*ee;s=rr+1
    def D(f):
        f=F(f)
        return f.diff(rr)+b*ee/jbase*f.diff(ee)+f.diff(ll)/jbase
    def integrate_coefficient(expr):
        """Integrate rational(r)*e^m by a finite real rational linear system."""
        num,den=sp.fraction(sp.cancel(expr))
        denparts=sp.Poly(den,e).terms()
        assert len(denparts)==1
        (denpower,),denpoly=denparts[0]
        rational=sp.Integer(0);angle=sp.Integer(0)
        for (num_power,),coef in sp.Poly(num,e).terms():
            m=num_power-denpower
            rn,rd=sp.fraction(sp.cancel(coef/denpoly))
            dp=sp.Poly(rd,r);N=0
            while dp.degree()>0:
                qp,rem=sp.div(dp,sp.Poly(1+r*r,r))
                assert rem.is_zero
                dp=qp;N+=1
            rhs=sp.expand(rn/dp.as_expr());k=N-1
            assert k>=0 and sp.degree(rhs,r)<=2*k
            coeffs=sp.symbols(f'c0:{2*k+1}')
            P=sum(coeffs[i]*r**i for i in range(len(coeffs)))
            angular=sp.Symbol('angular') if m==0 else sp.Integer(0)
            identity=(1+r*r)*sp.diff(P,r)+(m*b-2*k*r)*P+angular*(1+r*r)**k-rhs
            equations=sp.Poly(sp.expand(identity),r).all_coeffs()
            unknowns=list(coeffs)
            if m==0:
                equations.append(P.subs(r,0));unknowns.append(angular)
            matrix,vector=sp.linear_eq_to_matrix(equations,unknowns)
            solutions=list(sp.linsolve((matrix,vector),unknowns))
            assert len(solutions)==1
            answer=dict(zip(unknowns,solutions[0]))
            rational+=e**m*P.subs(answer)/(1+r*r)**k
            if m==0:angle+=answer[angular]
        return sp.cancel(rational),sp.cancel(angle)
    def integrate(f):
        pending=sp.Poly(f.as_expr(),L);primitive=F(0)
        for degree in range(pending.degree(),-1,-1):
            rat,angle=integrate_coefficient(pending.nth(degree))
            term=F.from_expr(rat*L**degree+angle*L**(degree+1)/sp.Integer(degree+1))
            primitive+=term
            pending=sp.Poly((F.from_expr(pending.as_expr())-D(term)).as_expr(),L)
            assert pending.nth(degree)==0
        assert not D(primitive)-f
        return primitive
    aa=[];rows=[]
    for n in range(1,order+1):
        j=[J]+[D(a) for a in aa];k=[(i+1)*a for i,a in enumerate(aa)]
        E=curl_coefficients(j,[D(a) for a in j],[D(D(a)) for a in j],
            k,[D(a) for a in k],[(i+1)*k[i+1] for i in range(len(k)-1)],s,F(3),F(1),n-1,F)[-1]
        source=E/(n*J**3);gamma=5*n-1
        integral=integrate(s**(gamma-1)*source/J**(n+1))
        base=integral.as_expr().subs({r:-1,e:1,L:0})
        a=J**(n+1)/s**gamma*(integral-F.from_expr(base))
        assert not s*D(a)+(gamma-(n+1)*s*D(J)/J)*a-source
        aa.append(a)
        rows.append(dict(n=n,coefficient=str(a.as_expr()),odeResidual=0,
                         containsLogGenerator=a.as_expr().has(L)))
        print('asymmetric exact order',n,'seconds',time.monotonic()-start,flush=True)
    j=[J]+[D(a) for a in aa];k=[(i+1)*a for i,a in enumerate(aa)]
    full=curl_coefficients(j,[D(a) for a in j],[D(D(a)) for a in j],k,[D(a) for a in k],
        [(i+1)*k[i+1] for i in range(len(k)-1)],s,F(3),F(1),order-1,F)
    assert all(not v for v in full)
    return dict(b='-4/5',rates=dict(lam=3,mu=1,d=1),
                coordinate='L=atan(r)+pi/4, e=exp(b*L)',
                rows=rows,fullCurlCoefficients=list(map(str,full)),seconds=time.monotonic()-start,
                scope='A rational b benchmark, not identification of the fitted b as exactly -4/5.')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--order',type=int,default=2)
    parser.add_argument('--skip-data',action='store_true')
    args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    if not args.skip_data:
        data=[data_audit(n) for n in ['L4_129x65','L4_257x129','L8_257x129_continued','tau20']]
        (OUT/'asymmetric_wall_data.json').write_text(json.dumps(data,indent=2)+'\n')
        print(json.dumps(data,indent=2),flush=True)
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('45 second symbolic bound')))
    signal.alarm(45)
    try:
        exact=symbolic(args.order)
        exact['wallCertificate']=wall_certificate()
        (OUT/'asymmetric_ring_exact.json').write_text(json.dumps(exact,indent=2)+'\n')
    finally:
        signal.alarm(0)
