"""A drift/skew resonance: exact wall identities and held-out geometry.

The candidate skew is fixed by v/mu, never fitted to the saved wall.
Only bounded one-dimensional quadratures and saved fields are used.
"""
import json
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.interpolate import RectBivariateSpline, PchipInterpolator
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
from asymmetric_elementary_completion import jacobian
from ridge_trace_equation import load
from ridge_generated_jets import WallLift

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'runs/drift_skew_closure_20260928'


def certificate():
    u,b,k,e=sp.symbols('u b k e', real=True)
    R=1+u*u;s=u-sp.Rational(3,4)*b
    ell=4+sp.Rational(9,4)*b*b;mu=1/(k*ell)
    D=lambda f:sp.factor(sp.diff(f,u)+b*e/R*sp.diff(f,e))
    J=R*e;K=k/e
    first=sp.factor(mu*s*D(K)+(4*mu-2*mu*s*D(J)/J)*K-1/J)
    assert first==0
    square_rate=sp.factor(2*J*D(K))
    assert square_rate==-2*k*b

    d,mm=sp.symbols('d mu', positive=True)
    h=sp.sqrt(1+b*b/4)
    dcanon=-5*b/(4*h)
    kcanon=h*h/(mm*ell)
    assert sp.simplify(kcanon-5/(4*mm*(5+dcanon*dcanon)))==0
    peak_first=sp.factor((-b/h)*kcanon)
    assert sp.simplify(peak_first-dcanon/(mm*(5+dcanon*dcanon)))==0
    # At the canonical peak u=-b/2, J=1 and e=1.
    canonical_K=h*h/(mm*ell*e)
    second_density_derivative=sp.factor(h*h*D(D(canonical_K)).subs({u:-b/2,e:1}))
    assert second_density_derivative==0

    # Constructive derivation: prescribe only the initial square-Jacobian rate.
    ss,lam0,mu0,B,N,Np=sp.symbols('s lambda mu B N Nprime')
    first_in_N=mu0*ss*B+2*(lam0-mu0)*N-2*mu0*ss*(Np-B)-1
    euler=2*(lam0-mu0)*N-2*mu0*ss*Np+3*mu0*ss*B-1
    assert sp.expand(first_in_N-euler)==0
    r,dd=sp.symbols('r d', real=True)
    P=1+sp.Rational(4,5)*dd*r+r*r
    kk=5/(4*mu0*(5+dd*dd));BB=sp.Rational(4,5)*dd*kk
    NN=kk*P
    euler_cubic=sp.factor(4*mu0*NN-2*mu0*(r+dd)*sp.diff(NN,r)
                           +3*mu0*(r+dd)*BB-1)
    assert euler_cubic==0
    assert sp.factor((sp.diff(NN,r)-BB)/NN-2*r/P)==0
    normalized_wall_gate=sp.factor(mu0*(r+dd)*BB
        +(4*mu0-4*mu0*(r+dd)*r/P)*NN-1)
    assert normalized_wall_gate==0
    higher=[]
    for n in [2,4,6]:
        kn=(2*n+1)/(2*mu0*(n*(2*n+1)+2*dd**2))
        Bn=4*dd*kn/(2*n+1)
        Cn=2*kn/(n*(n-1)*dd**(n-2))
        poly=1/(2*n*mu0)-3*Bn*ss/(2*(n-1))+Cn*ss**n
        euler_n=sp.factor(2*n*mu0*poly-2*mu0*ss*sp.diff(poly,ss)+3*mu0*ss*Bn-1)
        peak_n=[sp.factor(poly.subs(ss,dd)-kn),
                sp.factor(sp.diff(poly,ss).subs(ss,dd)-Bn),
                sp.factor(sp.diff(poly,ss,2).subs(ss,dd)/kn-2)]
        assert euler_n==0 and peak_n==[0,0,0]
        higher.append(dict(n=n,rateRatio=n+1,eulerResidual=str(euler_n),
                           normalizationResiduals=list(map(str,peak_n))))

    jp,jpp,kp,kpp=D(J),D(D(J)),D(K),D(D(K))
    E1=(mu*(1+K*K)*(s*jp*jp-J*jp-s*J*jpp)
        -2*mu*J*J*K*kp+2*mu*s*J*K*kp*jp
        +2*mu*s*J*J*K*kpp-3*mu*s*J*J*kp*kp+2*J*kp)
    rhs=sp.factor(E1/(2*mu*J**3))
    c0,c1=sp.symbols('c0 c1')
    gates=[]
    for power in [-1,-3]:
        force=sp.factor(sp.expand(rhs*e**3).coeff(e,power+3))
        if power==-3:force=sp.factor(force/k**2)
        f=(c0+c1*u)/R**2
        residual=sp.factor(s*sp.diff(f,u)+(9+s*(power*b-3*(2*u+b))/R)*f-force)
        equations=sp.Poly(sp.fraction(residual)[0],u).all_coeffs()
        gb=sp.groebner(equations,c0,c1,b)
        eliminated=[sp.factor(p.as_expr()) for p in gb.polys
                    if not p.as_expr().has(c0,c1)]
        gates.append(dict(exponentialPower=power,forcing=str(force),
                          eliminationPolynomials=list(map(str,eliminated))))
    return dict(firstLayerResidual=str(first),squareJacobianWallRate=str(square_rate),
                squareRateEulerIdentity=str(sp.expand(first_in_N-euler)),
                normalizedConstructiveGate=str(normalized_wall_gate),
                eulerCubicResidual=str(euler_cubic),
                higherDiscreteBranches=higher,
                canonicalDensityDrift=str(dcanon),canonicalK=str(kcanon),
                canonicalKSecondDerivative=str(second_density_derivative),
                secondLayerRationalCoefficientGates=gates,
                gateScope='The ansatz a2=e^-1 P1(u)/(1+u^2)^2 + k^2 e^-3 Q1(u)/(1+u^2)^2, where P1,Q1 have degree at most one; not a general second-layer exclusion.')


def audit(name,nodes=64):
    field,_=load(name)
    lam,mu,v,off=map(float,field['rates'])
    d=off/mu
    if abs(d)>=2.5:
        return dict(name=name,d=d,admissible=False)
    b=-4*d/np.sqrt(25-4*d*d)
    r=np.linspace(-.95,.7,171)
    inverse=PchipInterpolator(field['Q'][0],field['x'],extrapolate=False)
    x=inverse(r)
    assert np.all(np.isfinite(x))
    q=RectBivariateSpline(field['y'],field['x'],field['q'],kx=5,ky=5)
    Q=RectBivariateSpline(field['y'],field['x'],field['Q'],kx=5,ky=5)
    # Match the old wall audit's interpolation when comparing its metrics.
    q3=RectBivariateSpline(field['y'],field['x'],field['q'])
    measured=1/q3.ev(r*0,x)
    j=jacobian(r,b)
    t,w=leggauss(80);t=(t+1)/2;w=w/2
    X=r*np.sum(w*jacobian(r[:,None]*t,b),axis=1)
    train=(r>=-.5)&(r<=.35);hold=~train
    rel=lambda a,z,mask:float(np.linalg.norm((a-z)[mask])/np.linalg.norm(z[mask]))
    idx=int(np.argmax(field['q'][0]));xs=field['x']
    peak=brentq(lambda xx:q.ev(0,xx,dy=1),xs[idx-1],xs[idx+1])
    ap=float(q.ev(0,peak,dx=1))
    mp=-float(q.ev(0,peak,dx=1,dy=1))/float(q.ev(0,peak,dy=2))
    cp=float(Q.ev(0,peak,dx=1))+mp*float(q.ev(0,peak))
    k=5/(4*mu*(5+d*d));predA=-d/(mu*(5+d*d))
    lift=WallLift(b,3*mu,mu,off,nodes=nodes)
    generated=lift.ridge()
    assert abs(generated['a1'][0]-k)<1e-11
    assert abs(generated['a1'][1]+predA)<1e-11
    assert abs(generated['a1'][2])<1e-11
    return dict(name=name,admissible=True,sourceKind='finite-time trajectory' if name=='tau20' else 'conditional finite-box root',
        savedRates=dict(lam=lam,mu=mu,v=v,densityOffset=off),d=d,bFromResonance=float(b),
        imposedLambda=3*mu,lambdaRelativeChange=float((3*mu-lam)/lam),
        wallErrors=dict(jacobianCore=rel(j,measured,train),jacobianWings=rel(j,measured,hold),
                        inversePositionWings=rel(X,x,hold)),
        predictedSlopes=dict(A=predA,m=k,C=0.),
        measuredSlopes=dict(A=ap,m=mp,C=cp),
        generatedSecondOrder=generated,
        slopeRelativeErrors=dict(A=abs((predA-ap)/ap),m=abs((k-mp)/mp)),
        quadratureNodes=nodes,
        protocol='mu and density offset from saved source; impose lambda=3mu; b is fixed algebraically by d, with no wall, ridge, or velocity fit. Saved field is a comparison, not a certified whole-plane steady solution.')


def second_density_zero_scan():
    rows=[]
    for d in np.linspace(.1,2.3,12):
        b=-4*d/np.sqrt(25-4*d*d)
        values=[]
        for mu in [1.,2.]:
            w=WallLift(b,3*mu,mu,mu*d,nodes=48)
            values.append(float(w.second(0.,2)[2]))
        # a2''/2 = B0(d)+B2(d)/mu^2; C(y)=-a2''/2*y^2+...
        B2=(values[0]-values[1])/.75;B0=values[0]-B2
        ratio=-B2/B0 if B0 else float('nan')
        row=dict(d=float(d),b=float(b),B0=B0,B2=B2,
                 positiveMuSquared=ratio if ratio>0 else None)
        if ratio>0:
            mu=np.sqrt(ratio);w=WallLift(b,3*mu,mu,mu*d,nodes=80)
            row['independentResidual']=float(w.second(0.,2)[2])
        rows.append(row)
    return rows


def constant_density_special_root():
    """Select second/third jets, then test the unused fourth jet."""
    def parameters(d,nodes):
        b=-4*d/np.sqrt(25-4*d*d)
        vals=[WallLift(b,3*mu,mu,mu*d,nodes=nodes).second(0.,2)[2]
              for mu in [1.,2.]]
        B2=(vals[0]-vals[1])/.75;B0=vals[0]-B2
        return b,np.sqrt(-B2/B0)
    def third(d,nodes):
        b,mu=parameters(d,nodes)
        return -WallLift(b,3*mu,mu,mu*d,nodes=nodes).third(0.,2)[2]
    roots=[]
    for nodes in [24,32,48]:
        d=brentq(lambda d:third(d,nodes),1.3,1.5,xtol=5e-14)
        b,mu=parameters(d,nodes)
        roots.append(dict(nodes=nodes,d=float(d),b=float(b),mu=float(mu),
                          v=float(mu*d),thirdResidual=float(third(d,nodes))))
    fourths=[]
    for nodes in [16,20,24]:
        ridge=WallLift(b,3*mu,mu,mu*d,nodes=nodes).fourth_ridge()
        fourths.append(dict(nodes=nodes,ridge=ridge))
    return dict(roots=roots,fourthChecks=fourths,
                scope='A local numerical root in the stated bracket. The independent nonzero fourth coefficient rejects constant peak density at this candidate; no uniqueness or global exclusion claim.')


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    result=dict(symbolic=certificate())
    result['audits']=[audit(n) for n in ['L4_129x65','L4_257x129','L8_257x129_continued','tau20']]
    result['L4quadratureRefinement']=audit('L4_257x129',nodes=96)
    result['constantPeakDensitySecondGate']=second_density_zero_scan()
    (OUT/'resonant_skew_branch.json').write_text(json.dumps(result,indent=2)+'\n')
    special=constant_density_special_root()
    (OUT/'constant_density_special_root.json').write_text(json.dumps(special,indent=2)+'\n')
    source=result['audits'][1];mu=source['savedRates']['mu'];off=source['savedRates']['densityOffset']
    fourth=[dict(nodes=n,ridge=WallLift(source['bFromResonance'],3*mu,mu,off,nodes=n).fourth_ridge())
            for n in [20,28]]
    (OUT/'resonant_L4_fourth.json').write_text(json.dumps(fourth,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
