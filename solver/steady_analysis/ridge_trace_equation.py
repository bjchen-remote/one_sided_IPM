"""Exact ridge-trace reduction, checked against saved velocity fields.

This closes a velocity trace once A,m,C are supplied. It does NOT close the
two-dimensional profile or replace the centered half-plane Biot-Savart law.
"""
import json
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.interpolate import CubicSpline, RectBivariateSpline, make_interp_spline
from scipy.integrate import solve_ivp
from scipy.io import loadmat
from construct_profile import load_profile

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/ridge_parameter_20260928'


def symbolic_check():
    A,g,b,mp,kp,U,sigma,lam,mu=sp.symbols('A g b mp kp U sigma lambda mu')
    k=mp-b;D=1+k*k;nu=2*lam-mu
    ap=((1+mp*mp)*g+b*kp)/D
    force=nu+b*(2*(lam-mu)*k-A-mu*sigma*g/A)/D
    Up=ap*U+force
    # Recover partial_x V2 using the curl/divergence identities on the ridge.
    R=((2*mp-b)*Up+kp*U+mu*b-mu*sigma*g/A-A-2*lam*mp)/(1+mp*mp)
    gradient_transport=sp.factor(b*R-Up+nu+g*U)
    assert gradient_transport==0
    p=2*lam+mp*R-Up
    # V1=mu*sigma/A+k*U, C'=A*b, k'=kp.
    V1prime=mu*b-mu*sigma*g/A+kp*U+k*Up
    curl_trace=sp.factor(V1prime-(R+A+mp*p))
    assert curl_trace==0
    return dict(gradientTransportResidual=str(gradient_transport),
        curlTraceResidual=str(curl_trace),leadingCoefficient='1+(m_prime-C_prime/A)^2 >= 1',
        scope='Exact necessary ridge identities for a smooth monotone steady profile; not sufficient for a full profile.')


def peak(qfield,x,y):
    cs=CubicSpline(x,qfield(y,x)[0])
    roots=cs.derivative().roots(extrapolate=False)
    roots=roots[(roots>=x[0])&(roots<=x[-1])]
    xx=np.r_[x[0],roots,x[-1]]
    m=float(xx[np.argmax(cs(xx))])
    return m,float(cs(m))


def load(name):
    folder=ROOT/'runs/steady_20260928_tau20'
    p=load_profile(folder,name)
    if name=='tau20':
        d=loadmat(folder/'run_analysis.mat',simplify_cells=True)
        meta=d['metadata']
        psi=d['psi']/(meta['P']*meta['ell']**2)
    else:
        psi=np.load(folder/name/'state.npz')['psi']
    return p,RectBivariateSpline(p['y'],p['x'],psi)


def audit(name,stride,degree=3):
    p,psi=load(name)
    q=RectBivariateSpline(p['y'],p['x'],p['q'])
    Q=RectBivariateSpline(p['y'],p['x'],p['Q'])
    lam,mu,drift,offset=p['rates']
    ymax=min(2.,float(p['y'][-1]))
    native=p['y'][(p['y']>=0)&(p['y']<=ymax)]
    ys=np.unique(np.r_[native[::stride],native[-1],ymax])
    pairs=np.array([peak(q,p['x'],float(y)) for y in ys])
    ms,amps=pairs.T
    cc=Q.ev(ys,ms)
    M=make_interp_spline(ys,ms,k=degree)
    A=make_interp_spline(ys,amps,k=degree)
    C=make_interp_spline(ys,cc,k=degree)
    yy=np.linspace(0,ymax,401)
    mm=np.array([peak(q,p['x'],float(y))[0] for y in yy])
    actualU=psi.ev(yy,mm,dx=0,dy=1)+lam*yy
    originU1=float(-psi.ev(0,0,dx=1,dy=0))
    actualV1=-psi.ev(yy,mm,dx=1,dy=0)-originU1+lam*mm+drift

    results=[];predictions={}
    for choice in ['measured_density','quarter_density','constant_density']:
        def parameters(y):
            a,ap=A(y),A(y,1);mp=M(y,1);mpp=M(y,2)
            if choice=='measured_density':
                c,cp,cpp=C(y),C(y,1),C(y,2)
            elif choice=='quarter_density':
                c=C(0)+y*a/4
                cp=(a+y*ap)/4
                cpp=(2*ap+y*A(y,2))/4
            else:
                c,cp,cpp=C(0),0.,0.
            g=ap/a;b=cp/a;k=mp-b;kp=mpp-cpp/a+b*g
            sigma=c+offset/mu;D=1+k*k
            linear=((1+mp*mp)*g+b*kp)/D
            force=2*lam-mu+b*(2*(lam-mu)*k-a-mu*sigma*g/a)/D
            return a,g,b,k,kp,sigma,linear,force
        def rhs(y,u):
            *_,linear,force=parameters(y)
            return linear*u+force
        sol=solve_ivp(rhs,[0,ymax],[0.],rtol=3e-10,atol=3e-12,
                      dense_output=True,max_step=.02)
        if not sol.success: raise RuntimeError(sol.message)
        pred=sol.sol(yy)[0]
        a,g,b,k,kp,sigma,linear,force=parameters(yy)
        v1=mu*sigma/a+k*pred
        mp=M(yy,1);Dp=1+mp*mp
        r1=((2*mp-b)*linear+kp)/Dp
        r0=((2*mp-b)*force+mu*b-mu*sigma*g/a-a-2*lam*mp)/Dp
        p1=mp*r1-linear
        p0=2*lam+mp*r0-force
        r=r1*pred+r0
        # Differentiate the scalar coefficient curves, not the measured velocity.
        r1p=CubicSpline(yy,r1)(yy,1);r0p=CubicSpline(yy,r0)(yy,1)
        p1p=CubicSpline(yy,p1)(yy,1);p0p=CubicSpline(yy,p0)(yy,1)
        rp=(r1p+linear*r1)*pred+r0p+force*r1
        pp=(p1p+linear*p1)*pred+p0p+force*p1
        crossingVelocity=mu*sigma/a-b*pred
        curvature=-(2*r*A(yy,1)+a*(b*pp-(1+k*mp)*rp)/Dp)/crossingVelocity
        measuredCurvature=q.ev(yy,mm,dx=0,dy=2)
        e1=-b*measuredCurvature+2*A(yy,1)*r1+a*(b*(p1p+linear*p1)-(1+k*mp)*(r1p+linear*r1))/Dp
        e0=mu*sigma/a*measuredCurvature+2*A(yy,1)*r0+a*(b*(p0p+force*p1)-(1+k*mp)*(r0p+force*r1))/Dp
        geometricU=-e0/e1
        rows=[]
        for low,high in [(0,.25),(.25,1.),(1.,2.)]:
            mask=(yy>low)&(yy<=high)
            delta=pred[mask]-actualU[mask]
            rows.append(dict(heightBand=[low,high],verticalVelocityRelativeL2=float(np.linalg.norm(delta)/np.linalg.norm(actualU[mask])),
                verticalVelocityLinf=float(abs(delta).max()),
                horizontalVelocityRelativeL2=float(np.linalg.norm((v1-actualV1)[mask])/np.linalg.norm(actualV1[mask])),
                curvatureRelativeL2=float(np.linalg.norm((curvature-measuredCurvature)[mask])/np.linalg.norm(measuredCurvature[mask])),
                predictedNegativeCurvatureFraction=float(np.mean(curvature[mask]<0)),
                algebraicTraceRelativeL2=float(np.linalg.norm((geometricU-actualU)[mask])/np.linalg.norm(actualU[mask])),
                minimumAbsEliminationCoefficient=float(abs(e1[mask]).min())))
        a0,g0,b0,k0,kp0,sig0,_,force0=parameters(0.)
        wall_gate=2*(lam-mu)*k0-a0-mu*sig0*g0/a0
        results.append(dict(densityModel=choice,comparisons=rows,
            wallCompatibilityDefect=float(wall_gate),
            predictedWallUprime=float(force0),requiredWallUprime=float(2*lam-mu)))
        predictions[choice]=pred
        predictions[choice+'_predicted_curvature']=curvature
        predictions[choice+'_algebraic_U']=geometricU
    # Direct field diagnostics distinguish steady defects from interpolation effects.
    aq=q.ev(yy,mm);Qy=Q.ev(yy,mm,dx=1,dy=0);Qx=Q.ev(yy,mm,dx=0,dy=1)
    transport=actualV1*Qx+actualU*Qy-mu*Q.ev(yy,mm)-offset
    poisson=-psi.ev(yy,mm,dx=2,dy=0)-psi.ev(yy,mm,dx=0,dy=2)-aq
    v1x=-psi.ev(yy,mm,dx=1,dy=1)+lam
    v2x=psi.ev(yy,mm,dx=0,dy=2)
    qx=q.ev(yy,mm,dx=0,dy=1);qy=q.ev(yy,mm,dx=1,dy=0)
    gradient=actualV1*qx+actualU*qy+(v1x-mu)*aq+v2x*Qy
    result=dict(name=name,stride=stride,interpolationDegree=degree,nativeHeightNodes=len(ys),rates=p['rates'].tolist(),
        curves='A,m,C from saved ridge samples; interpolating splines; no fitted velocity',
        sourceKind='finite-time trajectory' if name=='tau20' else 'conditional finite-box root',
        directRidgeFieldDiagnostics=dict(transportRms=float(np.sqrt(np.mean(transport**2))),
            poissonRms=float(np.sqrt(np.mean(poisson**2))),
            differentiatedTransportRms=float(np.sqrt(np.mean(gradient**2))),
            qMinusDensityDerivativeRms=float(np.sqrt(np.mean((aq-Qx)**2)))),
        models=results)
    np.savez_compressed(OUT/f'trace_{name}_stride{stride}_degree{degree}.npz',y=yy,measuredU=actualU,
        measuredV1=actualV1,**predictions)
    print(name,'stride',stride,'degree',degree,'U relative errors',[r['comparisons'][1]['verticalVelocityRelativeL2'] for r in results],flush=True)
    return result


if __name__=='__main__':
    OUT.mkdir(exist_ok=True)
    exact=symbolic_check()
    rows=[audit(name,stride) for name,stride in [('L4_129x65',1),('L4_257x129',1),
        ('L4_257x129',2),('L8_257x129_continued',1),('tau20',1)]]
    rows += [audit('L4_257x129',stride,5) for stride in [1,2]]
    (OUT/'ridge_trace.json').write_text(json.dumps(dict(symbolic=exact,audits=rows),indent=2)+'\n')
