"""Test a finite density slaving law, with wall data and PDE holdouts.

C=C0+A*(chi*y+(p1*(t-1)+p2*(t-1)^2)/t), t=sqrt(A0/A).
This is a structural hypothesis, not a consequence of the steady equations.
Only p2 is fitted when p1 is fixed by the wall compatibility condition.
"""
import json
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.interpolate import RectBivariateSpline, make_interp_spline, CubicSpline
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from numpy.polynomial import Chebyshev
from ridge_trace_equation import load, peak

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/ridge_density_closure_20260928'


def error(pred,obs,mask):
    return float(np.linalg.norm((pred-obs)[mask])/np.linalg.norm(obs[mask]))


def symbolic():
    t,chi,p1,p2,Y,Y1,M1,lam,mu,v=sp.symbols('t chi p1 p2 Y Y1 M1 lambda mu v')
    P=(p1*(t-1)+p2*(t-1)**2)/t
    assert sp.factor(P-(p2*t+p1-2*p2+(p2-p1)/t))==0
    nu=2*lam-mu
    y1=2*(v*v+(lam-mu)*(nu-v*M1))/v
    c1=2*M1-nu/v
    seed=sp.factor(c1-chi*y1)
    # Derivatives of the composed density, used by the independent PDE audit.
    y=sp.symbols('y');a=sp.Function('A')(y)
    tt=1/sp.sqrt(a)
    poly=sp.lambdify(t,P,'sympy')
    c=a*(chi*y+poly(tt))
    Pt=sp.diff(P,t).subs(t,tt);Ptt=sp.diff(P,t,2).subs(t,tt)
    B=chi*y+poly(tt)-tt*Pt/2
    cp=sp.diff(a,y)*B+chi*a
    cpp=sp.diff(a,y,2)*B+2*chi*sp.diff(a,y)-sp.diff(a,y)**2*tt*(Pt-tt*Ptt)/(4*a)
    assert sp.simplify(sp.diff(c,y)-cp)==0
    assert sp.simplify(sp.diff(c,y,2)-cpp)==0
    return dict(law='C=C0+A*(chi*y+(p1*(t-1)+p2*(t-1)^2)/t), t=sqrt(A0/A)',
                poleFreeDomain='t>=1',p1FromIdealWall=str(seed),
                seedDefinitions='M1=dm/dt at t=1; y1=2*(v^2+(lambda-mu)*(nu-v*M1))/v; c1=2*M1-nu/v',
                derivativeChecks=[0,0],scope='Finite algebraic slaving candidate, not proven by PDE.')


def trajectory_fits():
    old=json.loads((ROOT/'runs/ridge_parameter_20260928/amplitude_parameter.json').read_text())
    wall=json.loads((ROOT/'runs/ridge_symbolic_20260928/checks.json').read_text())['data'][-1]
    lam,mu,v,_=wall['rates'];M=wall['empiricalRidgeWallSeed']['Mstar']
    Y1=wall['empiricalRidgeWallSeed']['predictedY_t']
    C1=10*M-(2*lam-mu)/v
    result=[]
    for f in old['frames']:
        vv=np.array([r['v'] for r in f['rows']]);t=np.sqrt(vv);z=t-1
        y=np.array([r['Y'] for r in f['rows']])/f['ell']
        a=1/vv
        c=(np.array([r['C'] for r in f['rows']])-f['C0'])/(f['A0']*f['ell'])
        train=(vv>1.1)&(vv<=4);w=1/(1+y)
        target=c/a-y/4
        rows=[]
        for denpow in [0,1,2]:
            basis=np.array([z/t**denpow,z*z/t**denpow]).T
            co=np.linalg.lstsq(basis[train]*w[train,None],target[train]*w[train],rcond=None)[0]
            pred=a*(y/4+basis@co)
            rows.append(dict(kind='two_density_fit_constants',denominatorPower=denpow,p1=float(co[0]),p2=float(co[1]),
                errors=[dict(vRange=[lo,hi],relativeL2=error(pred,c,(vv>lo)&(vv<=hi))) for lo,hi in [(1,1.1),(1.1,4),(4,8),(8,16)]]))
        if f['tau']==20:
            for denpow in [0,1,2]:
                p1=C1-Y1/4;b=z*z/t**denpow
                rhs=target-p1*z/t**denpow
                p2=np.dot(b[train]*w[train],rhs[train]*w[train])/np.dot(b[train]*w[train],b[train]*w[train])
                pred=a*(y/4+(p1*z+p2*z*z)/t**denpow)
                rows.append(dict(kind='wall_locked_one_density_constant',denominatorPower=denpow,p1=float(p1),p2=float(p2),
                    errors=[dict(vRange=[lo,hi],relativeL2=error(pred,c,(vv>lo)&(vv<=hi))) for lo,hi in [(1,1.1),(1.1,4),(4,8),(8,16)]]))
        result.append(dict(tau=f['tau'],models=rows,quarterErrors=[dict(vRange=[lo,hi],relativeL2=error(a*y/4,c,(vv>lo)&(vv<=hi))) for lo,hi in [(1.1,4),(4,8),(8,16)]]))
    return result


def pde_audit(name,degree=5,stride=1,smooth_degree=None):
    p,psi=load(name);lam,mu,drift,off=p['rates']
    interpolation=5 if smooth_degree is not None else 3
    q=RectBivariateSpline(p['y'],p['x'],p['q'],kx=interpolation,ky=interpolation)
    Q=RectBivariateSpline(p['y'],p['x'],p['Q'],kx=interpolation,ky=interpolation)
    ys=p['y'][(p['y']>=0)&(p['y']<=2.)][::stride]
    def take_peak(y):
        if smooth_degree is None:
            return peak(q,p['x'],float(y))
        j=int(np.argmax(q(y,p['x'])[0]))
        m=brentq(lambda x:q.ev(y,x,dx=0,dy=1),p['x'][j-1],p['x'][j+1],xtol=1e-13)
        return m,float(q.ev(y,m))
    ma=np.array([take_peak(float(y)) for y in ys]);ms,amps=ma.T
    def curve(values):
        if smooth_degree is None:
            return make_interp_spline(ys,values,k=degree)
        fit=Chebyshev.fit(ys,values,smooth_degree,domain=[0,ys[-1]])
        return lambda y,n=0:fit.deriv(n)(y)
    A=curve(amps);M=curve(ms);C=curve(Q.ev(ys,ms))
    yy=np.linspace(0, min(1.8,float(ys[-1])),601)
    a=A(yy);aa0=float(A(0));c0=float(C(0));aa1=float(A(0,1))
    Y1=-2*aa0/aa1
    k0=(aa0+(mu*c0+off)*aa1/aa0**2)/(2*(lam-mu))
    requiredCp0=aa0*(float(M(0,1))-k0)
    p1=(requiredCp0/aa0-.25)*Y1
    t=np.sqrt(aa0/a);z=t-1;b=z*z/t
    target=(C(yy)-c0)/a-yy/4-p1*z/t
    train=(yy>=.1)&(yy<=1.0)
    w=1/(1+yy)
    p2=np.dot(b[train]*w[train],target[train]*w[train])/np.dot(b[train]*w[train],b[train]*w[train])
    models=[]
    trueU=psi.ev(yy,M(yy),dx=0,dy=1)+lam*yy
    trueH=q.ev(yy,M(yy),dx=0,dy=2)
    out_arrays=dict(y=yy,measuredU=trueU,measuredCurvature=trueH,measuredC=C(yy))
    for model in ['measured','quarter','wall_corrected']:
        def density(y):
            a,ap,app=A(y),A(y,1),A(y,2)
            if model=='measured':
                return C(y),C(y,1),C(y,2)
            if model=='quarter':
                return c0+a*y/4,(a+ap*y)/4,(2*ap+app*y)/4
            t=np.sqrt(aa0/a)
            P=p2*t+p1-2*p2+(p2-p1)/t
            Pt=p2-(p2-p1)/t**2;Ptt=2*(p2-p1)/t**3
            B=y/4+P-t*Pt/2
            return c0+a*(y/4+P),ap*B+a/4,app*B+ap/2-ap*ap*t*(Pt-t*Ptt)/(4*a)
        def pars(y):
            a,ap=A(y),A(y,1);c,cp,cpp=density(y)
            g=ap/a;beta=cp/a;mp=M(y,1);k=mp-beta
            kp=M(y,2)-cpp/a+beta*g;S=(mu*c+off)/a
            r1=((mp+k)*g+kp)/(1+k*k)
            r0=(2*(lam-mu)*k-a-S*g)/(1+k*k)
            linear=g+beta*r1;force=2*lam-mu+beta*r0
            return a,g,beta,mp,k,S,r1,r0,linear,force
        def rhs(y,u):
            *_,linear,force=pars(y)
            return linear*u+force
        sol=solve_ivp(rhs,[0,yy[-1]],[0.],rtol=1e-10,atol=1e-12,dense_output=True,max_step=.01)
        assert sol.success
        u=sol.sol(yy)[0]
        a,g,beta,mp,k,S,r1,r0,linear,force=pars(yy)
        pc1=k*r1-g;pc0=mu+k*r0
        r=r1*u+r0
        rp=(CubicSpline(yy,r1)(yy,1)+linear*r1)*u+CubicSpline(yy,r0)(yy,1)+force*r1
        pp=(CubicSpline(yy,pc1)(yy,1)+linear*pc1)*u+CubicSpline(yy,pc0)(yy,1)+force*pc1
        crossing=S-beta*u
        h=-(2*r*A(yy,1)+a*(beta*pp-(1+k*mp)*rp)/(1+mp*mp))/crossing
        c=density(yy)[0]
        z0=pars(0.)
        wallGate=2*(lam-mu)*z0[4]-z0[0]-z0[5]*z0[1]
        rows=[]
        for lo,hi in [(.1,.5),(.5,1.),(1.,1.5),(1.5,float(yy[-1]))]:
            mask=(yy>lo)&(yy<=hi)
            rows.append(dict(heightBand=[lo,hi],densityRelativeL2=error(c-c0,C(yy)-c0,mask),
                verticalVelocityRelativeL2=error(u,trueU,mask),curvatureRelativeL2=error(h,trueH,mask),
                negativeCurvatureFraction=float(np.mean(h[mask]<0)),minimumCrossingVelocity=float(crossing[mask].min())))
        models.append(dict(model=model,wallCompatibilityDefect=float(wallGate),bands=rows))
        out_arrays[model+'_U']=u;out_arrays[model+'_C']=c;out_arrays[model+'_curvature']=h
    suffix=f'_smooth{smooth_degree}' if smooth_degree is not None else ''
    np.savez_compressed(OUT/f'{name}_degree{degree}_stride{stride}{suffix}.npz',**out_arrays)
    return dict(name=name,degree=degree,stride=stride,smoothDegree=smooth_degree,rates=p['rates'].tolist(),
        p1=float(p1),p2=float(p2),requiredCp0=float(requiredCp0),
        parameterSource='p1 from exact first wall compatibility for measured A,m; p2 from C only on y=.1..1; velocities and curvature held out.',
        sourceKind='finite-time trajectory' if name=='tau20' else 'conditional finite-box root',models=models)


if __name__=='__main__':
    OUT.mkdir(parents=True,exist_ok=True)
    result=dict(symbolic=symbolic(),trajectories=trajectory_fits(),
        pde=[pde_audit(n) for n in ['L4_257x129','L8_257x129_continued','tau20']])
    result['pde'] += [pde_audit('L4_257x129',degree=3),pde_audit('L4_257x129',degree=5,stride=2)]
    result['smoothedPde']=[pde_audit('L4_257x129',smooth_degree=n) for n in [10,14,18]]
    (OUT/'density_slaving.json').write_text(json.dumps(result,indent=2)+'\n')
    for row in result['pde']+result['smoothedPde']:
        print(row['name'],row['degree'],row['stride'],'smooth',row['smoothDegree'],'coefficients',row['p1'],row['p2'],flush=True)
        for m in row['models']:
            print(m['model'],'wall',m['wallCompatibilityDefect'],'band .5..1',m['bands'][1],flush=True)
