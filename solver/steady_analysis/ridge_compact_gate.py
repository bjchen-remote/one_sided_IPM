"""Apply the first two exact ridge equations to amplitude-parametric fits.

No model coefficient is tuned using the PDE or velocity data.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline
from ridge_trace_equation import load, peak

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/ridge_parameter_20260928'


def power_basis(v,p,n,density=False):
    def power(e):
        factor=1.
        for i in range(n): factor*=e-i
        return factor*v**(e-n)
    if density:
        b1=power(p-1)/(2**(p-1)-1)
        b2=-2*power(-1)
        if n==0:
            b1-=1/(2**(p-1)-1);b2+=2
    else:
        b1=power(p)/(2**p-1)
        b2=power(1)
        if n==0:
            b1-=1/(2**p-1);b2-=1
    return b1,b2


def audit(frame,powers):
    profile,psi=load('tau20')
    q=RectBivariateSpline(profile['y'],profile['x'],profile['q'])
    lam,mu,drift,offset=profile['rates']
    ell=frame['ell']
    # Use the exact saved normalization P, not a sampled wall maximum.
    from scipy.io import loadmat
    P=loadmat(ROOT/'runs/geometry_loop_20260928/frame_tau20.mat',simplify_cells=True)['P']
    factors={'Y':1/ell,'m':1/ell,'Wc':1/ell,'C':1/(P*ell)}
    def curve(key,v,n=0):
        power=powers[key]
        model=min(frame['models'],key=lambda m:abs(m['p']-power))['parts'][key]
        theta=model['theta']
        b1,b2=power_basis(v,power,n,key=='C')
        ans=frame['anchors'][key]*(theta*b1+(1-theta)*b2)
        if n==0 and key!='Y': ans+=frame[key+'0']
        return ans*factors[key]
    def pars(v):
        a=frame['A0']/P/v;av=-a/v
        yv=curve('Y',v,1);yvv=curve('Y',v,2)
        mp=curve('m',v,1)/yv
        mpv=curve('m',v,2)/yv-mp*yvv/yv
        b=curve('C',v,1)/(a*yv)
        bv=curve('C',v,2)/(a*yv)-b*(av/a+yvv/yv)
        k=mp-b;kp=(mpv-bv)/yv;g=av/(a*yv)
        sigma=curve('C',v)+offset/mu;D=1+k*k
        r1=((mp+k)*g+kp)/D
        r0=(2*(lam-mu)*k-a-mu*sigma*g/a)/D
        linear=g+b*r1;force=2*lam-mu+b*r0
        return dict(A=a,Yv=yv,mp=mp,b=b,k=k,kp=kp,g=g,sigma=sigma,
                    r1=r1,r0=r0,linear=linear,force=force)
    def rhs(v,u):
        z=pars(v)
        return z['Yv']*(z['linear']*u+z['force'])
    sol=solve_ivp(rhs,[1,8],[0.],rtol=1e-10,atol=1e-12,dense_output=True,max_step=.025)
    if not sol.success: raise RuntimeError(sol.message)
    v=np.geomspace(1.001,8,301);z=pars(v);u=sol.sol(v)[0]
    r=z['r1']*u+z['r0'];up=z['linear']*u+z['force']
    # Complex step gives derivatives of the explicit analytic coefficient functions.
    step=1e-20;zc=pars(v+1j*step)
    r1p=np.imag(zc['r1'])/step/z['Yv'];r0p=np.imag(zc['r0'])/step/z['Yv']
    gp=np.imag(zc['g'])/step/z['Yv']
    rp=(r1p+z['linear']*z['r1'])*u+r0p+z['force']*z['r1']
    pp=z['kp']*r+z['k']*rp-gp*u-z['g']*up
    T=mu*z['sigma']/z['A']-z['b']*u
    numerator=2*r*z['A']*z['g']+z['A']*(z['b']*pp-(1+z['k']*z['mp'])*rp)/(1+z['mp']**2)
    desiredCurvature=-2*z['A']/curve('Wc',v)**2
    curvatureDefect=T*desiredCurvature+numerator
    y=curve('Y',v);m=curve('m',v)
    measuredPeak=np.array([peak(q,profile['x'],float(h)) for h in y])
    realU=psi.ev(y,measuredPeak[:,0],dx=0,dy=1)+lam*y
    rows=[]
    for lo,hi in [(1,2),(2,4),(4,8)]:
        mask=(v>lo)&(v<=hi)
        scale=np.linalg.norm((T*desiredCurvature)[mask])+np.linalg.norm(numerator[mask])
        rows.append(dict(amplitudeCoordinate=[lo,hi],
            verticalVelocityRelativeL2=float(np.linalg.norm((u-realU)[mask])/np.linalg.norm(realU[mask])),
            curvatureIdentityRelativeDefect=float(np.linalg.norm(curvatureDefect[mask])/scale),
            minimumYprime=float(z['Yv'][mask].min()),minimumCurvatureWidth=float(curve('Wc',v)[mask].min())))
    first=pars(1.)
    wall=2*(lam-mu)*first['k']-first['A']-mu*first['sigma']*first['g']/first['A']
    return dict(powers=powers,wallCompatibilityDefect=float(wall),bands=rows,
        note='All coefficients inherited from trajectory geometry fits. A nonzero necessary-condition defect rejects these fixed coefficients as an exact steady skeleton, not the whole function family.')


if __name__=='__main__':
    data=json.loads((OUT/'amplitude_parameter.json').read_text())
    frame=data['frames'][-1]
    choices=[{key:p for key in ['Y','m','Wc','C']} for p in [1.5,5/3,2.]]
    choices.append(dict(Y=5/3,m=1.5,Wc=5/3,C=1.5))
    result=[audit(frame,choice) for choice in choices]
    (OUT/'compact_ridge_gates.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
