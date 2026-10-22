"""Correction-free mother-profile research: exact gates and saved-wall audits.

Only new outputs in runs/exact_master_20260928 are written. No PDE evolution,
no modification of saved fields, and no transverse correction modes are used.
The first two normal derivatives are generated from the steady equations,
then tested against the tangent/second tangent of ONE fixed mother family.
"""
from pathlib import Path
import json
import os

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import sympy as sp
from scipy.interpolate import CubicSpline, make_interp_spline
from scipy.integrate import cumulative_simpson, solve_ivp
from scipy.special import roots_jacobi
from geometry_measure import candidate, trajectory_frame

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'runs/exact_master_20260928'


def exact_checks():
    # Derive both wall equations directly by differentiating transport twice.
    x, y, lam, mu, v = sp.symbols('x y lambda mu v')
    Q, U, r, t = [sp.Function(n)(x) for n in ['Q', 'U', 'r', 't']]
    q = sp.diff(Q, x)
    V1 = U + q*y + (sp.diff(r,x)-sp.diff(U,x,2))*y*y/2
    V2 = (2*lam-sp.diff(U,x))*y - sp.diff(q,x)*y*y/2
    rho = Q+r*y+t*y*y/2
    E = V1*sp.diff(rho,x)+V2*sp.diff(rho,y)-mu*rho-v
    one = U*sp.diff(r,x)+(2*lam-mu-sp.diff(U,x))*r+q*q
    two = (U*sp.diff(t,x)+(4*lam-mu-2*sp.diff(U,x))*t
           -sp.diff(U,x,2)*q+3*q*sp.diff(r,x)-sp.diff(q,x)*r)
    assert sp.expand(sp.diff(E,y).subs(y,0)-one) == 0
    assert sp.expand(sp.diff(E,y,2).subs(y,0)-two) == 0

    # Generic geometry pullback: no special master or data used here.
    z = sp.symbols('z')
    A, W, k, c, m = [sp.Function(n)(y) for n in ['A','W','k','c','m']]
    f, g0, g1, h0, h1 = [sp.Function(n)(z) for n in ['f','g0','g1','h0','h1']]
    G = g0+k*g1
    H = h0+k*h1
    X = m+W*G
    rho = c+A*W*H
    rules = {sp.diff(h0,z):f*sp.diff(g0,z), sp.diff(h1,z):f*sp.diff(g1,z)}
    physical_y = lambda e: sp.diff(e,y)-sp.diff(X,y)/sp.diff(X,z)*sp.diff(e,z)
    r_actual = physical_y(rho).subs(rules)
    basis = [-A*f, A*(H-f*G), A*W*(h1-f*g1), W*H, sp.Integer(1)]
    theta = [m,W,k,A,c]
    r_expected = sum(b*sp.diff(p,y) for b,p in zip(basis,theta))
    assert sp.simplify(r_actual-r_expected) == 0
    T = sp.diff(X,y)
    curvature = (2*sp.diff(A,y)*sp.diff(W,y)*H
        +2*sp.diff(A*W,y)*sp.diff(k,y)*h1
        -2*sp.diff(A,y)*f*T+A*sp.diff(f,z)*T*T/sp.diff(X,z)
        -2*A*sp.diff(W,y)*sp.diff(k,y)*f*g1)
    t_expected = sum(b*sp.diff(p,y,2) for b,p in zip(basis,theta))+curvature
    t_actual = physical_y(r_actual).subs(rules)
    assert sp.simplify(t_actual-t_expected) == 0

    # The free master ODE is a rearrangement of the whole first wall equation.
    ff, fp, hh, H1, gg, T0, a, b, cc, kk, N = sp.symbols(
        'f fp h H1 g T a b c k N')
    B = cc+(a+b)*hh+kk*H1
    rr = B-T0*ff
    wall = N/ff*(a*ff-T0*fp/gg)+(2*(lam-mu)+N*fp/(gg*ff**2))*rr+ff**2
    compact = a*N+2*(lam-mu)*rr+ff**2+N*fp*(B-2*T0*ff)/(gg*ff**2)
    assert sp.factor(wall-compact) == 0
    return dict(wall_first='exact zero', wall_second='exact zero',
        tangent_identity='exact zero', second_tangent_identity='exact zero',
        free_master_equation='exact zero',
        restrictions='Finite normal-order identities only; global centered Biot-Savart selection is separate.')


def quadrature(power, n):
    z, w = roots_jacobi(n, 0, power-1)
    return (z+1)/2, w/(2**power)


class Wall:
    """Positive wall q; Q is its integral so Q_x=q is exact to spline accuracy."""
    def __init__(self, name, stride=1, density_nodes=8001):
        p = trajectory_frame(20) if name == 'tau20' else candidate(name)
        mask = abs(p['x']) <= (6 if name == 'tau20' else p['x'][-1])
        xx, qq = p['x'][mask][::stride], p['q'][0,mask][::stride]
        self.qx = CubicSpline(xx,qq,extrapolate=False)
        ant = self.qx.antiderivative()
        self.Qx = lambda x: ant(x)-ant(0.)
        self.lam, self.mu, self.v, _ = p['rates']
        self.lo, self.hi = max(float(xx[0]),-5.), min(float(xx[-1]),5.)
        dense = np.linspace(self.lo,self.hi,density_nodes)
        ss = self.Qx(dense)+self.v/self.mu
        assert ss[0] < 0 < ss[-1] and np.all(np.diff(ss)>0)
        self.qspline = make_interp_spline(ss,self.qx(dense),k=5)
        self.slo, self.shi = ss[0],ss[-1]
        self.xspline = make_interp_spline(ss,dense,k=5)
        roots = self.qx.derivative().roots(extrapolate=False)
        roots = roots[(roots>self.lo)&(roots<self.hi)]
        self.m = float(roots[np.argmax(self.qx(roots))])
        self.A = float(self.qx(self.m))
        self.c = float(self.Qx(self.m))
        half = CubicSpline(xx,qq-self.A/2).roots(extrapolate=False)
        self.left = self.m-max(half[half<self.m])
        self.right = min(half[half>self.m])-self.m
        self.name, self.stride = name,stride

    def r_rs(self, s, n):
        s = np.asarray(s)
        z,w = quadrature(2*(self.lam-self.mu)/self.mu,n)
        st = s[...,None]*z
        assert st.min() >= self.slo and st.max() <= self.shi
        qt = self.qspline(st)
        qst = self.qspline(st,1)
        I = (qt**3)@w
        Is = (3*qt**2*qst*z)@w
        q, qs = self.qspline(s),self.qspline(s,1)
        return -I/(self.mu*q), -Is/(self.mu*q)+I*qs/(self.mu*q*q)

    def forcing2(self,s,n):
        q,qs,qss = [self.qspline(s,i) for i in range(3)]
        r,rs = self.r_rs(s,n)
        return (-self.mu*q*qs-self.mu*s*q*qss+self.mu*s*qs*qs
                -3*q*q*rs+q*qs*r)

    def jets(self,s,n):
        s = np.asarray(s)
        r,_ = self.r_rs(s,n)
        z,w = quadrature((4*self.lam-3*self.mu)/self.mu,n)
        t = np.empty_like(s)
        for start in range(0,len(s),50):
            sl = slice(start,start+50)
            st = s[sl,None]*z
            q = self.qspline(st)
            R = self.forcing2(st,n)
            t[sl] = ((q*q*R)@w)/(self.mu*self.qspline(s[sl])**2)
        return r,t

    def geometry(self,x,beta):
        A,m = self.A,self.m
        if beta == 0:
            W=(self.left+self.right)/2
            k=0.
            xi=(x-m)/W
            G=xi; g=np.ones_like(x); Gk=np.zeros_like(x)
            H1=np.zeros_like(x)
        else:
            W=beta*(self.left+self.right)/(2*np.sinh(beta))
            k=(self.right-self.left)/(self.right+self.left)/np.tanh(beta/2)
            assert abs(k)<1
            theta=np.arctanh(k)
            xi=(np.arcsinh(np.sinh(theta)+beta*(x-m)/W*np.cosh(theta))-theta)/beta
            G=(np.sinh(beta*xi)+k*(np.cosh(beta*xi)-1))/beta
            g=np.cosh(beta*xi)+k*np.sinh(beta*xi)
            Gk=(np.cosh(beta*xi)-1)/beta
            # Integrate H1 using the saved wall, with no outside evaluation.
            zz=np.linspace(min(float(xi.min()),0),max(float(xi.max()),0),6001)
            GX=(np.sinh(beta*zz)+k*(np.cosh(beta*zz)-1))/beta
            hv=self.qx(m+W*GX)/A*np.sinh(beta*zz)
            hc=cumulative_simpson(hv,x=zz,initial=0.)
            hs=CubicSpline(zz,hc)
            H1=hs(xi)-hs(0.)
        f=self.qx(x)/A
        fp=self.qx(x,1)*W*g/A
        h=(self.Qx(x)-self.c)/(A*W)
        B=np.column_stack([-A*f,A*(h-f*G),A*W*(H1-f*Gk),W*h,np.ones_like(x)])
        return dict(A=A,m=m,W=W,k=k,c=self.c,xi=xi,G=G,g=g,Gk=Gk,
                    f=f,fp=fp,h=h,H1=H1,B=B)


def relative(a,b):
    return float(np.linalg.norm(a-b)/np.linalg.norm(b))


def fit_scaled(B,target,mask):
    scale=np.linalg.norm(B[mask],axis=0)
    valid=scale>1e-13
    value=np.zeros(B.shape[1])
    coef,_,rank,sv=np.linalg.lstsq(B[mask][:,valid]/scale[valid],target[mask],rcond=None)
    value[valid]=coef/scale[valid]
    return value,dict(rank=int(rank),scaledCondition=float(sv[0]/sv[-1]))


def tangent_audit(name,stride=1,n=64,beta=1.5):
    wall=Wall(name,stride)
    x=np.linspace(-2.,2.,801)
    s=wall.Qx(x)+wall.v/wall.mu
    r,t=wall.jets(s,n)
    g=wall.geometry(x,beta)
    B=g['B']; train=abs(x)<=1.2; hold=~train
    results=[]
    for region,mask in [('full',np.ones(len(x),bool)),('central_train',train)]:
        p,condition=fit_scaled(B,r,mask)
        mp,wp,kp,ap,cp=p
        T=mp+wp*g['G']+g['W']*kp*g['Gk']
        bend=(2*ap*wp*g['h']+2*(ap*g['W']+g['A']*wp)*kp*g['H1']
              -2*ap*g['f']*T+g['A']*g['fp']*T*T/(g['W']*g['g'])
              -2*g['A']*wp*kp*g['f']*g['Gk'])
        p2,condition2=fit_scaled(B,t-bend,mask)
        first=B@p; second=B@p2+bend
        results.append(dict(region=region,firstParameters=p.tolist(),secondParameters=p2.tolist(),
            firstRelativeL2=relative(first,r),secondRelativeL2=relative(second,t),
            firstHoldoutRelativeL2=relative(first[hold],r[hold]),
            secondHoldoutRelativeL2=relative(second[hold],t[hold]),
            firstAbsoluteMax=float(abs(first-r).max()),secondAbsoluteMax=float(abs(second-t).max()),
            condition=condition,condition2=condition2))
    base={k:float(g[k]) for k in ['m','A','W','k','c']}
    result=dict(name=name,stride=stride,nodes=n,beta=beta,base=base,
        rates=[float(wall.lam),float(wall.mu),float(wall.v)],
        stagnationX=float(wall.xspline(0)),xRange=[-2,2],trainingRange=[-1.2,1.2],
        source='Saved wall q; first/second normal jets generated by continuous steady wall equations.',
        fits=results,
        restriction='Second gate is diagnostic when first gate is nonzero; this is not a proof excluding the unknown limiting master.')
    return result,dict(x=x,r=r,t=t,B=B,wall_q=wall.qx(x),**{k:g[k] for k in ['xi','f','fp','h','H1']})


def free_master_local(audit):
    """Generate a NEW fixed mother from the exact first wall ODE, near its peak.

    No y-dependent shape modes. Wall geometry is rescaled to W=A=1 for this
    construction; rates are those of the selected conditional root.
    """
    lam,mu,v=audit['rates']
    base=audit['base']; W=base['W']; A=base['A']
    # Normalize spatial scale by W, q by A, density by A W; shift peak density.
    lam,mu,v=lam/A,mu/A,(mu*base['c']+v)/(A*A*W)
    k0=audit['base']['k']; beta=audit['beta']
    mp,wp,kp,ap,cp=audit['fits'][0]['firstParameters']
    a=ap*W/A; b=wp; kap=kp*W
    # f''(0) need not be -2 after the half-width normalization.
    wall=Wall(audit['name'],audit['stride'])
    curvature=-float(wall.qx(base['m'],2))*W*W/A
    d=2*(lam-mu); nu=2*lam-mu
    r0=-(1+a*v)/d
    # f'(0)=0 and f''(0)=-curvature determine m'(0),c'(0) exactly.
    m1=r0-nu*a/(curvature*v)
    c1=2*r0-nu*a/(curvature*v)

    def quantities(z,state):
        f,h,H1=state
        G=(np.sinh(beta*z)+k0*(np.cosh(beta*z)-1))/beta
        g=np.cosh(beta*z)+k0*np.sinh(beta*z)
        Gk=(np.cosh(beta*z)-1)/beta
        T=m1+b*G+kap*Gk
        B=c1+(a+b)*h+kap*H1
        N=mu*h+v
        return g,T,B,N

    def rhs(z,state):
        f,h,H1=state
        g,T,B,N=quantities(z,state)
        fp=-g*f*f*(a*N+d*(B-T*f)+f*f)/(N*(B-2*T*f))
        return [fp,f*g,f*np.sinh(beta*z)]

    def regularity(z,state):
        _,T,B,N=quantities(z,state)
        return abs(N*(B-2*T*state[0]))-1e-6
    regularity.terminal=True; regularity.direction=-1
    sols=[solve_ivp(rhs,[0,end],[1.,0.,0.],method='DOP853',rtol=2e-11,atol=2e-13,
                    dense_output=True,max_step=.003,events=regularity) for end in [-1.5,1.5]]
    z=np.linspace(max(sols[0].t[-1],-.8),min(sols[1].t[-1],1.2),601)
    values=np.column_stack([sols[0 if zz<0 else 1].sol(zz) for zz in z])
    generated=values[0]
    G=(np.sinh(beta*z)+k0*(np.cosh(beta*z)-1))/beta
    old=wall.qx(base['m']+W*G)/A
    residual=[]
    for zz,state in zip(z,values.T):
        f,h,H1=state;g,T,B,N=quantities(zz,state)
        fp=rhs(zz,state)[0]
        residual.append(a*N+d*(B-T*f)+f*f+N*fp*(B-2*T*f)/(g*f*f))
    return dict(parameters=dict(lam=lam,mu=mu,v=v,beta=beta,k0=k0,a=a,b=b,
                        kprime=kap,mprime=m1,cprime=c1,peakCurvature=curvature),
        maxFirstWallIdentityDefect=float(max(abs(np.array(residual)))),
        comparisonXiRange=[float(z[0]),float(z[-1])],
        measuredWallRelativeL2=relative(generated,old),minimumMother=float(generated.min()),
        integration=[dict(end=float(s.t[-1]),success=bool(s.success),message=s.message,
            terminatedByDenominatorEvent=bool(len(s.t_events[0]))) for s in sols],
        restriction='A local master generated by the first wall identity. No second-order or full-height closure is claimed.'),dict(
        xi=z,generatedMother=generated,savedWallMother=old,h=values[1],H1=values[2])


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    checks=exact_checks()
    print('symbolic identities verified',flush=True)
    rows=[]
    settings=[('L4_257x129',1,64,b) for b in [0.,1.25,1.5,1.75,2.]]
    settings += [('L4_257x129',1,128,1.5),('L4_257x129',2,64,1.5),
                 ('L4_129x65',1,64,1.5),('L8_257x129_continued',1,64,1.5),('tau20',1,64,1.5)]
    seed=None
    for name,stride,n,beta in settings:
        report,arrays=tangent_audit(name,stride,n,beta)
        rows.append(report)
        if name=='L4_257x129' and stride==1 and n==128:
            seed=report
            np.savez_compressed(OUT/'wall_tangent_arrays.npz',**arrays)
        first=report['fits'][0]
        print(name,stride,n,beta,'first',first['firstRelativeL2'],'second',first['secondRelativeL2'],flush=True)
    construction,arrays=free_master_local(seed)
    np.savez_compressed(OUT/'generated_local_master.npz',**arrays)
    result=dict(symbolic=checks,audits=rows,localConstruction=construction,
        parameterOrder=['m_prime','W_prime','k_prime','A_prime','c_prime'])
    (OUT/'checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(construction,indent=2),flush=True)


if __name__=='__main__':
    main()
