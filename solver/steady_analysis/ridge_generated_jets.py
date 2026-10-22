"""Generate local ridge orders from a finite wall family, not ridge fits.

All density derivatives are analytic Taylor arithmetic; one-dimensional
Gauss-Jacobi integrals impose regularity at Q=-v/mu. No velocity or ridge data
enters the generation. Comparisons use saved conditional roots separately.
"""
import json
from pathlib import Path
import time
import numpy as np
import sympy as sp
from scipy.special import roots_jacobi
from scipy.interpolate import RectBivariateSpline, PchipInterpolator, make_interp_spline
from scipy.optimize import brentq
from numpy.polynomial import Polynomial, Chebyshev
from ridge_trace_equation import load
from normal_lift_algebra import curl_coefficients

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/ridge_generated_20260928'


def mul(a,b):
    n=min(len(a),len(b))
    return np.stack([sum(a[k]*b[i-k] for k in range(i+1)) for i in range(n)])


def inv(a):
    out=[1/a[0]]
    for i in range(1,len(a)):
        out.append(-sum(a[k]*out[i-k] for k in range(1,i+1))/a[0])
    return np.stack(out)


def power(a,n):
    if n<0:
        return power(inv(a),-n)
    out=np.zeros_like(a);out[0]=1
    for _ in range(n):
        out=mul(out,a)
    return out


def derivative(a,n=1):
    out=a
    for _ in range(n):
        out=np.stack([(k+1)*out[k+1] for k in range(len(out)-1)])
    return out


class Jet:
    def __init__(self,a):self.a=np.asarray(a)
    def other(self,b):
        if isinstance(b,Jet):return b.a
        c=np.zeros_like(self.a);c[0]=b;return c
    def __add__(self,b):return Jet(self.a+self.other(b))
    __radd__=__add__
    def __neg__(self):return Jet(-self.a)
    def __sub__(self,b):return self+-b
    def __rsub__(self,b):return -self+b
    def __mul__(self,b):return Jet(mul(self.a,self.other(b)))
    __rmul__=__mul__
    def __truediv__(self,b):return Jet(mul(self.a,inv(self.other(b))))


def compose(a,b,N):
    out=np.zeros(N+1);base=np.zeros(N+1);base[0]=1
    for val in a[:N+1]:
        out+=val*base;base=mul(base,b)
    return out


def ridge_series(coefficients,N):
    """The r=0 wall maximum is normalized; coefficients[n] is a_n(r)."""
    c=np.zeros(N+1)
    def curve(deriv):
        ans=np.zeros(N+1)
        for n,a in enumerate(coefficients):
            part=compose(derivative(a,deriv),c,N)
            ans[n:]+=part[:N+1-n]
        return ans
    for n in range(1,N+1):
        c[n]=-curve(2)[n]/(6*coefficients[0][3])
    j=curve(1);m=curve(0)
    t=np.zeros(N+1);t[0]=1
    for n in range(1,N+1):t[n]=(j[n]-sum(t[k]*t[n-k] for k in range(1,n)))/2
    yy=np.zeros(N+1)
    for n in range(1,N+1):yy[n]=((1. if n==1 else 0.)-compose(t,yy,N)[n])/t[1]
    mt=compose(m,yy,N)
    return dict(C_ySeries=c.tolist(),A_ySeries=inv(j).tolist(),m_ySeries=m.tolist(),
        t_ySeries=t.tolist(),m_tSeries=mt.tolist(),C_tSeries=compose(c,yy,N).tolist(),
        Y_tSeries=yy.tolist(),quadraticToLinear=float(mt[2]/mt[1]),
        cubicToLinear=float(mt[3]/mt[1]) if N>=3 else None,
        twoPowerCubicGate=float((mt[1]-2*mt[2]+3*mt[3])/mt[1]) if N>=3 else None)


class WallLift:
    def __init__(self,b,lam,mu,v,nodes=96):
        self.b=b;self.lam=lam;self.mu=mu;self.v=v;self.d=v/mu
        self.nodes=nodes;self.quadrature={}

    def rule(self,gamma):
        if gamma not in self.quadrature:
            xx,ww=roots_jacobi(self.nodes,0,gamma-1)
            self.quadrature[gamma]=((xx+1)/2,ww/2**gamma)
        return self.quadrature[gamma]

    def J(self,r,order):
        r=np.asarray(r);b=self.b;u0=-b/2;h=np.sqrt(1+u0*u0);u=u0+h*r
        den=np.zeros((order+1,)+r.shape);num=np.zeros_like(den)
        den[0]=1+u*u;num[0]=h*(2*u+b)
        if order>=1:den[1]=2*u*h;num[1]=2*h*h
        if order>=2:den[2]=h*h
        g=mul(num,inv(den))
        out=[(1+u*u)/(1+u0*u0)*np.exp(b*(np.arctan(u)-np.arctan(u0)))]
        for n in range(order):
            out.append(sum(g[k]*out[n-k] for k in range(n+1))/(n+1))
        return np.stack(out)

    def first(self,r,order):
        r=np.asarray(r);t,w=self.rule(2*self.lam/self.mu-2)
        inner=-self.d+(r[...,None]+self.d)*t
        ji=self.J(inner,order)
        reciprocal=power(ji,-3)
        for k in range(order+1):reciprocal[k]*=t**k
        integral=np.sum(reciprocal*w,axis=-1)
        return mul(power(self.J(r,order),2),integral)/self.mu

    def E1(self,r,order):
        r=np.asarray(r);mu=self.mu;lam=self.lam
        jj=self.J(r,order+2);kk=self.first(r,order+2)
        j=jj[:order+1];jp=derivative(jj)[:order+1];jpp=derivative(jj,2)
        k=kk[:order+1];kp=derivative(kk)[:order+1];kpp=derivative(kk,2)
        s=np.zeros_like(j);s[0]=r+self.d
        if order>=1:s[1]=1
        one=np.zeros_like(j);one[0]=1
        return (mu*mul(one+mul(k,k),mul(s,mul(jp,jp))-mul(j,jp)-mul(s,mul(j,jpp)))
            +(4*mu-2*lam)*mul(mul(j,j),mul(k,kp))
            +2*mu*mul(s,mul(j,mul(k,mul(kp,jp))))
            +2*mu*mul(s,mul(mul(j,j),mul(k,kpp)))
            -3*mu*mul(s,mul(mul(j,j),mul(kp,kp)))+2*mul(j,kp))

    def second(self,r,order):
        r=np.asarray(r);t,w=self.rule(4*self.lam/self.mu-3)
        inner=-self.d+(r[...,None]+self.d)*t
        integrand=mul(self.E1(inner,order),power(self.J(inner,order),-6))
        for k in range(order+1):integrand[k]*=t**k
        integral=np.sum(integrand*w,axis=-1)
        return mul(power(self.J(r,order),3),integral)/(2*self.mu)

    def E2(self,r,order):
        # E(y) through y^2 with a3 set to zero. Exact polynomial in local jets.
        r=np.asarray(r)
        J=self.J(r,order+2);a1=self.first(r,order+3);a2=self.second(r,order+2)
        take=lambda a,n=0:Jet(derivative(a,n)[:order+1])
        zero=Jet(np.zeros_like(J[:order+1]))
        s=np.zeros_like(J[:order+1]);s[0]=r+self.d
        if order>=1:s[1]=1
        return curl_coefficients([take(J),take(a1,1),take(a2,1)],
            [take(J,1),take(a1,2),take(a2,2)],[take(J,2),take(a1,3),zero],
            [take(a1),2*take(a2)],[take(a1,1),2*take(a2,1)],[2*take(a2)],
            Jet(s),self.lam,self.mu,2)[2].a

    def third(self,r,order):
        r=np.asarray(r);t,w=self.rule(6*self.lam/self.mu-4)
        inner=-self.d+(r[...,None]+self.d)*t
        integrand=mul(self.E2(inner,order),power(self.J(inner,order),-7))
        for k in range(order+1):integrand[k]*=t**k
        return mul(power(self.J(r,order),4),np.sum(integrand*w,axis=-1))/(3*self.mu)

    def third_ridge(self):
        j=self.J(0.,4);a0=np.r_[0.,j/np.arange(1,6)]
        return ridge_series([a0,self.first(0.,4),self.second(0.,3),self.third(0.,2)],3)

    def E3(self,r,order):
        r=np.asarray(r)
        J=self.J(r,order+2)
        aa=[self.first(r,order+3),self.second(r,order+3),self.third(r,order+2)]
        take=lambda a,n=0:Jet(derivative(a,n)[:order+1])
        zero=Jet(np.zeros_like(J[:order+1]))
        s=np.zeros_like(J[:order+1]);s[0]=r+self.d
        if order>=1:s[1]=1
        return curl_coefficients([take(J)]+[take(a,1) for a in aa],
            [take(J,1)]+[take(a,2) for a in aa],
            [take(J,2)]+[take(a,3) for a in aa[:-1]]+[zero],
            [(i+1)*take(a) for i,a in enumerate(aa)],
            [(i+1)*take(a,1) for i,a in enumerate(aa)],
            [(i+1)*(i+2)*take(a) for i,a in enumerate(aa[1:])],
            Jet(s),self.lam,self.mu,3)[3].a

    def fourth(self,r,order):
        r=np.asarray(r);t,w=self.rule(8*self.lam/self.mu-5)
        inner=-self.d+(r[...,None]+self.d)*t
        integrand=mul(self.E3(inner,order),power(self.J(inner,order),-8))
        for k in range(order+1):integrand[k]*=t**k
        return mul(power(self.J(r,order),5),np.sum(integrand*w,axis=-1))/(4*self.mu)

    def fourth_ridge(self):
        j=self.J(0.,5);a0=np.r_[0.,j/np.arange(1,7)]
        return ridge_series([a0,self.first(0.,5),self.second(0.,4),self.third(0.,3),self.fourth(0.,2)],4)

    def ridge(self):
        j=self.J(0.,3);a1=self.first(0.,3);a2=self.second(0.,2)
        c1=-a1[2]
        c2=-(3*j[3]*c1*c1+6*a1[3]*c1+2*a2[2])/2
        p=a1[1];q=a2[1]-c1*c1
        t1=p/2;t2=q/2-p*p/8
        m1=a1[0]+c1;m2=c2+p*c1+a2[0]
        mt=m1/t1;mt2=(m2-mt*t2)/t1**2
        return dict(wallJ=list(map(float,j)),a1=list(map(float,a1)),a2=list(map(float,a2)),
            C_ySeries=[0.,float(c1),float(c2)],A_ySeries=[1.,float(-p),float(p*p-q)],
            m_ySeries=[0.,float(m1),float(m2)],t_ySeries=[1.,float(t1),float(t2)],
            m_tSeries=[0.,float(mt),float(mt2)],quadraticToLinear=float(mt2/mt),
            scope='Ordinary Taylor coefficients. All generated from the wall Jacobian and supplied rates; no ridge fitting.')


def symbolic():
    J,Jp,Jpp,K,Kp,Kpp,s,lam,mu=sp.symbols('J Jp Jpp K Kp Kpp s lam mu')
    direct=curl_coefficients([J,Kp],[Jp,Kpp],[Jpp,0],[K],[Kp],[],s,lam,mu,1)[1]
    compact=(mu*(1+K*K)*(s*Jp**2-J*Jp-s*J*Jpp)+(4*mu-2*lam)*J**2*K*Kp
             +2*mu*s*J*K*Kp*Jp+2*mu*s*J**2*K*Kpp-3*mu*s*J**2*Kp**2+2*J*Kp)
    assert sp.expand(direct-compact)==0
    return dict(secondOrderSourceResidual=0,E1=str(compact),empiricalCubicQuadraticToLinear='13/5')


def measured(name):
    p,_=load(name)
    q=RectBivariateSpline(p['y'],p['x'],p['q'],kx=5,ky=5)
    Q=RectBivariateSpline(p['y'],p['x'],p['Q'],kx=5,ky=5)
    ys=p['y'][(p['y']>=0)&(p['y']<=2.)]
    mm=[];aa=[];cc=[]
    for y in ys:
        j=int(np.argmax(q(y,p['x'])[0]))
        m=brentq(lambda x:q.ev(y,x,dx=0,dy=1),p['x'][j-1],p['x'][j+1],xtol=1e-13)
        mm.append(m);aa.append(q.ev(y,m));cc.append(Q.ev(y,m))
    mm=np.array(mm);aa=np.array(aa);cc=np.array(cc);tt=np.sqrt(aa[0]/aa)
    fits=[]
    for ymax in [.5,1.,2.]:
        mask=ys<=ymax
        for degree in [3,4,5,6]:
            poly=Polynomial.fit(tt[mask]-1,mm[mask]-mm[0],degree).convert()
            fy=Polynomial.fit(ys[mask],mm[mask]-mm[0],degree).convert()
            fa=Polynomial.fit(ys[mask],aa[mask],degree).convert()
            fc=Polynomial.fit(ys[mask],cc[mask]-cc[0],degree).convert()
            fits.append(dict(ymax=ymax,degree=degree,m_tSeries=poly.coef[:4].tolist(),
                quadraticToLinear=float(poly.coef[2]/poly.coef[1]),m_ySeries=fy.coef[:3].tolist(),
                A_ySeries=fa.coef[:3].tolist(),C_ySeries=fc.coef[:3].tolist()))
    np.savez_compressed(OUT/f'{name}_ridge.npz',y=ys,m=mm,A=aa,C=cc,t=tt)
    return dict(name=name,localPolynomialFits=fits)


class DiagnosticWall(WallLift):
    """A flexible measured-wall control, NOT a proposed symbolically closed family."""
    def __init__(self,name,degree,inverse_degree=5,nodes=96):
        p,_=load(name)
        lam,mu,v,_=p['rates'];super().__init__(0.,lam,mu,v,nodes)
        r=np.linspace(-1.,.7,301)
        inverse=(PchipInterpolator(p['Q'][0],p['x']) if inverse_degree==3
                 else make_interp_spline(p['Q'][0],p['x'],k=5))
        q=RectBivariateSpline(p['y'],p['x'],p['q'],kx=5,ky=5)
        measuredJ=1/q.ev(r*0,inverse(r))
        # Lock the three normalization jets, then fit only measured wall values.
        basis=np.polynomial.chebyshev.chebvander((2*r+.3)/1.7,degree)
        co=np.linalg.lstsq(r[:,None]**3*basis,np.log(measuredJ)-r*r,rcond=None)[0]
        self.logJ=Polynomial([0,0,1])+Polynomial([0,0,0,1])*Chebyshev(co,domain=[-1.,.7]).convert(kind=Polynomial)
        self.fitError=float(np.linalg.norm(np.exp(self.logJ(r))-measuredJ)/np.linalg.norm(measuredJ))

    def J(self,r,order):
        r=np.asarray(r)
        from math import factorial
        f=[self.logJ.deriv(k)(r)/factorial(k) for k in range(order+1)]
        out=[np.exp(f[0])]
        for n in range(order):
            out.append(sum((k+1)*f[k+1]*out[n-k] for k in range(n+1))/(n+1))
        return np.stack(out)


def main():
    start=time.monotonic();OUT.mkdir(parents=True,exist_ok=True)
    fits=json.loads((ROOT/'runs/ridge_density_closure_20260928/asymmetric_wall_data.json').read_text())
    rows=[]
    for row in fits:
        if row['name']=='tau20':continue
        lam,mu,v,_=row['rates'];b=row['fit']['b']
        entry=dict(name=row['name'],b=b,rates=row['rates'],generated=[])
        for resonance in [False,True]:
            for nodes in [48,96,192]:
                model=WallLift(b,3*mu if resonance else lam,mu,v,nodes)
                val=model.ridge();val.update(nodes=nodes,enforcedResonance=resonance)
                entry['generated'].append(val)
                print(row['name'],'resonance',resonance,'nodes',nodes,'m_t',val['m_tSeries'],
                      'ratio',val['quadraticToLinear'],flush=True)
        entry['measured']=measured(row['name']);rows.append(entry)
        entry['thirdOrder']=[]
        if row['name']=='L4_257x129':
            for resonance in [False,True]:
                for nodes in [24,32,48]:
                    lift=WallLift(b,3*mu if resonance else lam,mu,v,nodes)
                    third=lift.third_ridge()
                    second=lift.ridge()
                    assert np.max(abs(np.array(third['m_tSeries'][:3])-second['m_tSeries']))<1e-10
                    entry['thirdOrder'].append(dict(nodes=nodes,enforcedResonance=resonance,**third))
                    print('third ridge',resonance,nodes,third['m_tSeries'],flush=True)
    exact=symbolic()
    # Independent exact expression check for first normal order at rational b.
    r=sp.Symbol('r')
    rr=np.array([-.5,0.,.5])
    # The old exact benchmark is unnormalized; use its normalized b=0 counterpart.
    old=json.loads((ROOT/'runs/ridge_symbolic_20260928/elementary_lift.json').read_text())
    L=sp.Symbol('L');expr=sp.sympify(old['rows'][0]['coefficient'],locals={'r':r,'L':L})
    expected=sp.lambdify((r,L),expr,'numpy')(rr,np.arctan(rr)+np.pi/4)
    actual=WallLift(0.,3.,1.,1.,192).first(rr,0)[0]
    assert np.max(np.abs(actual-expected))<1e-11
    exact['symmetricFirstLiftMaxError']=float(np.max(np.abs(actual-expected)))
    expr2=sp.sympify(old['rows'][1]['coefficient'],locals={'r':r,'L':L})
    exact['symmetricSecondLiftDerivativeErrors']=[]
    second=WallLift(0.,3.,1.,1.,192).second(rr,2)
    from math import factorial
    for n in range(3):
        want=sp.lambdify((r,L),expr2,'numpy')(rr,np.arctan(rr)+np.pi/4)
        err=float(np.max(np.abs(factorial(n)*second[n]-want)))
        assert err<2e-9
        exact['symmetricSecondLiftDerivativeErrors'].append(err)
        expr2=sp.diff(expr2,r)+sp.diff(expr2,L)/(1+r*r)
    controls=[]
    for name in ['L4_257x129','L8_257x129_continued']:
        for degree in [6,10,14]:
            for inverse_degree in [3,5]:
                diag=DiagnosticWall(name,degree,inverse_degree)
                prediction=diag.ridge()
                controls.append(dict(name=name,polynomialDegree=degree,inverseDegree=inverse_degree,
                    wallRelativeL2=diag.fitError,prediction=prediction,
                    scope='Measured-wall diagnostic with more coefficients; not a new closure family.'))
                print('wall control',name,degree,inverse_degree,diag.fitError,
                    prediction['m_tSeries'],flush=True)
    result=dict(symbolic=exact,rows=rows,seconds=time.monotonic()-start,
        measuredWallControls=controls,
        scope='Finite second/third-order local predictions, not a summed solution or global closure.')
    (OUT/'generated_ridge_jets.json').write_text(json.dumps(result,indent=2)+'\n')
    print('seconds',result['seconds'],flush=True)


if __name__=='__main__':main()
