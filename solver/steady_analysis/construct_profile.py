"""Construct and audit a curved, nonsymmetric hodograph using full profiles.

The wall Jacobian is positive by construction. The first two normal jets
come from the steady PDE. A single positive Padé parameter is fitted to
two-dimensional level sets, with taller rows reserved for prediction.
This is a finite-core construction, not a verified global steady state.
"""
from pathlib import Path
import json
import sys

import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
from scipy.io import loadmat, savemat
from scipy.interpolate import PchipInterpolator, RectBivariateSpline, CubicSpline
from scipy.optimize import brentq, least_squares, minimize_scalar
from scipy.integrate import cumulative_trapezoid

from fixed_point import fd_matrix


def json_save(path, data):
    path.write_text(json.dumps(data, indent=2)+'\n')


def wall_functions():
    r,beta,a,b,c,s=sp.symbols('r beta a b c s',real=True)
    ell=beta*sp.log(1+r*r/beta)
    ell+=a*r**3/(1+s*r*r)**sp.Rational(3,2)
    ell+=b*r**4/(1+s*r*r)**2+c*r**5/(1+s*r*r)**sp.Rational(5,2)
    J=sp.exp(ell)
    values=[J,J*sp.diff(ell,r),J*(sp.diff(ell,r,2)+sp.diff(ell,r)**2)]
    return sp.lambdify((r,beta,a,b,c,s),values,'numpy')


WALL=wall_functions()
tn,wn=leggauss(48);T=(tn+1)/2;W=wn/2


class Wall:
    def __init__(self, rates, parameters):
        self.lam,self.mu,self.v,_=rates
        self.beta=(self.lam/self.mu-1)/2
        self.d=self.v/self.mu
        self.parameters=np.asarray(parameters)

    def jacobian(self,r):
        return WALL(np.asarray(r),self.beta,*self.parameters)

    def F(self,r):
        r=np.asarray(r)
        return r*np.sum(W*self.jacobian(r[...,None]*T)[0],axis=-1)

    def first_jet(self,r):
        r=np.asarray(r);z=r+self.d
        J,Jp,Jpp=self.jacobian(r)
        jt,jtp,jtpp=self.jacobian(z[...,None]*T-self.d)
        gamma=2*(self.lam-self.mu)/self.mu
        h=np.sum(W*T**(gamma-1)/jt**3,axis=-1)
        hp=np.sum(W*(-3)*T**gamma*jtp/jt**4,axis=-1)
        hpp=np.sum(W*T**(gamma+1)*(12*jtp*jtp/jt**5-3*jtpp/jt**4),axis=-1)
        A=J*J*h/self.mu
        Ap=(2*J*Jp*h+J*J*hp)/self.mu
        App=(2*(Jp*Jp+J*Jpp)*h+4*J*Jp*hp+J*J*hpp)/self.mu
        return A,Ap,App

    def second_jet(self,r):
        r=np.asarray(r);z=r+self.d;J=self.jacobian(r)[0]
        rt=z[...,None]*T-self.d;zt=rt+self.d
        j,jp,jpp=self.jacobian(rt);A,Ap,App=self.first_jet(rt)
        mu,lam=self.mu,self.lam
        e1=mu*(1+A*A)*(zt*(j*jpp-jp*jp)+j*jp)/j**3
        e1+=4*mu*zt*A*Ap*jp/j**2-2*mu*zt*A*App/j
        e1+=(2*mu-4*lam)*A*Ap/j+Ap/j**2
        gamma=4*lam/mu-3
        return J**3/mu*np.sum(W*T**(gamma-1)*(-e1/2)/j**3,axis=-1)

    def jets(self,r):
        r=np.asarray(r)
        B=np.concatenate([self.second_jet(part) for part in np.array_split(r,max(1,len(r)//64))])
        return self.F(r),self.first_jet(r)[0],B


def load_profile(folder,name):
    if name=='tau20':
        d=loadmat(folder/'run_analysis.mat',simplify_cells=True)
        m=d['metadata'];a,p,ell=m['a'],m['P'],m['ell']
        rho=RectBivariateSpline(d['yFull'],d['xFull'],d['rho'])
        return dict(x=(d['xFull']-a)/ell,y=d['yFull']/ell,
                    Q=(d['rho']-rho(0,a)[0,0])/(p*ell),
                    q=d['rhoX']/p,rates=np.array(m['innerRates']))
    d=np.load(folder/name/'state.npz')
    return dict(x=d['x'],y=d['y'],Q=d['Q'],q=d['qx'],rates=d['rates'])


def inverse_rows(profile,r,y):
    x=profile['x']
    interp=RectBivariateSpline(profile['y'],x,profile['Q'])
    xq=np.linspace(max(x[0],-10),min(x[-1],12),10001)
    rows=[]
    for h in y:
        qq=interp(h,xq)[0]
        if np.any(np.diff(qq)<=0) or min(r)<qq[0] or max(r)>qq[-1]:
            raise ValueError(f'Inverse outside monotone saved region at y={h}')
        rows.append(PchipInterpolator(qq,xq)(r))
    return np.array(rows)


def fit_wall(profile):
    q=CubicSpline(profile['x'],profile['q'][0])
    r=np.linspace(-1.3,1.3,301)
    x=inverse_rows(profile,r,[0])[0]
    logj=np.log(1/q(x))
    def residual(par):
        wall=Wall(profile['rates'],par)
        J=wall.jacobian(r)[0]
        return np.r_[np.log(J)-logj,1e-4*np.asarray(par[:3])]
    fits=[least_squares(residual,start,bounds=([-8,-8,-8,.05],[8,8,8,20]),
                       max_nfev=400,ftol=1e-11,xtol=1e-11,gtol=1e-11)
          for start in ([-.5,0,0,1],[-1,0,0,.2],[0,0,0,4])]
    fit=min(fits,key=lambda a:np.linalg.norm(a.fun))
    wall=Wall(profile['rates'],fit.x)
    metrics=dict(parameters=fit.x.tolist(),beta=wall.beta,
        logJacobianRms=float(np.sqrt(np.mean(residual(fit.x)[:len(r)]**2))),
        inverseWallLinf=float(np.max(abs(wall.F(r)-x))),
        inverseWallRelativeL2=float(np.linalg.norm(wall.F(r)-x)/np.linalg.norm(x)),
        minimumJacobianOverMinus10To10=float(np.min(wall.jacobian(np.linspace(-10,10,4001))[0])))
    return wall,metrics


def x_trial(F,A,B,y,kappa):
    y=np.asarray(y)[:,None]
    return F[None,:]+y*A[None,:]+B[None,:]*y*y/(1+kappa*y)


def phi_trial(wall,r,y,kappa):
    F,A,B=wall.jets(r)
    h=np.asarray(y)[:,None]
    X=x_trial(F,A,B,y,kappa)
    if kappa<1e-7:
        H=h**3/3
    else:
        z=kappa*h
        # Avoid cancellation near the physical wall.
        H=(z*z/2-z+np.log1p(z))/kappa**3
        small=abs(z)<1e-3
        H[small]=(h**3/3-kappa*h**4/4+kappa**2*h**5/5-kappa**3*h**6/6)[small]
    I=h*F[None,:]+h*h*A[None,:]/2+H*B[None,:]
    D=fd_matrix(r,1,9)
    Iz=I@D.T
    phi=wall.lam*(2*I-h*X)+wall.v*h-wall.mu*(r[None,:]+wall.d)*Iz
    return X,phi


def residual_audit(wall,kappa,nr=321,ny=161):
    r=np.linspace(-1.25,.5,nr);y=np.linspace(0,2,ny)
    X,phi=phi_trial(wall,r,y,kappa)
    Dr=fd_matrix(r,1,9);Dy=fd_matrix(y,1,9)
    J=X@Dr.T;K=Dy@X
    pz=phi@Dr.T;py=Dy@phi
    fluxr=(1+K*K)/J*pz-K*py
    fluxy=J*py-K*pz
    E=-(fluxr@Dr.T+Dy@fluxy)/J-1/J
    transport=-py+wall.lam*(X-y[:,None]*K)+wall.v-wall.mu*(r[None,:]+wall.d)*J
    interior_r=(r>=r[0]+.05*np.ptp(r))&(r<=r[-1]-.05*np.ptp(r))
    interior_y=(y>=y[0]+.05*np.ptp(y))&(y<=y[-1]-.05*np.ptp(y))
    core=E[np.ix_(interior_y,interior_r)]
    return dict(nr=nr,ny=ny,minimumJ=float(J.min()),
        poissonResidualInteriorLinf=float(np.max(abs(core))),
        poissonResidualInteriorRms=float(np.sqrt(np.mean(core**2))),
        poissonResidualWallLinf=float(np.max(abs(E[0,4:-4]))),
        transportResidualInteriorLinf=float(np.max(abs(transport[4:-4,4:-4]))))


def construct(profile,out,name):
    wall,wall_metrics=fit_wall(profile)
    r=np.linspace(-1.25,.5,151)
    ymax=2 if name=='L4_257x129' else 4
    y=np.linspace(0,ymax,161)
    observed=inverse_rows(profile,r,y)
    F,A,B=wall.jets(r)
    train=y<=1
    def cost(kappa):
        return np.mean((x_trial(F,A,B,y[train],kappa)-observed[train])**2)
    fit=minimize_scalar(cost,bounds=(0,20),method='bounded',
                        options=dict(xatol=1e-8))
    kappa=float(fit.x)
    prediction=x_trial(F,A,B,y,kappa)
    linear=x_trial(F,A,B*0,y,0)
    taylor=x_trial(F,A,B,y,0)
    displacement=observed-observed[0]
    sv=np.linalg.svd(displacement,compute_uv=False)
    def errors(prediction,mask):
        delta=prediction[mask]-observed[mask]
        return dict(Linf=float(np.max(abs(delta))),rms=float(np.sqrt(np.mean(delta**2))),
                    relativeDisplacementL2=float(np.linalg.norm(delta)/np.linalg.norm(displacement[mask])))
    metrics=dict(name=name,rates=profile['rates'].tolist(),wall=wall_metrics,kappa=kappa,
        source='saved trajectory' if name=='tau20' else 'conditional frozen-exterior root',
        domain=dict(densityMin=float(r[0]),densityMax=float(r[-1]),heightMax=float(y[-1]),trainingHeightMax=1),
        displacementSecondSingularRatio=float(sv[1]/sv[0]),
        displacementRank2RelativeL2=float(np.linalg.norm(sv[2:])/np.linalg.norm(sv)),
        firstNormalJet=errors(linear,y>0),secondOrderTaylor=errors(taylor,y>0),
        curvedTraining=errors(prediction,train),
        curvedHeldout=errors(prediction,~train),
        residualAudits=[residual_audit(wall,kappa,161,81),residual_audit(wall,kappa,321,161)])
    np.savez_compressed(out/f'{name}.npz',r=r,y=y,observed=observed,
                        prediction=prediction,linear=linear,taylor=taylor,
                        F=F,A=A,B=B,parameters=wall.parameters,kappa=kappa,
                        rates=profile['rates'])
    print(name,json.dumps(metrics),flush=True)
    return metrics


def plots(folder,names):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(14,4.5),constrained_layout=True)
    for name,ax in zip(names,axes):
        d=np.load(folder/f'{name}.npz')
        for i in np.linspace(0,len(d['r'])-1,7).astype(int):
            line,=ax.plot(d['observed'][:,i],d['y'],lw=1.8)
            ax.plot(d['prediction'][:,i],d['y'],'--',color=line.get_color(),lw=1.1)
        ax.set(xlabel='x',ylabel='y',title=name)
        ax.grid(alpha=.2)
    axes[0].plot([],[],'k-',label='saved density contours')
    axes[0].plot([],[],'k--',label='transport-compatible construction')
    axes[0].legend(fontsize=8)
    fig.savefig(folder/'curved_level_sets.png',dpi=180)
    plt.close(fig)


def third_forcing():
    """Exact y^2 Poisson coefficient before adding the third normal jet."""
    z,y,lam,mu,v=sp.symbols('z y lam mu v')
    F=sp.Function('F')(z);A=sp.Function('A')(z);B=sp.Function('B')(z)
    X=F+A*y+B*y*y;J=sp.diff(X,z);K=sp.diff(X,y)
    phi=(lam*F+v-mu*z*sp.diff(F,z))*y-mu*z*sp.diff(A,z)*y*y/2
    phi+=(-lam*B-mu*z*sp.diff(B,z))*y**3/3
    dx=lambda expr:sp.diff(expr,z)/J
    dy=lambda expr:sp.diff(expr,y)-K*sp.diff(expr,z)/J
    E=-dx(dx(phi))-dy(dy(phi))-1/J
    coefficient=sp.factor(sp.diff(E,y,2).subs(y,0)/2)
    args=[z,lam,mu]+[sp.diff(F,z,i) for i in range(1,4)]
    args+=[sp.diff(A,z,i) for i in range(4)]+[sp.diff(B,z,i) for i in range(3)]
    assert v not in coefficient.free_symbols
    plain=sp.symbols('p0:'+str(len(args)))
    return sp.lambdify(plain,coefficient.xreplace(dict(zip(args,plain))),'numpy',cse=True)


class JetExpansion:
    def __init__(self,wall,degree,forcing):
        from numpy.polynomial import Chebyshev
        self.wall=wall
        # Include every segment from the transport stagnation density to the
        # audit interval. These are spectral interpolants, not exact formulas.
        self.domain=[min(-1.4,-wall.d-.15),.7]
        self.A=Chebyshev.interpolate(lambda r:wall.first_jet(r)[0],degree,self.domain)
        self.B=Chebyshev.interpolate(lambda r:wall.second_jet(r),degree,self.domain)
        self.forcing=forcing
        self.C=Chebyshev.interpolate(self.third,degree,self.domain)

    def third(self,r):
        w=self.wall;r=np.asarray(r);rt=(r[...,None]+w.d)*T-w.d
        z=rt+w.d;j,jp,jpp=w.jacobian(rt)
        e2=self.forcing(z,w.lam,w.mu,j,jp,jpp,
            *[self.A.deriv(n)(rt) for n in range(4)],
            *[self.B.deriv(n)(rt) for n in range(3)])
        gamma=6*w.lam/w.mu-4
        return w.jacobian(r)[0]**4/w.mu*np.sum(W*T**(gamma-1)*(-e2/3)/j**4,axis=-1)


def flux_audit(wall,r,y,X,I):
    Dr=fd_matrix(r,1,9);Dy=fd_matrix(y,1,9);h=y[:,None]
    phi=wall.lam*(2*I-h*X)+wall.v*h-wall.mu*(r[None,:]+wall.d)*(I@Dr.T)
    J=X@Dr.T;K=Dy@X;pz=phi@Dr.T;py=Dy@phi
    E=-(((1+K*K)/J*pz-K*py)@Dr.T+Dy@(J*py-K*pz))/J-1/J
    transport=-py+wall.lam*(X-h*K)+wall.v-wall.mu*(r[None,:]+wall.d)*J
    interior_r=(r>=r[0]+.05*np.ptp(r))&(r<=r[-1]-.05*np.ptp(r))
    interior_y=(y>=y[0]+.05*np.ptp(y))&(y<=y[-1]-.05*np.ptp(y))
    core=E[np.ix_(interior_y,interior_r)]
    return dict(minimumJ=float(J.min()),poissonLinf=float(abs(core).max()),
        poissonRms=float(np.sqrt(np.mean(core*core))),
        auditInteriorFraction=[.05,.95],
        wallPoissonLinf=float(abs(E[0,4:-4]).max()),
        transportLinf=float(abs(transport[4:-4,4:-4]).max()))


def variable_trial(wall,r,y,coefficients):
    """Positive density-dependent denominator; exact first and second jets."""
    F,A,B=wall.jets(r);h=y[:,None]
    k=np.exp(np.polynomial.chebyshev.chebval((r+.375)/.875,coefficients))[None,:]
    z=k*h
    H=(z*z/2-z+np.log1p(z))/k**3
    H=np.where(abs(z)<1e-3,h**3/3-k*h**4/4+k*k*h**5/5-k**3*h**6/6,H)
    X=F[None,:]+h*A[None,:]+h*h*B[None,:]/(1+k*h)
    I=h*F[None,:]+h*h*A[None,:]/2+H*B[None,:]
    return X,I


def refine(source,out):
    forcing=third_forcing();results=[]
    for name in ['L4_257x129','L8_257x129_continued']:
        profile=load_profile(source,name)
        saved=np.load(out/f'{name}.npz')
        wall=Wall(saved['rates'],saved['parameters'])
        r=np.linspace(-1.25,.5,161);y=np.linspace(0,2,81)
        observed=inverse_rows(profile,r,y)
        jets=[JetExpansion(wall,n,forcing) for n in (60,90)]
        A,B,C=jets[-1].A(r),jets[-1].B(r),jets[-1].C(r)
        k=-C/B;singular=abs(B)<1e-5
        third=dict(maximumCoefficientChange=float(abs(jets[0].C(r)-C).max()),
            coefficientMaximum=float(abs(C).max()),
            forcedKappaMin=float(k[~singular].min()),forcedKappaMax=float(k[~singular].max()),
            nonpositiveKappaFraction=float(np.mean(k[~singular]<=0)),
            polesBeforeHeight2Fraction=float(np.mean((k<-.5)&(~singular))))
        train=y<=1
        start=np.r_[np.log(float(saved['kappa'])),0.,0.,0.]
        fit=least_squares(lambda c:(variable_trial(wall,r,y[train],c)[0]-observed[train]).ravel(),
            start,bounds=([-4,-3,-3,-3],[4,3,3,3]),max_nfev=150,
            ftol=1e-10,xtol=1e-10,gtol=1e-10)
        X,I=variable_trial(wall,r,y,fit.x)
        def rel(mask):
            return float(np.linalg.norm((X-observed)[mask])/np.linalg.norm((observed-observed[0])[mask]))
        audit=flux_audit(wall,r,y,X,I)
        rf=np.linspace(-1.25,.5,321);yf=np.linspace(0,2,161)
        xf,If=variable_trial(wall,rf,yf,fit.x)
        entry=dict(name=name,thirdJet=third,variableDenominator=dict(coefficients=fit.x.tolist(),
            trainingRelativeDisplacement=rel(train),heldoutRelativeDisplacement=rel(~train),
            audits=[audit,flux_audit(wall,rf,yf,xf,If)]))
        # Check actual local improvement of the third jet on a narrow strip.
        local=[]
        for ymax in (.025,.05,.1,.2):
            yl=np.linspace(0,ymax,81);h=yl[:,None];F=wall.F(r)[None,:]
            pol=F+h*A+h*h*B+h**3*C
            integral=h*F+h*h*A/2+h**3*B/3+h**4*C/4
            check=flux_audit(wall,r,yl,pol,integral)
            check['heightMax']=ymax;local.append(check)
        entry['thirdOrderLocalAudits']=local
        results.append(entry);json_save(out/'refinement.json',results)
        np.savez_compressed(out/f'{name}_refined.npz',r=r,y=y,observed=observed,
                            prediction=X,A=A,B=B,C=C,forcedKappa=k,coefficients=fit.x)
        print(json.dumps(entry),flush=True)


def global_structure(source,out):
    """Measure the native full source, not a frozen-exterior Newton root."""
    from fourier_wall_probe import polar_strain
    profile=load_profile(source,'tau20')
    q=RectBivariateSpline(profile['y'],profile['x'],profile['q'])
    density=RectBivariateSpline(profile['y'],profile['x'],profile['Q'])
    alpha=profile['rates'][1]/profile['rates'][0]
    heights=np.r_[0,.25,.5,2.**np.arange(0,15)]
    ridge=[]
    for y in heights:
        fit=minimize_scalar(lambda x:-q(y,x)[0,0],
            bounds=(-max(2,y),max(2,3*y)),method='bounded')
        ridge.append([float(y),float(fit.x),float(-fit.fun),float(density(y,fit.x)[0,0])])
    ridge=np.asarray(ridge)
    radii=2.**np.arange(0,14)
    angles=np.linspace(0,np.pi,361)
    xx=radii[:,None]*np.cos(angles);yy=radii[:,None]*np.sin(angles)
    rays=q.ev(yy.ravel(),xx.ravel()).reshape(xx.shape)
    powers=np.log2(rays[1:]/rays[:-1])
    paired=lambda xx,yy:q.ev(yy,xx)-q.ev(yy,-xx)
    lx=min(-profile['x'][0],profile['x'][-1]);ly=profile['y'][-1]
    strains=[polar_strain(paired,lx,ly,order=n,
                         shell_edges=[0]+[2.**k for k in range(-2,19)]) for n in (48,96)]
    fractions=[]
    for radius in [1,2,4,8,16,32,64,128,256,1024,4096]:
        interior=sum(s['contribution'] for s in strains[-1]['shells'] if s['rmax']<=radius)
        fractions.append(dict(radius=radius,interiorStrain=interior,
            exteriorFraction=(strains[-1]['total']-interior)/strains[-1]['total']))
    nodes,weights=leggauss(96)
    theta=(nodes+1)*np.pi/2;weights=weights*np.pi/2
    angular_moment=[]
    for row in fractions:
        radius=row['radius']
        d2=2/np.pi*np.dot(weights,q.ev(radius*np.sin(theta),radius*np.cos(theta))*np.sin(2*theta))
        tail=strains[-1]['total']-row['interiorStrain']
        angular_moment.append(dict(radius=radius,sine2Coefficient=float(d2),
            homogeneousTailPrediction=float(-d2/(2*(1-alpha))),actualSavedTailStrain=tail,
            effectiveExponentFromTailMoment=float(1+d2/(2*tail))))

    # Test a translated and dilated wall shape. Fit only q below y=1;
    # both taller rows and the density primitive are independent checks.
    collapse=[]
    for name in ['tau20','L4_257x129','L8_257x129_continued']:
        p=load_profile(source,name);a=p['rates'][1]/p['rates'][0];d=p['rates'][2]/p['rates'][1]
        qq=RectBivariateSpline(p['y'],p['x'],p['q'])
        RR=RectBivariateSpline(p['y'],p['x'],p['Q'])
        qw=CubicSpline(p['x'],p['q'][0]);F=CubicSpline(p['x'],p['Q'][0])
        x=np.linspace(-2,3,151);y=np.linspace(0,2,81);X,Y=np.meshgrid(x,y)
        observed=qq(y,x);observedR=RR(y,x);train=Y<=1
        def prediction(par):
            tilt,kappa=par;H=1+kappa*Y;xi=(X-tilt*Y)/H
            return H**(a-1)*qw(xi),H**a*(F(xi)+d)-d
        fit=least_squares(lambda par:(prediction(par)[0]-observed)[train],
                          [1.4,1],bounds=([0,.001],[5,10]))
        predicted,predictedR=prediction(fit.x)
        collapse.append(dict(name=name,alpha=float(a),tilt=float(fit.x[0]),
            dilation=float(fit.x[1]),qTrainingRelativeL2=float(np.linalg.norm((predicted-observed)[train])/np.linalg.norm(observed[train])),
            qHeldoutRelativeL2=float(np.linalg.norm((predicted-observed)[~train])/np.linalg.norm(observed[~train])),
            densityRelativeL2=float(np.linalg.norm(predictedR-observedR)/np.linalg.norm(observedR))))

    # A common angular amplitude would be unchanged after this normalization.
    angular_changes=[]
    for i in range(len(radii)-1):
        base=radii[i]**(1-alpha)*rays[i]
        other=radii[i+1]**(1-alpha)*rays[i+1]
        angular_changes.append(dict(radius=float(radii[i]),
            doubledRadiusNormalizedQRelativeChange=float(np.linalg.norm(other-base)/np.linalg.norm(base))))
    result=dict(source='tau=20 finite-time native trajectory, inner coordinates',
        alphaFromInstantaneousRates=float(alpha),ridgeColumns=['y','xPeak','qPeak','QatPeak'],
        ridge=ridge.tolist(),strainByQuadratureOrder={'48':strains[0],'96':strains[1]},
        strainOrderChange=float(abs(strains[0]['total']-strains[1]['total'])),
        referenceMuMinusLambda=float(profile['rates'][1]-profile['rates'][0]),
        exteriorStrainFractions=fractions,homotheticTrials=collapse,angularScaleChanges=angular_changes,
        strainSelectingAngularMoment=angular_moment,
        localRidgeTiltFromInterpolatedDerivatives=float(-q.ev(0,0,dx=1,dy=1)/q.ev(0,0,dx=0,dy=2)),
        caveats=['Quadrature refinement does not bound interpolation or original discretization error.',
                 'Exterior fractions refer to the saved finite domain; its unsaved infinite tail is not included.',
                 'A q fit alone omits the independent height-dependent density primitive.',
                 'Neither the native field nor the conditional candidates verify an infinite-domain steady state.'])
    json_save(out/'global_structure.json',result)
    np.savez_compressed(out/'global_structure.npz',ridge=ridge,radii=radii,angles=angles,
                        rayQ=rays,rayExponents=powers)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(12,8.2),constrained_layout=True)
    xp=np.linspace(-4,12,321);yp=np.linspace(0,8,241)
    field=q(yp,xp);ax=axes[0,0]
    image=ax.pcolormesh(xp,yp,field,shading='auto',cmap='viridis',vmin=0,vmax=1)
    ax.contour(xp,yp,field,levels=[.1,.2,.3,.5,.7],colors='white',linewidths=.5,alpha=.6)
    inside=ridge[:,0]<=8
    ax.plot(ridge[inside,1],ridge[inside,0],'r.-',lw=1.3,label='gradient ridge')
    ax.set(xlabel='x',ylabel='y',title='Native q = Q_x at tau = 20');ax.legend(fontsize=8)
    fig.colorbar(image,ax=ax,label='q')
    ax=axes[0,1];law=collapse[0];xi=np.linspace(-2,2,401)
    for y in [0,.5,1,2,4]:
        H=1+law['dilation']*y
        values=q(y,law['tilt']*y+H*xi)[0]*H**(1-alpha)
        ax.plot(xi,values,label=f'y={y:g}')
    ax.set(xlabel='(x - s y) / (1 + k y)',ylabel='(1 + k y)^(1-alpha) q',
           title='Tilt + dilation captures most of the core shape')
    ax.legend(fontsize=8);ax.grid(alpha=.2)
    ax=axes[1,0]
    ax.semilogx([f['radius'] for f in fractions],[100*f['exteriorFraction'] for f in fractions],'o-')
    ax.set(xlabel='Radius from wall peak',ylabel='Strain from larger radii (%)',
           title='The outer source matters at the center');ax.grid(alpha=.2)
    ax=axes[1,1]
    for radius in [8,32,128,512]:
        i=np.flatnonzero(radii==radius)[0]
        ax.plot(angles*180/np.pi,radius**(1-alpha)*rays[i],label=f'r={radius}')
    ax.set(xlabel='Polar angle (degrees)',ylabel='r^(1-alpha) q',
           title='One homogeneous far-field shape is not yet established')
    ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.savefig(out/'global_structure.png',dpi=170);plt.close(fig)
    print(json.dumps(dict(strain=strains[-1]['total'],fractions=fractions,collapse=collapse)),flush=True)


def construction_identities(out):
    """Verify the three wall selections and first far-field correction."""
    x=sp.symbols('x');lam,mu,v=sp.symbols('lam mu v')
    w=sp.Function('w')(x);q=sp.Function('q')(x)
    E=(w+lam*x+v)*sp.diff(q,x)-(mu-lam-sp.diff(w,x))*q
    q3,w3=sp.symbols('q3 w3')
    at_peak={q.subs(x,0):1,sp.diff(q,x).subs(x,0):0,
        sp.diff(q,x,2).subs(x,0):-2,sp.diff(q,x,3).subs(x,0):q3,
        w.subs(x,0):0,sp.diff(w,x).subs(x,0):mu-lam,
        sp.diff(w,x,2).subs(x,0):2*v,sp.diff(w,x,3).subs(x,0):4*mu-v*q3}
    wall_checks=[sp.simplify(sp.diff(E,x,n).subs(x,0).subs(at_peak)) for n in range(3)]
    assert wall_checks==[0,0,0]
    r,theta,alpha=sp.symbols('r theta alpha',positive=True)
    A=sp.Function('A')(theta);B=sp.Function('B')(theta)
    zeta=r**alpha*A;psi=r**(alpha+1)*B
    radial=-sp.diff(psi,theta)/r;angular=sp.diff(psi,r)
    nonlinear=sp.simplify(radial*sp.diff(zeta,r)+angular/r*sp.diff(zeta,theta))
    N=-alpha*sp.diff(B,theta)*A+(alpha+1)*B*sp.diff(A,theta)
    assert sp.simplify(nonlinear-r**(2*alpha-1)*N)==0
    correction=r**(2*alpha-1)*N/(lam*(1-alpha))
    assert sp.simplify(lam*r*sp.diff(correction,r)-lam*alpha*correction+nonlinear)==0
    Dx=lambda f:sp.cos(theta)*sp.diff(f,r)-sp.sin(theta)*sp.diff(f,theta)/r
    lap=lambda f:sp.diff(f,r,2)+sp.diff(f,r)/r+sp.diff(f,theta,2)/r**2
    D=alpha*sp.cos(theta)*A-sp.sin(theta)*sp.diff(A,theta)
    assert sp.simplify((-lap(psi)-Dx(zeta))/r**(alpha-1)+sp.diff(B,theta,2)+(alpha+1)**2*B+D)==0
    json_save(out/'construction_identities.json',dict(wallSelection='three exact symbolic zeros',
        farFieldPoisson='exact symbolic zero',farFieldFirstCorrection='exact symbolic zero',
        claimsExcluded=['far-field remainder estimate','series convergence','global steady existence']))


def ridge_row(x,q):
    """Global horizontal maximum, including all cubic critical points/edges."""
    spline=CubicSpline(x,q)
    roots=spline.derivative().roots(extrapolate=False)
    roots=roots[np.isfinite(roots)]
    points=np.r_[x[0],roots[(roots>=x[0])&(roots<=x[-1])],x[-1]]
    peak=float(points[np.argmax(spline(points))])
    i=int(np.argmax(q));parabolic=float(x[i]);dx=np.nan
    if 0<i<len(x)-1:
        local=x[i-1:i+2]-x[i]
        poly=np.polyfit(local,q[i-1:i+2],2)
        vertex=-poly[1]/(2*poly[0]) if poly[0]<0 else np.nan
        if np.isfinite(vertex) and local[0]<=vertex<=local[-1]:
            parabolic=float(x[i]+vertex)
        dx=float(max(x[i]-x[i-1],x[i+1]-x[i]))
    return [peak,float(spline(peak)),parabolic,float(x[i]),dx,
            float(peak==x[0] or peak==x[-1])]


def ridge_fit(source,out):
    """Fit the native tau=20 gradient ridge in main's saved X,Y coordinates."""
    data=loadmat(source/'run_analysis.mat',simplify_cells=True)
    x,y,q=data['xFull'],data['yFull'],data['rhoX']
    extracted=np.array([ridge_row(x,row) for row in q])
    origin=float(extracted[0,0]);position=extracted[:,0]-origin
    alternative=extracted[:,2]-extracted[0,2]
    field=RectBivariateSpline(y,x,q)

    def curve(heights):
        return np.array([ridge_row(x,field(h,x)[0])[0] for h in heights])

    def quadrature_weights(heights,logarithmic):
        z=np.log(heights) if logarithmic else heights
        widths=np.r_[np.diff(z)[0]/2,(z[2:]-z[:-2])/2,np.diff(z)[-1]/2]
        return widths/widths.sum()

    def prediction(kind,z,p):
        if kind=='line':return p[0]*z
        if kind=='cubic':return p[0]*z+p[1]*z*z+p[2]*z**3
        if kind=='power':return p[0]*z**p[1]
        if kind=='bent_power':return p[0]*z*(1+p[1]*z)/(1+p[2]*z)**p[3]
        raise ValueError(kind)

    specs=[
        ('core_line','line',0,4,[1.3],([-5],[5]),False),
        ('near_wall_cubic','cubic',0,1,[1.2,.2,-.1],([-10]*3,[10]*3),False),
        ('intermediate_power','power',16,1024,[1.45,.94],([0,.1],[10,2]),True),
        ('wide_power','power',4,10000,[1.4,.94],([0,.1],[10,2]),True),
        ('whole_bent_power','bent_power',float(y[1]),10000,[1.15,3.2,2.4,1.06],
         ([.5,.001,.001,.8],[2,100,100,1.5]),True)]
    models=[]
    for name,kind,low,high,start,bounds,relative in specs:
        use=(y>=low)&(y<=high)&(extracted[:,5]==0)
        heights=y[use];observed=position[use]
        weights=quadrature_weights(heights,relative)
        denominator=observed if relative else np.ones_like(observed)
        def solve(mask,start):
            fit=least_squares(lambda p:np.sqrt(weights[mask])*(prediction(kind,heights[mask],p)-observed[mask])/denominator[mask],
                start,bounds=bounds,max_nfev=2000,ftol=1e-12,xtol=1e-12,gtol=1e-12)
            if not fit.success:raise RuntimeError(f'{name}: {fit.message}')
            return fit
        fit=solve(np.ones(len(heights),dtype=bool),start)
        # Held-out native rows check interpolation/fit stability only, not an
        # independent trajectory or a convergence estimate for the solver.
        train=np.arange(len(heights))%2==0
        held=solve(train,fit.x)
        error=prediction(kind,heights,fit.x)-observed
        test_error=prediction(kind,heights[~train],held.x)-observed[~train]
        dense=np.geomspace(low,high,501) if relative else np.linspace(low,high,501)
        dense_observed=curve(dense)-origin
        dense_error=prediction(kind,dense,fit.x)-dense_observed
        dense_weights=quadrature_weights(dense,relative)
        row=dict(name=name,kind=kind,coefficients=fit.x.tolist(),requestedRange=[low,high],
            nativeRange=[float(heights[0]),float(heights[-1])],nativeRows=int(len(heights)),
            objective='relative position error, uniform in log(height)' if relative else 'absolute position error, uniform in height',
            nativeRmsAbsolute=float(np.sqrt(np.dot(weights,error**2))),nativeMaxAbsolute=float(abs(error).max()),
            denseRmsAbsolute=float(np.sqrt(np.dot(dense_weights,dense_error**2))),denseMaxAbsolute=float(abs(dense_error).max()),
            heldoutCoefficients=held.x.tolist(),heldoutMaxAbsolute=float(abs(test_error).max()))
        if relative:
            row.update(nativeRmsRelative=float(np.sqrt(np.dot(weights,(error/observed)**2))),
                nativeMaxRelative=float(abs(error/observed).max()),
                denseRmsRelative=float(np.sqrt(np.dot(dense_weights,(dense_error/dense_observed)**2))),
                denseMaxRelative=float(abs(dense_error/dense_observed).max()),
                heldoutMaxRelative=float(abs(test_error/observed[~train]).max()))
        if kind=='bent_power':row['formalFarPower']=float(2-fit.x[3])
        models.append(row)

    samples=np.array([0,.1,.25,.5,1,2,4,8,16,32,64,128,256,1024,4096,10000])
    positions=curve(samples)
    sensitivity=[]
    for low,high in [(0,1),(0,4),(16,1024),(4,10000)]:
        use=(y>=low)&(y<=high)
        shift=position[use]-alternative[use]
        entry=dict(range=[low,high],maximumCubicVersusParabolicDisplacement=float(abs(shift).max()),
            maximumDifferenceInLocalGridCells=float(np.nanmax(abs(extracted[use,0]-extracted[use,2])/extracted[use,4])))
        if low>0:entry['maximumRelativeDisplacementDifference']=float(np.max(abs(shift/position[use])))
        sensitivity.append(entry)

    # Refit the intermediate power using the alternate peak extraction.
    use=(y>=16)&(y<=1024);yy=y[use];aa=alternative[use]
    ww=quadrature_weights(yy,True)
    alt=least_squares(lambda p:np.sqrt(ww)*(p[0]*yy**p[1]/aa-1),[1.45,.94],
                      bounds=([0,.1],[10,2]),ftol=1e-12,xtol=1e-12,gtol=1e-12)
    result=dict(source=data['metadata']['sourceResultFile'],tau=float(data['metadata']['tau']),
        coordinates='X=x1, Y=x2 from the saved main run; not divided by ell and not physical-space coordinates',
        definition='m(Y)=argmax_X rhoX(X,Y)-argmax_X rhoX(X,0)',wallPeakX=origin,
        originalWallPeakFromParabola=float(extracted[0,2]),ell=float(data['metadata']['ell']),
        fitFunctions={'line':'m=a Y','cubic':'m=a Y+b Y^2+c Y^3','power':'m=A Y^p',
                      'bent_power':'m=a Y (1+b Y)/(1+c Y)^d'},
        models=models,sampleColumns=['x2','x1Peak','peakDisplacementFromWall'],
        samples=np.column_stack([samples,positions,positions-origin]).tolist(),
        extractionSensitivity=sensitivity,alternateIntermediatePower=alt.x.tolist(),
        caveats=['Fits describe this finite-time saved field, not a proved limiting ridge.',
                 'The bent-power formula is empirical on the fitted interval; its infinite-height power is unverified.',
                 'Cubic/parabolic differences measure extraction sensitivity, not full numerical error.',
                 'Dense checks interpolate the same saved field; held-out rows are not independent simulations.'])
    json_save(out/'ridge_fit.json',result)
    savemat(out/'ridge_fit.mat',dict(x2=y,x1Peak=extracted[:,0],ridgeDisplacement=position,
        rhoXPeak=extracted[:,1],x1PeakParabolic=extracted[:,2],x1PeakGrid=extracted[:,3],
        localGridDx=extracted[:,4],boundaryMaximum=extracted[:,5],wallPeakX=origin,
        tau=float(data['metadata']['tau']),models=models,samples=np.column_stack([samples,positions,positions-origin])),
        long_field_names=True,do_compression=True)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(15,4.6),constrained_layout=True)
    core=np.linspace(0,4,301);core_positions=curve(core)-origin
    ax=axes[0];ax.plot(core,core_positions,'k',lw=2,label='extracted ridge')
    line=models[0];near=models[1];whole=models[4];power=models[2]
    ax.plot(core,prediction('line',core,line['coefficients']),'--',label='line: 0 <= x2 <= 4')
    ax.plot(core,prediction('bent_power',core,whole['coefficients']),':',lw=2,label='bent power: full fit')
    use=y<=4;ax.plot(y[use],position[use],'k.',ms=3)
    ax.set(xlabel=r'$x_2$',ylabel=r'$x_1^*(x_2)-x_1^*(0)$',title='Core peak position');ax.legend(fontsize=8);ax.grid(alpha=.2)
    ax=axes[1];use=(y>0)&(y<=10000)
    ax.semilogx(y[use],position[use]/y[use],'k.',ms=3,label='native rows')
    z=np.geomspace(y[1],10000,701)
    ax.semilogx(z,prediction('bent_power',z,whole['coefficients'])/z,'r',lw=1.4,label='bent-power fit')
    ax.set(xlabel=r'$x_2$',ylabel=r'$[x_1^*(x_2)-x_1^*(0)]/x_2$',title='Slope from the wall peak');ax.legend(fontsize=8);ax.grid(alpha=.2)
    ax=axes[2];z=np.geomspace(16,1024,501);observed=curve(z)-origin
    ax.semilogx(z,100*(prediction('power',z,power['coefficients'])/observed-1),label='power-fit relative error')
    use=(y>=16)&(y<=1024)
    ax.semilogx(y[use],100*(alternative[use]/position[use]-1),'.',ms=3,label='parabolic vs cubic peak')
    ax.axhline(0,color='k',lw=.6)
    ax.set(xlabel=r'$x_2$',ylabel='Difference (%)',title=rf"$m={power['coefficients'][0]:.4f}\,x_2^{{{power['coefficients'][1]:.5f}}}$")
    ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.savefig(out/'ridge_fit.png',dpi=180);plt.close(fig)
    print(json.dumps(result),flush=True)


def geometric_transform(source,out):
    """Straighten the ridge and normalize both half widths with simple constants.

    The map uses beta=3/2 exactly. Its height-dependent quantities are measured
    from the data; decimal regression coefficients are deliberately not built
    into the coordinate system. The rational master profile is a diagnostic on
    |xi|<=3/2, not an asserted global tail.
    """
    data=loadmat(source/'run_analysis.mat',simplify_cells=True)
    x=data['xFull'];all_y=data['yFull'];all_q=data['rhoX'];all_rho=data['rho']
    reliable=all_y<=10000
    y=all_y[reliable];q=all_q[reliable];rho=all_rho[reliable]
    qfield=RectBivariateSpline(y,x,q);rhofield=RectBivariateSpline(y,x,rho)
    rho_origin=float(rhofield(0,0)[0,0])
    beta=sp.Rational(3,2)
    beta_float=float(beta)

    def G(xi,k):
        return (np.sinh(beta_float*xi)+k*(np.cosh(beta_float*xi)-1))/beta_float

    def Gxi(xi,k):
        return np.cosh(beta_float*xi)+k*np.sinh(beta_float*xi)

    def measure(height):
        row=qfield(height,x)[0]
        spline=CubicSpline(x,row)
        roots=spline.derivative().roots(extrapolate=False)
        roots=roots[(roots>x[0])&(roots<x[-1])]
        peak=float(roots[np.argmax(spline(roots))]);amplitude=float(spline(peak))
        half=CubicSpline(x,row-amplitude/2).roots(extrapolate=False)
        left=half[half<peak];right=half[half>peak]
        if len(left)==0 or len(right)==0:
            raise ValueError(f'Half maximum not resolved at x2={height:g}')
        left_width=peak-float(left[-1]);right_width=float(right[0])-peak
        delta=(right_width-left_width)/(right_width+left_width)
        skew=delta/np.tanh(beta_float/2)
        width=(right_width+left_width)*beta_float/(2*np.sinh(beta_float))
        if abs(skew)>=1:
            raise ValueError(f'Geometric map is not onto at x2={height:g}: k={skew:g}')
        center=float(rhofield(height,peak)[0,0])-rho_origin
        return [peak,amplitude,width,skew,center,left_width,right_width]

    measured=np.array([measure(height) for height in y])
    wall_peak=float(measured[0,0])
    measured[:,0]-=wall_peak
    # Logarithmic sampling prevents the many outer rows from outweighing the
    # resolved near-wall core. Even samples define the master; odd samples audit it.
    heights=np.r_[0,np.geomspace(max(float(y[1]),.02),10000,140)]
    sample=np.array([measure(height) for height in heights])
    sample[:,0]-=wall_peak
    xi=np.linspace(-1.5,1.5,401);zero=np.argmin(abs(xi))
    gradient=[];density=[]
    for height,(peak,amplitude,width,skew,center,_,_) in zip(heights,sample):
        physical_x=peak+wall_peak+width*G(xi,skew)
        gradient.append(qfield.ev(np.full_like(xi,height),physical_x)/amplitude)
        density.append((rhofield.ev(np.full_like(xi,height),physical_x)
                        -(center+rho_origin))/(amplitude*width))
    gradient=np.asarray(gradient);density=np.asarray(density)
    train=np.arange(len(heights))%2==0
    master=gradient[train].mean(axis=0)
    density_prediction=[]
    for skew in sample[:,3]:
        derivative=master*Gxi(xi,skew)
        primitive=cumulative_trapezoid(derivative,xi,initial=0)
        density_prediction.append(primitive-primitive[zero])
    density_prediction=np.asarray(density_prediction)

    # Simple integer-coefficient candidate: positive for every real xi, but
    # only compared with data on the measured central interval.
    denominator=(xi**6-3*xi**5+3*xi**4+3*xi**3+28*xi**2+32)
    rational=32/denominator
    zeta=sp.symbols('zeta',real=True)
    exact_denominator=sp.Poly(zeta**6-3*zeta**5+3*zeta**4+3*zeta**3+28*zeta**2+32,zeta)
    assert exact_denominator.count_roots(-sp.oo,sp.oo)==0

    def relative_l2(value,reference):
        return float(np.linalg.norm(value-reference)/np.linalg.norm(reference))

    collapse=[]
    for high in (2,8,64,256,1024,10000):
        use=heights<=high
        collapse.append(dict(heightMax=high,
            gradientToLearnedMasterRelativeL2=relative_l2(gradient[use],np.broadcast_to(master,gradient[use].shape)),
            gradientToIntegerRationalRelativeL2=relative_l2(gradient[use],np.broadcast_to(rational,gradient[use].shape)),
            densityReconstructionRelativeL2=relative_l2(density[use],density_prediction[use]),
            maximumNormalizedGradientErrorToIntegerRational=float(np.max(abs(gradient[use]-rational)))))

    heldout=dict(
        gradientToLearnedMasterRelativeL2=relative_l2(gradient[~train],np.broadcast_to(master,gradient[~train].shape)),
        densityReconstructionRelativeL2=relative_l2(density[~train],density_prediction[~train]),
        maximumNormalizedGradientError=float(np.max(abs(gradient[~train]-master))),
        maximumNormalizedDensityError=float(np.max(abs(density[~train]-density_prediction[~train]))))

    # Test only small rational exponents. An intercept is fitted because the
    # normalization is dimensional; the exponent itself remains fixed.
    exponent_checks=[]
    for low,high in ((.25,8),(1,64),(4,256),(16,1024)):
        use=(y>=low)&(y<=high);log_width=np.log(measured[use,2]);log_amplitude=np.log(measured[use,1])
        alternatives=[]
        for numerator,denominator_alpha in ((1,3),(3,8),(2,5)):
            alpha=numerator/denominator_alpha
            residual=log_amplitude-(alpha-1)*log_width
            residual-=residual.mean()
            alternatives.append(dict(alpha=f'{numerator}/{denominator_alpha}',
                logarithmicRms=float(np.sqrt(np.mean(residual**2))),
                maximumRelativeDeviation=float(np.max(abs(np.exp(residual)-1)))))
        exponent_checks.append(dict(range=[low,high],alternatives=alternatives))

    # Medium-height ratios are reported as tests of simple fractions, not used
    # to overwrite the measured geometry.
    simple=[]
    for height in (32,64,128,256,512,1024,2048,4096):
        row=np.array(measure(height));row[0]-=wall_peak
        peak,amplitude,width,skew,center,_,_=row
        simple.append(dict(x2=height,ridgeOverX2=float(peak/height),
            twiceWidthOverX2=float(2*width/height),
            skewOverThreeQuarters=float(skew/(3/4)),
            peakOverThreeQuartersX2PowerMinusThreeFifths=float(amplitude/((3/4)*height**(-3/5))),
            fourCenterDensityOverPeakTimesX2=float(4*center/(amplitude*height))))

    time_convergence=[]
    history_path=out/'geometry_history.mat'
    if history_path.exists():
        history=loadmat(history_path,simplify_cells=True)
        for time,xh,yh,qh,rhoh in zip(history['tau'],history['x'],history['y'],history['rhoX'],history['rho']):
            qh_field=RectBivariateSpline(yh,xh,qh);rhoh_field=RectBivariateSpline(yh,xh,rhoh)
            rhoh_origin=float(rhoh_field(0,0)[0,0])

            def history_measure(height):
                row=qh_field(height,xh)[0];spline=CubicSpline(xh,row)
                roots=spline.derivative().roots(extrapolate=False)
                roots=roots[(roots>xh[0])&(roots<xh[-1])]
                peak=float(roots[np.argmax(spline(roots))]);amplitude=float(spline(peak))
                half=CubicSpline(xh,row-amplitude/2).roots(extrapolate=False)
                left_width=peak-float(half[half<peak][-1]);right_width=float(half[half>peak][0])-peak
                delta=(right_width-left_width)/(right_width+left_width)
                skew=delta/np.tanh(beta_float/2)
                width=(right_width+left_width)*beta_float/(2*np.sinh(beta_float))
                center=float(rhoh_field(height,peak)[0,0])-rhoh_origin
                return [peak,amplitude,width,skew,center]

            wall=history_measure(0);fit_heights=np.geomspace(4,256,60)
            rows=[];geometry=[]
            for height in fit_heights:
                peak,amplitude,width,skew,center=history_measure(height)
                physical_x=peak+width*G(xi,skew)
                rows.append(qh_field.ev(np.full_like(xi,height),physical_x)/amplitude)
                geometry.append([peak-wall[0],amplitude,width,skew,center])
            rows=np.asarray(rows);geometry=np.asarray(geometry)
            laws=[]
            for numerator,denominator_alpha in ((3,8),(2,5)):
                alpha=numerator/denominator_alpha
                residual=np.log(geometry[:,1])-(alpha-1)*np.log(geometry[:,2])
                residual-=residual.mean()
                laws.append(dict(alpha=f'{numerator}/{denominator_alpha}',
                    logarithmicRms=float(np.sqrt(np.mean(residual**2)))))
            ratios=[]
            for height in (32,64,128,256,512,1024):
                peak,amplitude,width,skew,center=history_measure(height)
                ratios.append(dict(x2=height,ridgeOverX2=float((peak-wall[0])/height),
                    twiceWidthOverX2=float(2*width/height),skew=float(skew),
                    fourCenterDensityOverPeakTimesX2=float(4*center/(amplitude*height))))
            time_convergence.append(dict(tau=float(time),heightRange=[4,256],
                integerRationalRelativeL2=relative_l2(rows,np.broadcast_to(rational,rows.shape)),
                exponentSensitivity=laws,ratios=ratios))

    result=dict(source=data['metadata']['sourceResultFile'],tau=float(data['metadata']['tau']),
        coordinates='saved rescaled X=x1, Y=x2 coordinates from main',
        map=dict(
            forward='x1=m(x2)+W(x2)*G_k(xi)',
            G='(2/3)*(sinh(3*xi/2)+k*(cosh(3*xi/2)-1))',
            inverse='xi=(2/3)*(asinh(sinh(theta)+(3/2)*z*cosh(theta))-theta), z=(x1-m)/W, theta=atanh(k)',
            definitions=['m is the horizontal rho_x1 maximum',
                         'a is rho_x1(m,x2)',
                         'left/right widths end at rho_x1=a/2',
                         'k=((right-left)/(right+left))/tanh(3/4)',
                         'W=3*(right+left)/(4*sinh(3/2))',
                         'eta=log(W/W(0))']),
        transformedRepresentation=[
            'rho_x1(x1,x2) approximately equals a(x2)*f(xi)',
            'rho-rho(m,x2) approximately equals a*W*integral_0^xi f(s)*(cosh(3s/2)+k*sinh(3s/2)) ds'],
        integerRationalMaster='f0(xi)=32/(xi^6-3xi^5+3xi^4+3xi^3+28xi^2+32)',
        integerRationalProperties=['f0(0)=1','f0(-1)=f0(1)=1/2','denominator is positive for all real xi'],
        collapseByHeight=collapse,heldoutAlternatingLogHeights=heldout,
        rationalExponentSensitivity=exponent_checks,simpleFractionDiagnostics=simple,
        timeConvergence=time_convergence,
        transformedPDE=dict(
            geometry='x=X(xi,eta), y=Y(eta), J=X_xi, K=X_eta, H=Y_eta',
            derivatives='partial_x=J^(-1) partial_xi; partial_y=H^(-1)*(partial_eta-(K/J) partial_xi)',
            poisson='-1/(JH)*[partial_xi(((K^2+H^2)/(JH))*psi_xi-(K/H)*psi_eta)+partial_eta(-(K/H)*psi_xi+(J/H)*psi_eta)]=Q_xi/J',
            transport='(psi_xi*Q_eta-psi_eta*Q_xi)/(JH)+lambda*((X-YK/H)*Q_xi/J+(Y/H)*Q_eta)+v*Q_xi/J=mu*Q+v'),
        cautions=['The integer rational is fitted only on |xi|<=3/2 and does not specify the physical far tail.',
                  'The simple fractions are candidate structure, not exact constants extracted from an unconverged field.',
                  'Rows beyond x2=1024 are used only as trend checks; they do not select the analytic constants.',
                  'Alternating heights test representation stability in the same saved field, not convergence of the PDE trajectory.'])
    json_save(out/'geometry_transform.json',result)
    savemat(out/'geometry_transform.mat',dict(x2=y,ridge=measured[:,0],peak=measured[:,1],
        width=measured[:,2],skew=measured[:,3],centerDensity=measured[:,4],
        leftHalfWidth=measured[:,5],rightHalfWidth=measured[:,6],sampleHeights=heights,
        xi=xi,transformedGradient=gradient,learnedMaster=master,integerRationalMaster=rational,
        transformedDensity=density,densityPrediction=density_prediction),do_compression=True)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    ax=axes[0,0]
    for height in (0,.5,2,8,32,128,512,1024):
        i=np.argmin(abs(heights-height))
        ax.plot(xi,gradient[i],lw=1,label=f'x2={height:g}')
    ax.plot(xi,rational,'k--',lw=2,label='integer rational')
    ax.set(xlabel='xi',ylabel='rho_x1 / a',title='Profiles after the inverse geometry map')
    ax.legend(fontsize=7,ncol=2);ax.grid(alpha=.2)

    ax=axes[0,1];positive=(y>=32)&(y<=4096)
    ax.semilogx(y[positive],measured[positive,0]/y[positive],label='m/x2; target 1')
    ax.semilogx(y[positive],2*measured[positive,2]/y[positive],label='2W/x2; target 1')
    ax.semilogx(y[positive],4*measured[positive,4]/(measured[positive,1]*y[positive]),label='4c/(a x2); target 1')
    ax.axhline(1,color='k',lw=.7);ax.axvspan(1024,4096,color='gray',alpha=.12,label='more tail-sensitive')
    ax.set(xlabel='x2',ylabel='dimensionless ratio',ylim=(.7,1.3),title='Simple fractions on the medium-height band')
    ax.legend(fontsize=8);ax.grid(alpha=.2)

    ax=axes[1,0]
    if time_convergence:
        times=[entry['tau'] for entry in time_convergence]
        ax.plot(times,[100*entry['integerRationalRelativeL2'] for entry in time_convergence],'o-',label='integer rational profile')
        for alpha in ('3/8','2/5'):
            ax.plot(times,[100*next(item['logarithmicRms'] for item in entry['exponentSensitivity'] if item['alpha']==alpha)
                           for entry in time_convergence],'o--',label=f'amplitude law alpha={alpha}')
    ax.set(xlabel='tau',ylabel='RMS mismatch (%)',title='Profile converges; exponent selection still drifts')
    ax.legend(fontsize=8);ax.grid(alpha=.2)

    ax=axes[1,1]
    if time_convergence:
        for entry in time_convergence:
            ax.semilogx([r['x2'] for r in entry['ratios']],[r['skew'] for r in entry['ratios']],marker='.',label=f"tau={entry['tau']:.0f}")
    ax.axhline(1/2,color='k',ls=':',lw=.8);ax.axhline(3/4,color='k',ls='--',lw=.8)
    ax.set(xlabel='x2',ylabel='k',title='Residual skew remains a modulation')
    ax.legend(fontsize=7,ncol=2);ax.grid(alpha=.2)
    fig.savefig(out/'geometry_transform.png',dpi=180);plt.close(fig)
    print(json.dumps(dict(heldout=heldout,collapse=collapse,simple=simple)),flush=True)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path);parser.add_argument('out',type=Path)
    parser.add_argument('--phase',choices=['base','refine','structure','ridge','geometry','all'],default='all')
    args=parser.parse_args();source=args.source;out=args.out
    out.mkdir(parents=True,exist_ok=True)
    names=['tau20','L4_257x129','L8_257x129_continued']
    if args.phase in ['base','all']:
        metrics=[]
        for name in names:
            metrics.append(construct(load_profile(source,name),out,name))
            json_save(out/'construction.json',metrics)
        plots(out,names)
    if args.phase in ['refine','all']:
        refine(source,out)
    if args.phase in ['structure','all']:
        global_structure(source,out)
        construction_identities(out)
    if args.phase in ['ridge','all']:
        ridge_fit(source,out)
    if args.phase in ['geometry','all']:
        geometric_transform(source,out)
