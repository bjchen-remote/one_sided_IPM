"""Fresh, positive low-parameter model competition on original IPM data.

No old fitted constants are imported. Fits and held-out residuals are numerical
diagnostics, not certificates of exact closure or new half-plane solutions.
Run with numpy/scipy/sympy available; all outputs are in a new result directory.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.io import loadmat
from scipy.optimize import least_squares
from numpy.polynomial.legendre import leggauss
from fixed_point import fd_matrix

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'runs/fresh_positive_search_20260928'
BASE = ROOT / 'runs/steady_20260928_tau20'
GN, GW = leggauss(64)
NAMES = ('student', 'log_skew', 'algebraic_warp', 'atan_skew')


def sources():
    path = BASE/'run_analysis.mat'
    d = loadmat(path, simplify_cells=True)
    meta = d['metadata']; a, P, ell = (float(meta[k]) for k in ('a', 'P', 'ell'))
    origin = float(RectBivariateSpline(d['yFull'], d['xFull'], d['rho'])(0, a)[0, 0])
    yield dict(name='tau20', kind='native_finite_time', path=path,
               x=(d['xFull']-a)/ell, y=d['yFull']/ell,
               Q=(d['rho']-origin)/(P*ell), q=d['rhoX']/P,
               u1=(d['u1']-float(meta['uAtPeak']))/(P*ell),
               u2=d['u2']/(P*ell), rates=np.asarray(meta['innerRates']))
    for name in ('L4_129x65', 'L4_257x129', 'L8_257x129_continued'):
        path = BASE/name/'state.npz'
        with np.load(path) as z:
            d = {k:z[k].copy() for k in z.files}
        u1 = -fd_matrix(d['y'], 1, 5)@d['psi']
        u1 -= float(RectBivariateSpline(d['y'],d['x'],u1)(0,0)[0,0])
        yield dict(name=name, kind='conditional_frozen_exterior_root', path=path,
                   x=d['x'], y=d['y'], Q=d['Q'], q=d['qx'], u1=u1,
                   u2=d['psi']@fd_matrix(d['x'],1,5).T, rates=d['rates'])


def f(z, p, k, family):
    z = np.asarray(z)
    if family == 'algebraic_warp':
        w = z+k
        xi = (w-k*np.sqrt(1-k*k+w*w))/(1-k*k)
        return (1+xi*xi)**(-p)
    if family == 'atan_skew':
        t = k+np.sqrt(1+k*k)*z
        return np.exp(p*(np.log1p(k*k)-np.log1p(t*t))
                      +2*p*k*(np.arctan(t)-np.arctan(k)))
    logf = -p*np.log1p(z*z)
    if family == 'log_skew':
        logf += 2*p*k*(np.arcsinh(z)-z/np.sqrt(1+z*z))
    return np.exp(logf)


def primitive(z, p, k, family):
    z = np.asarray(z)
    return z/2*np.sum(f(z[...,None]*(GN+1)/2,p,k,family)*GW,axis=-1)


def wall_model(x, pars, family):
    p = pars[0]; k = pars[1] if len(pars)>1 else 0.
    W = np.sqrt(p)
    return f(x/W,p,k,family), W*primitive(x/W,p,k,family)


def score(qh,Qh,q,Q,mask):
    eq = float(np.linalg.norm((qh-q)[mask])/np.linalg.norm(q[mask]))
    eQ = float(np.linalg.norm((Qh-Q)[mask])/np.linalg.norm(Q[mask]))
    return dict(q_relative_l2=eq,Q_relative_l2=eQ,joint=np.hypot(eq,eQ).item(),
                q_linf=float(np.max(abs(qh[mask]-q[mask]))))


def fit_wall(source):
    x = np.linspace(-2.5,2.5,401)
    q = RectBivariateSpline(source['y'],source['x'],source['q'])(0,x)[0]
    Q = RectBivariateSpline(source['y'],source['x'],source['Q'])(0,x)[0]
    core = abs(x)<=1.25
    train = core & (np.arange(x.size)%2==0)
    hold = core & ~train
    outer = ~core
    rows = []
    for family in NAMES:
        ks = [0.] if family=='student' else [-.4,0.,.4]
        n = 1 if family=='student' else 2
        kbound = 3. if family=='atan_skew' else .97
        low = [.04]+([-kbound] if n==2 else [])
        high = [5.]+([kbound] if n==2 else [])
        def residual(pars):
            qh,Qh = wall_model(x[train],pars,family)
            return np.r_[(qh-q[train])/np.sqrt(np.mean(q[train]**2)),
                         (Qh-Q[train])/np.sqrt(np.mean(Q[train]**2))]
        fits = [least_squares(residual,[p]+([k] if n==2 else []),
                             bounds=(low,high),max_nfev=200,
                             ftol=1e-11,xtol=1e-11,gtol=1e-11)
                for p in [.15,.4,1.] for k in ks]
        best = min(fits,key=lambda z:np.linalg.norm(z.fun))
        qh,Qh = wall_model(x,best.x,family)
        rows.append(dict(family=family,shape_parameter_count=n,
                         parameters=best.x.tolist(),train=score(qh,Qh,q,Q,train),
                         interleaved_holdout=score(qh,Qh,q,Q,hold),
                         spatial_holdout=score(qh,Qh,q,Q,outer),
                         active_bounds=best.active_mask.tolist()))
    return rows


def geometry_values(x, pars, p, k, family):
    m,logA,logW,c = pars
    A,W = np.exp(logA),np.exp(logW); z=(x-m)/W
    ff=f(z,p,k,family); h=primitive(z,p,k,family)
    return A*ff,c+A*W*h,np.column_stack([-A*ff,W*h,A*(h-z*ff),np.ones_like(x)])


def project(M,R,train):
    norms = np.linalg.norm(M[train],axis=0)
    z,_,rank,sing = np.linalg.lstsq(M[train]/norms,R[train],rcond=None)
    coeff=z/norms; res=R-M@coeff; hold=~train
    return dict(coefficients=coeff.tolist(),rank=int(rank),
                condition=float(sing[0]/sing[-1]),
                train_relative_l2=float(np.linalg.norm(res[train])/np.linalg.norm(R[train])),
                holdout_relative_l2=float(np.linalg.norm(res[hold])/np.linalg.norm(R[hold])),
                holdout_linf=float(np.max(abs(res[hold])))),res


def transverse_audit(source, selected):
    family=selected['family']; p=selected['parameters'][0]
    k=selected['parameters'][1] if len(selected['parameters'])>1 else 0.
    interp={key:RectBivariateSpline(source['y'],source['x'],source[key])
            for key in ('q','Q','u1','u2')}
    lam,mu,v=map(float,source['rates'][:3])
    x=np.linspace(-1.5,1.5,241); train=np.arange(x.size)%2==0
    heights=[.025,.05,.1,.15,.2,.3,.4,.5,.6,.7,.8]
    rows=[]; all_res=[]
    for y in heights:
        q,Q,u1,u2=(interp[key](y,x)[0] for key in ('q','Q','u1','u2'))
        Qy=interp['Q'](y,x,dx=1)[0]
        def fun(pars):
            qh,Qh,_=geometry_values(x[train],pars,p,k,family)
            return np.r_[(qh-q[train])/np.sqrt(np.mean(q[train]**2)),
                         (Qh-Q[train])/max(np.sqrt(np.mean(Q[train]**2)),.1)]
        m0=x[np.argmax(q)]
        fit=least_squares(fun,[m0,np.log(max(q)),np.log(np.sqrt(p)*(1+y)),np.interp(m0,x,Q)],
                          bounds=([-5,-5,-5,-10],[5,3,4,10]),max_nfev=250)
        qh,Qh,B=geometry_values(x,fit.x,p,k,family)
        V1=u1+lam*x+v; V2=u2+lam*y
        R=mu*Qh+v-V1*qh; M=V2[:,None]*B
        proj,rr=project(M,R,train)
        Qyh=B@np.asarray(proj['coefficients'])
        scale=np.sqrt(np.sum((V1*qh)**2+(V2*Qyh)**2+(mu*Qh)**2+v*v))
        baseline=V1*q+V2*Qy-mu*Q-v
        basescale=np.sqrt(np.sum((V1*q)**2+(V2*Qy)**2+(mu*Q)**2+v*v))
        rows.append(dict(y=y,parameters_m_logA_logW_c=fit.x.tolist(),
                         shape_holdout=score(qh,Qh,q,Q,~train),projection=proj,
                         residual_term_relative_l2=float(np.linalg.norm(rr)/scale),
                         measured_baseline_term_relative_l2=float(np.linalg.norm(baseline)/basescale)))
        all_res.append(rr)
    # Projection chooses a derivative independently at each height. It must
    # also agree with the derivative of the fitted geometry to be a solution.
    from scipy.interpolate import CubicSpline
    parameters=np.array([row['parameters_m_logA_logW_c'] for row in rows])
    dparameters=CubicSpline(heights,parameters).derivative()(heights)
    for row,pars,dpars in zip(rows,parameters,dparameters):
        y=row['y']; qh,Qh,B=geometry_values(x,pars,p,k,family)
        deriv=dpars.copy(); deriv[1]*=np.exp(pars[1]); deriv[2]*=np.exp(pars[2])
        V1=interp['u1'](y,x)[0]+lam*x+v; V2=interp['u2'](y,x)[0]+lam*y
        qyh=B@deriv; residual=V1*qh+V2*qyh-mu*Qh-v
        scale=np.sqrt(np.sum((V1*qh)**2+(V2*qyh)**2+(mu*Qh)**2+v*v))
        row['fitted_geometry_derivative']=deriv.tolist()
        row['actual_geometry_transport_term_relative_l2']=float(np.linalg.norm(residual)/scale)
        row['projected_derivative_mismatch_relative_l2']=float(np.linalg.norm(
            deriv-np.asarray(row['projection']['coefficients']))/np.linalg.norm(deriv))
    # Necessary first-normal wall condition, independent of a velocity fit.
    q,Q=wall_model(x,selected['parameters'],family)
    _,_,B=geometry_values(x,[0.,0.,.5*np.log(p),0.],p,k,family)
    V1=(mu*Q+v)/q
    # Analytic primitive derivatives and complex-step evaluation avoid a
    # finite-difference/spline floor in this necessary wall diagnostic.
    z=x/np.sqrt(p); eps=1e-25
    fz=np.imag(f(z+1j*eps,p,k,family))/eps
    V1x=mu-(mu*Q+v)*fz/(np.sqrt(p)*q*q)
    Bx=np.column_stack([-fz/np.sqrt(p),q,-z*fz/np.sqrt(p),np.zeros_like(x)])
    C=V1[:,None]*Bx+(2*lam-mu-V1x)[:,None]*B
    wallproj,wallres=project(C,-q*q,train)
    extra=[]; extra_x=[]
    for label,par in [('vary_k','k'),('vary_k_and_p','p')]:
        eps=1e-25
        hder=np.imag(primitive(z,p+1j*eps if par=='p' else p,
                              k+1j*eps if par=='k' else k,family))/eps
        fder=np.imag(f(z,p+1j*eps if par=='p' else p,
                      k+1j*eps if par=='k' else k,family))/eps
        extra.append(np.sqrt(p)*hder)
        extra_x.append(fder)
        BB=np.column_stack([B]+extra)
        CC=V1[:,None]*np.column_stack([Bx]+extra_x)+(2*lam-mu-V1x)[:,None]*BB
        pp,_=project(CC,-q*q,train)
        wallproj[label]=pp
    return dict(heights=rows,wall_first_normal=wallproj,
                fixed_shape_parameters=[p,k],geometry_function_count=4,
                velocity_scope='Saved measured velocity held fixed: necessary diagnostic only; candidate Biot-Savart velocity is not computed.'), dict(x=x,heights=np.asarray(heights),residuals=np.asarray(all_res),wall_residual=wallres)


def resonant_exponent_probe(source):
    """Necessary tail-power resonances from a separate symbolic calculation.

    This is a new constrained experiment, never a redefinition of the original
    unconstrained data competition. The tighter even-power gate allows only
    1/6 and 1/4.  The older coarse-gate 1/3 is retained as a labelled control.
    """
    x=np.linspace(-2.5,2.5,401)
    q=RectBivariateSpline(source['y'],source['x'],source['q'])(0,x)[0]
    Q=RectBivariateSpline(source['y'],source['x'],source['Q'])(0,x)[0]
    train=(abs(x)<=1.25)&(np.arange(x.size)%2==0)
    hold=(abs(x)<=1.25)&~train; outer=abs(x)>1.25
    rows=[]
    for p in [1/6,1/4,1/3]:
        def residual(kk):
            qh,Qh=wall_model(x[train],[p,kk[0]],'algebraic_warp')
            return np.r_[(qh-q[train])/np.sqrt(np.mean(q[train]**2)),
                         (Qh-Q[train])/np.sqrt(np.mean(Q[train]**2))]
        fit=least_squares(residual,[.3],bounds=([-.97],[.97]),max_nfev=100)
        qh,Qh=wall_model(x,[p,fit.x[0]],'algebraic_warp')
        selected=dict(family='algebraic_warp',parameters=[p,float(fit.x[0])])
        audit,_=transverse_audit(source,selected)
        rows.append(dict(p=p,k=float(fit.x[0]),
                         exponent_status=('necessary_tail_resonance_only' if p<.3
                                          else 'coarse_gate_control_excluded_by_even_power_gate'),
                         train=score(qh,Qh,q,Q,train),
                         interleaved_holdout=score(qh,Qh,q,Q,hold),
                         spatial_holdout=score(qh,Qh,q,Q,outer),
                         wall_first_normal=audit['wall_first_normal']))
    return rows


def temporal_probe(selected):
    from scipy.interpolate import CubicSpline
    d=loadmat(BASE/'run_analysis.mat',simplify_cells=True)
    x=np.linspace(-2.5,2.5,401); rows=[]
    for target in [10.,15.,18.,20.]:
        i=int(np.argmin(abs(d['tau']-target)))
        a,P,_,ell=map(float,d['geometry'][i])
        cs=CubicSpline(d['xFull'],d['wallQ'][i]); anti=cs.antiderivative()
        q=cs(a+ell*x)/P
        Q=(anti(a+ell*x)-anti(a))/(P*ell)
        train=(abs(x)<=1.25)&(np.arange(x.size)%2==0)
        hold=(abs(x)<=1.25)&~train; outer=abs(x)>1.25
        def residual(pars):
            qh,Qh=wall_model(x[train],pars,'algebraic_warp')
            return np.r_[(qh-q[train])/np.sqrt(np.mean(q[train]**2)),
                         (Qh-Q[train])/np.sqrt(np.mean(Q[train]**2))]
        fit=least_squares(residual,selected['parameters'],bounds=([.04,-.97],[5,.97]))
        qh,Qh=wall_model(x,fit.x,'algebraic_warp')
        qlock,Qlock=wall_model(x,selected['parameters'],'algebraic_warp')
        rows.append(dict(tau=float(d['tau'][i]),refit_parameters=fit.x.tolist(),
                         refit_core=score(qh,Qh,q,Q,hold),refit_outer=score(qh,Qh,q,Q,outer),
                         locked_tau20_core=score(qlock,Qlock,q,Q,hold),
                         locked_tau20_outer=score(qlock,Qlock,q,Q,outer)))
    return dict(density_construction='Q is the integral of each saved wall-gradient spline, with Q(0)=0. It is not an independently stored earlier density.',records=rows)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    data=list(sources()); result=dict(protocol={
        'selection':'Fit native tau20 wall on |x|<=1.25, alternating samples. Compare one/two shape parameters without rational snapping. Spatial validation on 1.25<|x|<=2.5.',
        'normalization':'q(0)=1, q_x(0)=0, q_xx(0)=-2, Q(0)=0 imposed analytically.',
        'uncertainty':'Interpolated held-out points are not independent experiments; L4 roots share exterior provenance. No validated infinite-time or half-plane limit is available.',
        'exact_closure':'No numerical error threshold is accepted as a proof of R_perp identically zero.',
        'positivity':'Every competing f is globally positive; log_skew and algebraic_warp use |k|<1. Scale and amplitude are exponentials.',
        'coordinates':'Physical normalized x, not a previously chosen hyperbolic coordinate.'},sources={})
    for source in data:
        rows=fit_wall(source)
        result['sources'][source['name']]=dict(kind=source['kind'],
            input=str(source['path'].relative_to(ROOT)),sha256=hashlib.sha256(source['path'].read_bytes()).hexdigest(),
            rates=source['rates'].tolist(),wall_models=rows)
        print(source['name'],[(r['family'],r['parameters'],r['spatial_holdout']['joint']) for r in rows],flush=True)
    native=result['sources']['tau20']['wall_models']
    # Native core holdout decides; outer region remains a genuine validation.
    selected=min(native,key=lambda r:r['interleaved_holdout']['joint'])
    result['selected_by_native_core_holdout']=selected
    for source in data:
        audit,arrays=transverse_audit(source,selected)
        result['sources'][source['name']]['native_selected_shape_audit']=audit
        np.savez_compressed(OUT/(source['name']+'_diagnostics.npz'),**arrays)
    result['native_symbolic_resonance_followup']=resonant_exponent_probe(data[0])
    result['temporal_validation']=temporal_probe(selected)
    result['status']='Positive low-parameter candidates tested; no exact full closure certificate.'
    (OUT/'results.json').write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    print('selected',selected['family'],selected['parameters'],flush=True)


if __name__=='__main__':
    main()
