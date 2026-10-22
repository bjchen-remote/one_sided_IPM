"""Solve a finite reduced global model with centered BS built into its basis.

Unknowns are coefficients of a few angular shapes with radial functions.
Every velocity is recomputed by the half-plane radial Green transform.
No exterior field or boundary velocity is borrowed from a saved PDE root.
Finite angular/radial truncation and the unprojected transport residual
are reported separately. This is a numerical trial space, not exact closure.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import trapezoid
from scipy.interpolate import CubicSpline, RectBivariateSpline
from scipy.io import loadmat, savemat
from scipy.optimize import least_squares

from geometry_bs_reduce import ROOT, OUT, Seed


def radial_filters(d, log_r, alpha):
    """Exact integration for piecewise-linear d(log r); stable even for high n."""
    modes=np.arange(1,d.shape[1]+1)[:,None]
    step=log_r[1]-log_r[0]
    def causal(rate, values, initial):
        decay=np.exp(-rate*step)
        integral=-np.expm1(-rate*step)/rate
        wnew=(step-integral)/(rate*step);wold=integral-wnew
        out=np.empty_like(values);out[0]=initial
        for i in range(1,len(values)):
            out[i]=decay*out[i-1]+wold*values[i-1]+wnew*values[i]
        return out
    lower=causal(modes+alpha+1,d,np.zeros_like(d[0]))
    upper=np.empty_like(d)
    upper[:,1:]=causal(modes[1:]-alpha-1,d[::-1,1:],d[-1,1:]/(modes[1:]-alpha-1))[::-1]
    upper[:,0:1]=causal(np.array([[alpha]]),d[:,0:1],np.zeros_like(d[0,0:1]))
    b=(lower+upper)/(2*modes)
    bt=(-(modes+alpha+1)*lower+(modes-alpha-1)*upper)/(2*modes)
    b[:,0]=(lower[:,0]-upper[:,0])/2
    bt[:,0]=(-(alpha+2)*lower[:,0]+alpha*upper[:,0])/2
    return b,bt


class ReducedModel:
    def __init__(self,alpha,modes=32,radial_nodes=2401,scales=(.5,1.,2.,4.,8.,16.)):
        self.alpha=alpha;self.modes=modes;self.radial_nodes=radial_nodes
        self.seed=Seed(alpha,3.)
        self.descriptors=[('base',0,0,1.),('tail',0,0,1.),('offset',0,0,1.)]
        monomials=((0,1),(1,1),(0,2),(2,1),(1,2),(0,3),(4,0),(5,0))
        self.descriptors += [('local',i,j,L) for L in scales for i,j in monomials]
        self.logr=np.linspace(-16.,16.,radial_nodes)
        r=np.exp(self.logr)
        z,w=leggauss(max(256,4*modes));theta=(z+1)*np.pi/2;w=w*np.pi/2
        sine=np.sin(theta[:,None]*np.arange(1,modes+1))
        d=np.empty((len(r),modes,len(self.descriptors)))
        for start in range(0,len(r),64):
            rr=r[start:start+64]
            q,qr,qt=self.density(rr,theta)
            source=np.cos(theta)[None,:,None]*qr-np.sin(theta)[None,:,None]*qt/rr[:,None,None]
            d[start:start+64]=np.einsum('rtb,tn,t->rnb',source,sine,w,optimize=True)*(2/np.pi)*rr[:,None,None]**(1-alpha)
        del q,qr,qt,source
        b,bt=radial_filters(d,self.logr,alpha)
        self.bs_b=CubicSpline(self.logr,b,axis=0)
        self.bs_bt=CubicSpline(self.logr,bt,axis=0)

    def density(self,r,theta):
        r=np.asarray(r);theta=np.asarray(theta)
        cos,sin=np.cos(theta),np.sin(theta)
        aa,ap,_=self.seed.values(theta)
        arrays=[[],[],[]]
        for kind,i,j,L in self.descriptors:
            if kind=='base':
                e=np.exp(-r*r)
                R=r*e;Rp=(1-2*r*r)*e
                T=2/3*r**3*e;Tp=(2*r*r-4/3*r**4)*e
                value=R[:,None]*cos+T[:,None]*cos**3
                dr=Rp[:,None]*cos+Tp[:,None]*cos**3
                dt=-R[:,None]*sin-3*T[:,None]*cos**2*sin
            else:
                if kind in ('tail','offset'):
                    degree=self.alpha if kind=='tail' else 0.
                    R=r**degree*(r*r/(1+r*r))**2
                    Rp=R/r*(degree+4/(1+r*r))
                    phi,phip=(aa,ap) if kind=='tail' else (np.ones_like(theta),np.zeros_like(theta))
                else:
                    n=i+j;s=r/L
                    R=s**n*np.exp(-s*s);Rp=R/r*(n-2*s*s)
                    phi=cos**i*sin**j
                    phip=np.zeros_like(theta)
                    if i:phip-=i*cos**(i-1)*sin**(j+1)
                    if j:phip+=j*cos**(i+1)*sin**(j-1)
                value=R[:,None]*phi;dr=Rp[:,None]*phi;dt=R[:,None]*phip
            for group,v in zip(arrays,(value,dr,dt)):group.append(v)
        return tuple(np.stack(group,axis=-1) for group in arrays)

    def sample(self,r,theta):
        r=np.asarray(r);theta=np.asarray(theta)
        Q,Pr,Pt=self.density(r,theta)
        cs,ss=np.cos(theta),np.sin(theta)
        Qx=cs[None,:,None]*Pr-ss[None,:,None]*Pt/r[:,None,None]
        Qy=ss[None,:,None]*Pr+cs[None,:,None]*Pt/r[:,None,None]
        b=self.bs_b(np.log(r));bt=self.bs_bt(np.log(r))
        n=np.arange(1,self.modes+1)
        ur=-r[:,None,None]**self.alpha*np.einsum('rnb,tn,n->rtb',b,np.cos(theta[:,None]*n),n,optimize=True)
        ut=r[:,None,None]**self.alpha*np.einsum('rnb,tn->rtb',(self.alpha+1)*b+bt,np.sin(theta[:,None]*n),optimize=True)
        u1=cs[None,:,None]*ur-ss[None,:,None]*ut
        u2=ss[None,:,None]*ur+cs[None,:,None]*ut
        return dict(Q=Q.reshape(-1,Q.shape[-1]),Qx=Qx.reshape(-1,Q.shape[-1]),
                    Qy=Qy.reshape(-1,Q.shape[-1]),u1=u1.reshape(-1,Q.shape[-1]),
                    u2=u2.reshape(-1,Q.shape[-1]),
                    dilation=(r[:,None,None]*Pr-self.alpha*Q).reshape(-1,Q.shape[-1]),
                    x=(r[:,None]*cs).ravel(),y=(r[:,None]*ss).ravel(),
                    r=np.broadcast_to(r[:,None],(len(r),len(theta))).ravel(),
                    theta=np.broadcast_to(theta[None,:],(len(r),len(theta))).ravel())


def evaluate(c,lam,v,s):
    vals={k:s[k]@c for k in ('Q','Qx','Qy','u1','u2','dilation')}
    residual=vals['u1']*vals['Qx']+vals['u2']*vals['Qy']+lam*vals['dilation']+v*(vals['Qx']-1)
    return residual,vals


def audit(c,lam,v,s,alpha):
    res,vals=evaluate(c,lam,v,s)
    scale=(1+s['r'])**alpha
    bins=[]
    for lo,hi in ((0.,1.),(1.,4.),(4.,16.),(16.,64.),(64.,256.)):
        mask=(s['r']>=lo)&(s['r']<hi)
        if not mask.any():continue
        bins.append(dict(radius=[lo,hi],scaled_transport_rms=float(np.sqrt(np.mean((res[mask]/scale[mask])**2))),
                         scaled_transport_linf=float(np.max(abs(res[mask]/scale[mask]))),
                         min_scaled_q=float(np.min(vals['Qx'][mask]*(1+s['r'][mask])**(1-alpha)))))
    return dict(scaled_transport_rms=float(np.sqrt(np.mean((res/scale)**2))),
                scaled_transport_linf=float(np.max(abs(res/scale))),
                min_scaled_q=float(np.min(vals['Qx']*(1+s['r'])**(1-alpha))),
                negative_q_fraction=float(np.mean(vals['Qx']<0)),radial_bins=bins)


def initial_guess(model,s):
    data=loadmat(ROOT/'runs/geometry_loop_20260928/frame_tau20.mat',simplify_cells=True)
    field=RectBivariateSpline(data['y'],data['x'],data['Q'])
    target=field.ev(s['y'],s['x'])
    weight=(1+s['r'])**(-model.alpha)
    matrix=s['Q'][:,1:]*weight[:,None]
    coef=np.linalg.lstsq(matrix,(target-s['Q'][:,0])*weight,rcond=1e-7)[0]
    coef[0]=max(coef[0],.05)
    return np.r_[coef,.60,.175],float(np.sqrt(np.mean((matrix@coef-(target-s['Q'][:,0])*weight)**2)))


def run(alpha,modes,radial_nodes,max_nfev,out,positive_weight=100.,tag='positive',resume=None):
    model=ReducedModel(alpha,modes,radial_nodes)
    radii=np.geomspace(.04,32.,45)
    theta=np.linspace(0,np.pi,49)
    s=model.sample(radii,theta)
    fitted,fiterr=initial_guess(model,s)
    # Enforce the far constant balance exactly: Q_infinity contains -v/mu.
    start=np.r_[fitted[0],fitted[2:-2],fitted[-2:]]
    if resume:
        previous=json.loads(Path(resume).read_text())
        pc=np.asarray(previous['coefficients'])
        assert abs(previous['alpha']-alpha)<1e-14 and len(pc)==len(model.descriptors)
        start=np.r_[pc[1],pc[3:],previous['rates']['lambda_'],previous['rates']['v']]
    def coefficients(z):
        lam,v=z[-2:]
        return np.r_[1.,z[0],-v/(alpha*lam),z[1:-2]]
    active=np.r_[1,np.arange(3,len(model.descriptors))]
    weight=(1+s['r'])**(-alpha)*np.where(s['r']<=8,1.,.3)
    pr=np.geomspace(.025,128.,89);pt=np.linspace(0,np.pi,97)
    _,dr,dt=model.density(pr,pt)
    positive_q=(np.cos(pt)[None,:,None]*dr-np.sin(pt)[None,:,None]*dt/pr[:,None,None]).reshape(-1,len(model.descriptors))
    qweight=np.broadcast_to((1+pr[:,None])**(1-alpha),(len(pr),len(pt))).ravel()
    def fun(z, jac=False):
        c=coefficients(z);lam,v=z[-2:]
        res,a=evaluate(c,lam,v,s)
        scaled_q=qweight*(positive_q@c)-1e-4
        neg=scaled_q<0
        if not jac:
            return np.r_[weight*res,positive_weight*np.minimum(scaled_q,0)]
        mat=(s['u1']*a['Qx'][:,None]+s['u2']*a['Qy'][:,None]
             +s['Qx']*(a['u1']+v)[:,None]+s['Qy']*a['u2'][:,None]+lam*s['dilation'])
        dc_lam=v/(alpha*lam**2);dc_v=-1/(alpha*lam)
        main=np.column_stack([mat[:,active],a['dilation']+mat[:,2]*dc_lam,
                              a['Qx']-1+mat[:,2]*dc_v])*weight[:,None]
        positive=np.column_stack([positive_q[:,active],positive_q[:,2]*dc_lam,
                                  positive_q[:,2]*dc_v])*(positive_weight*qweight*neg)[:,None]
        return np.vstack([main,positive])
    # Independent directional finite difference checks the nonlinear BS response.
    rng=np.random.default_rng(927);direction=rng.normal(size=len(start));direction/=np.linalg.norm(direction)
    step=1e-5
    exact=fun(start,True)@direction
    numeric=(fun(start+step*direction)-fun(start-step*direction))/(2*step)
    jacerr=float(np.linalg.norm(exact-numeric)/np.linalg.norm(exact))
    lower=np.full_like(start,-np.inf);upper=np.full_like(start,np.inf)
    lower[0]=.01;upper[0]=10
    lower[-2]=.05;upper[-2]=2.;lower[-1]=-2.;upper[-1]=2.
    start=np.minimum(np.maximum(start,lower+1e-8),upper-1e-8)
    initial_audit=audit(coefficients(start),start[-2],start[-1],s,alpha)
    sol=least_squares(fun,start,jac=lambda z:fun(z,True),bounds=(lower,upper),
                      x_scale='jac',max_nfev=max_nfev,ftol=2e-9,xtol=2e-9,gtol=2e-8)
    c=coefficients(sol.x);lam,v=sol.x[-2:]
    held_r=np.geomspace(.045,200.,59);held_theta=np.r_[0.,(np.arange(64)+.5)*np.pi/64,np.pi]
    held=model.sample(held_r,held_theta)
    result=dict(alpha=alpha,unknown_count=len(sol.x),angular_bs_modes=modes,radial_bs_nodes=radial_nodes,
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                resumed_from=str(resume) if resume else None,
                positivity_penalty_weight=positive_weight,positivity_sample_count=len(positive_q),
                initial_data_fit_scaled_rms=fiterr,initial_pde_audit=initial_audit,
                optimizer=dict(success=bool(sol.success),message=str(sol.message),nfev=sol.nfev,
                               cost=float(sol.cost),optimality=float(sol.optimality),
                               active_bounds=np.flatnonzero(sol.active_mask).tolist()),
                rates=dict(lambda_=float(lam),mu=float(alpha*lam),v=float(v)),
                jacobian_directional_relative_error=jacerr,
                tail_density_offset_defect=float(alpha*lam*c[2]+v),
                training=audit(c,lam,v,s,alpha),heldout=audit(c,lam,v,held,alpha),
                coefficients=c.tolist(),basis=[list(d) for d in model.descriptors])
    refined=ReducedModel(alpha,2*modes,2*radial_nodes-1)
    fine=refined.sample(held_r,held_theta)
    result['refined_bs_heldout']=audit(c,lam,v,fine,alpha)
    _,old=evaluate(c,lam,v,held);_,new=evaluate(c,lam,v,fine)
    result['bs_refinement_relative_l2']={k:float(np.linalg.norm(old[k]-new[k])/np.linalg.norm(new[k])) for k in ('u1','u2')}
    result['status']='finite_reduced_least_squares_trial; global root not certified'
    file=out/f'global_{tag}_alpha_{alpha:g}.json'
    file.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    payload=dict(r=held_r,theta=held_theta,coefficients=c,rates=np.array([lam,alpha*lam,v]),
                 **{k:value.reshape(len(held_r),len(held_theta)) for k,value in new.items()})
    np.savez_compressed(out/f'global_{tag}_alpha_{alpha:g}.npz',**payload)
    savemat(out/f'global_{tag}_alpha_{alpha:g}.mat',payload,do_compression=True)
    print(json.dumps({k:result[k] for k in ('alpha','optimizer','rates','training','heldout','refined_bs_heldout','bs_refinement_relative_l2')}),flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--alpha',type=float,default=.4)
    parser.add_argument('--modes',type=int,default=24);parser.add_argument('--radial-nodes',type=int,default=1601)
    parser.add_argument('--max-nfev',type=int,default=150);parser.add_argument('--out',type=Path,default=OUT)
    parser.add_argument('--positive-weight',type=float,default=100.);parser.add_argument('--tag',default='positive')
    parser.add_argument('--resume',type=Path)
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    run(args.alpha,args.modes,args.radial_nodes,args.max_nfev,args.out,args.positive_weight,args.tag,args.resume)


if __name__=='__main__':main()
