"""Finite-parameter shooting for a fixed, PDE-generated mother; no corrections.

The first normal wall identity generates f. Three initial geometry rates are
selected by the SECOND normal identity, whose five accelerations are eliminated
linearly. A finite residual is not a symbolic closure certificate.
"""
from pathlib import Path
import json
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

OUT=Path(__file__).resolve().parents[2]/'runs/exact_master_20260928'


def build(parameters,config,z):
    a,b,kap=parameters
    lam,mu,v,beta,k,C=[config[k] for k in ['lam','mu','v','beta','k0','peakCurvature']]
    d=2*(lam-mu);nu=2*lam-mu
    r0=-(1+a*v)/d
    m=r0-nu*a/(C*v);c=2*r0-nu*a/(C*v)

    def geom(zz,state):
        f,h,H1=state
        sh,ch=np.sinh(beta*zz),np.cosh(beta*zz)
        G=(sh+k*(ch-1))/beta;g=ch+k*sh;Gk=(ch-1)/beta
        T=m+b*G+kap*Gk;B=c+(a+b)*h+kap*H1;N=mu*h+v
        return sh,ch,G,g,Gk,T,B,N

    def rhs(zz,state):
        f,h,H1=state
        sh,ch,G,g,Gk,T,B,N=geom(zz,state)
        fp=-g*f*f*(a*N+d*(B-T*f)+f*f)/(N*(B-2*T*f))
        return np.array([fp,f*g,f*sh])

    sols=[solve_ivp(rhs,[0,end],[1.,0.,0.],rtol=2e-10,atol=2e-12,
                    method='DOP853',dense_output=True,max_step=.035)
          for end in [float(z.min()),float(z.max())]]
    if not all(sol.success for sol in sols):
        raise ValueError('The mother ODE encountered a degeneracy before covering the requested interval')
    states=np.column_stack([sols[0 if zz<0 else 1].sol(zz) for zz in z])
    f,h,H1=states
    sh,ch,G,g,Gk,T,B,N=geom(z,states)
    deriv=rhs(z,states);fp=deriv[0]
    # Analytic directional derivative of the rational/trigonometric ODE.
    # Complex step avoids finite difference cancellation and adds no shape modes.
    eps=1e-28
    fpp=np.imag(rhs(z+1j*eps,states+1j*eps*deriv)[0])/eps
    gp=beta*(sh+k*ch);Tp=b*g+kap*sh
    qx=fp/g;qxx=(fpp*g-fp*gp)/(g**3)
    U=N/f;Ux=mu-N*qx/(f*f)
    Uxx=-mu*qx/f-N*qxx/(f*f)+2*N*qx*qx/(f**3)
    r=B-T*f;rx=a*f-T*qx
    basis=np.column_stack([-f,h-f*G,H1-f*Gk,h,np.ones_like(z)])
    basis_xi=np.column_stack([-fp,-fp*G,-fp*Gk,f*g,np.zeros_like(z)])
    bend=2*a*b*h+2*(a+b)*kap*H1-2*a*f*T+fp*T*T/g-2*b*kap*f*Gk
    bend_xi=(2*a*b*f*g+2*(a+b)*kap*f*sh-2*a*(fp*T+f*Tp)
              +fpp*T*T/g+2*fp*T*Tp/g-fp*T*T*gp/(g*g)
              -2*b*kap*(fp*Gk+f*sh))
    coefficient=4*lam-mu-2*Ux
    matrix=U[:,None]/g[:,None]*basis_xi+coefficient[:,None]*basis
    forcing=Uxx*f-3*f*rx+qx*r
    target=forcing-U/g*bend_xi-coefficient*bend
    first=U*rx+(2*lam-mu-Ux)*r+f*f
    return dict(f=f,h=h,H1=H1,matrix=matrix,target=target,forcing=forcing,
                first=first,minDen=float(np.min(abs(N*(B-2*T*f)))),
                mprime=float(m),cprime=float(c))


def eliminate(state,mask=None):
    mat=state['matrix'];rhs=state['target']
    if mask is None:mask=np.ones(len(rhs),bool)
    scales=np.linalg.norm(mat[mask],axis=0)
    normalized=mat[mask]/scales
    value,_,rank,sv=np.linalg.lstsq(normalized,rhs[mask],rcond=None)
    accel=value/scales
    return accel,mat@accel-rhs,dict(rank=int(rank),condition=float(sv[0]/sv[-1]))


def shoot(prior,train_range,hold_range,free_rates=False):
    z=np.linspace(*train_range,161)
    initial=np.array([prior[k] for k in ['a','b','kprime']])
    if free_rates:
        initial=np.r_[initial,[prior[k] for k in ['lam','mu','v']]]
    def evaluate(par,zz):
        config=dict(prior)
        if free_rates:
            config.update(zip(['lam','mu','v'],par[3:]))
        return build(par[:3],config,zz)
    reference=evaluate(initial,z)
    _,res0,_=eliminate(reference)
    normalization=float(np.linalg.norm(reference['forcing'])/np.sqrt(len(z)))
    calls=0
    def objective(par):
        nonlocal calls
        calls+=1
        try:
            state=evaluate(par,z)
            _,res,_=eliminate(state)
            return res/normalization
        except (ValueError,FloatingPointError):
            return np.ones(len(z))*1e3
    lower=[-1.4,.3,-.5];upper=[-.35,2.8,1.0]
    if free_rates:
        lower += [.38,.13,.10];upper += [.68,.25,.26]
    result=least_squares(objective,initial,bounds=(lower,upper),
                         max_nfev=160,ftol=3e-12,xtol=3e-12,gtol=3e-12,
                         diff_step=2e-5)
    solved=evaluate(result.x,z)
    accel,res,cond=eliminate(solved)
    # Freeze all fitted constants and score a denser, wider interval.
    holdz=np.linspace(*hold_range,601)
    hold_error=None
    try:
        hold=evaluate(result.x,holdz)
    except ValueError as exc:
        hold_error=str(exc)
        holdz=np.linspace(*train_range,601)
        hold=evaluate(result.x,holdz)
    holdres=hold['matrix']@accel-hold['target']
    outside=(holdz<z.min())|(holdz>z.max())
    changed=evaluate(initial,holdz)
    data=dict(config=prior,freeRates=free_rates,initial=initial.tolist(),optimized=result.x.tolist(),
        firstGeometryRates=[solved['mprime'],float(result.x[1]),float(result.x[2]),
                            float(result.x[0]),solved['cprime']],
        secondGeometryRates=accel.tolist(),trainingXiRange=[float(z.min()),float(z.max())],
        requestedHoldoutXiRange=hold_range,
        evaluatedXiRange=[float(holdz.min()),float(holdz.max())],
        holdoutFailure=hold_error,
        initialSecondRms=float(np.sqrt(np.mean(res0**2))),
        optimizedSecondRms=float(np.sqrt(np.mean(res**2))),
        trainingSecondRelativeToForcing=float(np.linalg.norm(res)/np.linalg.norm(solved['forcing'])),
        holdoutSecondRelativeToForcing=(float(np.linalg.norm(holdres[outside])/np.linalg.norm(hold['forcing'][outside]))
            if outside.any() else None),
        holdoutSecondMax=float(abs(holdres[outside]).max()) if outside.any() else None,
        firstIdentityMax=float(abs(hold['first']).max()),minimumMother=float(hold['f'].min()),
        minimumAbsOdeDenominator=hold['minDen'],
        optimizer=dict(success=bool(result.success),message=result.message,nfev=int(result.nfev),calls=calls),
        condition=cond,
        restriction='No correction modes. First wall equation exact by ODE. Nonzero second residual remains; no full-height, all-order or global closure claim.')
    return data,dict(xi=holdz,mother=hold['f'],initialMother=changed['f'],
                    firstResidual=hold['first'],secondResidual=holdres,forcing=hold['forcing'])


def main():
    prior=json.loads((OUT/'checks.json').read_text())['localConstruction']['parameters']
    rows=[]
    for label,train,hold,free in [
        ('narrow_fixed_rates',[-.65,.85],[-.9,1.15],False),
        ('wide_fixed_rates',[-.85,1.1],[-.95,1.25],False),
        ('wide_free_rates',[-.85,1.1],[-.95,1.25],True)]:
        data,arrays=shoot(prior,train,hold,free)
        data['case']=label
        rows.append(data)
        np.savez_compressed(OUT/f'{label}.npz',**arrays)
        print(json.dumps(data,indent=2),flush=True)
    (OUT/'shooting.json').write_text(json.dumps(rows,indent=2)+'\n')


if __name__=='__main__':main()
