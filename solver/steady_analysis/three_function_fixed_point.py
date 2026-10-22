"""Design checks for a residual-exact three-function fixed-point iteration.

This does NOT solve the IPM profile. It verifies the algebraic reduction,
linearized iteration stability, and a consistent/inconsistent manufactured
three-block problem. Saved fields are only read for a width-slaving audit.
"""
from pathlib import Path
import json
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
import sympy as sp
from geometry_measure import candidate, trajectory_frame, interp_fields, peak_measure

OUT=Path(__file__).resolve().parents[2]/'runs/three_function_fixed_point_20260928'


def step(z,dual,residual,jacobian,metric,tau=1.,gamma=1.):
    """One preconditioned primal/dual step; dual stores defects, not q modes.

    Residual and Jacobian must already use identical, fixed positive weights.
    The optional finite rate parameters belong in z, not in the dual state.
    """
    R=residual(z);J=jacobian(z)
    delta=np.linalg.solve(metric/tau+gamma*J.T@J,-J.T@(dual+gamma*R))
    new=z+delta
    defect=residual(new)
    return new,dual+gamma*defect,dict(stepNorm=float(np.linalg.norm(delta)),
        residualNorm=float(np.linalg.norm(defect)),
        dualIncrementNorm=float(gamma*np.linalg.norm(defect)))


def symbolic():
    # Exact modal stability in the whitened primal metric, sigma=gamma.
    tau,gamma,s=sp.symbols('tau gamma s',positive=True)
    c=tau*gamma*s*s
    T=sp.Matrix([[1/(1+c),-tau*s/(1+c)],
                 [gamma*s/(1+c),1-c/(1+c)]])
    det=sp.factor(T.det());tr=sp.factor(sp.trace(T))
    assert sp.factor(det-1/(1+c))==0
    assert sp.factor(tr-2/(1+c))==0
    assert sp.factor(tr*tr-4*det+4*c/(1+c)**2)==0

    # Linearization of q=A f((x-m)/W), W=A^(-p), with a=log A.
    eps,x,m,a,p,dm,da=sp.symbols('eps x m a p dm da')
    f=sp.Function('f');df=sp.Function('df')
    xi=(x-m)*sp.exp(p*a)
    xie=(x-m-eps*dm)*sp.exp(p*(a+eps*da))
    qe=sp.exp(a+eps*da)*(f(xie)+eps*df(xie))
    # Keep derivatives in a separate dummy variable for SymPy substitution.
    zz=sp.symbols('zz')
    fp=sp.diff(f(zz),zz).subs(zz,xi)
    expected=sp.exp(a)*(df(xi)-sp.exp(p*a)*fp*dm+(f(xi)+p*xi*fp)*da)
    assert sp.simplify(sp.diff(qe,eps).subs(eps,0)-expected)==0

    # C elimination, including the wall-desingularized covariance identity.
    y,vr,rr,meanR,meanB,mu=sp.symbols('y var_b cov_b_r meanR meanB mu',nonzero=True)
    zold=-(y*y*rr)/(y*y*vr);znew=-rr/vr
    assert sp.simplify(zold-znew)==0
    C=(meanR+y*meanB*znew)/mu
    assert sp.simplify(meanR+y*meanB*znew-mu*C)==0

    # Peak transport gives log-amplitude and moving-peak equations.
    A,C,Cp,V1,V2,V1x,V2x,v,lam,ap,mp=sp.symbols(
        'A C Cp V1 V2 V1x V2x v lambda ap mp',nonzero=True)
    k=(V1-(mu*C+v)/A)/V2
    rhsA=(mu-V1x+V2x*k)/V2
    rhsM=Cp/A+k
    transport=V1*A+V2*(Cp-A*rhsM)-mu*C-v
    differentiated=V2*A*rhsA+(V1x-mu)*A+V2x*(Cp-A*rhsM)
    assert sp.factor(transport)==0 and sp.factor(differentiated)==0
    # All three finite rates can also be selected by normalized wall jets.
    xx,ux,uxx,uxxx,qxxx=sp.symbols('xx ux uxx uxxx qxxx')
    qq=1-xx**2+qxxx*xx**3/6
    QQ=xx-xx**3/3+qxxx*xx**4/24
    uu=ux*xx+uxx*xx**2/2+uxxx*xx**3/6
    defect=(uu+lam*xx+v)*qq-mu*QQ-v
    selected_v=uxx/2
    selected_mu=(uxxx+selected_v*qxxx)/4
    selected_lam=selected_mu-ux
    rate_errors=[sp.factor(sp.diff(defect,xx,j).subs(xx,0).subs(
        {lam:selected_lam,mu:selected_mu,v:selected_v},simultaneous=True)) for j in [1,2,3]]
    assert rate_errors==[0,0,0]
    return dict(shapeLinearization='exact zero',peakTransport='exact zero',
        peakGradientTransport='exact zero',wallCovarianceCancellation='exact zero',
        normalizedRateSelection=[str(e) for e in rate_errors],
        iterationTrace=str(tr),iterationDeterminant=str(det),
        activeModeSpectralRadius='1/sqrt(1+tau*gamma*sigma**2)',
        scope='Linearized active modes only; not a nonlinear IPM contraction proof.')


def manufactured():
    rng=np.random.default_rng(728)
    U=np.linalg.qr(rng.normal(size=(18,18)))[0]
    V=np.linalg.qr(rng.normal(size=(12,12)))[0]
    singular=np.geomspace(.3,3.,12)
    K=U[:,:12]@np.diag(singular)@V.T
    truth=rng.normal(size=12)*.15
    transform=lambda z:z+.1*np.sin(z)
    target=K@transform(truth)
    jac=lambda z:K*(1+.1*np.cos(z))[None,:]
    # This well-conditioned reference metric isolates the fixed-point mechanism.
    metric=jac(truth).T@jac(truth)
    rows=[]
    for label,extra in [('consistent',np.zeros(18)),('inconsistent',.08*U[:,12])]:
        residual=lambda z:K@transform(z)-target-extra
        z=np.zeros(12);dual=np.zeros(18);history=[]
        for iteration in range(100):
            z,dual,info=step(z,dual,residual,jac,metric)
            history.append(dict(iteration=iteration+1,**info,dualNorm=float(np.linalg.norm(dual))))
        result=dict(case=label,primalError=float(np.linalg.norm(z-truth)),
            finalResidual=float(np.linalg.norm(residual(z))),
            projectedGradient=float(np.linalg.norm(jac(z).T@residual(z))),
            dualNorm=float(np.linalg.norm(dual)),history=history)
        rows.append(result)
    assert rows[0]['primalError']<1e-10 and rows[0]['finalResidual']<1e-10
    assert abs(rows[1]['finalResidual']-.08)<1e-10
    assert rows[1]['projectedGradient']<1e-10
    assert rows[1]['dualNorm']>7.99
    return dict(description='Manufactured 12-variable, 18-equation problem; 3 groups of 4 coefficients, not an IPM run.',
        consistent=rows[0],inconsistent=rows[1],
        conclusion='An incompatible full residual is not accepted as a fixed point even when the primal least-squares gradient vanishes.')


def width_audit():
    result=[]
    for name in ['L4_257x129','L8_257x129_continued','tau20']:
        p=trajectory_frame(20) if name=='tau20' else candidate(name)
        q,Q=interp_fields(p)
        lam,mu,*_=p['rates'];alpha=float(mu/lam);power=1/(1-alpha)
        rows=[]
        for y in [0.,.125,.25,.5,1.]:
            geo=peak_measure(p,y,q,Q)
            rows.append(dict(y=y,A=geo['A'],m=geo['m'],curvatureWidth=geo['w']))
        A0=rows[0]['A'];W0=rows[0]['curvatureWidth']
        for row in rows:
            predicted=W0*(A0/row['A'])**power
            row['predictedWidth']=float(predicted)
            row['widthLawRelativeDefect']=float(row['curvatureWidth']/predicted-1)
        result.append(dict(source=name,alpha=alpha,power=power,rows=rows,
            caveat='Saved conditional roots/finite-time trajectory only. Curvature width uses cubic interpolation; no claim of converged tail law.'))
    return result


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    checks=symbolic();toy=manufactured();width=width_audit()
    (OUT/'checks.json').write_text(json.dumps(dict(symbolic=checks,manufactured=toy,widthAudit=width),indent=2)+'\n')
    print(json.dumps(dict(symbolic=checks,
        manufactured={name:{k:v for k,v in toy[name].items() if k!='history'} for name in ['consistent','inconsistent']},
        widthAudit=width),indent=2))


if __name__=='__main__':main()
