"""Inspect a saved trajectory and attempt a conditional centered steady root.

fixed_point.py is copied without modifications from Research/IPM/fixed_point_lab.
The source snapshot exterior is retained; no reflection symmetry is imposed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.io import loadmat, savemat
from scipy.interpolate import CubicSpline, RectBivariateSpline
from fixed_point import ProfileProblem, SolverOptions, solve


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False)+'\n')


def inspect(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d = loadmat(folder/'run_analysis.mat', simplify_cells=True)
    x, tau, geom, wall = (d[k] for k in ('xFull','tau','geometry','wallQ'))
    radii = np.array([0.25,0.5,1.,2.])
    asym = np.zeros((len(tau),len(radii)))
    for i in range(len(tau)):
        a,p,_,ell = geom[i]
        f = CubicSpline(x,wall[i])
        for j,r in enumerate(radii):
            s = np.linspace(0,r,501)
            asym[i,j] = np.max(abs(f(a+s)-f(a-s)))/p
    indices = sorted(set(np.argmin(abs(tau-t)) for t in (0,5,10,15,20)))
    rows = [dict(tau=float(tau[i]),a=float(geom[i,0]),P=float(geom[i,1]),
                 ell=float(geom[i,3]),
                 asymmetry={str(r):float(asym[i,j]) for j,r in enumerate(radii)})
            for i in indices]
    # Compare right/left values at identical distances from the tracked peak.
    a,p,_,ell = geom[-1]
    f = CubicSpline(x,wall[-1])
    paired = {str(r):[float(f(a-r)),float(f(a+r))] for r in radii}
    u = CubicSpline(x,d['u1'][0])
    q = f(a)
    # Independent local continuous wall-identity audit (not a discrete gauge).
    meta = d['metadata']
    speed = float(u(a)+meta['c_l']*a+meta['c_r'])
    derivative_audit = dict(wallSpeedAtPeak=speed,
        strainAtPeak=float(u(a,1)),
        cwMinusCl=float(meta['c_omega']-meta['c_l']),
        peakQxSpline=float(f(a,1)))
    core = (abs(d['xFull'][None,:]-a)<=2)&(d['yFull'][:,None]<=2)
    center_rate = float(CubicSpline(x,d['rhoRate'][0])(a))
    residuals = dict(fullRhoRate=float(np.max(abs(d['rhoRate']))),
        centeredCoreRhoRate=float(np.max(abs(d['rhoRate'][core]-center_rate))),
        fullRhoXRate=float(np.max(abs(d['rhoXRate']))),
        centerDensityRate=center_rate)
    summary = dict(symmetryDefinition='max_{0<=s<=r}|R_X(a+s,0)-R_X(a-s,0)|/P',
                   radiusUnits='saved rescaled X coordinate', rows=rows,
                   finalLeftRight=paired,derivativeAudit=derivative_audit,
                   residuals=residuals)
    dump(folder/'trajectory_analysis.json',summary)
    print(json.dumps(summary,indent=2),flush=True)
    fig,ax=plt.subplots(1,3,figsize=(14.4,4.2),constrained_layout=True)
    s=np.linspace(0,2.5,501)
    ax[0].plot(s,f(a-s)/p,label='left: q(a-s)/P',lw=2)
    ax[0].plot(s,f(a+s)/p,label='right: q(a+s)/P',lw=2)
    ax[0].set(xlabel='Distance from wall peak s',ylabel='Normalized wall gradient',title='Final wall profile, tau = 20')
    ax[0].legend()
    for j,r in enumerate(radii):
        ax[1].semilogy(tau,asym[:,j],label=f'r = {r:g}')
    ax[1].set(xlabel='Rescaled time tau',ylabel='Left/right difference / peak',title='Finite-width symmetry error')
    ax[1].legend()
    for i in indices[1:]:
        a0,p0,_,e0=geom[i]
        z=np.linspace(-3,3,601)
        ax[2].plot(z,CubicSpline(x,wall[i])(a0+e0*z)/p0,label=f'tau = {tau[i]:.1f}')
    ax[2].set(xlabel='Inner coordinate (X-a)/ell',ylabel='q/P',title='Centered and width-normalized profiles')
    ax[2].legend()
    for aplot in ax: aplot.grid(alpha=.25)
    fig.savefig(folder/'wall_symmetry.png',dpi=180)
    plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    s=np.linspace(0,2,501)
    for i in indices[1:]:
        a,p,_,ell=geom[i];f=CubicSpline(x,wall[i])
        plus=f(a+s)/p;minus=f(a-s)/p
        axs[0].plot(s,(plus+minus)/2,label=f'tau={tau[i]:.0f}')
        axs[1].plot(s,(plus-minus)/2,label=f'tau={tau[i]:.0f}')
    for ax,title in zip(axs,['Even part of the wall gradient','Odd part of the wall gradient']):
        ax.set(xlabel='Rescaled distance s from peak',ylabel='Fraction of peak P',title=title)
        ax.grid(alpha=.25);ax.legend()
    fig.savefig(folder/'wall_parity.png',dpi=180);plt.close(fig)
    a,p,_,ell=geom[-1]
    q=RectBivariateSpline(d['yFull'],d['xFull'],d['rhoX'])
    by_height={str(y):float(np.max(abs(q.ev(y,a+s)-q.ev(y,a-s)))/(2*p)) for y in [0,.25,.5,1,2]}
    def parts(i):
        a,p,_,_=geom[i];f=CubicSpline(x,wall[i])
        return np.array([(f(a+s)+f(a-s))/(2*p),(f(a+s)-f(a-s))/(2*p)])
    terminal_parts=parts(-1)
    changes={str(int(round(tau[i]))):np.max(abs(parts(i)-terminal_parts),axis=1).tolist()
             for i in indices[1:-1]}
    dump(folder/'parity_analysis.json',dict(
        definition='max_{0<=s<=2}|R_X(a+s,Y)-R_X(a-s,Y)|/(2P)',
        oddAmplitudeByHeight=by_height,
        wallChangeToTau20EvenThenOdd=changes))


def newton(seed,out,iterations=12,stencil=5,restart=None):
    if out.exists():
        raise ValueError('Refusing to overwrite an experiment')
    d = loadmat(seed,simplify_cells=True)
    x,y,q,boundary,rates = (np.asarray(d[k]) for k in ('x','y','Q','psi_boundary','rates'))
    if q.size>35000 or not 0<=iterations<=60:
        raise ValueError('This pilot is bounded to 35,000 field nodes and 60 iterations')
    problem = ProfileProblem(x,y,boundary,stencil=stencil,green_reference=q)
    initial_q=q
    if restart is not None:
        with np.load(restart) as previous:
            assert np.array_equal(x,previous['x']) and np.array_equal(y,previous['y'])
            initial_q=previous['Q'];rates=previous['rates']
    normalized = problem.normalize_seed(initial_q)
    z0 = np.r_[normalized,rates]
    r,c=problem.evaluate(z0)
    # Test the exact nonlocal Jacobian against a central directional difference.
    rng=np.random.default_rng(230928)
    direction=rng.standard_normal(z0.size)
    direction/=np.linalg.norm(direction)
    epsilon=1e-5
    difference=(problem.evaluate(z0+epsilon*direction)[0]-problem.evaluate(z0-epsilon*direction)[0])/(2*epsilon)
    exact=problem.jvp(direction,c)
    jacobian_error=float(np.linalg.norm(difference-exact)/np.linalg.norm(exact))
    assert jacobian_error<1e-5, jacobian_error
    out.mkdir()
    initial=problem.metrics(z0,c)
    config=dict(seed=str(seed.resolve()),seed_sha256=hashlib.sha256(seed.read_bytes()).hexdigest(),
                stencil=stencil,closure='frozen_source_exterior',
                restart=str(restart) if restart else None,
                initial=initial,normalizationCorrection=float(np.max(abs(normalized-initial_q.ravel()))),
                jacobianDirectionalRelativeError=jacobian_error)
    dump(out/'config.json',config)
    start=time.time()
    with (out/'iterations.jsonl').open('w') as stream:
        def callback(row):
            stream.write(json.dumps(row)+'\n');stream.flush()
            print(f"iter={row['iteration']} residual={row['residual_inf']:.6e} "
                  f"rates=({row['cl']:.7g},{row['cw']:.7g},{row['drift']:.7g}) "
                  f"inflow={row['inflow_boundary_count']}",flush=True)
        result=solve(problem,z0,SolverOptions(max_iterations=iterations),callback)
    z=result.pop('z')
    _,c=problem.evaluate(z)
    result['seconds']=time.time()-start
    result['interpretation']='Conditional finite-box solution with frozen exterior, not a verified half-plane steady state.'
    dump(out/'report.json',result)
    np.savez_compressed(out/'state.npz',x=x,y=y,Q=z[:problem.N].reshape(q.shape),
                        rates=z[problem.N:],psi=c['psi'].reshape(q.shape),
                        residual=c['raw'].reshape(q.shape),
                        qx=c['qx'].reshape(q.shape),seedQ=q)
    print(json.dumps({k:v for k,v in result.items() if k!='history'},indent=2),flush=True)


def seed_grid(folder,L,nx,ny):
    """Export another box directly from the native final data, not a root."""
    d=loadmat(folder/'run_analysis.mat',simplify_cells=True)
    meta=d['metadata']; a,p,ell=(meta[k] for k in ('a','P','ell'))
    x=np.linspace(-L,L,nx);y=np.linspace(0,L,ny)
    rho=RectBivariateSpline(d['yFull'],d['xFull'],d['rho'])
    psi=RectBivariateSpline(d['yFull'],d['xFull'],d['psi']+meta['uAtPeak']*d['yFull'][:,None])
    q=(rho(ell*y,a+ell*x)-rho(0,a)[0,0])/(p*ell)
    psiq=psi(ell*y,a+ell*x)/(p*ell**2)
    boundary=np.zeros_like(q);boundary[-1]=psiq[-1];boundary[:,[0,-1]]=psiq[:,[0,-1]];boundary[0]=0
    out=folder/f'seed_L{L:g}_{nx}x{ny}.mat'
    assert not out.exists()
    savemat(out,dict(x=x,y=y,Q=q,psi_boundary=boundary,rates=meta['innerRates'],source='native run_analysis.mat via bicubic spline'))
    print(out)


def compare(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    native=loadmat(folder/'run_analysis.mat',simplify_cells=True)
    names=['L4_129x65','L4_257x129','L6_193x97','L8_257x129_continued']
    labels=['L=4, h=1/16','L=4, h=1/32','L=6, h=1/16','L=8, h=1/16']
    xx=np.linspace(-2,2,201);yy=np.linspace(0,2,101)
    profiles={};rows=[];candidates={}
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True)
    a,p,_,ell=native['geometry'][-1]
    wall_seed=CubicSpline(native['xFull'],native['wallQ'][-1])
    z=np.linspace(-3,3,601)
    axes[0].plot(z,wall_seed(a+ell*z)/p,'k--',label='tau=20 seed',lw=2)
    for name,label in zip(names,labels):
        state=np.load(folder/name/'state.npz')
        report=json.loads((folder/name/'report.json').read_text())
        m=report['metrics']
        candidates[name]=dict(x=state['x'],y=state['y'],Q=state['Q'],
            Qxi=state['qx'],psi=state['psi'],rates=state['rates'],
            residual=state['residual'],status='conditional finite-box root',
            closure='native finite-time exterior frozen',infiniteDomainVerified=False)
        profile=RectBivariateSpline(state['y'],state['x'],state['qx'])
        profiles[name]=profile(yy,xx)
        wall=CubicSpline(state['x'],state['qx'][0])
        s=np.linspace(0,2,501)
        rows.append(dict(name=name,cl=m['cl'],cw=m['cw'],v=m['drift'],alpha=m['alpha'],
            residual=m['residual_inf'],inflow=m['inflow_boundary_count'],status=report['status'],
            innerRadius2Asymmetry=float(np.max(abs(wall(s)-wall(-s))))))
        axes[0].plot(z,wall(z),label=label)
        hist=report['history']
        axes[1].semilogy([v['iteration'] for v in hist],[v['residual_inf'] for v in hist],label=label)
    base=profiles[names[0]]
    for row in rows:
        delta=profiles[row['name']]-base
        row['coreQxRelativeL2ToL4']=float(np.linalg.norm(delta)/np.linalg.norm(base))
        row['coreQxLinfToL4']=float(np.max(abs(delta)))
    # Re-evaluate the fine-grid root with a different differential stencil.
    fine=np.load(folder/'L4_257x129'/'state.npz')
    seed=loadmat(folder/'seed_L4_257x129.mat',simplify_cells=True)
    alternate=ProfileProblem(fine['x'],fine['y'],seed['psi_boundary'],
                             stencil=7,green_reference=seed['Q'])
    alternate_metrics=alternate.metrics(np.r_[fine['Q'].ravel(),fine['rates']])
    dump(folder/'steady_comparison.json',dict(rows=rows,
        coreDomain='|xi|<=2, 0<=eta<=2, sampled uniformly',
        fineRootReevaluatedWithStencil7=alternate_metrics))
    savemat(folder/'steady_candidates.mat',dict(candidates=candidates,
        normalization='Q(0)=0,Qxi(0)=1,Qxixi(0)=0,Qxixixi(0)=-2',
        rateConvention='rates=(lambda,mu,v,b); b=v; centered frame and ds/dtau=P',
        pde='(U-U(0)+lambda*z+v*e1).grad Q - mu*Q - b = 0',
        sourceRun='run_20260927_133548_651, tau=20',
        infiniteDomainVerified=False),do_compression=True)
    print(json.dumps(dict(rows=rows,fineRootReevaluatedWithStencil7=alternate_metrics),indent=2))
    axes[0].set(xlabel='Inner horizontal coordinate xi',ylabel='Q_xi(xi,0)',title='Conditional steady profiles')
    axes[1].set(xlabel='Newton iteration (L=8: continuation)',ylabel='Discrete steady residual infinity norm',title='Nonlinear convergence on each grid')
    for ax in axes:ax.legend(fontsize=9);ax.grid(alpha=.25)
    fig.savefig(folder/'steady_candidates.png',dpi=180)
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('inspect');a.add_argument('folder',type=Path)
    a=sub.add_parser('compare');a.add_argument('folder',type=Path)
    a=sub.add_parser('solve');a.add_argument('seed',type=Path);a.add_argument('out',type=Path)
    a.add_argument('--iterations',type=int,default=12);a.add_argument('--stencil',type=int,default=5)
    a.add_argument('--restart',type=Path)
    a=sub.add_parser('seed');a.add_argument('folder',type=Path);a.add_argument('L',type=float);a.add_argument('nx',type=int);a.add_argument('ny',type=int)
    a=p.parse_args()
    if a.command=='inspect':inspect(a.folder)
    elif a.command=='compare':compare(a.folder)
    elif a.command=='seed':seed_grid(a.folder,a.L,a.nx,a.ny)
    else:newton(a.seed,a.out,a.iterations,a.stencil,a.restart)
