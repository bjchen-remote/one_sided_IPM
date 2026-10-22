"""Small actual radial ODE pilot for data-learned streamfunction/scalar modes."""
from pathlib import Path
import json, os
os.environ.setdefault("OPENBLAS_NUM_THREADS","2")
os.environ.setdefault("OMP_NUM_THREADS","2")
os.environ.setdefault("MPLCONFIGDIR","/private/tmp/ipm_polar_mpl")
import numpy as np
from scipy.interpolate import RectBivariateSpline,CubicSpline
from scipy.integrate import solve_ivp
from streamfunction_sources import load_native

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/streamfunction_polar_closure_20260928'


def rms(a): return float(np.sqrt(np.mean(np.asarray(a)**2)))
def relative(a,b): return float(np.linalg.norm(a-b)/np.linalg.norm(b))


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    p=load_native(20);ell=p['ell'];lam,mu,v=map(float,p['rates'][:3])
    fs={k:RectBivariateSpline(p['y'],p['x'],p[k]) for k in ('psi','Q','q','u1','u2')}
    theta=np.linspace(0,np.pi,321);ct=np.cos(theta);st=np.sin(theta)
    coll=np.arange(len(theta))%2==0
    def truth(r):
        x=r*ct;y=r*st
        psi=fs['psi'].ev(y,x);Q=fs['Q'].ev(y,x)
        px=fs['psi'].ev(y,x,dy=1);py=fs['psi'].ev(y,x,dx=1)
        Qx=fs['Q'].ev(y,x,dy=1);Qy=fs['Q'].ev(y,x,dx=1)
        qr=ct*Qx+st*Qy;qt=r*(-st*Qx+ct*Qy)
        ar=ct*px+st*py
        u1=fs['u1'].ev(y,x);u2=fs['u2'].ev(y,x)
        ur=ct*u1+st*u2;ut=-st*u1+ct*u2
        lap=fs['psi'].ev(y,x,dy=2)+fs['psi'].ev(y,x,dx=2)
        P=lap+Qx
        T=(u1+lam*x+v)*Qx+(u2+lam*y)*Qy-mu*Q-v
        pscale=np.sqrt(np.mean(lap**2+Qx**2))
        # Use the same polar term decomposition as the model diagnostics;
        # separate Cartesian advection terms give a different denominator.
        terms=np.array([(lam*r+v*ct+ur)*qr,(-v*st+ut)*qt/r,-mu*Q,np.full_like(Q,-v)])
        tscale=np.sqrt(np.mean(np.sum(terms**2,axis=0)))
        return dict(psi=psi,Q=Q,q=fs['q'].ev(y,x),ar=ar,qr=qr,qt=qt,
                    velocity=np.array([ur,ut]),P=P,T=T,pscale=pscale,tscale=tscale)
    records=[]
    for band,trainR,start,end in [('middle',2.**np.arange(2,8.01,.5),16.,32.),('outer',2.**np.arange(8,12.01,.5),512.,1024.)]:
        observed=[truth(R/ell) for R in trainR]
        psi=np.array([o['psi'] for o in observed]);Q=np.array([o['Q'] for o in observed])
        normalized=psi/np.sqrt(np.mean(psi**2,axis=1))[:,None]
        _,_,pv=np.linalg.svd(normalized,full_matrices=False)
        # Data modes retain nonzero second derivatives at the wall. Values only
        # are pinned to zero; no sine truncation or natural spline is imposed.
        phi=pv[:3]*np.sqrt(len(theta));phi[:,[0,-1]]=0
        psp=CubicSpline(theta,phi,axis=1)
        ph,ph1,ph2=[psp(theta,nu=i).T for i in (0,1,2)]
        residual=Q-Q.mean(axis=1)[:,None]
        residual/=np.sqrt(np.mean(Q**2,axis=1))[:,None]
        _,_,gv=np.linalg.svd(residual,full_matrices=False)
        gg=np.vstack([np.ones(len(theta)),gv[:3]*np.sqrt(len(theta))])
        gsp=CubicSpline(theta,gg,axis=1)
        G,G1=gsp(theta).T,gsp(theta,nu=1).T
        r0=start/ell;r1=end/ell;sampleR=np.geomspace(start,end,17)
        initial=truth(r0)
        sourceRows=[]
        for R in sampleR:
            o=truth(R/ell)
            sourceRows.append(dict(nativeR=float(R),poissonComponentRelativeRms=rms(o['P'])/o['pscale'],
                transportComponentRelativeRms=rms(o['T'])/o['tscale'],
                poissonRawRms=rms(o['P']),transportRawRms=rms(o['T'])))
        for m in (2,3):
            Phi,Dphi,DDphi=ph[:,:m],ph1[:,:m],ph2[:,:m]
            a0=np.linalg.lstsq(Phi,initial['psi'],rcond=None)[0]
            d0=np.linalg.lstsq(Phi,initial['ar'],rcond=None)[0]
            b0=np.linalg.lstsq(G,initial['Q'],rcond=None)[0]
            state0=np.r_[a0,d0,b0]
            def solve_high(r,state,weight):
                a,d,b=state[:m],state[m:2*m],state[2*m:]
                psi=Phi@a;ar=Phi@d;qt=G1@b;Q=G@b
                ur=-(Dphi@a)/r;ut=ar
                Vr=lam*r+v*ct+ur;Vt=-v*st+ut
                knownP=ar/r+(DDphi@a)/r**2-st*qt/r
                knownT=Vt*qt/r-mu*Q-v
                # State-dependent term scales set an explicitly declared inner
                # product. They need no data outside the initial radius.
                sp=max(np.sqrt(np.mean((psi/r**2)**2+(ar/r)**2+(Q/r)**2)),1e-12)
                ts=max(np.sqrt(np.mean((mu*Q)**2+v*v+(Vr*Q/r)**2+(Vt*qt/r)**2)),1e-12)
                MP=np.column_stack([Phi,ct[:,None]*G])
                MT=np.column_stack([np.zeros_like(Phi),Vr[:,None]*G])
                M=np.vstack([MP[coll]/sp,weight*MT[coll]/ts])
                target=np.r_[-knownP[coll]/sp,-weight*knownT[coll]/ts]
                cs=np.linalg.norm(M,axis=0)
                z,_,rank,sing=np.linalg.lstsq(M/cs,target,rcond=None);z=z/cs
                P=MP@z+knownP;T=MT@z+knownT
                acc,bprime=z[:m],z[m:]
                lap=Phi@acc+ar/r+DDphi@a/r**2
                Qx=ct*(G@bprime)-st*qt/r
                qr=G@bprime
                pden=np.sqrt(np.mean(lap**2+Qx**2))
                td=np.array([Vr*qr,Vt*qt/r,-mu*Q,np.full_like(Q,-v)])
                tden=np.sqrt(np.mean(np.sum(td**2,axis=0)))
                return acc,bprime,dict(P=P,T=T,pscale=sp,tscale=ts,pden=pden,tden=tden,
                    psi=psi,Q=Q,velocity=np.array([ur,ut]),q=-lap,qx=Qx,
                    cond=float(sing[0]/sing[-1]),rank=int(rank))
            for weight in (.5,1.,2.):
                def rhs(s,state):
                    r=np.exp(s);acc,bprime,_=solve_high(r,state,weight)
                    return r*np.r_[state[m:2*m],acc,bprime]
                solution=solve_ivp(rhs,(np.log(r0),np.log(r1)),state0,method='DOP853',
                    t_eval=np.log(sampleR/ell),rtol=2e-6,atol=1e-9,max_step=np.log(2)/8)
                rows=[]
                for s,state in zip(solution.t,solution.y.T):
                    r=float(np.exp(s));_,_,z=solve_high(r,state,weight);o=truth(r)
                    rows.append(dict(nativeR=r*ell,psiRelativeL2=relative(z['psi'],o['psi']),
                        velocityRelativeL2=relative(z['velocity'],o['velocity']),rhoRelativeL2=relative(z['Q'],o['Q']),
                        qFromPsiRelativeL2=relative(z['q'],o['q']),qFromRhoRelativeL2=relative(z['qx'],o['q']),
                        poisson=dict(rawRms=rms(z['P']),balancingScale=z['pscale'],componentRelativeRms=rms(z['P'])/max(z['pden'],1e-14),
                                     heldoutRawRms=rms(z['P'][~coll]),wallValues=z['P'][[0,-1]].tolist()),
                        transport=dict(rawRms=rms(z['T']),balancingScale=z['tscale'],componentRelativeRms=rms(z['T'])/max(z['tden'],1e-14),
                                       heldoutRawRms=rms(z['T'][~coll]),wallValues=z['T'][[0,-1]].tolist()),
                        columnNormalizedCondition=z['cond'],rank=z['rank'],minimumQFromPsi=float(min(z['q'])),
                        minimumQFromRho=float(min(z['qx'])),
                        minimumSavedQ=float(min(o['q']))))
                records.append(dict(band=band,psiModes=m,rhoModes=4,transportEquationWeight=weight,
                    nativeInterval=[start,end],innerInterval=[r0,r1],trainingNativeRadii=trainR.tolist(),
                    initialState=state0.tolist(),success=bool(solution.success),message=solution.message,nfev=solution.nfev,
                    finalState=solution.y[:,-1].tolist(),rows=rows,sourceResidualBaseline=sourceRows))
                print(band,m,weight,'ok',solution.success,'nfev',solution.nfev,'end', {k:rows[-1][k] for k in ('psiRelativeL2','velocityRelativeL2','rhoRelativeL2')},
                    'P',rows[-1]['poisson']['componentRelativeRms'],'T',rows[-1]['transport']['componentRelativeRms'],flush=True)
    result=dict(title='Actual short projected radial ODE pilot',source=p['source'],tau=20,
        innerRates=dict(lambda_=lam,mu=mu,v=v),nativeToInnerEll=ell,
        state='a_j, d_j=a_j_prime, b_l for centered inner psi=sum a Phi and original scalar Q=sum b G.',
        highestDerivativeSolve='At each r solve full-angle collocation least squares for (a_second,b_prime); state is then integrated in log r using d a/d log r=r*d, d d/d log r=r*a_second, d b/d log r=r*b_prime.',
        protocol=dict(psiModes='Data SVD of row-normalized psi on fixed training radii. Compare 2 and 3 modes; endpoint values are zero, second derivatives free.',
            rhoModes='Original centered scalar Q, constant mode plus three learned angular modes of row-demeaned Q. No rho_x shape is used.',
            initialization='Project psi, its radial derivative, and Q only at the starting radius. No source values are injected after this point.',
            angularGrid='321 uniform theta points including both walls. Alternate points collocate; all 321 and held-out points are audited.',
            weighting='Poisson divided by sqrt(mean[(psi/r^2)^2+(psi_r/r)^2+(Q/r)^2]); transport divided by sqrt(mean[(mu Q)^2+v^2+(Vr Q/r)^2+(Vtheta Qtheta/r)^2]), then multiplied by weight 0.5,1,2.',
            residualNorm='Full equation residual evaluated on all angles, divided by RMS square-sum of its component terms; raw residual and holdout RMS also recorded. Velocity combines both polar components.',
            sourceBaseline='Same inner steady equations evaluated on saved u,Q and spline derivatives, including original finite-time transport residual.',
            caution='Projection is an approximate model with a nonzero angular remainder. A successful short integration is not a proof of exact closure, steady validity, or large-radius stability.'),runs=records)
    (OUT/'radial_pilot.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    write_report(result);plot(result)


def write_report(result):
    allrows=[z for run in result['runs'] for z in run['rows']]
    minpsi=min(z['minimumQFromPsi'] for z in allrows)
    minrho=min(z['minimumQFromRho'] for z in allrows)
    maxcond=max(z['columnNormalizedCondition'] for z in allrows)
    lines=['# 流函数少角模态：实际短径向 ODE 试跑','',
        '由单一起始半径上的 ψ、ψ_r、原始标量 Q 初始化，之后实际积分投影径向 ODE，没有把后续数据重新注入系数。比较 ψ 的 2/3 个数据角模态，Q 固定为常数模态加三个数据角模态。两个半径段分别独立初始化。','',
        f"12 组短积分均成功完成；列归一线性系统的最大条件数为 {maxcond:.3f}。所有试跑、所有记录角点上的 `−Δψ` 最小值为 {minpsi:.8g}，`Q_x` 最小值为 {minrho:.8g}。两种 q 的正性分别核查，因为 Poisson 余量非零时二者不相同。这个检查只覆盖所列有限半径与角点，不能宣称全域正性。",'',
        '| 区间 | ψ 模态 | 输运权重 | 终点 ψ 误差 | 终点速度误差 | 终点 Q 误差 | 终点 Poisson 残差 | 终点输运残差 |',
        '|---|---:|---:|---:|---:|---:|---:|---:|']
    for run in result['runs']:
        z=run['rows'][-1]
        lines.append(f"| {run['band']} | {run['psiModes']} | {run['transportEquationWeight']:g} | {100*z['psiRelativeL2']:.3f}% | {100*z['velocityRelativeL2']:.3f}% | {100*z['rhoRelativeL2']:.3f}% | {100*z['poisson']['componentRelativeRms']:.3f}% | {100*z['transport']['componentRelativeRms']:.3f}% |")
    lines += ['', '误差单位为内坐标归一化后的相对 L2；速度合并两个极分量。PDE 残差使用完整方程在所有角点求值，按各项平方和 RMS 归一，包含未用于 collocation 的角点。它不是只看被投影到有限维空间后的残差。','',
        '## 与原始数据的比较','',
        '| 区间 | 源场终点 Poisson 残差 | 源场终点输运残差 |','|---|---:|---:|']
    for band in ('middle','outer'):
        run=next(z for z in result['runs'] if z['band']==band);s=run['sourceResidualBaseline'][-1]
        lines.append(f"| {band} | {100*s['poissonComponentRelativeRms']:.4f}% | {100*s['transportComponentRelativeRms']:.3f}% |")
    lines += ['', '中心化内坐标半径为 `r=R/ell`，本次 `ell=0.6009309468198614`；率为 `lambda=.6073592379, mu=.2219323977, v=.1747320183`。原始 Q 方程保留右端常数 v。Psi 角模态只强制端点值为零，不强制二阶导数为零；未换成有限 sine 模态。','',
        '权重 0.5/1/2 是同一个短区间模型的必要敏感性对照。每步解 `(a_second,b_prime)` 的列归一条件数、逐半径完整残差、q 误差、壁面残差、初始化与终值均保存在 `radial_pilot.json`。有限时间源场自身并不满足零输运残差，不能以零作为数据真实性基线。','',
        '![短径向 ODE 误差与完整残差](radial_pilot.png)']
    (OUT/'RADIAL_PILOT.md').write_text('\n'.join(lines)+'\n')


def plot(result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
    for row,band in enumerate(('middle','outer')):
        for m,color in ((2,'#bc6054'),(3,'#397fa0')):
            run=next(z for z in result['runs'] if z['band']==band and z['psiModes']==m and z['transportEquationWeight']==1.)
            R=[z['nativeR'] for z in run['rows']]
            for key,style,label in [('psiRelativeL2','-','psi'),('velocityRelativeL2','--','u'),('rhoRelativeL2',':','Q')]:
                axs[row,0].semilogy(R,[100*z[key] for z in run['rows']],style,color=color,label=f'm={m}, {label}')
            for key,style in [('poisson','-'),('transport','--')]:
                axs[row,1].plot(R,[100*z[key]['componentRelativeRms'] for z in run['rows']],style,color=color,label=f'm={m}, {key}')
        for key,style in [('poissonComponentRelativeRms','-'),('transportComponentRelativeRms','--')]:
            axs[row,1].plot(R,[100*z[key] for z in run['sourceResidualBaseline']],style,color='gray',alpha=.6,label='source '+key.split('Component')[0])
        axs[row,0].set(title=f'{band}: actual radial prediction',xlabel='Native radius R',ylabel='Relative L2 error (%)')
        axs[row,1].set(title=f'{band}: full angular equation residuals',xlabel='Native radius R',ylabel='Component-scaled RMS (%)')
        for ax in axs[row]:ax.legend(fontsize=8,ncol=2);ax.grid(alpha=.2)
    fig.suptitle('Projected radial ODE | fixed angular modes | one-radius initialization | weight=1')
    fig.savefig(OUT/'radial_pilot.png',dpi=160);plt.close(fig)


if __name__=='__main__':main()
