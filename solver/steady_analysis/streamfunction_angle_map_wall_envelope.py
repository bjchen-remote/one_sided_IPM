"""Fixed-angle mother with a wall-curvature-controlled smooth envelope.

This checks the predeclared 3-function and 5-function geometry only. All
parameters are observable calibrations; there is no internal-shape fitting.
"""
import os,json
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/ipm_polar_mpl')
import numpy as np
from scipy.interpolate import CubicSpline,make_interp_spline
from scipy.optimize import brentq
from streamfunction_sources import load_native
from streamfunction_angle_map import source,mapping,WINDOWS

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/streamfunction_polar_structure_20260929'

def main():
    old=np.load(ROOT/'runs/streamfunction_fresh_20260928/polar_geometry.npz')
    R,theta=old['nativeRadii'],old['thetaRadians']
    raw={}
    for tau in (10,15,20):
        z=source(tau,R,theta);p=load_native(tau)
        u=CubicSpline(p['x'],p['u1'][0]);q=CubicSpline(p['x'],p['q'][0])
        z['savedSlope0']=-R*p['P']*p['ell']*u(R/p['ell'])
        z['savedSlopePi']=R*p['P']*p['ell']*u(-R/p['ell'])
        z['second0']=-R*R*p['P']*q(R/p['ell'])
        z['secondPi']=-R*R*p['P']*q(-R/p['ell'])
        z['L0']=z['second0']/z['savedSlope0']
        z['Lpi']=z['secondPi']/z['savedSlopePi']
        raw[tau]=z
    ir=int(np.argmin(abs(R-32)))
    ref=raw[20];norm=float(np.sqrt(np.mean(ref['psi'][ir]**2)))
    mother=make_interp_spline(np.r_[0,theta,np.pi],np.r_[0,ref['psi'][ir]/norm,0],k=5,
        bc_type=([(1,ref['savedSlope0'][ir]/norm),(2,ref['second0'][ir]/norm)],
                 [(1,ref['savedSlopePi'][ir]/norm),(2,ref['secondPi'][ir]/norm)]))
    npos=np.flatnonzero(ref['psi'][ir]>0)
    target=float(brentq(mother,theta[npos[-1]],theta[npos[-1]+1]))
    ref['zero'][ir]=target
    l0=float(mother(0,nu=2)/mother(0,nu=1));lp=float(mother(np.pi,nu=2)/mother(np.pi,nu=1))
    arrays=dict(nativeRadii=R,theta=theta,referenceMother=mother(theta),
        referenceMotherDerivative1=mother(theta,nu=1),referenceMotherDerivative2=mother(theta,nu=2),
        referenceMotherSplineKnots=mother.t,referenceMotherSplineCoefficients=mother.c,
        referenceMotherSplineDegree=np.array(mother.k),referenceMotherWallRatios=np.array([l0,lp]),
        referenceMotherWallSlopes=np.array([mother(0,nu=1),mother(np.pi,nu=1)]),
        referenceMotherWallSeconds=np.array([mother(0,nu=2),mother(np.pi,nu=2)]),
        referenceZeroAngleRadians=np.array(target))
    runs=[];allfits={}
    for tau in (10,15,20):
        z=raw[tau]
        arrays[f'tau{tau}_sourcePsi']=z['psi']
        arrays[f'tau{tau}_sourceL0']=z['L0']
        arrays[f'tau{tau}_sourceLpi']=z['Lpi']
        for count in (3,5):
            rows=[];predictions=[]
            for i,RR in enumerate(R):
                th,th1,th2,s,_,_=mapping(theta,z['zero'][i],target,'half_angle_mobius')
                d0=z['L0'][i]/s-l0;dp=s*z['Lpi'][i]-lp
                eta=float((d0-dp)/4) if count==5 else 0.
                zeta=float((d0+dp)/4) if count==5 else 0.
                h=eta*np.sin(th)+zeta*np.sin(th)*np.cos(th)
                base=mother(th)*np.exp(h)
                pos=theta<z['zero'][i];neg=~pos;obs=z['psi'][i]
                pnorm=np.linalg.norm(obs[pos]);nnorm=np.linalg.norm(obs[neg])
                def eq(k):
                    zz=base*np.exp(k*np.cos(th))
                    return np.log(np.linalg.norm(zz[pos])/np.linalg.norm(zz[neg]))-np.log(pnorm/nnorm)
                k=float(brentq(eq,-40,40));unscaled=base*np.exp(k*np.cos(th))
                A=float(nnorm/np.linalg.norm(unscaled[neg]));pred=A*unscaled
                L0pred=s*(l0+2*(eta+zeta));Lppred=(lp+2*(-eta+zeta))/s
                slope0pred=A*np.exp(k)*float(mother(0,nu=1))*s
                slopepipred=A*np.exp(-k)*float(mother(np.pi,nu=1))/s
                second0pred=slope0pred*L0pred;secondpipred=slopepipred*Lppred
                row=dict(nativeR=float(RR),s=float(s),A=A,kappa=k,eta=eta,zeta=zeta,
                    positiveShapeRelativeL2=float(np.linalg.norm((pred-obs)[pos])/pnorm),
                    negativeShapeRelativeL2=float(np.linalg.norm((pred-obs)[neg])/nnorm),
                    sourceL0=float(z['L0'][i]),sourceLpi=float(z['Lpi'][i]),
                    predictedL0=float(L0pred),predictedLpi=float(Lppred),
                    sourceLProduct=float(z['L0'][i]*z['Lpi'][i]),
                    predictedLProduct=float(L0pred*Lppred),
                    wallSlope0RelativeError=float(slope0pred/z['savedSlope0'][i]-1),
                    wallSlopePiRelativeError=float(slopepipred/z['savedSlopePi'][i]-1),
                    wallSecond0RelativeError=float(second0pred/z['second0'][i]-1),
                    wallSecondPiRelativeError=float(secondpipred/z['secondPi'][i]-1),
                    positiveNormRatio=float(np.linalg.norm(pred[pos])/pnorm),
                    negativeNormRatio=float(np.linalg.norm(pred[neg])/nnorm))
                rows.append(row);predictions.append(pred)
            windows={}
            for name,(lo,hi) in WINDOWS.items():
                rr=[r for r in rows if lo<=r['nativeR']<=hi and not(tau==20 and r['nativeR']==32)]
                windows[name]=dict(nativeRadiusRange=[lo,hi],count=len(rr),
                    positiveShapeRms=float(np.sqrt(np.mean([r['positiveShapeRelativeL2']**2 for r in rr]))),
                    negativeShapeRms=float(np.sqrt(np.mean([r['negativeShapeRelativeL2']**2 for r in rr]))),
                    etaRange=[min(r['eta'] for r in rr),max(r['eta'] for r in rr)],
                    zetaRange=[min(r['zeta'] for r in rr),max(r['zeta'] for r in rr)],
                    sourceLProductRange=[min(r['sourceLProduct'] for r in rr),max(r['sourceLProduct'] for r in rr)])
            runs.append(dict(tau=tau,functions=count,rows=rows,windows=windows))
            arrays[f'tau{tau}_{count}function_prediction']=np.array(predictions)
            for key in ('A','s','kappa','eta','zeta'):
                arrays[f'tau{tau}_{count}function_{key}']=np.array([r[key] for r in rows])
            allfits[tau,count]=rows
    result=dict(candidate='Psi=A exp(kappa cosTheta+eta sinTheta+zeta sinTheta cosTheta) Phi(Theta), Theta=2atan(s tan(theta/2))',
        reference=dict(tau=20,nativeR=32,zeroAngleRadians=target,L0=l0,Lpi=lp,Lproduct=l0*lp,
            mother='Fixed measured psi with quintic interpolation, endpoint first derivative from saved u and endpoint second derivative from saved q. No curvature is invented by interpolator; not an analytic mother.'),
        calibration=dict(s='observed zero angle',A_kappa='two observed lobe norms',
            eta='(L0/s - Fsecond0/Fprime0 - s*Lpi + FsecondPi/FprimePi)/4',
            zeta='(L0/s - Fsecond0/Fprime0 + s*Lpi - FsecondPi/FprimePi)/4',
            scope='Conditional shape transfer only. All five radial functions recalibrated from observations, mother frozen. Reference row omitted from errors.'),
        endpointCaution='Matching L=psi_theta_theta/psi_theta fixes wall curvature-to-slope ratios, not the slope or q separately. Actual wall first and second derivative errors are separately recorded.',
        derivative=dict(Hprime='-kappa sinTheta+eta cosTheta+zeta cos(2Theta)',
            Hsecond='-kappa cosTheta-eta sinTheta-2zeta sin(2Theta)',
            first='A exp(H)(Fprime+Hprime F)Theta_theta',
            second='A exp(H)[(Fsecond+2Hprime Fprime+(Hsecond+Hprime^2)F)Theta_theta^2+(Fprime+Hprime F)Theta_theta_theta]'),
        runs=runs)
    (OUT/'angle_map_wall_envelope.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(OUT/'angle_map_wall_envelope.npz',**arrays)
    plot(R,theta,raw,allfits,arrays)
    report(result)
    for run in runs:
        print(run['tau'],run['functions'],{name:(w['positiveShapeRms'],w['negativeShapeRms']) for name,w in run['windows'].items()})

def plot(R,theta,raw,fits,arrays):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,3,figsize=(15,8),constrained_layout=True)
    for j,tau in enumerate((20,15,10)):
        for count,color in [(3,'#426b91'),(5,'#b25e35')]:
            rows=fits[tau,count]
            for key,style,lobe in [('positiveShapeRelativeL2','-','positive'),('negativeShapeRelativeL2','--','negative')]:
                err=np.array([r[key] for r in rows]);err[(R==32)&(tau==20)]=np.nan
                axs[0,j].loglog(R,err,style,color=color,label=f'{count} functions, {lobe}')
        axs[0,j].set(title=f'tau={tau}: interior shape with wall ratios calibrated',xlabel='Native R',ylabel='Relative lobe L2 discrepancy')
        axs[0,j].legend(fontsize=8);axs[0,j].grid(alpha=.2)
    axs[1,0].semilogx(R,[r['sourceLProduct'] for r in fits[20,3]],label='source')
    axs[1,0].axhline(fits[20,3][0]['predictedLProduct'],color='#b25e35',linestyle='--',label='3-function invariant')
    axs[1,0].set(title='tau20: forbidden constant wall-ratio product',xlabel='Native R',ylabel='L0 Lpi');axs[1,0].legend()
    for tau in (20,15,10):
        rows=fits[tau,5]
        axs[1,1].semilogx(R,[r['eta'] for r in rows],label=f'tau={tau}, eta')
        axs[1,1].semilogx(R,[r['zeta'] for r in rows],'--',label=f'tau={tau}, zeta')
    axs[1,1].set(title='Geometrically fixed envelope corrections',xlabel='Native R',ylabel='eta, zeta');axs[1,1].legend(fontsize=8)
    rows=fits[20,5]
    for key,style,name in [('wallSecond0RelativeError','-','positive-side wall q'),('wallSecondPiRelativeError','--','negative-side wall q')]:
        axs[1,2].semilogx(R,[r[key] for r in rows],style,label=name)
    axs[1,2].set(title='tau20: ratio matching does not fix q amplitude',xlabel='Native R',ylabel='Predicted q / source q - 1');axs[1,2].legend(fontsize=8)
    for ax in axs[1]:ax.grid(alpha=.2)
    fig.savefig(OUT/'angle_map_wall_envelope.png',dpi=165);plt.close(fig)

def report(j):
    lines=[r'''# 角映射的壁面二阶约束与五函数对照

令 `L0=(Ψθθ/Ψθ)(R,0)`、`Lπ=(Ψθθ/Ψθ)(R,π)`。半角映射两端二阶角导数为零，基础候选 `Ψ=A exp(κcosΘ) F(Θ)` 因而强制

\[
L_0=s\frac{F''(0)}{F'(0)},\qquad L_\pi=s^{-1}\frac{F''(\pi)}{F'(\pi)},\qquad
L_0L_\pi=\text{常数}.
\]

实测不满足该不变量。小流函数形状差异不能据此推广成壁面 q 的骨架。

限定的补充候选是

\[
\Psi=Ae^{H(\Theta)}F(\Theta),\quad H=\kappa\cos\Theta+\eta\sin\Theta+\zeta\sin\Theta\cos\Theta.
\]

令 `d0=L0/s−F″(0)/F′(0)`、`dπ=sLπ−F″(π)/F′(π)`，则

\[
\eta=(d_0-d_\pi)/4,\qquad\zeta=(d_0+d_\pi)/4.
\]

这些来自实测壁面曲率/斜率比；仍用零角定 s、两瓣范数定 A、κ，未优化内部形状。壁面导数用保存 u 与 q 计算：`Ψθ(0)=−R u1(R,0)`、`Ψθ(π)=R u1(−R,0)`、`Ψθθ=−R²q`。固定母形的评价插值同时固定其两端一、二阶导数；3 与 5 函数共用同一母形。

| τ | 函数数 | 窗口 | 正瓣差异 | 负瓣差异 |
|---:|---:|---|---:|---:|''']
    for run in j['runs']:
        for name in ('near_wall','middle','outer'):
            w=run['windows'][name]
            lines.append(f"| {run['tau']} | {run['functions']} | {name} | {w['positiveShapeRms']:.6g} | {w['negativeShapeRms']:.6g} |")
    lines+=[r'''
壁面曲率形状确实是需要保留的信息，但这一最小补充也有内部代价。τ20 中区、远区的正瓣形状均改善；负瓣中区基本不变，远区却变差。近壁正瓣明显改善，负瓣仍没有收拢。τ15 的两瓣通常都改善；τ10 中区负瓣仍保留明显差异。因此五函数候选比三函数消除了一个错误的壁面不变量，却仍不能认作跨时刻固定母形已经成立。

这里匹配的是曲率与斜率的比值，**并没有同时匹配壁面速度和 q 的幅度**。JSON 另存一、二阶端点导数误差，不把比值命中当作 q 验证。

对所有有限参数，指数权重严格正且光滑，半角映射严格单调；因此保留两个符号瓣与同一条零支。角导数为

\[
H'=-\kappa\sin\Theta+\eta\cos\Theta+\zeta\cos2\Theta,\qquad
H''=-\kappa\cos\Theta-\eta\sin\Theta-2\zeta\sin2\Theta,
\]
\[
\Psi_\theta=Ae^H(F'+H'F)\Theta_\theta,
\]
\[
\Psi_{\theta\theta}=Ae^H\{[F''+2H'F'+(H''+H'^2)F]\Theta_\theta^2+(F'+H'F)\Theta_{\theta\theta}\}.
\]

母形仍是实测经验函数；五个径向函数均由观测重新标定。这是条件形状检验，不是径向预测、时间预测或 PDE 闭合。没有继续增加自由形变来消除剩余差异。

![壁面约束对内部两瓣的代价](/Users/cbjin/Documents/one_sided_IPM/runs/streamfunction_polar_structure_20260929/angle_map_wall_envelope.png)

数据为 angle_map_wall_envelope.json/.npz；NPZ 包含每 R 的 A、s、κ、η、ζ、原始 Ψ 矩阵、母形的一二阶角导数及完整 BSpline 的 knots、coefficients、degree，可复现同一母形的导数。脚本为 solver/steady_analysis/streamfunction_angle_map_wall_envelope.py。
''']
    (OUT/'ANGLE_MAP_WALL_ENVELOPE.md').write_text('\n'.join(lines)+'\n')

if __name__=='__main__':main()
