"""Small, geometry-calibrated angular-map comparison on direct psi data.

The fixed mother is one observed tau20 half-circle. Radial functions are
measured from zero angle, two lobe norms, and optionally endpoint slope ratio.
No free angular spline deformation, radial law, or PDE evolution is fitted.
"""
from pathlib import Path
import json,os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/ipm_polar_mpl')
import numpy as np
from scipy.interpolate import CubicSpline,RectBivariateSpline
from scipy.optimize import brentq
from streamfunction_sources import load_native

ROOT=Path(__file__).resolve().parents[2]
OLD=ROOT/'runs/streamfunction_fresh_20260928'
OUT=ROOT/'runs/streamfunction_polar_structure_20260929'
FAMILIES=('half_angle_mobius','interval_projective','half_angle_plus_slope_warp')
WINDOWS={'near_wall':(.0625,.5),'transition':(.5,4.),'middle':(4.,256.),'outer':(256.,4096.)}
EPS_BOUND=.45


def source(tau,R,theta):
    p=load_native(tau)
    f=RectBivariateSpline(p['y'],p['x'],p['psi'])
    RR,TT=np.meshgrid(R,theta,indexing='ij')
    x=RR*np.cos(TT)/p['ell'];y=RR*np.sin(TT)/p['ell']
    if x.min()<p['x'][0] or x.max()>p['x'][-1] or y.max()>p['y'][-1]:raise ValueError('No extrapolation')
    psi=p['P']*p['ell']**2*f.ev(y,x)
    left=R*p['P']*p['ell']*f.ev(np.zeros_like(R),R/p['ell'],dx=1)
    right=-R*p['P']*p['ell']*f.ev(np.zeros_like(R),-R/p['ell'],dx=1)
    zero=[]
    for row,l,r in zip(psi,left,right):
        cs=CubicSpline(np.r_[0,theta,np.pi],np.r_[0,row,0],bc_type=((1,l),(1,r)))
        roots=cs.roots(extrapolate=False)
        roots=roots[(roots>1e-5)&(roots<np.pi-1e-5)]
        if len(roots)!=1:raise ValueError('Expected exactly one observed interior nodal angle')
        zero.append(float(roots[0]))
    return dict(psi=psi,zero=np.array(zero),slope0=left,slopepi=right,source=p['source'])


def mapping(theta,zero,target,family,epsilon=0.):
    if family=='interval_projective':
        u=zero/np.pi;u0=target/np.pi;s=u0*(1-u)/(u*(1-u0))
        D=1+(s-1)*theta/np.pi
        q=np.pi*s*(theta/np.pi)/D
        q1=s/D**2;q2=-2*s*(s-1)/(np.pi*D**3)
        end0=s;endpi=1/s
    else:
        s=np.tan(target/2)/np.tan(zero/2)
        D=np.cos(theta/2)**2+s*s*np.sin(theta/2)**2
        q=2*np.arctan2(s*np.sin(theta/2),np.cos(theta/2))
        q1=s/D;q2=-s*(s*s-1)*np.sin(theta)/(2*D*D)
        end0=s;endpi=1/s
    if family=='half_angle_plus_slope_warp':
        c0=np.cos(target)
        w=np.sin(q)*(np.cos(q)-c0)
        w1=np.cos(2*q)-c0*np.cos(q)
        w2=-2*np.sin(2*q)+c0*np.sin(q)
        q,q1,q2=q+epsilon*w,(1+epsilon*w1)*q1,epsilon*w2*q1*q1+(1+epsilon*w1)*q2
        end0*=1+epsilon*(1-c0);endpi*=1+epsilon*(1+c0)
    assert np.all(q1>0)
    return q,q1,q2,float(s),float(end0),float(endpi)


def fit_geometry(obs,zero,slope_ratio,mother,reference_zero,theta,family):
    positive=theta<zero;negative=theta>zero
    plusnorm=np.linalg.norm(obs[positive]);minusnorm=np.linalg.norm(obs[negative])
    logratio=np.log(plusnorm/minusnorm)
    ref_slope_ratio=abs(mother(0,1)/mother(np.pi,1))
    def compute(epsilon):
        q,q1,q2,s,e0,epi=mapping(theta,zero,reference_zero,family,epsilon)
        base=mother(q)
        def ratio_equation(k):
            w=base*np.exp(k*np.cos(q))
            return np.log(np.linalg.norm(w[positive])/np.linalg.norm(w[negative]))-logratio
        k=brentq(ratio_equation,-12,12,xtol=1e-12)
        w=base*np.exp(k*np.cos(q));A=minusnorm/np.linalg.norm(w[negative])
        pred=A*w
        endratio=2*k+np.log(ref_slope_ratio*e0/epi)-np.log(slope_ratio)
        return dict(theta=q,dtheta=q1,ddtheta=q2,s=s,kappa=k,A=A,pred=pred,
                    epsilon=epsilon,endpointLogRatioResidual=float(endratio),end0=e0,endpi=epi)
    if family=='half_angle_plus_slope_warp':
        # One scalar geometric condition, not minimization of shape error.
        grid=np.linspace(-EPS_BOUND,EPS_BOUND,13)
        f=[compute(e)['endpointLogRatioResidual'] for e in grid]
        roots=[]
        for i in range(len(grid)-1):
            if f[i]==0:roots.append(float(grid[i]))
            if f[i]*f[i+1]<0:
                roots.append(float(brentq(lambda e:compute(e)['endpointLogRatioResidual'],grid[i],grid[i+1],xtol=1e-11)))
        if f[-1]==0:roots.append(float(grid[-1]))
        if not roots:return dict(available=False,reason='Endpoint slope ratio has no root in the preset safe monotonic epsilon interval',endpointResidualBounds=[f[0],f[-1]])
        epsilon=min(roots,key=abs)
        fit=compute(epsilon);fit['endpointRootCount']=len(roots)
    else:
        fit=compute(0.);fit['endpointRootCount']=None
    fit.update(available=True,
        positiveShapeRelativeL2=float(np.linalg.norm((fit['pred']-obs)[positive])/plusnorm),
        negativeShapeRelativeL2=float(np.linalg.norm((fit['pred']-obs)[negative])/minusnorm),
        positiveNormRatio=float(np.linalg.norm(fit['pred'][positive])/plusnorm),
        negativeNormRatio=float(np.linalg.norm(fit['pred'][negative])/minusnorm),
        observedLobeNormRatio=float(np.exp(logratio)),
        mapMinimumAngularDerivative=float(fit['dtheta'].min()))
    return fit


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    z=np.load(OLD/'polar_geometry.npz');R=z['nativeRadii'];theta=z['thetaRadians']
    raw={tau:source(tau,R,theta) for tau in (10,15,20)}
    iref=int(np.argmin(abs(R-32.)))
    ref=raw[20]['psi'][iref];scale=np.sqrt(np.mean(ref**2))
    mother=CubicSpline(np.r_[0,theta,np.pi],np.r_[0,ref/scale,0],
        bc_type=((1,raw[20]['slope0'][iref]/scale),(1,raw[20]['slopepi'][iref]/scale)))
    target=raw[20]['zero'][iref]
    results=[];allfits={};arrays=dict(nativeRadii=R,theta=theta,referenceMother=ref/scale)
    for tau in (10,15,20):
        p=raw[tau]
        for family in FAMILIES:
            fits=[];rows=[]
            for i,RR in enumerate(R):
                f=fit_geometry(p['psi'][i],p['zero'][i],abs(p['slope0'][i]/p['slopepi'][i]),mother,target,theta,family)
                fits.append(f)
                rows.append(dict(nativeR=float(RR),observedZeroAngleDegrees=float(p['zero'][i]*180/np.pi),
                    **{k:float(v) if isinstance(v,np.floating) else v for k,v in f.items() if k not in ('theta','dtheta','ddtheta','pred')}))
            windows={}
            for name,(lo,hi) in WINDOWS.items():
                requested=[r for r in rows if lo<=r['nativeR']<=hi]
                valid=[r for r in requested if r['available']]
                evaluated=[r for r in valid if not (tau==20 and r['nativeR']==float(R[iref]))]
                windows[name]=dict(nativeRadiusRange=[lo,hi],requestedRows=len(requested),availableRows=len(valid),
                    evaluatedRows=len(evaluated),referenceExcluded=len(valid)!=len(evaluated),
                    positiveShapeRms=float(np.sqrt(np.mean([r['positiveShapeRelativeL2']**2 for r in evaluated]))) if evaluated else None,
                    negativeShapeRms=float(np.sqrt(np.mean([r['negativeShapeRelativeL2']**2 for r in evaluated]))) if evaluated else None,
                    epsilonRange=[min(r['epsilon'] for r in valid),max(r['epsilon'] for r in valid)] if valid else None)
            results.append(dict(tau=tau,family=family,rows=rows,windows=windows))
            allfits[tau,family]=fits
            arrays[f'tau{tau}_{family}_mappedTheta']=np.array([f.get('theta',np.full_like(theta,np.nan)) for f in fits])
            arrays[f'tau{tau}_{family}_prediction']=np.array([f.get('pred',np.full_like(theta,np.nan)) for f in fits])
    result=dict(title='Geometry-calibrated smooth angular maps with independent lobe amplitudes',
        reference=dict(tau=20,nativeR=float(R[iref]),zeroAngleRadians=float(target),zeroAngleDegrees=float(target*180/np.pi),
            choice='R=32 is the geometric midpoint of the pre-existing middle band 4..256; not chosen by model error.',
            mother='One measured centred native psi half-circle, normalized by its angular RMS. Cubic interpolation is only a fixed empirical evaluation device, with measured endpoint derivatives. It is not a claimed analytic solution.',
            evaluation='The self-reference tau20 R32 row is excluded from window shape scores and the shape-discrepancy plot; at tau10/15 R32 remains a genuine transfer observation.'),
        candidate='Psi(R,theta)=A(R)*exp(kappa(R)*cos(Theta(R,theta)))*Phi_star(Theta(R,theta))',
        calibration=dict(s='Determined exactly by observed interior zero angle and fixed reference zero angle.',
            A_kappa='Two observed lobe L2 norms determine A and kappa. This is two independent amplitude functions with a smooth envelope, not a single global amplitude or a piecewise amplitude jump.',
            epsilon='Only in the third predeclared family: determined by the observed ratio of the two endpoint angular derivatives. No angular shape error is minimized.',
            safeEpsilonInterval=[-EPS_BOUND,EPS_BOUND],
            transfer='The same tau20 R32 mother is frozen at all radii and at tau10/15; measured zero, lobe norms and endpoint ratio are recalibrated. This is conditional shape transfer, not time or radial prediction.'),
        formulas=dict(halfAngle='Theta=2 atan(s tan(theta/2)); s=tan(Theta0/2)/tan(theta0/2)',
            intervalProjective='Theta=pi*s*(theta/pi)/[1+(s-1)*(theta/pi)]; s=(Theta0/pi)*(1-theta0/pi)/[(theta0/pi)*(1-Theta0/pi)]',
            smoothCorrection='phi=M_s(theta); Theta=phi+epsilon*sin(phi)*(cos(phi)-cos(Theta0)); preserves endpoints and the reference node.',
            amplitudeUniqueness='d/dkappa log(norm_positive/norm_negative)=mean_positive(cos Theta)-mean_negative(cos Theta)>0, because the two lobes are ordered by Theta and separated by Theta0.',
            angularDerivative1='Psi_theta=A exp(kappa cos Theta) (Phi_prime-kappa sin Theta Phi)*Theta_theta',
            angularDerivative2='Psi_thetatheta=A exp(kappa cos Theta)*[(Phi_second-2kappa sin Theta Phi_prime+(kappa^2 sin^2 Theta-kappa cos Theta)Phi)*Theta_theta^2+(Phi_prime-kappa sin Theta Phi)*Theta_thetatheta]'),
        runs=results,limitations=['The mother is empirical, so this is a concrete geometric organization, not symbolic closure.',
            'Matching both lobe norms and the zero does not certify either lobe shape; separate lobe errors are reported.',
            'The map is angularly smooth and strictly monotone in the stated interval, but radial smoothness and origin compatibility require additional equations for its radial parameters.',
            'The extra slope warp is a geometric sensitivity test, not a guarantee of better shape matching.',
            'No claim of global velocity/Poisson/transport closure is made; no PDE was evolved.'])
    (OUT/'angle_map.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(OUT/'angle_map.npz',**arrays)
    plot(R,theta,raw,allfits,mother,target)
    report(result)
    for run in results:
        print(run['tau'],run['family'], {b:run['windows'][b] for b in ('middle','outer')})


def plot(R,theta,raw,fits,mother,target):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,3,figsize=(15,8),constrained_layout=True)
    selected=[.0625,.5,4,16,64,256,1024,4096]
    for i,RR in enumerate(selected):
        k=np.argmin(abs(R-RR));f=fits[20,'half_angle_mobius'][k]
        mapped=f['theta'];profile=raw[20]['psi'][k]/(f['A']*np.exp(f['kappa']*np.cos(mapped)))
        color=plt.cm.viridis(i/(len(selected)-1))
        for ax,mask in [(axs[0,0],mapped<=target),(axs[0,1],mapped>=target)]:
            ax.plot(mapped[mask]*180/np.pi,profile[mask],color=color,label=f'R={RR:g}')
    for ax,low,high in [(axs[0,0],0,target),(axs[0,1],target,np.pi)]:
        q=np.linspace(low,high,301);ax.plot(q*180/np.pi,mother(q),'k--',lw=1.5,label='fixed mother')
        ax.set(xlabel='Mapped angle Theta (degrees)',ylabel='Psi / [A exp(kappa cos Theta)]');ax.grid(alpha=.15)
    axs[0,0].set_title('Positive lobe after half-angle map, tau=20');axs[0,0].legend(fontsize=7,ncol=2)
    axs[0,1].set_title('Negative lobe after the same map')
    colors={10:'#a7612b',15:'#40957e',20:'#315889'}
    for tau in (10,15,20):
        ff=fits[tau,'half_angle_mobius']
        axs[0,2].semilogx(R,[f['s'] for f in ff],color=colors[tau],label=f'tau={tau}')
        axs[1,0].semilogx(R,[f['kappa'] for f in ff],color=colors[tau],label=f'tau={tau}')
        ef=fits[tau,'half_angle_plus_slope_warp']
        axs[1,1].semilogx(R,[f.get('epsilon',np.nan) for f in ef],color=colors[tau],label=f'tau={tau}')
        for field,style,lobe in [('positiveShapeRelativeL2','-','positive'),('negativeShapeRelativeL2','--','negative')]:
            errors=np.array([f[field] for f in ff])
            if tau==20:errors[R==32.]=np.nan
            axs[1,2].loglog(R,errors,style,color=colors[tau],label=f'tau={tau}, {lobe}')
    axs[0,2].set(title='Zero-angle map parameter',xlabel='Native radius R',ylabel='s(R)');axs[0,2].legend()
    axs[1,0].set(title='Independent lobe contrast envelope',xlabel='Native radius R',ylabel='kappa(R)');axs[1,0].legend()
    axs[1,1].set(title='Optional endpoint-slope correction',xlabel='Native radius R',ylabel='epsilon(R)');axs[1,1].legend()
    axs[1,2].set(title='Frozen mother: shape check (self-reference omitted)',xlabel='Native radius R',ylabel='Relative lobe L2 discrepancy');axs[1,2].legend(fontsize=7,ncol=2)
    for ax in axs.flat:ax.grid(alpha=.15)
    fig.suptitle('Smooth angular geometry + two independent lobe amplitudes | conditional observations',fontsize=14)
    fig.savefig(OUT/'angle_map.png',dpi=165);plt.close(fig)


def report(j):
    lines=[r'''# 用光滑角映射组织直接流函数数据

候选不是固定 rank 的角展开，而是

\[
\boxed{\Psi(R,\theta)=A(R)e^{\kappa(R)\cos\Theta(R,\theta)}\Phi_*(\Theta(R,\theta)).}
\]

固定母形 Φ* 取 τ20、原生 R=32 的实测半圆剖面；这个半径是既有中区 4…256 的几何中点，没有按误差调选。Φ* 目前仍是经验函数，不是已经找到的符号解。

## 1. 为什么保留两个幅度函数

若在正、负瓣直接乘两个不同常数，零线两侧的一阶角导数一般跳跃。这里用全角光滑正权重 `A exp(κ cosΘ)`，由实测正、负瓣的两个 L2 幅度分别确定 A、κ；不把小正瓣压进一个由负瓣主导的整体幅度。

对已知角映射，预测的正负瓣幅度比随 κ 严格单调：其对数导数是两瓣加权平均 cosΘ 的差，而正瓣的全部 Θ 小于负瓣的 Θ。因此两个幅度条件唯一决定 κ 和 A。匹配两瓣幅度并不等于匹配瓣形，后文分别检查两瓣形状。

## 2. 两个最简端点固定映射

令 θ₀(R) 为实测内部零角，Θ₀ 为固定母形零角。比较两个事先确定的映射：

\[
\text{半角 Möbius:}\quad
\Theta=2\arctan\!\left(s\tan\frac\theta2\right),\qquad
s=\frac{\tan(\Theta_0/2)}{\tan(\theta_0/2)}>0.
\]

记 `D=cos²(θ/2)+s²sin²(θ/2)`，则

\[
\Theta_\theta=\frac{s}{D}>0,\qquad
\Theta_{\theta\theta}=-\frac{s(s^2-1)\sin\theta}{2D^2}.
\]

普通区间 projective 对照为 `t=θ/π`、

\[
\Theta=\frac{\pi st}{1+(s-1)t},\quad
s=\frac{(\Theta_0/\pi)(1-\theta_0/\pi)}{(\theta_0/\pi)(1-\Theta_0/\pi)},
\]
\[
\Theta_\theta=\frac{s}{[1+(s-1)t]^2},\qquad
\Theta_{\theta\theta}=-\frac{2s(s-1)}{\pi[1+(s-1)t]^3}.
\]

两者均严格单调，固定 0、π，精确对应实测零角。它们的完整流函数候选仍保留 A、κ 两个独立幅度信息，绝不是只用零角解释形状。

## 3. 仅一个额外几何修正的对照

为检查端点斜率是否带来必要信息，另预定

\[
\varphi=M_s(\theta),\quad
\Theta=\varphi+\varepsilon\sin\varphi(\cos\varphi-\cos\Theta_0).
\]

它固定两端与 Θ₀。令

\[
w'=\cos2\varphi-\cos\Theta_0\cos\varphi,\quad
w''=-2\sin2\varphi+\cos\Theta_0\sin\varphi.
\]

则 `Θ_θ=(1+εw′)M_s′`、`Θ_θθ=εw″(M_s′)²+(1+εw′)M_s″`。固定 `|ε|≤0.45` 给出 `1+εw′≥0.1`，故单调性有明确保证。

ε 由**实测两端角导数之比**确定；A、κ 仍由两瓣幅度确定。没有通过最小化瓣形误差来选 ε，也没有添加任意角样条形变。若规定范围内无解，结果记为不可用。

## 4. 实测支持与限制

以下是中区 R=4…256、各半径等权的两瓣形状核查。τ20 的参考半径 R=32 不计入统计；τ10/15 的 R=32 仍是实际跨时刻检验。每个半径的零角与两瓣幅度都来自实测；迁移到 τ10/15 时只冻结 Φ*，不会把这些几何量误称为预测。

| 时刻 | 映射 | 正瓣 L2 差异 | 负瓣 L2 差异 |
|---:|---|---:|---:|''']
    for run in j['runs']:
        w=run['windows']['middle']
        lines.append(f"| {run['tau']} | {run['family']} | {w['positiveShapeRms']:.5f} | {w['negativeShapeRms']:.5f} |")
    lines += ['\nτ20 半角 Möbius 候选的窗口差别：\n',
        '| 窗口 | 原生半径 | 正瓣 L2 差异 | 负瓣 L2 差异 |',
        '|---|---|---:|---:|']
    primary=next(r for r in j['runs'] if r['tau']==20 and r['family']=='half_angle_mobius')
    for name in ('near_wall','transition','middle','outer'):
        w=primary['windows'][name];lo,hi=w['nativeRadiusRange']
        lines.append(f"| {name} | {lo:g}…{hi:g} | {w['positiveShapeRms']:.5f} | {w['negativeShapeRms']:.5f} |")
    lines += [r'''
半角 Möbius 加光滑幅度对比，在 τ20 中区对正瓣已相当有效；负瓣仍留有随半径变化的形状。普通 projective 在当前时刻接近，却在较早时刻更不稳定。端点斜率修正在 τ20 中区只要求很小 ε；较早时刻需要更大的 ε，且不能同时保证两瓣形状都改善。因此数据支持的是一个具体、平滑、低参数的组织方式，尚不支持“所有时刻共享完全固定母形”。

近壁两瓣差异显著大于中、远区，因此该中区母形不能直接提升为覆盖平滑原点的固定母形。中、远、近壁窗口分别保存在 JSON；这里没有为消除剩余差异继续加高阶形变。

## 5. 角导数可直接写出

记 `E=A exp(κ cosΘ)`，则

\[
\Psi_\theta=E(\Phi_*'-\kappa\sin\Theta\,\Phi_*)\Theta_\theta,
\]

\[
\begin{aligned}
\Psi_{\theta\theta}=E\big[&
(\Phi_*''-2\kappa\sin\Theta\,\Phi_*'
+(\kappa^2\sin^2\Theta-\kappa\cos\Theta)\Phi_*)\Theta_\theta^2\\
&+(\Phi_*'-\kappa\sin\Theta\,\Phi_*)\Theta_{\theta\theta}\big].
\end{aligned}
\]

因此该形式没有在零线上偷偷引入速度跳跃，角导数也能显式计算。不过 A、κ、s（以及可选 ε）的径向演化尚未由方程决定；原点正则性也需另行匹配。当前只是实测几何支持的候选公式，不能宣称已经降维闭合。

进一步壁面核查发现，基础三函数形式强制 `(Ψθθ/Ψθ)(0)·(Ψθθ/Ψθ)(π)` 与半径无关，实测并不满足。故它只能作为流函数值的组织方式，不能原样接受为速度/q骨架。独立报告 `ANGLE_MAP_WALL_ENVELOPE.md` 已比较加入 `ηsinΘ+ζsinΘcosΘ` 的五函数包络：它解除这个错误不变量、改善正瓣，但在部分窗口会恶化负瓣内部形状，仍需保留剩余结构。

![角映射与冻结母形的实测核查](/Users/cbjin/Documents/one_sided_IPM/runs/streamfunction_polar_structure_20260929/angle_map.png)

数据：angle_map.json/.npz。复现脚本：solver/steady_analysis/streamfunction_angle_map.py。未作 PDE 时间推进；固定母形的三次插值仅用于评价已有实测剖面，不是径向任意角形变。
''']
    (OUT/'ANGLE_MAP.md').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':main()
