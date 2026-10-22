"""One additional observed negative-peak constraint, no model search.

H=kappa*cosTheta+eta*sinTheta+zeta*sinTheta*cosTheta+xi*sinTheta**2.
The peak equation eliminates xi; the two lobe norms determine kappa and A.
"""
import os,json
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/ipm_polar_mpl')
import numpy as np
from scipy.interpolate import CubicSpline,BSpline
from scipy.optimize import brentq
from streamfunction_sources import load_native
from streamfunction_angle_map import WINDOWS

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/streamfunction_polar_structure_20260929'

def main():
    baseline=json.loads((OUT/'angle_map_wall_envelope.json').read_text())
    archive=np.load(OUT/'angle_map_wall_envelope.npz')
    arrays={k:archive[k] for k in archive.files}
    R,theta=arrays['nativeRadii'],arrays['theta']
    mother=BSpline(arrays['referenceMotherSplineKnots'],arrays['referenceMotherSplineCoefficients'],int(arrays['referenceMotherSplineDegree']))
    target=float(arrays['referenceZeroAngleRadians'])
    out=[]
    for tau in (10,15,20):
        old=next(r for r in baseline['runs'] if r['tau']==tau and r['functions']==5)
        p=load_native(tau);u=CubicSpline(p['x'],p['u1'][0]);q=CubicSpline(p['x'],p['q'][0])
        slope0=-R*p['P']*p['ell']*u(R/p['ell']);slopepi=R*p['P']*p['ell']*u(-R/p['ell'])
        sec0=-R*R*p['P']*q(R/p['ell']);secpi=-R*R*p['P']*q(-R/p['ell'])
        rows=[];predictions=[]
        for i,RR in enumerate(R):
            prior=old['rows'][i];s=prior['s'];eta=prior['eta'];zeta=prior['zeta']
            zero=2*np.arctan(np.tan(target/2)/s)
            obs=arrays[f'tau{tau}_sourcePsi'][i]
            cs=CubicSpline(np.r_[0,theta,np.pi],np.r_[0,obs,0],bc_type=((1,slope0[i]),(1,slopepi[i])))
            roots=cs.derivative().roots(extrapolate=False)
            minima=[float(t) for t in roots if zero<t<np.pi and cs(t)<0 and cs(t,2)>0]
            if not minima:raise ValueError(f'No observed negative minimum at tau={tau}, R={RR}')
            theta_minus=min(minima,key=lambda t:float(cs(t)))
            theta_mapped=2*np.arctan(s*np.tan(theta_minus/2))
            den=2*np.sin(theta_mapped)*np.cos(theta_mapped)
            if abs(den)<.05:raise ValueError('Predetermined peak elimination is ill-conditioned')
            th=2*np.arctan(s*np.tan(theta/2));S=np.sin(th);C=np.cos(th)
            aa=np.sin(theta_mapped)/den
            bb=(-eta*np.cos(theta_mapped)-zeta*np.cos(2*theta_mapped)-float(mother(theta_mapped,nu=1)/mother(theta_mapped)))/den
            pos=theta<zero;neg=~pos;pnorm=np.linalg.norm(obs[pos]);nnorm=np.linalg.norm(obs[neg])
            def evalk(k):
                xi=aa*k+bb
                return mother(th)*np.exp(k*C+eta*S+zeta*S*C+xi*S*S)
            def eq(k):
                w=evalk(k)
                return np.log(np.linalg.norm(w[pos])/np.linalg.norm(w[neg]))-np.log(pnorm/nnorm)
            grid=np.linspace(-16,16,129);vals=np.array([eq(k) for k in grid]);kroots=[]
            for j in range(len(grid)-1):
                if vals[j]==0:kroots.append(float(grid[j]))
                if vals[j]*vals[j+1]<0:kroots.append(float(brentq(eq,grid[j],grid[j+1],xtol=1e-12)))
            if vals[-1]==0:kroots.append(float(grid[-1]))
            if not kroots:raise ValueError(f'No kappa root in the predetermined interval at tau={tau}, R={RR}')
            k=min(kroots,key=lambda x:abs(x-prior['kappa']));xi=aa*k+bb
            w=evalk(k);A=float(nnorm/np.linalg.norm(w[neg]));pred=A*w
            def critical(t):
                hp=-k*np.sin(t)+eta*np.cos(t)+zeta*np.cos(2*t)+2*xi*np.sin(t)*np.cos(t)
                return mother(t,nu=1)+hp*mother(t)
            dense=np.linspace(target+1e-7,np.pi-1e-7,2401);dval=critical(dense);crit=[]
            for j in range(len(dense)-1):
                if dval[j]*dval[j+1]<0:crit.append(float(brentq(critical,dense[j],dense[j+1],xtol=1e-12)))
            wall0=A*np.exp(k)*float(mother(0,nu=1))*s
            wallpi=A*np.exp(-k)*float(mother(np.pi,nu=1))/s
            L0=prior['predictedL0'];Lpi=prior['predictedLpi']
            cancellation0=abs(prior['sourceL0']/s)/(abs(baseline['reference']['L0'])+2*abs(eta)+2*abs(zeta))
            cancellationpi=abs(s*prior['sourceLpi'])/(abs(baseline['reference']['Lpi'])+2*abs(eta)+2*abs(zeta))
            hp=-k*np.sin(theta_mapped)+eta*np.cos(theta_mapped)+zeta*np.cos(2*theta_mapped)+2*xi*np.sin(theta_mapped)*np.cos(theta_mapped)
            hpp=-k*np.cos(theta_mapped)-eta*np.sin(theta_mapped)-2*zeta*np.sin(2*theta_mapped)+2*xi*np.cos(2*theta_mapped)
            peak_second=float(mother(theta_mapped,nu=2)+2*hp*mother(theta_mapped,nu=1)+(hpp+hp*hp)*mother(theta_mapped))
            row=dict(nativeR=float(RR),A=A,s=s,kappa=float(k),eta=eta,zeta=zeta,xi=float(xi),
                positiveShapeRelativeL2=float(np.linalg.norm((pred-obs)[pos])/pnorm),
                negativeShapeRelativeL2=float(np.linalg.norm((pred-obs)[neg])/nnorm),
                positiveNormRatio=float(np.linalg.norm(pred[pos])/pnorm),negativeNormRatio=float(np.linalg.norm(pred[neg])/nnorm),
                sourceNegativePeakAngle=float(theta_minus),mappedNegativePeakAngle=float(theta_mapped),
                sourceNegativeMinimaCount=len(minima),predictedNegativeCriticalPointCount=len(crit),
                predictedNegativeCriticalAngles=crit,peakIsStrictMinimum=peak_second>0,
                negativePeakIsUniqueMinimum=len(crit)==1 and peak_second>0,
                normalizedPeakDerivative=float(critical(theta_mapped)),
                peakEliminationDenominator=float(den),kappaRootCount=len(kroots),
                wallSlope0RelativeError=float(wall0/slope0[i]-1),wallSlopePiRelativeError=float(wallpi/slopepi[i]-1),
                wallSecond0RelativeError=float(wall0*L0/sec0[i]-1),wallSecondPiRelativeError=float(wallpi*Lpi/secpi[i]-1),
                sourceWallQ0=float(-sec0[i]/RR**2),sourceWallQpi=float(-secpi[i]/RR**2),
                evaluatedWallQ0=float(-wall0*L0/RR**2),evaluatedWallQpi=float(-wallpi*Lpi/RR**2),
                wallRatio0CancellationIndicator=float(cancellation0),wallRatioPiCancellationIndicator=float(cancellationpi),
                wallQ0RelativeComparisonResolved=bool(cancellation0>1e-12),wallQpiRelativeComparisonResolved=bool(cancellationpi>1e-12),
                predictedL0=L0,predictedLpi=Lpi)
            rows.append(row);predictions.append(pred)
        windows={}
        for name,(lo,hi) in WINDOWS.items():
            rr=[r for r in rows if lo<=r['nativeR']<=hi and not(tau==20 and r['nativeR']==32)]
            windows[name]=dict(nativeRadiusRange=[lo,hi],count=len(rr),
                **{key+'Rms':float(np.sqrt(np.mean([r[key]**2 for r in rr]))) for key in
                    ('positiveShapeRelativeL2','negativeShapeRelativeL2','wallSlope0RelativeError','wallSlopePiRelativeError','wallSecond0RelativeError','wallSecondPiRelativeError')},
                xiRange=[min(r['xi'] for r in rr),max(r['xi'] for r in rr)],
                uniqueNegativePeakAll=all(r['negativePeakIsUniqueMinimum'] for r in rr),
                uniqueKappaRootAll=all(r['kappaRootCount']==1 for r in rr))
        out.append(dict(tau=tau,functions=6,rows=rows,windows=windows))
        arrays[f'tau{tau}_6function_prediction']=np.array(predictions)
        for key in ('A','s','kappa','eta','zeta','xi','sourceNegativePeakAngle'):
            arrays[f'tau{tau}_6function_{key}']=np.array([r[key] for r in rows])
    result=dict(candidate='Psi=A exp(kappa C+eta S+zeta SC+xi S^2)F(Theta); Theta=2atan(s tan(theta/2))',
        reference=baseline['reference'],
        calibration='s,eta,zeta remain as in the five-function model. Measured negative minimum angle eliminates xi linearly in kappa; lobe norm ratio fixes kappa, and one norm fixes A. No shape minimization.',
        scalarKappaSearch=dict(interval=[-16,16],gridPoints=129,selection='nearest five-function kappa if multiple roots; root count reported'),
        peakMeasurement='Cubic interpolation of observed angular psi, with endpoint slopes taken from saved u; deepest observed negative minimum selected. Critical points separately counted.',
        derivative=dict(Hprime='-kappa S+eta C+zeta cos(2Theta)+2xi SC',Hsecond='-kappa C-eta S-2zeta sin(2Theta)+2xi cos(2Theta)'),
        limitations=['This matches an observed negative peak angle, not an independently predicted location.',
            'Unchanged endpoint curvature/slope ratios do not force correct endpoint speed or q amplitudes.',
            'At the tau10 far left wall, source Lpi reaches 1e-20 and subtraction of order-one terms cannot resolve it in float64. Raw wall q relative errors are retained but flagged unresolved there; velocity amplitude remains a meaningful comparison.',
            'No interior q or PDE closure claim is made; full radial derivatives require a separate audit.',
            'A root-count check on a fine angular mesh is empirical, not an analytic proof of uniqueness.'],runs=out)
    (OUT/'angle_map_lobe_peak.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(OUT/'angle_map_lobe_peak.npz',**arrays)
    report(result,baseline);plot(result,baseline,R)
    for r in out:print(r['tau'],r['windows'])

def report(j,baseline):
    text=[r'''# 只增加负峰角约束的六函数候选

\[
\Psi=Ae^{\kappa C+\eta S+\zeta SC+\xi S^2}F(\Theta),\quad C=\cos\Theta,\ S=\sin\Theta.
\]

新增项及其一阶角导数在两端均为零，故保留五函数模型的壁面曲率/斜率比。以实测负峰 θ₋ 定 Θ₋，驻点条件给出

\[
\xi=\frac{\kappa\sin\Theta_- -\eta\cos\Theta_- -\zeta\cos2\Theta_- -F'(\Theta_-)/F(\Theta_-)}{2\sin\Theta_-\cos\Theta_-}.
\]

s、η、ζ 完全沿用前一模型的观测几何；代回两瓣范数比后只有 κ 的一个标量方程，再用任一瓣的范数定 A。没有优化内部形状，也没有增加更多模式。母形及其完整样条与五函数模型完全相同。

| τ | 窗口 | 函数数 | 正瓣差异 | 负瓣差异 | 右壁幅度差异 | 左壁幅度差异 |
|---:|---|---:|---:|---:|---:|---:|''']
    for r in j['runs']:
        old=next(x for x in baseline['runs'] if x['tau']==r['tau'] and x['functions']==5)
        for name,(lo,hi) in WINDOWS.items():
            if name=='transition':continue
            for run in (old,r):
                rr=[a for a in run['rows'] if lo<=a['nativeR']<=hi and not(run['tau']==20 and a['nativeR']==32)]
                vals=[np.sqrt(np.mean([a[k]**2 for a in rr])) for k in ('positiveShapeRelativeL2','negativeShapeRelativeL2','wallSlope0RelativeError','wallSlopePiRelativeError')]
                text.append(f"| {r['tau']} | {name} | {run['functions']} | "+' | '.join(f'{x:.6g}' for x in vals)+' |')
    unique=all(row['negativePeakIsUniqueMinimum'] for r in j['runs'] for row in r['rows'])
    roots=all(row['kappaRootCount']==1 for r in j['runs'] for row in r['rows'])
    text+=[f'\n全部样本负极值唯一且确为极小：{unique}；预定 κ 区间内根数全部为一：{roots}。角向极值通过细网格括根核查，不是解析唯一性证明。',r'''
壁面表中误差为速度相对幅度差异的径向 RMS。精确算术下，由于曲率/斜率比已经固定，同一侧壁面速度与 q 的相对幅度误差相同。τ10 最远左壁源 Lπ 达到 10⁻²⁰ 量级，但按 `F″/F′+2(−η+ζ)` 计算要做 O(1) 数相减，双精度已不能相对分辨；原始巨大 q 相对误差仍在 JSON 保留，并用 cancellation indicator 和 resolved 标志说明，不能把该处的浮点商当作模型的有效相对精度。每个时刻仍重测全部几何参数，仅母形冻结；参考 τ20、R=32 不计入形状统计。

τ20 的中区两瓣形状均改善；远区负瓣改善但正瓣略退。中、远区两壁速度幅度更接近源值，近壁负侧仍有明显偏差。跨时刻的负瓣通常改善，正瓣却不全都改善；负峰位置确实提供了额外信息，但一个新增几何量仍没有抹平所有差异。

这一步不能凭更小的流函数误差认定 q 保正。径向参数导数和完整 Poisson 算子需另作留出审核，尤其已有五函数候选在远区内角产生负 q 的失败。

![五函数与六函数比较](/Users/cbjin/Documents/one_sided_IPM/runs/streamfunction_polar_structure_20260929/angle_map_lobe_peak.png)

NPZ 保存六个参数、预测矩阵、原始矩阵、固定母形样条以及先前 3/5 函数资料，键名为 `tau20_6function_A/s/kappa/eta/zeta/xi`，各时刻同理。
''']
    (OUT/'ANGLE_MAP_LOBE_PEAK.md').write_text('\n'.join(text)+'\n')

def plot(j,baseline,R):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,3,figsize=(15,8),constrained_layout=True)
    for col,tau in enumerate((20,15,10)):
        for count,color in [(5,'#b57542'),(6,'#376f80')]:
            run=next(r for r in (baseline['runs'] if count==5 else j['runs']) if r['tau']==tau and r['functions']==count)
            for key,style,name in [('positiveShapeRelativeL2','-','positive'),('negativeShapeRelativeL2','--','negative')]:
                err=np.array([r[key] for r in run['rows']]);err[(R==32)&(tau==20)]=np.nan
                axs[0,col].loglog(R,err,style,color=color,label=f'{count}, {name}')
            for key,style,name in [('wallSlope0RelativeError','-','right wall'),('wallSlopePiRelativeError','--','left wall')]:
                axs[1,col].semilogx(R,[r[key] for r in run['rows']],style,color=color,label=f'{count}, {name}')
        axs[0,col].set(title=f'tau={tau}: lobe shape',xlabel='Native R',ylabel='Relative lobe L2 discrepancy')
        axs[1,col].set(title=f'tau={tau}: wall velocity amplitude',xlabel='Native R',ylabel='Predicted slope / source slope - 1')
    for ax in axs.flat:ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.savefig(OUT/'angle_map_lobe_peak.png',dpi=165);plt.close(fig)

if __name__=='__main__':main()
