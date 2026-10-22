"""Observe radial/lobe structure directly from saved centred psi.

No prescribed radial family, fixed rank, fitting, or PDE evolution.
"""
from pathlib import Path
import json,os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/ipm_polar_mpl')
import numpy as np
from scipy.interpolate import RectBivariateSpline,CubicSpline
from streamfunction_sources import load_native

ROOT=Path(__file__).resolve().parents[2]
OLD=ROOT/'runs/streamfunction_fresh_20260928'
OUT=ROOT/'runs/streamfunction_polar_observation_20260929'
WINDOWS={'near_wall':(.0625,.5),'transition':(.5,4.),'middle':(4.,256.),'outer':(256.,4096.)}
ANGLES=np.array([15.,30.,60.,90.,120.,150.])


def field(p,R,theta):
    RR,TT=np.meshgrid(np.asarray(R),np.asarray(theta),indexing='ij')
    x=RR*np.cos(TT)/p['ell'];y=RR*np.sin(TT)/p['ell']
    if x.min()<p['x'][0] or x.max()>p['x'][-1] or y.min()<0 or y.max()>p['y'][-1]:
        raise ValueError('No extrapolation')
    f=RectBivariateSpline(p['y'],p['x'],p['psi'])
    psi=p['P']*p['ell']**2*f.ev(y,x)
    dt=p['P']*p['ell']*RR*(np.cos(TT)*f.ev(y,x,dy=1)+np.sin(TT)*f.ev(y,x,dx=1))
    return psi,dt


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    metadata=json.loads((OLD/'polar_geometry.json').read_text())
    arrays=np.load(OLD/'polar_geometry.npz')
    R=arrays['nativeRadii'];theta=arrays['thetaRadians'];t=np.log(R)
    data={};frames={}
    for tau in (10,15,20):
        p=load_native(tau)
        psi,dt=field(p,R,theta)
        previous=arrays[f'wall_saddle_tau{tau}_psi']
        if np.max(abs(psi-previous))>1e-9*max(abs(psi).max(),1.):
            raise AssertionError('Saved polar matrix differs from source')
        positive=np.maximum(psi,0);negative=np.minimum(psi,0)
        A=np.sqrt(np.mean(psi**2,axis=1));Ap=np.sqrt(np.mean(positive**2,axis=1));An=np.sqrt(np.mean(negative**2,axis=1))
        pp=np.sum(positive*dt,axis=1)/np.sum(positive**2,axis=1)
        pn=np.sum(negative*dt,axis=1)/np.sum(negative**2,axis=1)
        degree=np.sum(psi*dt,axis=1)/np.sum(psi**2,axis=1)
        phi=psi/A[:,None];phi_t=(dt-degree[:,None]*psi)/A[:,None]
        ratio=Ap/An
        zero=np.array(metadata['centers']['wall_saddle']['radialDiagnostics'][str(tau)]['zeroAngleDegrees'])
        anglepsi,angledt=field(p,R,ANGLES*np.pi/180)
        valid=abs(anglepsi)>.025*A[:,None]
        angledegree=np.divide(angledt,anglepsi,out=np.full_like(anglepsi,np.nan),where=valid)
        rows=[]
        for rr in (.0625,.25,1,4,16,64,256,1024,4096):
            i=int(np.argmin(abs(R-rr)))
            rows.append(dict(nativeR=float(R[i]),positiveLobeAngularL2=float(Ap[i]),negativeLobeAngularL2=float(An[i]),
                positiveNegativeNormRatio=float(ratio[i]),degreePositive=float(pp[i]),degreeNegative=float(pn[i]),
                degreeWholeCircle=float(degree[i]),zeroAngleDegrees=float(zero[i]),
                positiveMaximumAngleDegrees=float(theta[np.argmax(psi[i])]*180/np.pi),
                negativeMinimumAngleDegrees=float(theta[np.argmin(psi[i])]*180/np.pi)))
        frames[str(tau)]=dict(source=p['source'],selectedRadii=rows,
            sampledRatioMinimum=dict(nativeR=float(R[np.argmin(ratio)]),value=float(ratio.min())),
            sampledZeroAngleMinimum=dict(nativeR=float(R[np.argmin(zero)]),degrees=float(zero.min())),
            sampledZeroAngleMaximum=dict(nativeR=float(R[np.argmax(zero)]),degrees=float(zero.max())),
            positiveDegreeRange=[float(pp.min()),float(pp.max())],negativeDegreeRange=[float(pn.min()),float(pn.max())])
        data[tau]=dict(psi=psi,dt=dt,A=A,Ap=Ap,An=An,pp=pp,pn=pn,p=degree,phi=phi,phi_t=phi_t,
                       ratio=ratio,zero=zero,anglepsi=anglepsi,angledegree=angledegree)
    d=data[20];window=[];directions=[]
    for name,(lo,hi) in WINDOWS.items():
        mask=(R>=lo)&(R<=hi);mat=d['phi'][mask];centered=mat-mat.mean(axis=0)
        _,sv,vv=np.linalg.svd(centered,full_matrices=False)
        direction=vv[0]*np.sign((mat[-1]-mat[0])@vv[0]);directions.append(direction)
        window.append(dict(name=name,nativeRadiusRange=[lo,hi],centeredAngularSingularStrengthRatios=(sv[:8]/sv[0]).tolist(),
            directionSign='Oriented toward phi at the final radius minus phi at the first radius.'))
    dirs=np.array(directions)
    directed=np.arccos(np.clip(dirs@dirs.T,-1,1))*180/np.pi
    axisangles=np.arccos(np.clip(abs(dirs@dirs.T),0,1))*180/np.pi
    centered=d['phi']-d['phi'].mean(axis=0)
    spectrum=np.linalg.svd(centered,compute_uv=False)
    # A chosen shared nodal angle has two radial occurrences: inspect the
    # actual source at those radii, not a fitted angular family.
    roots=CubicSpline(t,d['zero']-36.,extrapolate=False).roots(extrapolate=False)
    rootR=np.exp(roots)
    paired,paired_dt=field(load_native(20),rootR,theta)
    pairnorm=np.sqrt(np.mean(paired**2,axis=1))
    pairratio=np.sqrt(np.mean(np.maximum(paired,0)**2,axis=1)/np.mean(np.minimum(paired,0)**2,axis=1))
    repeated=dict(chosenZeroAngleDegrees=36.,nativeRadii=rootR.tolist(),positiveNegativeNormRatios=pairratio.tolist(),
        method='Cubic interpolation of measured zero angle versus log R locates the two radii; psi is then reevaluated directly from the saved source. No mother shape is fitted.',
        conclusion='The same nodal angle occurs with distinct lobe contrast; nodal angle alone cannot parameterize all measured angular shapes globally.')
    result=dict(title='Radial structure observed directly from centred streamfunction',
        sourceMatrices=str(OLD/'polar_geometry.npz'),nativeRadii=R.tolist(),
        protocol=dict(coordinates='Native centred half-circles X=a+R cos(theta),Y=R sin(theta),Psi=P ell^2 psi_inner.',
            lobeNorms='Aplus=sqrt(mean max(Psi,0)^2), Aminus=sqrt(mean min(Psi,0)^2), using the whole half-circle measure for both. Moving zero boundaries contribute no boundary term because Psi vanishes there.',
            localDegrees='pplus=<Psi_plus,R Psi_R>/<Psi_plus,Psi_plus>, analogously pminus; whole p is weighted by the two lobe energies and is not their equal average.',
            angleSamples=401,fixedAnglesDegrees=ANGLES.tolist(),pointwiseDegreeMask='Only display p_theta where abs(Psi(theta))>0.025 A(R); near a nodal angle the quotient is not interpreted.',
            angularDrift='Observe phi=Psi/A; subtract each radial-window mean before SVD to inspect variation, not the shared large background. No rank is selected and no fixed angular model is fitted.',
            noPdeEvolution=True,noPrescribedRadialFamily=True,noPrescribedRank=True),
        frames=frames,repeatedNodalAngle=repeated,
        windowDrift=dict(windows=window,directedAnglesDegrees=directed.tolist(),unorientedAxisAnglesDegrees=axisangles.tolist(),
            allRadiiCenteredAngularSingularStrengthRatios=(spectrum[:10]/spectrum[0]).tolist(),
            qualification='Different dominant window directions do not prove that no larger common angular subspace exists. A single radial trajectory is trivially a one-parameter curve; the nontrivial issues are monotone geometric parametrization, a common curve across times, and fixed versus rotating deformation directions.'),
        limitations=['The positive lobe is much smaller in L2 norm; whole-circle summaries largely describe the negative lobe.',
            'All turning radii are observations of a finite saved field, not exact constants or asymptotic exponents.',
            'The far positive lobe and nodal ray drift substantially across tau10/15/20; no universal far-field profile is inferred.',
            'A curved path in angular-function space can still lie on a nonlinear one-parameter family; these plots neither construct nor exclude such a family.',
            'No error-percentage fit contest is used to choose a symbolic form.'])
    (OUT/'radial_observation.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    save=dict(nativeRadii=R,theta=theta,directions=dirs,sameZeroAngleRadii=rootR,sameZeroAngleProfiles=paired/pairnorm[:,None])
    for tau,dd in data.items():
        for key,val in dd.items():save[f'tau{tau}_{key}']=val
    np.savez_compressed(OUT/'radial_observation.npz',**save)
    plot(R,theta,data,window,dirs,spectrum,rootR,paired/pairnorm[:,None])
    report(result)
    print(json.dumps(dict(tau20=frames['20'],repeatedNodalAngle=repeated,
        windowAxisAnglesDegrees=axisangles.tolist(),allRadiusStrengths=(spectrum[:5]/spectrum[0]).tolist()),indent=2))


def plot(R,theta,data,window,dirs,spectrum,pairedR,paired):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    colors={10:'#a7612b',15:'#40957e',20:'#315889'}
    d=data[20]
    fig,axs=plt.subplots(2,3,figsize=(15,8),constrained_layout=True)
    axs[0,0].loglog(R,d['Ap'],label='positive lobe',color='#b44f44')
    axs[0,0].loglog(R,d['An'],label='negative lobe',color='#315889')
    axs[0,0].set(title='Separate radial amplitudes, tau=20',xlabel='Native radius R',ylabel='Angular L2 norm');axs[0,0].legend()
    axs[0,1].semilogx(R,d['pp'],label='positive',color='#b44f44');axs[0,1].semilogx(R,d['pn'],label='negative',color='#315889')
    axs[0,1].semilogx(R,d['p'],':',color='black',label='whole half-circle')
    axs[0,1].set(title='Observed local radial degrees',xlabel='Native radius R',ylabel='d log amplitude / d log R');axs[0,1].legend()
    for key,pkey,color,label in [('Ap','pp','#b44f44','positive'),('An','pn','#315889','negative')]:
        axs[0,2].semilogx(d[key],d[pkey],color=color,label=label)
        for rr in (1,16,256,4096):
            i=np.argmin(abs(R-rr));axs[0,2].plot(d[key][i],d[pkey][i],'o',color=color,ms=3)
            axs[0,2].annotate(f'R={rr}',(d[key][i],d[pkey][i]),xytext=(3,3),textcoords='offset points',fontsize=7)
    axs[0,2].set(title='Degree versus its own lobe amplitude',xlabel='Angular L2 amplitude',ylabel='Local radial degree');axs[0,2].legend()
    for tau,dd in data.items():
        axs[1,0].semilogx(R,dd['ratio'],color=colors[tau],label=f'tau={tau}')
        axs[1,2].plot(dd['zero'],dd['ratio'],color=colors[tau],label=f'tau={tau}')
        for rr in (1,16,256,4096):
            i=np.argmin(abs(R-rr));axs[1,2].plot(dd['zero'][i],dd['ratio'][i],'o',color=colors[tau],ms=3)
            if tau==20:axs[1,2].annotate(f'R={rr}',(dd['zero'][i],dd['ratio'][i]),xytext=(3,3),textcoords='offset points',fontsize=7)
    axs[1,0].set(title='Positive / negative lobe contrast',xlabel='Native radius R',ylabel='Aplus / Aminus');axs[1,0].legend()
    for i,ang in enumerate(ANGLES):axs[1,1].semilogx(R,d['angledegree'][:,i],label=f'{ang:g} deg')
    axs[1,1].set(title='Fixed-direction radial degrees, tau=20',xlabel='Native radius R',ylabel='R Psi_R / Psi');axs[1,1].legend(ncol=2,fontsize=8)
    axs[1,2].set(title='Nodal angle and contrast need not move together',xlabel='Nodal angle (degrees)',ylabel='Aplus / Aminus');axs[1,2].legend()
    for ax in axs.flat:ax.grid(alpha=.15)
    fig.suptitle('Direct radial observations | no fitted radial family or chosen rank',fontsize=15)
    fig.savefig(OUT/'radial_amplitudes.png',dpi=165);plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(15,4.5),constrained_layout=True)
    for i,w in enumerate(window):
        axs[0].plot(theta*180/np.pi,dirs[i]*np.sqrt(len(theta)),label=w['name'])
        sv=w['centeredAngularSingularStrengthRatios'];axs[1].semilogy(np.arange(1,len(sv)+1),sv,'o-',label=w['name'])
    axs[0].set(title='Observed main drift direction in each window',xlabel='Angle (degrees)',ylabel='Unit angular norm direction');axs[0].legend(fontsize=8)
    axs[1].semilogy(np.arange(1,9),spectrum[:8]/spectrum[0],'ko--',label='all radii')
    axs[1].set(title='Spectrum of angular variation after centering',xlabel='Singular direction index',ylabel='Singular strength / first strength');axs[1].legend(fontsize=8)
    for R0,row in zip(pairedR,paired):axs[2].plot(theta*180/np.pi,row,label=f'R={R0:.3g}')
    axs[2].axvline(36,color='.6',ls=':',lw=1)
    axs[2].set(title='Same nodal angle, distinct angular shapes',xlabel='Angle (degrees)',ylabel='Psi / angular L2 norm');axs[2].legend()
    for ax in axs:ax.grid(alpha=.15)
    fig.savefig(OUT/'angular_drift_windows.png',dpi=165);plt.close(fig)


def report(j):
    selected=j['frames']['20']['selectedRadii'];w=j['windowDrift'];pair=j['repeatedNodalAngle']
    text=['# 先观察极坐标径向结构：不先选幂律或角基','',
        '本次只观察已保存的中心化流函数。半圆矩阵来自原有 polar_geometry.npz，径向导数重新从原始 ψ 样条求出并核对矩阵一致。没有拟合径向母式，没有先选 rank 2，也没有用拟合误差百分比筛选公式。','',
        '## 1. 整体幅度主要在描述负瓣','',
        '把正、负部分各自在整个半圆测度上取 L2 范数，记 A₊、A₋。它们不能简单平均：整体 p(R)=d log||ψ||/d log R 按两瓣能量加权。精确观测恒等式是 `p=(A₊²p₊+A₋²p₋)/(A₊²+A₋²)`、`d log(A₊/A₋)/d log R=p₊−p₋`。τ20 的 A₊/A₋ 仅在约 0.034 到 0.062 之间，因此整体指数几乎跟随负瓣，正瓣的径向变化容易被遮住。','',
        '| 原生 R | p₊ | p₋ | 整体 p | A₊/A₋ | 零支角度 |','|---:|---:|---:|---:|---:|---:|']
    for r in selected:
        text.append(f"| {r['nativeR']:g} | {r['degreePositive']:.4f} | {r['degreeNegative']:.4f} | {r['degreeWholeCircle']:.4f} | {r['positiveNegativeNormRatio']:.5f} | {r['zeroAngleDegrees']:.3f}° |")
    text += ['', '负瓣的局部指数持续下降；正瓣先下降、再上升、到所测最外层又下降。R=64…256 时正瓣指数在上升，负瓣指数仍在下降。不能用“一个逐渐趋于常数的全场指数”替代这两个不同过程。两瓣幅度本身在 τ20 仍增长，变化的是增长率。','',
        '## 2. 零支角和瓣强度的转折不同步','',
        'τ20 的零支角在采样 R≈6.73 附近最小，而正负瓣比在 R≈19.0 附近最小。这是两个不同的可见径向标志；不是同一个缩放量的重复测量。','',
        f"更直接地，零支角同为 36° 的两个半径约为 R={pair['nativeRadii'][0]:.5f} 与 R={pair['nativeRadii'][1]:.5f}；重新采样原始 ψ 后，两处 A₊/A₋ 分别为 {pair['positiveNegativeNormRatios'][0]:.6f} 与 {pair['positiveNegativeNormRatios'][1]:.6f}。因此零支角本身不是全径向范围内的唯一角形坐标。",'',
        '固定方向上的指数也随角度而变。图中直接列出多个角度的 Rψ_R/ψ；靠近零支时该商不稳定，故只在 |ψ| 大于整体角向幅度 0.025 倍处显示，不用零线附近的大商猜指数。','',
        '## 3. 每个窗口内的形状变化简单，但方向会转弯','',
        '先取 φ=ψ/||ψ||，再减去各径向窗口自身的平均角形，才观察“变化方向”。如果直接分解原矩阵，强大的共同负瓣背景会盖住这些方向差异。这里保留完整奇异强度序列，仅作为观察工具，不据此选定模型秩。','',
        '每个单独窗口有一个最强变化方向，第二方向相对强度约为 0.071、0.170、0.131、0.204；合并所有半径后，第二方向相对第一方向的强度升到 0.681。这与一条在角函数空间中明显弯曲的径向轨迹相符，不像沿同一个固定角向修正方向一直移动。','',
        f"中区与外区主变化轴夹角约 {w['unorientedAxisAnglesDegrees'][2][3]:.1f}°，近壁与外区约 {w['unorientedAxisAnglesDegrees'][0][3]:.1f}°。这是直接观察到的窗口转向；它不证明不存在更大的公共角空间，也不排除一个非线性的单参数角形族。",'',
        '任何单一时刻的径向路径都可以形式上用 log R 参数化，因此“画出一条曲线”本身不是降维证据。真正待检验的是：是否存在不依赖径向标签的几何参数，能覆盖两次零支角相同但瓣强度不同的情况；以及三个时刻是否落在同一关系上。当前相图有回折，各时刻轨迹也未合为一条显然固定的曲线。','',
        '## 4. 对下一步符号假设的要求','',
        '可以从这些事实筛掉过强的假设，而不急着指定新的基：','',
        '- 生成式应允许正、负瓣的径向增长率不同，并允许 p₊−p₋ 换号。',
        '- 零支位置与瓣强度不能未经证明就绑定成同一个单值函数。',
        '- 局部窗口的少量角向变化，不等于整个半平面共享一个固定有限角空间。',
        '- 若寻找非线性单参数族，应先说明该参数如何穿过相图回折，再研究它是否由方程选定。','',
        'τ10 的外层正瓣和零支还在明显变化，τ15/20 也未重合，所以不把当前 R=4096 处的指数或角度当作渐近常数。所有半径均处于保存域内；没有向无穷远外推。','',
        '![径向幅度与两瓣局部指数](radial_amplitudes.png)','',
        '![角向变化方向与窗口比较](angular_drift_windows.png)','',
        '详细观测见 radial_observation.json/.npz；复现脚本为 solver/steady_analysis/streamfunction_radial_observation.py。']
    (OUT/'RADIAL_OBSERVATION.md').write_text('\n'.join(text)+'\n')


if __name__=='__main__':main()
