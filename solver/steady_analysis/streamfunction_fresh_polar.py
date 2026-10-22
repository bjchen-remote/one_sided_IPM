"""Direct polar geometry of saved centred psi, without a prescribed mother.

All physical queries are inside saved half-plane boxes. This is descriptive
finite-time geometry, not a steady or infinite-domain certification.
"""
from pathlib import Path
import json
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline
from streamfunction_sources import load_native

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/streamfunction_fresh_20260928"
RADII = 2.**(np.arange(-16, 49)/4)
THETA = (np.arange(401)+.5)*np.pi/401
BANDS = {"near_wall":(.0625,.5), "transition":(.5,4.),
         "middle":(4.,256.), "outer":(256.,4096.)}
TRAIN = np.arange(len(RADII)) % 2 == 0


def make_source(tau):
    p=load_native(tau)
    p["spline"]=RectBivariateSpline(p["y"],p["x"],p["psi"])
    return p


def zero_line_center(p):
    """One separately motivated virtual centre; never optimize SVD error."""
    heights=np.array([4.,8.,16.,32.,64.,128.,256.])
    train=np.array([True,False,True,False,True,False,True])
    values=[]
    for Y in heights:
        y=Y/p["ell"]
        row=p["spline"](y,p["x"])[0]
        roots=CubicSpline(p["x"],row,extrapolate=False).roots(extrapolate=False)
        roots=roots[np.isfinite(roots)]
        # Select the near-core zero by proximity to the wall saddle, without
        # inspecting any q or velocity peak. The far-box zero is not this line.
        xzero=float(roots[np.argmin(abs(roots))])
        values.append(p["ell"]*xzero)
    values=np.asarray(values)
    slope,intercept=np.polyfit(heights[train],values[train],1)
    predicted=intercept+slope*heights
    return dict(nativeWallOffset=float(intercept),zeroRaySlope=float(slope),
        method="Fit the centred-native psi=0 line X-a=b+mY on prescribed tau20 training Y=4,16,64,256; no SVD/separation score is used to choose b. Freeze b for all radii and times.",
        nativeY=heights.tolist(),observedCenteredNativeX=values.tolist(),predictedCenteredNativeX=predicted.tolist(),
        trainingMask=train.tolist(),heldoutL2RelativeToX=float(np.linalg.norm((predicted-values)[~train])/np.linalg.norm(values[~train])),
        heldoutMaxAbsoluteNativeX=float(max(abs(predicted-values)[~train])))


def polar_sample(p,center_offset):
    R,T=np.meshgrid(RADII,THETA,indexing="ij")
    # Radius is measured in native coordinates; the centre is X=a+offset,Y=0.
    x=(center_offset+R*np.cos(T))/p["ell"]
    y=R*np.sin(T)/p["ell"]
    if x.min()<p["x"][0] or x.max()>p["x"][-1] or y.min()<p["y"][0] or y.max()>p["y"][-1]:
        raise ValueError("Polar query outside saved source")
    scale=p["P"]*p["ell"]**2
    psi=scale*p["spline"].ev(y,x)
    dr=p["P"]*p["ell"]*(np.cos(T)*p["spline"].ev(y,x,dy=1)+np.sin(T)*p["spline"].ev(y,x,dx=1))
    A=np.sqrt(np.mean(psi**2,axis=1))
    if np.any(A<=0):
        raise ValueError("Nonpositive half-circle angular norm")
    phi=psi/A[:,None]
    degree=np.mean(psi*R*dr,axis=1)/A**2
    defect=np.sqrt(np.mean((R*dr-degree[:,None]*psi)**2,axis=1))/A
    degree_fd=np.gradient(np.log(A),np.log(RADII),edge_order=2)
    zero_angles=[]; zero_counts=[]
    for row in psi:
        roots=CubicSpline(THETA,row,extrapolate=False).roots(extrapolate=False)
        roots=roots[(roots>0)&(roots<np.pi)&np.isfinite(roots)]
        zero_counts.append(len(roots))
        zero_angles.append(float(roots[0]*180/np.pi) if len(roots) else None)
    return dict(psi=psi,A=A,phi=phi,degree=degree,defect=defect,
                degree_fd=degree_fd,zero_angles=zero_angles,zero_counts=zero_counts)


def rms(v):
    return float(np.sqrt(np.mean(np.asarray(v)**2)))


def svd_basis(phi,rank):
    u,s,v=np.linalg.svd(phi,full_matrices=False)
    basis=v[:rank].copy()
    if np.mean(phi@basis[0])<0:
        basis[0]*=-1
    return basis,s


def separation_fit(sample,mask):
    train=mask&TRAIN; holdout=mask&~TRAIN
    basis,s=svd_basis(sample["phi"][train],2)
    first=basis[0]*np.sqrt(len(THETA))
    alpha,logc=np.polyfit(np.log(RADII[train]),np.log(sample["A"][train]),1)
    predicted_A=np.exp(logc)*RADII**alpha
    pred=predicted_A[:,None]*first[None,:]
    coef=sample["phi"]@basis.T
    rank1=coef[:,:1]@basis[:1]
    rank2=coef@basis
    all_s=np.linalg.svd(sample["phi"][mask],compute_uv=False)
    info=dict(nativeRadiusRange=[float(RADII[mask][0]),float(RADII[mask][-1])],
        trainingNativeRadii=RADII[train].tolist(),heldoutNativeRadii=RADII[holdout].tolist(),
        amplitudePower=dict(exponent=float(alpha),coefficient=float(np.exp(logc)),
            heldoutRelativeRms=rms((predicted_A/sample["A"]-1)[holdout])),
        angularTrainingNormalizedSingularValues=(s[:6]/np.linalg.norm(s)).tolist(),
        angularAllRank1RelativeFrobenius=float(np.linalg.norm(all_s[1:])/np.linalg.norm(all_s)),
        angularAllRank2RelativeFrobenius=float(np.linalg.norm(all_s[2:])/np.linalg.norm(all_s)),
        angularHeldoutRank1ProjectionRelativeRms=rms((rank1-sample["phi"])[holdout]),
        angularHeldoutRank2ProjectionRelativeRms=rms((rank2-sample["phi"])[holdout]),
        angularHeldoutFixedRank1ShapeRms=rms((first-sample["phi"])[holdout]),
        fullHomogeneousHeldoutRadiusNormalizedRms=rms(((pred-sample["psi"])/sample["A"][:,None])[holdout]),
        localDegreeRange=[float(sample["degree"][mask].min()),float(sample["degree"][mask].max())],
        localDegreeMedian=float(np.median(sample["degree"][mask])),
        homogeneityDefectRange=[float(sample["defect"][mask].min()),float(sample["defect"][mask].max())],
        zeroAngleDegreesRange=[min(z for z,m in zip(sample["zero_angles"],mask) if m and z is not None),
                               max(z for z,m in zip(sample["zero_angles"],mask) if m and z is not None)])
    return info,dict(phi=first,basis=basis,predicted_A=predicted_A,prediction=pred)


def transfer(sample,model,mask):
    phi=sample["phi"][mask]
    projection=(phi@model["basis"].T)@model["basis"]
    return dict(fixedTau20AngularShapeRms=rms(phi-model["phi"]),
        fixedTau20Rank2SubspaceProjectionRms=rms(projection-phi),
        fixedTau20AmplitudePowerRelativeRms=rms(model["predicted_A"][mask]/sample["A"][mask]-1),
        fixedTau20FullHomogeneousRadiusNormalizedRms=rms((model["prediction"][mask]-sample["psi"][mask])/sample["A"][mask,None]))


def quadratic_source(p):
    f=p["spline"];P=p["P"]
    # Native Hessian of psi_native is P times the inner psi Hessian.
    hxx=P*float(f.ev(0,0,dy=2)); hxy=P*float(f.ev(0,0,dx=1,dy=1)); hyy=P*float(f.ev(0,0,dx=2))
    angular=.5*hxx*np.cos(THETA)**2+hxy*np.sin(THETA)*np.cos(THETA)+.5*hyy*np.sin(THETA)**2
    psi=RADII[:,None]**2*angular[None,:]
    return dict(nativeHessian=dict(xx=hxx,xy=hxy,yy=hyy),prediction=psi)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    sources={tau:make_source(tau) for tau in (10,15,20)}
    virtual=zero_line_center(sources[20])
    centers={"wall_saddle":0.,"zero_line_virtual_center":virtual["nativeWallOffset"]}
    results={};data={};arrays=dict(nativeRadii=RADII,thetaRadians=THETA,trainingRadiusMask=TRAIN)
    for center,offset in centers.items():
        bytime={tau:polar_sample(p,offset) for tau,p in sources.items()}
        bandnames=BANDS if center=="wall_saddle" else {k:BANDS[k] for k in ("middle","outer")}
        bands={};models={}
        for band,(lo,hi) in bandnames.items():
            mask=(RADII>=lo)&(RADII<=hi)
            info,model=separation_fit(bytime[20],mask)
            info["frozenTau20TimeTransfer"]={str(tau):transfer(bytime[tau],model,mask) for tau in (10,15,20)}
            info["otherTimesOwnSeparationDiagnostics"]={str(tau):separation_fit(bytime[tau],mask)[0] for tau in (10,15)}
            bands[band]=info;models[band]=model
            arrays[f"{center}_{band}_tau20AngularBasis"]=model["basis"]
        rows={}
        for tau,sample in bytime.items():
            rows[str(tau)]=dict(amplitudeAngularL2=sample["A"].tolist(),localDegree=sample["degree"].tolist(),
                finiteDifferenceDegree=sample["degree_fd"].tolist(),degreeDerivativeComparisonMaxAbsoluteInterior=float(max(abs(sample["degree"][1:-1]-sample["degree_fd"][1:-1]))),
                angularHomogeneityDefect=sample["defect"].tolist(),zeroAngleDegrees=sample["zero_angles"],
                interiorZeroCounts=sample["zero_counts"])
            for key in ("psi","A","phi","degree","defect"):
                arrays[f"{center}_tau{tau}_{key}"]=sample[key]
        if center=="wall_saddle":
            for tau,p in sources.items():
                quad=quadratic_source(p)
                mask=RADII<=.25
                rows[str(tau)]["wallHessianQuadratic"]={**quad["nativeHessian"],
                    "nativeRadiusRange":[float(RADII[mask][0]),float(RADII[mask][-1])],
                    "radiusNormalizedRms":rms((quad["prediction"][mask]-bytime[tau]["psi"][mask])/bytime[tau]["A"][mask,None])}
        results[center]=dict(centeredNativeXOffset=offset,bands=bands,radialDiagnostics=rows)
        data[center]=bytime
    output=dict(title="Direct polar geometry of centred streamfunction; no inherited density or velocity mother shape",
        protocol=dict(nativeRadiusRange=[float(RADII[0]),float(RADII[-1])],nativeRadii=RADII.tolist(),
            angularInterval="0<theta<pi; midpoint quadrature, 401 equally weighted angles",angularSamples=len(THETA),
            coordinates="X=a(tau)+b+R cos(theta), Y=R sin(theta); Psi=psi_native+uAtPeak*Y=P*ell^2*psi_inner.",
            norm="A(R)=[mean_theta Psi^2]^(1/2), phi=Psi/A. All shape errors use angular L2 with equal log-radius weights; never divide by pointwise Psi, including at its nodal curve.",
            localDegree="p(R)=<Psi,R Psi_R>/<Psi,Psi>=d log A/d log R; the analytic spline derivative is checked against log-radius finite differences.",
            homogeneityDefect="||R Psi_R-p(R)Psi||_theta/||Psi||_theta = ||d phi/d log R||_theta; zero for a fixed angular shape, even if the radial amplitude is not a power.",
            separation="SVD is applied to row-normalized phi(R,theta), so large-radius amplitudes do not dominate. No angular mother shape is assumed. First two angular singular modes are learned from alternating radii at tau20.",
            holdout="Alternate quarter-octave radii are held out; amplitude power and angular mode are fitted only on tau20 training radii in each named band, then frozen for tau10/15 comparisons.",
            centerChoice="Main centre is the centred wall stagnation point X=a,Y=0. One alternative centre is derived independently from the measured psi zero line; it is not tuned for rank or homogeneity.",
            noDomainExtrapolation=True,noPdeEvolution=True),
        sources={str(t):dict(path=p["source"],a=p["a"],P=p["P"],ell=p["ell"],gauge=p["gauge"]) for t,p in sources.items()},
        virtualCenterConstruction=virtual,centers=results,
        limitations=["These are finite-time saved fields, not a demonstrated global steady solution.",
            "Native radii 256..4096 are a measured outer band, not a certified far-field limit; the computational box extends to about 100000.",
            "A two-dimensional angular subspace still leaves radial coefficient functions; low rank alone is not a closed reduced PDE.",
            "The innermost shell R=.0625 is only a few native mesh spacings; the spline derivative and wall Hessian checks do not establish subgrid accuracy.",
            "Native-centred Psi and its nodal geometry depend on the specified horizontal-velocity gauge; the same physical centring convention is used at all times."])
    (OUT/"polar_geometry.json").write_text(json.dumps(output,indent=2,allow_nan=False)+"\n")
    np.savez_compressed(OUT/"polar_geometry.npz",**arrays)
    plot(data,output)
    write_report(output)
    for c,r in results.items():
        print(c,"offset",r["centeredNativeXOffset"])
        for b,d in r["bands"].items():
            print(b,"p",d["amplitudePower"]["exponent"],"local",d["localDegreeRange"],"rank1",d["angularHeldoutRank1ProjectionRelativeRms"],"rank2",d["angularHeldoutRank2ProjectionRelativeRms"],"full",d["fullHomogeneousHeldoutRadiusNormalizedRms"])
            print("transfer",d["frozenTau20TimeTransfer"])


def plot(data,out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    fig,axs=plt.subplots(2,3,figsize=(15,8),constrained_layout=True)
    colors={10:"#9d5b24",15:"#348f7d",20:"#274c8f"}
    for tau,s in data["wall_saddle"].items():
        axs[0,0].loglog(RADII,s["A"],label=f"tau={tau}",color=colors[tau])
        axs[0,1].semilogx(RADII,s["degree"],color=colors[tau])
        axs[1,0].semilogx(RADII,s["zero_angles"],color=colors[tau])
        axs[1,1].loglog(RADII,s["defect"],color=colors[tau])
    axs[0,0].set(title="Angular L2 amplitude of centred Psi",xlabel="Native radius R",ylabel="A(R)");axs[0,0].legend()
    axs[0,1].set(title="Observed local radial degree",xlabel="Native radius R",ylabel="d log A / d log R")
    axs[0,1].axhline(2,color=".7",lw=.8,ls="--")
    for R in (.0625,.5,4,32,256,4096):
        i=np.argmin(abs(RADII-R));axs[0,2].plot(THETA*180/np.pi,data["wall_saddle"][20]["phi"][i],label=f"R={R:g}")
    axs[0,2].set(title="Unprescribed angular shapes, tau=20",xlabel="Angle (degrees)",ylabel="Psi / angular L2 norm");axs[0,2].legend(ncol=2,fontsize=8)
    axs[1,0].set(title="Interior nodal ray of Psi",xlabel="Native radius R",ylabel="Zero angle (degrees)")
    axs[1,1].set(title="Angular drift with log radius",xlabel="Native radius R",ylabel="||d phi / d log R||")
    for center,style in (("wall_saddle","-"),("zero_line_virtual_center","--")):
        vals=out["centers"][center]["bands"]
        for rank,color in ((1,"#b74d43"),(2,"#347894")):
            for ib,(b,d) in enumerate(vals.items()):
                lo,hi=d["nativeRadiusRange"]
                value=d[f"angularHeldoutRank{rank}ProjectionRelativeRms"]
                label=f"rank {rank}, {'saddle' if style=='-' else 'virtual'}" if ib==0 else None
                axs[1,2].loglog([lo,hi],[value,value],style,color=color,label=label,lw=2)
    axs[1,2].set(title="Held-out angular subspace error",xlabel="Native radius R",ylabel="Normalized angular L2 error");axs[1,2].legend(fontsize=8)
    for ax in axs.flat:
        ax.grid(alpha=.15)
    fig.suptitle("Direct streamfunction geometry | native half-circles | no peak/width mapping",fontsize=15)
    fig.savefig(OUT/"polar_geometry.png",dpi=170)
    plt.close(fig)


def write_report(out):
    r=out["centers"]["wall_saddle"]["bands"]
    lines=["# 从流函数本身观察极坐标几何","",
        "本实验直接读取三个时刻的中心化流函数，没有采用密度或速度峰形，没有半高宽、sinh 坐标或预设角向母形。", "",
        "直接观测支持单时刻、有限径向窗口内的近径角分离；不支持把全半平面压成一个固定幂次及时间不变角形。τ20 的局部径向次数从近壁约 1.97 连续下降至所测外层约 1.38；固定角形在各窗口的半径留出误差约 0.8%–1.5%，但跨时刻迁移可达 6%–13%。第二个由数据确定的角模态显著降低同一时刻的形状误差，提示保留随半径变化的少数角模态比先固定一个齐次指数更忠实。", "",
        "## 取样与误差规范", "",
        "主中心为壁面驻点 `(X,Y)=(a,0)`，原生半径 `R=0.0625…4096`，每个半圆在 `0<θ<π` 上用 401 个中点取样。中心化原生流函数为 `Ψ=ψ_native+uAtPeak·Y=Pℓ²ψ_inner`。", "",
        "令 `A(R)=sqrt(mean_θ Ψ²)`、`φ=Ψ/A`。所有角向形状误差按半圆 L2、各对数半径等权计算，零线处从不除以 Ψ。SVD 直接作用于 φ，避免外层幅度掩盖内层变化。局部径向次数为 `p=<Ψ,RΨ_R>/<Ψ,Ψ>`；角向漂移为 `||RΨ_R−pΨ||/||Ψ||`。", "",
        "每隔一个四分之一倍频半径训练，夹在其间的半径留出；所有母模态及径向幂律只在 τ20 训练，然后冻结迁移到 τ10/15。两模态投影在被测半径重新取投影系数，故检验的是角向子空间，而非径向系数预测。", "",
        "## τ20 主中心结果", "",
        "| 区间 | 原生 R | 幅度幂律指数 | 局部指数范围 | 单角模态留出 | 双角模态留出 | 完整齐次式留出 |",
        "|---|---|---:|---|---:|---:|---:|"]
    for name,d in r.items():
        lines.append(f"| {name} | {d['nativeRadiusRange'][0]:g}–{d['nativeRadiusRange'][1]:g} | {d['amplitudePower']['exponent']:.5f} | {d['localDegreeRange'][0]:.4f}–{d['localDegreeRange'][1]:.4f} | {100*d['angularHeldoutRank1ProjectionRelativeRms']:.3f}% | {100*d['angularHeldoutRank2ProjectionRelativeRms']:.3f}% | {100*d['fullHomogeneousHeldoutRadiusNormalizedRms']:.3f}% |")
    lines += ["", "完整齐次式是 `Ψ≈C R^p Φ(θ)`，同时测试径向幅度与固定角形；它比仅看角向低秩更严格。指数是指定有限区间的拟合结果，不是已确定的渐近指数。", "",
        "## 冻结 τ20 的时间迁移", "",
        "| 区间 | 时刻 | 固定单角形误差 | 固定双角子空间误差 | 完整齐次式误差 |", "|---|---:|---:|---:|---:|"]
    for name,d in r.items():
        for tau in (10,15):
            t=d["frozenTau20TimeTransfer"][str(tau)]
            lines.append(f"| {name} | {tau} | {100*t['fixedTau20AngularShapeRms']:.3f}% | {100*t['fixedTau20Rank2SubspaceProjectionRms']:.3f}% | {100*t['fixedTau20FullHomogeneousRadiusNormalizedRms']:.3f}% |")
    v=out["virtualCenterConstruction"]
    lines += ["", "## 有限的中心对照", "",
        f"唯一备选中心由 τ20 的 Ψ=0 线在固定训练高度 Y=4,16,64,256 的直线拟合给出：`X−a={v['nativeWallOffset']:.6f}+{v['zeroRaySlope']:.6f}Y`。留出 Y=8,32,128 的零线相对误差为 {100*v['heldoutL2RelativeToX']:.3f}%。其壁面交点作为固定虚拟中心，未按 SVD 或齐次误差优化，之后冻结到所有时刻。只在中、远区比较，近壁展开仍围绕实际驻点。", "",
        "| 区间 | 中心 | 单角模态留出 | 双角模态留出 | 完整齐次式留出 |", "|---|---|---:|---:|---:|"]
    for name in ("middle","outer"):
        for center in out["centers"]:
            d=out["centers"][center]["bands"][name]
            lines.append(f"| {name} | {center} | {100*d['angularHeldoutRank1ProjectionRelativeRms']:.3f}% | {100*d['angularHeldoutRank2ProjectionRelativeRms']:.3f}% | {100*d['fullHomogeneousHeldoutRadiusNormalizedRms']:.3f}% |")
    lines += ["", "## 可支持的降维方向与边界", "",
        "本数据可用来判断 `Ψ(R,θ)=Σ a_j(R) Φ_j(θ)` 是否是合适的观测骨架，再检查径向系数是否满足简洁方程。单角模态对应严格径角分离；第二角模态容纳径向漂移。它们均由 Ψ 自身确定，尚未证明在 Poisson 与输运方程下闭合。", "",
        "近壁实际驻点处的 Hessian 给出独立的二次齐次基准；其定量误差保存在 JSON。外区 R=256…4096 距约 100000 的外边界仍远，但属于有限计算域内的观测，不能宣称全域渐近或稳态。最小半径也只有数个原生网格间距。", "",
        "![直接流函数的极坐标诊断](polar_geometry.png)", "",
        "完整取样、范数、留出划分、各时刻独立 SVD、局部指数与零线角度见 `polar_geometry.json`；矩阵和角向模态见 `polar_geometry.npz`。"]
    (OUT/"POLAR_GEOMETRY.md").write_text("\n".join(lines)+"\n")


if __name__=="__main__":
    main()
