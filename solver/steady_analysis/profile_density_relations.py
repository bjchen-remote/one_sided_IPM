"""Density/gradient geometry on saved frames, with explicit time/range holdouts.

Only reads saved fields.  No PDE evolution, closure, or tail extrapolation.
Run with PYTHONPATH=/private/tmp/ipm_steady_deps and at most two BLAS threads.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_simpson, quad
from scipy.interpolate import CubicSpline, RectBivariateSpline
from scipy.io import loadmat
from scipy.optimize import brentq

from geometry_measure import BETA, interp_fields, peak_measure, trajectory_frame

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/profile_geometry_20260928"
NATIVE_HEIGHTS = np.array([.25, .5, 1., 2., 4., 8., 16., 32., 64., 128., 256.])
INNER_HEIGHTS = np.array([0., .25, .5, .75, 1., 1.25, 1.5, 1.75, 2.])
DENSITY_LEVELS = np.linspace(-.4, .75, 47)
GRADIENT_LEVELS = [.5, .75, .9]
XI = np.linspace(-1.5, 1.5, 301)


def f0(x):
    return 32. / (x**6 - 3*x**5 + 3*x**4 + 3*x**3 + 28*x**2 + 32)


def G(xi, k):
    return (np.sinh(BETA*xi) + k*(np.cosh(BETA*xi)-1)) / BETA


def primitive_model():
    # A core-only primitive: never use the rational f0 outside [-1.5,1.5].
    xi = np.linspace(-1.5, 1.5, 6001)
    iz = len(xi)//2
    a = cumulative_simpson(f0(xi)*np.cosh(BETA*xi), x=xi, initial=0)
    b = cumulative_simpson(f0(xi)*np.sinh(BETA*xi), x=xi, initial=0)
    return CubicSpline(xi, a-a[iz]), CubicSpline(xi, b-b[iz])


H0, H1 = primitive_model()


def model_at_density(sigma, k):
    fn = lambda xi: float(H0(xi) + k*H1(xi) - sigma)
    if fn(-1.5)*fn(1.5) > 0:
        return None
    xi = brentq(fn, -1.5, 1.5, xtol=1e-13)
    return [float(G(xi, k)), float(1/f0(xi)), float(xi)]


def rms(x):
    x = np.asarray(x)
    return float(np.sqrt(np.mean(x*x)))


def relative_rms(delta, target):
    return rms(delta)/rms(target) if rms(target) else None


def measure_row(p, qf, Qf, coarse, y, track, height, offset):
    assert p["y"][0] <= y <= p["y"][-1]
    r = peak_measure(p, float(y), qf, Qf)
    m, A, c = r["m"], r["A"], r["c"]
    W, k = r["halfWidthGeometry"]["W"], r["halfWidthGeometry"]["skew"]
    assert W > 0 and abs(k) < 1, (track, height, W, k)
    xcore = m + W*G(XI, k)
    assert xcore[0] >= p["x"][0] and xcore[-1] <= p["x"][-1]
    qcore = qf.ev(np.full(len(XI), y), xcore)/A
    hcore = (Qf.ev(np.full(len(XI), y), xcore)-c)/(A*W)
    assert np.min(np.diff(hcore)) > 0, (track, height, "density nonmonotone")
    # Qf has its first spline coordinate y, second x: dy=1 differentiates x.
    Qxcore = Qf.ev(np.full(len(XI), y), xcore, dy=1)/A
    qc, Qc = coarse
    qc_core = qc.ev(np.full(len(XI), y), xcore)/A
    hc_core = (Qc.ev(np.full(len(XI), y), xcore)-Qc(y, m)[0, 0])/(A*W)
    hmodel = H0(XI) + k*H1(XI)

    def inverse(sigma):
        assert hcore[0] <= sigma <= hcore[-1], (track, height, sigma, hcore[[0,-1]])
        xx = brentq(lambda x: float((Qf(y, x)[0, 0]-c)/(A*W)-sigma),
                    float(xcore[0]), float(xcore[-1]), xtol=1e-11)
        return (xx-m)/W, A/float(qf(y, xx)[0, 0])

    Z, J = np.array([inverse(s) for s in DENSITY_LEVELS]).T
    Zmodel, Jmodel, _ = np.array([model_at_density(s, k) for s in DENSITY_LEVELS]).T
    csq = CubicSpline(p["x"], qf(y, p["x"])[0])
    gradient = []
    for level in GRADIENT_LEVELS:
        roots = CubicSpline(p["x"], csq(p["x"])-A*level).roots(extrapolate=False)
        lr, rr = roots[roots < m], roots[roots > m]
        for side, xx in [("left", float(lr[-1])), ("right", float(rr[0]))]:
            assert xcore[0] <= xx <= xcore[-1]
            xi0 = brentq(lambda t: float(f0(t)-level), -1.5, 0.) if side == "left" else brentq(lambda t: float(f0(t)-level), 0., 1.5)
            observed = float((Qf(y, xx)[0, 0]-c)/(A*W))
            predicted = float(H0(xi0)+k*H1(xi0))
            gradient.append(dict(level=level, side=side, densityInterval=observed,
                                 f0DensityInterval=predicted, error=predicted-observed))

    sigma = c/(A*W)
    zero = None
    if hcore[0] <= -sigma <= hcore[-1]:
        z0, j0 = inverse(-sigma)
        pred = model_at_density(-sigma, k)
        zero = dict(xInner=float(m+W*z0), offsetFromPeakOverW=float(-z0),
                    qZeroOverA=float(1/j0), linearOffsetRatio=float(-z0/sigma) if abs(sigma)>1e-10 else None,
                    modelOffsetFromPeakOverW=float(-pred[0]), modelQZeroOverA=float(1/pred[1]),
                    xErrorOverW=float(pred[0]-z0), qErrorOverA=float(1/pred[1]-1/j0))
    return dict(track=track, height=float(height), yInner=float(y), m=m, A=A, W=W, k=k,
                cPeakBaseline=c, cOriginBaseline=float(c+offset), cOverAW=float(sigma),
                cOriginOverAW=float((c+offset)/(A*W)),
                fourCPeakOverYA=float(4*c/(y*A)) if y else None,
                fourCOriginOverYA=float(4*(c+offset)/(y*A)) if y else None,
                twoWOverY=float(2*W/y) if y else None,
                fractionIdentityError=float(4*c/(y*A)-2*sigma*2*W/y) if y else None,
                densitySupport=[float(hcore[0]),float(hcore[-1])],
                inverseZ=Z.tolist(), inverseAJ=J.tolist(), f0InverseZ=Zmodel.tolist(), f0InverseAJ=Jmodel.tolist(),
                gradientDensityIntervals=gradient, zeroDensity=zero,
                coreErrors=dict(f0GradientRelative=relative_rms(f0(XI)-qcore,qcore),
                                f0PrimitiveRelative=relative_rms(hmodel-hcore,hcore),
                                f0InversePositionOverW=rms(Zmodel-Z),
                                f0InverseAJRelative=relative_rms(Jmodel-J,J),
                                savedQxVsQRelative=relative_rms(Qxcore-qcore,qcore),
                                decimatedGradientRelative=relative_rms(qc_core-qcore,qcore),
                                decimatedPrimitiveRelative=relative_rms(hc_core-hcore,hcore)))


def select(rows, track, lo, hi):
    return [r for r in rows if r["track"] == track and lo <= r["height"] <= hi]


def pooled_errors(rows, sharedZ, sharedJ, fixedZ, fixedJ):
    Z = np.array([r["inverseZ"] for r in rows])
    J = np.array([r["inverseAJ"] for r in rows])
    Zf = np.array([r["f0InverseZ"] for r in rows])
    Jf = np.array([r["f0InverseAJ"] for r in rows])
    hd = np.array([v["error"] for r in rows for v in r["gradientDensityIntervals"]])
    h = np.array([v["densityInterval"] for r in rows for v in r["gradientDensityIntervals"]])
    return dict(count=len(rows), f0InversePositionOverW=rms(Zf-Z),
                f0InverseAJRelative=relative_rms(Jf-J,J),
                translationScaleInversePositionOverW=rms(sharedZ-Z),
                translationScaleInverseAJRelative=relative_rms(sharedJ-J,J),
                sameF0FixedKPositionOverW=rms(fixedZ-Z),
                sameF0FixedKAJRelative=relative_rms(fixedJ-J,J),
                f0GradientLevelDensityAbsolute=rms(hd),
                f0GradientLevelDensityRelative=relative_rms(hd,h),
                meanCPeakOverAW=float(np.mean([r["cOverAW"] for r in rows])),
                halfFractionRelativeRMS=rms(np.array([r["cOverAW"] for r in rows])/.5-1),
                meanQZeroOverA=float(np.mean([r["zeroDensity"]["qZeroOverA"] for r in rows])),
                zeroPredictionXOverW=rms([r["zeroDensity"]["xErrorOverW"] for r in rows]),
                zeroPredictionQOverA=rms([r["zeroDensity"]["qErrorOverA"] for r in rows]),
                derivativeConsistencyRMS=float(np.sqrt(np.mean([r["coreErrors"]["savedQxVsQRelative"]**2 for r in rows]))),
                interpolationSensitivityRMS=float(np.sqrt(np.mean([r["coreErrors"]["decimatedPrimitiveRelative"]**2 for r in rows]))))


def run():
    frames=[]
    for tau in [10,15,20]:
        p=trajectory_frame(tau)
        d=loadmat(ROOT/f"runs/geometry_loop_20260928/frame_tau{tau}.mat",simplify_cells=True)
        ell, a, P = float(d["ell"]), float(d["a"]), float(d["P"])
        qf,Qf=interp_fields(p)
        coarse=(RectBivariateSpline(p["y"][::2],p["x"][::2],p["q"][::2,::2]),
                RectBivariateSpline(p["y"][::2],p["x"][::2],p["Q"][::2,::2]))
        wall=peak_measure(p,0.,qf,Qf)
        offset=float(CubicSpline(p["x"],qf(0.,p["x"])[0]).integrate(-a/ell,0.))
        rows=[measure_row(p,qf,Qf,coarse,Y/ell,"native",Y,offset) for Y in NATIVE_HEIGHTS]
        rows += [measure_row(p,qf,Qf,coarse,y,"inner",y,offset) for y in INNER_HEIGHTS]
        frames.append(dict(tauLabel=tau,tau=float(p["tau"]),ell=ell,a=a,P=P,
                           wallOriginOffsetNormalized=offset,wallOriginOffsetNativeDensity=P*ell*offset,
                           reestimatedWallPeak=dict(mInner=wall["m"],QAtPeak=wall["c"],
                                                   note="Spline peak versus saved three-node centering; no re-zeroing is performed."),
                           rows=rows))

    # No optimization of rational exponents or small fractions.  The two
    # learned comparators each use only tau10 native Y=4,8,16,32.
    train=select(frames[0]["rows"],"native",4,32)
    sharedZ=np.mean([r["inverseZ"] for r in train],axis=0)
    sharedJ=np.mean([r["inverseAJ"] for r in train],axis=0)
    fixedK=float(np.mean([r["k"] for r in train]))
    fixedZ,fixedJ,_=np.array([model_at_density(s,fixedK) for s in DENSITY_LEVELS]).T
    affine=[]
    ks=np.array([r["k"] for r in train]);design=np.c_[np.ones(len(ks)),ks]
    for level in GRADIENT_LEVELS:
        for side in ["left","right"]:
            obs=np.array([next(v["densityInterval"] for v in r["gradientDensityIntervals"] if v["level"]==level and v["side"]==side) for r in train])
            coef=np.linalg.lstsq(design,obs,rcond=None)[0]
            xi=brentq(lambda t:float(f0(t)-level),-1.5,0.) if side=="left" else brentq(lambda t:float(f0(t)-level),0.,1.5)
            affine.append(dict(level=level,side=side,fitIntercept=float(coef[0]),fitSlope=float(coef[1]),
                               f0Intercept=float(H0(xi)),f0Slope=float(H1(xi))))
    groups=[]
    domains=[("native",.25,2,"native_near"),("native",4,32,"native_train_range"),
             ("native",64,256,"native_range_holdout"),("inner",0,1,"inner_near"),
             ("inner",1.25,2,"inner_range_holdout")]
    for frame in frames:
        for track,lo,hi,name in domains:
            rows=select(frame["rows"],track,lo,hi)
            err=pooled_errors(rows,sharedZ,sharedJ,fixedZ,fixedJ)
            aff_err=[];hobs=[]
            for r in rows:
                for v, fit in zip(r["gradientDensityIntervals"],affine):
                    aff_err.append(fit["fitIntercept"]+fit["fitSlope"]*r["k"]-v["densityInterval"])
                    hobs.append(v["densityInterval"])
            err["learnedAffineGradientDensityRelative"]=relative_rms(aff_err,hobs)
            groups.append(dict(tau=frame["tauLabel"],domain=name,**err))

    # Same-coordinate temporal drift is observational, without convergence
    # assumptions.  Density levels are normalized by each row's own AW.
    temporal=[]
    for track,heights in [("native",NATIVE_HEIGHTS),("inner",INNER_HEIGHTS)]:
        for height in heights:
            r=[next(r for r in f["rows"] if r["track"]==track and r["height"]==height) for f in frames]
            for i,j in [(0,1),(1,2),(0,2)]:
                temporal.append(dict(track=track,height=float(height),tauFrom=frames[i]["tauLabel"],tauTo=frames[j]["tauLabel"],
                                     inversePositionRMSOverW=rms(np.array(r[j]["inverseZ"])-r[i]["inverseZ"]),
                                     inverseAJRelative=relative_rms(np.array(r[j]["inverseAJ"])-r[i]["inverseAJ"],r[j]["inverseAJ"]),
                                     cOverAWChange=float(r[j]["cOverAW"]-r[i]["cOverAW"]),
                                     qZeroOverAChange=float(r[j]["zeroDensity"]["qZeroOverA"]-r[i]["zeroDensity"]["qZeroOverA"])))

    # An explicit physical-space density-level geometry comparator.  Each
    # line is fitted on Y=4..32; higher Y and all later times are held out.
    zero_lines=[]
    tau10_line=None
    for frame in frames:
        rows=select(frame["rows"],"native",4,32)
        yy=np.array([r["height"] for r in rows])
        xx=np.array([frame["a"]+frame["ell"]*r["zeroDensity"]["xInner"] for r in rows])
        coef=np.linalg.lstsq(np.c_[yy,np.ones(len(yy))],xx,rcond=None)[0]
        if tau10_line is None:
            tau10_line=coef.copy()
        checks=[]
        for lo,hi,name in [(4,32,"train_range"),(64,256,"range_holdout")]:
            test=select(frame["rows"],"native",lo,hi)
            def line_err(coeff):
                return rms([(coeff[0]*r["height"]+coeff[1]-frame["a"]-frame["ell"]*r["zeroDensity"]["xInner"])/(frame["ell"]*r["W"]) for r in test])
            checks.append(dict(domain=name,currentFrameLinePositionOverW=line_err(coef),
                               tau10LinePositionOverW=line_err(tau10_line)))
        zero_lines.append(dict(tau=frame["tauLabel"],slope=float(coef[0]),interceptNative=float(coef[1]),checks=checks))

    result=dict(title="Density relations and hodograph geometry on saved finite-time frames",
                scope="Saved-frame geometry only; no evolution, PDE closure, or asymptotic extrapolation.",
                densityReference="Q=(rho-rho(a,0))/(P ell). Density inverse/zero curve uses this wall-peak reference. c_origin adds the independently recovered rho(a,0)-rho(0,0) offset.",
                coordinateRestoration="X=a+ell*x, Y=ell*y, rho-rho(a,0)=P*ell*Q, rho_X=P*q",
                tracks={"native":"Same native Y across time, each y_inner=Y/ell(tau).",
                        "inner":"Same inner y across time, native Y=ell(tau)*y; a different physical comparison."},
                predeclaredFractionCandidates={"quarter":"c_origin/(Y*A_native)=1/4", "half":"c_peak/(A_inner*W_inner)=1/2", "localLinearZero":"(m-x_zero)*A/c_peak=1"},
                exactRatioIdentity="4*c_peak/(y*A) = (2*c_peak/(A*W))*(2*W/y). These three ratios cannot count as independent coincidences.",
                support="All measured and model density inversions use only |xi|<=1.5 and the original stored x,y domain. No f0 continuation beyond this core.",
                densityLevels=DENSITY_LEVELS.tolist(),gradientLevels=GRADIENT_LEVELS,
                masterModel="f0=32/(xi^6-3xi^5+3xi^4+3xi^3+28xi^2+32); H=H0+k*H1, H0=int_0^xi f0 cosh(3t/2)dt, H1=int_0^xi f0 sinh(3t/2)dt; inverse solves H=sigma. No fitted constants.",
                learnedComparators=dict(training="tau10 native Y=4,8,16,32 only; all tau15/20 and native Y64..256 are held out.",
                                        translationScaleZ=sharedZ.tolist(),translationScaleAJ=sharedJ.tolist(),
                                        affineGradientDensity=affine,designCondition=float(np.linalg.cond(design)),
                                        trainingKRange=[float(ks.min()),float(ks.max())],
                                        f0FixedSkewAblation=dict(k=fixedK,inverseZ=fixedZ.tolist(),inverseAJ=fixedJ.tolist(),
                                                                note="Same f0 as the measured-k predictor; isolates skew adaptation without changing the master.")),
                errorDefinitions={"position":"RMS(delta X/W)=RMS(delta Z), absolute in measured width units.",
                                  "AJ":"L2(delta AJ)/L2(observed AJ), AJ=A/q=partial_sigma Z.",
                                  "gradientDensity":"RMS error in H=(Q-c)/(AW), divided by RMS observed H for relative column.",
                                  "decimation":"Remove every other stored x and y node, interpolate back at identical physical sample points; sensitivity only, not an independent grid-convergence proof.",
                                  "derivative":"L2(partial_x Q-q)/L2(q) on 301 xi samples; checks saved/interpolated primitive consistency."},
                maxRatioIdentityError=max(abs(r["fractionIdentityError"]) for f in frames for r in f["rows"] if r["fractionIdentityError"] is not None),
                zeroDensityStraightLine=dict(training="Each-frame native Y=4,8,16,32; X_zero=slope*Y+intercept. Tau10 line additionally transfers across time.",fits=zero_lines),
                frames=frames,groupErrors=groups,temporalDrift=temporal)
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"density_relations.json").write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    write_report(result)
    print(json.dumps({"output":str(OUT),"groupErrors":groups,"affineCoefficients":affine},ensure_ascii=False,indent=2))


def write_report(d):
    pct=lambda x:f"{100*x:.4g}%"
    lines=["# 密度等值线与梯度等值线：三帧几何检索", "",
           "本报告仅分析已存 τ10、τ15、τ20 数据。远场尚未收敛，时间差是单帧形变与收敛趋势的观测，不排除极限几何结构。没有求解方程、拟合分数幂或外推尾部。", "",
           "最明确的新证据是密度坐标的交叉检验：native Y4..32 的 f0+k 逆位置误差从 τ10 的 0.00366W 降至 τ20 的 0.000266W，AJ 相对误差从 1.549% 降至 0.183%；τ20 native Y64..256 分别为 0.000303W 与 0.146%。密度零线的 q/A 则仍有显著时间漂移，因此不支持把它当成新的固定分数。", "",
           "## 基准、范围与预先指定的候选", "",
           r"保存坐标满足 $X=a+\ell x,Y=\ell y,\rho-\rho(a,0)=P\ell Q,\rho_X=Pq$。密度逆映射使用壁峰基准 $Q=0$；用户的 quarter 使用原点基准 $\rho(0,0)$，对应 $c_{\rm origin}=Q(m,y)+\int_{-a/\ell}^{0}q(x,0)dx$。JSON 同时保存两者，不混用。", "",
           "native 轨固定 Y；inner 轨固定 y，所以不同 τ 对应不同原生高度。内部计算复用已有峰、半高宽与偏斜测量，仅作为归一化，不将它们重复列为新发现。所有样本、密度逆解和整数母剖面预测限制在 |ξ|≤1.5 与原保存域内。", "",
           r"少量预指定检验为 $c_{\rm origin}/(YA)=1/4$、$c_{\rm peak}/(AW)=1/2$、峰处局部线性预测 $(m-x_0)A/c_{\rm peak}=1$；没有枚举分数或优化小数后寻找巧合。", "",
           "| τ | ρ(a,0)−ρ(0,0) 除以 Pℓ | spline 重测壁峰的 Q |", "|---:|---:|---:|"]
    for f in d["frames"]:
        lines.append(f"| {f['tauLabel']} | {f['wallOriginOffsetNormalized']:.3e} | {f['reestimatedWallPeak']['QAtPeak']:.3e} |")
    lines += ["", "后列 spline 峰相对原三点峰的微小差异是取峰算法/插值灵敏度，未重新平移基准。壁面 y=0 不作 quarter 除法。", "",
              "## 密度坐标下的逆剖面", "",
              r"定义 $\sigma=(Q-c)/(AW)$，$Z(\sigma,y)=[X(c+AW\sigma,y)-m]/W$，$AJ=A/q$。连续数据中 $\partial_\sigma Z=AJ$；这里对保存的 Q 与 q 分别插值，并另报一致性误差。共有密度带为 −0.4≤σ≤0.75，共 47 点。", "",
              "仅平移、尺度归一化的参照，是 τ10、native Y=4,8,16,32 的平均 Z 和 AJ。该参照没有偏斜自由度。第二个预测没有本轮拟合：使用已报告的整数核心 f0，积分得 H0、H1，以实测 k 解 H0(ξ)+kH1(ξ)=σ，然后 Z=Gk(ξ)、AJ=1/f0(ξ)。因此它同时接受 Q 与 q 的检验。", "",
              r"具体模型为 $f_0=32/(\xi^6-3\xi^5+3\xi^4+3\xi^3+28\xi^2+32)$，$G_k=[\sinh(3\xi/2)+k(\cosh(3\xi/2)-1)]/(3/2)$；$H_0=\int_0^\xi f_0(t)\cosh(3t/2)dt$，$H_1=\int_0^\xi f_0(t)\sinh(3t/2)dt$。W、k 由已有左右半高交点归一化；不从密度等值线重拟合它们。", "",
              "| τ | 域 | 仅平移缩放：δX/W RMS | 加 k 与 f0：δX/W RMS | 仅平移缩放：AJ 相对误差 | 加 k 与 f0：AJ 相对误差 |", "|---:|---|---:|---:|---:|---:|"]
    labels={"native_near":"native Y=.25..2", "native_train_range":"native Y=4..32", "native_range_holdout":"native Y=64..256", "inner_near":"inner y=0..1", "inner_range_holdout":"inner y=1.25..2"}
    for g in d["groupErrors"]:
        lines.append(f"| {g['tau']} | {labels[g['domain']]} | {g['translationScaleInversePositionOverW']:.5f} | {g['f0InversePositionOverW']:.5f} | {pct(g['translationScaleInverseAJRelative'])} | {pct(g['f0InverseAJRelative'])} |")
    lines += ["", "位置误差是绝对宽度单位，不能与相对百分数混读。τ10 native Y4..32 是平均逆剖面的训练域；τ15/20 全部以及 native Y64..256 是留出。inner 全域也是该参照的范围迁移检验。整数 f0 已由其他研究从 τ20 数据提出，所以其 τ20 数字是独立物理量检查，不能冒称无历史泄漏的时间预测；τ10、15 是逆时间验证。", "",
              f"为隔离偏斜自由度，另将同一个 f0 的 k 固定为 τ10 四训练高度的平均 {d['learnedComparators']['f0FixedSkewAblation']['k']:.6f}，不改变母剖面。这样不会把 f0 与 τ10 平均母线本身的差异误归因于 k。", "",
              "| τ | native 域 | 固定 k：δX/W RMS | 实测 k：δX/W RMS | 固定 k：AJ 相对误差 | 实测 k：AJ 相对误差 |", "|---:|---|---:|---:|---:|---:|"]
    for g in d["groupErrors"]:
        if g["domain"] in ["native_train_range","native_range_holdout"]:
            lines.append(f"| {g['tau']} | {labels[g['domain']]} | {g['sameF0FixedKPositionOverW']:.6f} | {g['f0InversePositionOverW']:.6f} | {pct(g['sameF0FixedKAJRelative'])} | {pct(g['f0InverseAJRelative'])} |")
    lines += ["", "当前 τ20 的 k 调整明显改善 AJ；τ10 的母剖面误差较大，加入 k 不保证每个指标改善，不能据此声称精确单主曲线。", "",
              "## 梯度等值线上的密度", "",
              r"在预指定 $q/A=r\in\{1/2,3/4,9/10\}$ 的左右交点直接测 $H_r^{\pm}=(Q(x_r^{\pm},y)-c)/(AW)$。若单一 q 母剖面与 sinh 几何同时成立，则每个 H 对 k 为仿射：$H_r^\pm=H_0(\xi_r^\pm)+kH_1(\xi_r^\pm)$。这里不是再次比较半高宽，而是检验梯度等值线上密度的位置。", "",
              "| r | 侧 | τ10 训练截距 | 训练 k 斜率 | f0 积分截距 | f0 积分 k 斜率 |", "|---:|---|---:|---:|---:|---:|"]
    for a in d["learnedComparators"]["affineGradientDensity"]:
        lines.append(f"| {a['level']:g} | {a['side']} | {a['fitIntercept']:.7f} | {a['fitSlope']:.7f} | {a['f0Intercept']:.7f} | {a['f0Slope']:.7f} |")
    kr=d["learnedComparators"]["trainingKRange"]
    lines += ["",f"仿射训练只含四个高度，k∈[{kr[0]:.5f},{kr[1]:.5f}]，未缩放设计矩阵条件数 {d['learnedComparators']['designCondition']:.2f}。窄 k 区间会放大斜率不确定性，必须看留出。积分常数仅为数值求积结果，未识别成新的简单分数。", "",
              "| τ | 域 | τ10 仿射 H 留出相对误差 | f0 积分 H 相对误差 | f0 积分 H 绝对 RMS |", "|---:|---|---:|---:|---:|"]
    for g in d["groupErrors"]:
        lines.append(f"| {g['tau']} | {labels[g['domain']]} | {pct(g['learnedAffineGradientDensityRelative'])} | {pct(g['f0GradientLevelDensityRelative'])} | {g['f0GradientLevelDensityAbsolute']:.6f} |")
    lines += ["", "## 密度零线与峰：不可当常数的关系", "",
              "x0 为壁峰基准 Q=0 的密度等值线；它不是 q 的半高线。下表只列代表性 native 高度，完整 near/inner 值在 JSON。", "",
              "| τ | Y | c_peak/(AW) | q(x0)/A | (m−x0)/W | (m−x0)A/c_peak | f0 对 x0 的误差 / W | f0 对 q(x0) 的误差 / A |", "|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for f in d["frames"]:
        for r in f["rows"]:
            if r["track"]=="native" and r["height"] in [4.,16.,64.,256.]:
                z=r["zeroDensity"]
                lines.append(f"| {f['tauLabel']} | {r['height']:g} | {r['cOverAW']:.6f} | {z['qZeroOverA']:.6f} | {z['offsetFromPeakOverW']:.6f} | {z['linearOffsetRatio']:.6f} | {z['xErrorOverW']:.6f} | {z['qErrorOverA']:.6f} |")
    lines += ["", "σ_zero=−c_peak/(AW)，所以 x0 和 q(x0)/A 本来就是密度幅度比与偏斜 k 的共同函数；不能仅因一个局部数值接近小分数就认作新的常数。峰处局部线性近似给 (m−x0)A/c=1，表中的系统偏差量化曲率作用。", "",
              "另在原生平面检验这条固定密度零线能否是单一直线 X0=sY+b。各帧仅 Y=4,8,16,32 拟合，误差除以同高度 native 宽度 ℓW；没有将 s 识别成分数。", "",
              "| τ | slope s | native 截距 b | 同帧训练 δX/W | 同帧 Y64..256 留出 | τ10 直线对该帧 Y4..32 | τ10 直线对该帧 Y64..256 |", "|---:|---:|---:|---:|---:|---:|---:|"]
    for f in d["zeroDensityStraightLine"]["fits"]:
        a,b=f["checks"]
        lines.append(f"| {f['tau']} | {f['slope']:.6f} | {f['interceptNative']:.6f} | {a['currentFrameLinePositionOverW']:.6f} | {b['currentFrameLinePositionOverW']:.6f} | {a['tau10LinePositionOverW']:.6f} | {b['tau10LinePositionOverW']:.6f} |")
    lines += ["", "当前密度零线的方向随时间变化，范围留出还有可辨曲率；单一直线/固定旋转角不足以精确描述这些帧。这一观测不否定极限密度线可能变直。", "",
              "## quarter 与 half 的依赖，以及时间漂移", "",
              r"确切恒等式 $4c/(yA)=[2c/(AW)]\,[2W/y]$。若先已观察 quarter 与 W≈Y/2，则 half 不再是独立巧合。反过来，两者在同一有限域的偏离也可能补偿。quarter 的完整基础表由其他任务报告，本处只给此依赖与密度零线新信息。", "",
              f"逐行恒等式的最大数值差为 {d['maxRatioIdentityError']:.3e}。例如 τ20 native Y=4 的原点 quarter 比为 1.06775，但 c_peak/(AW)=0.28108、2W/Y=1.89939；quarter 在此较准不意味着 half 较准。", "",
              "| τ | 域 | mean c_peak/(AW) | 相对 1/2 的 RMS 偏差 | mean q(x0)/A |", "|---:|---|---:|---:|---:|"]
    for g in d["groupErrors"]:
        lines.append(f"| {g['tau']} | {labels[g['domain']]} | {g['meanCPeakOverAW']:.6f} | {pct(g['halfFractionRelativeRMS'])} | {g['meanQZeroOverA']:.6f} |")
    lines += ["", "| 轨 | 高度 | 时段 | 归一逆位置漂移 RMS | AJ 漂移相对后帧 | c/(AW) 改变量 | q(x0)/A 改变量 |", "|---|---:|---|---:|---:|---:|---:|"]
    for t in d["temporalDrift"]:
        if t["tauFrom"]==10 and t["tauTo"]==20 and ((t["track"]=="native" and t["height"] in [4.,16.,64.,256.]) or (t["track"]=="inner" and t["height"] in [.25,1.,2.])):
            lines.append(f"| {t['track']} | {t['height']:g} | 10→20 | {t['inversePositionRMSOverW']:.6f} | {pct(t['inverseAJRelative'])} | {t['cOverAWChange']:+.6f} | {t['qZeroOverAChange']:+.6f} |")
    lines += ["", "所有时段 10→15、15→20、10→20 均保存在 JSON。远场的变化没有被当成极限结构的反证；half 与 q_zero 当前不够时间稳定，不能从单帧挑常数。", "",
              "## 插值与原函数一致性", "",
              "| τ | 域 | Q_x 对 q 相对 RMS | 隔点抽稀原函数灵敏度 |", "|---:|---|---:|---:|"]
    for g in d["groupErrors"]:
        lines.append(f"| {g['tau']} | {labels[g['domain']]} | {pct(g['derivativeConsistencyRMS'])} | {pct(g['interpolationSensitivityRMS'])} |")
    lines += ["", "抽稀同时删除每隔一个 x、y 节点，再在完全相同物理点比较；只作插值敏感性检查，不冒充独立网格收敛误差。它没有重新拟合峰参数，因此也不包含全部取峰不确定性。当前最有信息的低维关系，是将实测 k 带入密度原函数/逆映射后，是否同时减小 Q、AJ 与零线误差；漂亮的单一密度逆曲线及 half 常数都应按时间和范围留出判断。", "",
              "复现：`OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 PYTHONPATH=/private/tmp/ipm_steady_deps /Users/cbjin/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 solver/steady_analysis/profile_density_relations.py`。"]
    (OUT/"DENSITY_RELATIONS.md").write_text("\n".join(lines)+"\n")


if __name__ == "__main__":
    run()
