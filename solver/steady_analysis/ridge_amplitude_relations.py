"""Joint ridge/amplitude relations; reuse Luna's measured height curves.

Saved-field data analysis only. No PDE closure or profile reconstruction.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator, RectBivariateSpline
from scipy.io import loadmat
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/ridge_amplitude_20260928"
SOURCE = OUT / "amplitude_detail.json"
PS = [("1/2", 1/2), ("3/5", 3/5), ("5/8", 5/8), ("2/3", 2/3)]
DOMAINS = [("near", .025, 2.), ("middle", 4., 64.), ("far", 128., 1024.)]


def rms(x):
    x = np.asarray(x, dtype=float)
    return float(np.sqrt(np.mean(x*x))) if x.size else None


def rr(delta, obs):
    return rms(delta)/rms(obs) if len(obs) else None


def mask(h, lo, hi):
    return (h >= lo) & (h <= hi)


def load_curves():
    source = json.loads(SOURCE.read_text())
    all_frames = []
    for f in source["timeFrames"]:
        tau = int(round(float(f["tau"])))
        rows = f["profileRows"]
        Y = np.array([r["Ynative"] for r in rows], dtype=float)
        m = np.array([r["mCenteredNative"] for r in rows], dtype=float)
        A = np.array([r["A_native"] for r in rows], dtype=float)
        order = np.argsort(Y); Y,m,A=Y[order],m[order],A[order]
        assert Y[0] == 0 and np.all(np.diff(Y)>0)
        A0 = float(A[0]); B = A/A0
        mat = loadmat(ROOT/f"runs/geometry_loop_20260928/frame_tau{tau}.mat", simplify_cells=True)
        ell,a,P = (float(mat[k]) for k in ["ell","a","P"])
        qf = RectBivariateSpline(mat["y"],mat["x"],mat["q"])
        coarse = RectBivariateSpline(mat["y"][::2],mat["x"][::2],mat["q"][::2,::2])
        mp = PchipInterpolator(Y,m); bp = PchipInterpolator(Y,B)
        # Exact same requested height grid, interpreted as inner y.  This
        # reparametrizes already measured curves and does not re-extract peaks.
        tracks = []
        for track in ["native","inner"]:
            h = Y.copy()
            YY = h if track == "native" else ell*h
            valid = (YY >= Y[0]) & (YY <= Y[-1])
            h,YY = h[valid],YY[valid]
            mm = m[valid] if track == "native" else mp(YY)
            BB = B[valid] if track == "native" else bp(YY)
            xinner = mm/ell; yinner = YY/ell
            assert yinner[-1] <= mat["y"][-1]
            qxx = qf.ev(yinner,xinner,dy=2)
            qxy = qf.ev(yinner,xinner,dx=1,dy=1)
            qy = qf.ev(yinner,xinner,dx=1)
            qq = qf.ev(yinner,xinner)
            deriv = -qxy/qxx
            # dA_native/dY=P*q_y/ell. This is the derivative of a peak
            # envelope when q_x=0; below we report ridge-stationarity error.
            dB = P*qy/(ell*A0)
            pY = -YY*dB/BB
            pm = -mm*dB/(BB*deriv)
            cxx = coarse.ev(yinner,xinner,dy=2)
            cxy = coarse.ev(yinner,xinner,dx=1,dy=1)
            cy = coarse.ev(yinner,xinner,dx=1)
            dc = -cxy/cxx
            pc = -mm*(P*cy/(ell*A0))/(BB*dc)
            mcs = CubicSpline(YY,mm)
            bcs = CubicSpline(YY,BB)
            d1 = mcs(YY,1)
            p1 = -mm*bcs(YY,1)/(BB*d1)
            coordm = mm if track == "native" else mm/ell
            tracks.append(dict(track=track,height=h,Ynative=YY,m=coordm,
                               mCenteredNative=mm,B=BB,mPrime=deriv,
                               pY=pY,pm=pm,mPrimeCurve=d1,pmCurve=p1,
                               mPrimeDecimated=dc,pmDecimated=pc,
                               amplitudeDecimatedRelative=coarse.ev(yinner,xinner)/qq-1,
                               amplitudeCurveConsistencyRelative=BB*A0/(P*qq)-1,
                               ridgeStationarity=np.abs(qf.ev(yinner,xinner,dy=1))*np.maximum(1,yinner)/qq,
                               wallMPrime=float(deriv[0]),wallM=float(coordm[0])))
        H=next(r["Ynative"] for r in f["thresholdCrossings"] if r["fraction"]==.5)
        all_frames.append(dict(tau=tau,a=a,P=P,ell=ell,A0Native=A0,
                               halfHeightNative=float(H),halfHeightInner=float(H/ell),tracks=tracks))
    return sorted(all_frames,key=lambda f:f["tau"])


def fit_affine(t, coord, p):
    use = mask(t["height"],4,64)
    x = t["height"][use] if coord=="height" else t["m"][use]
    z = t["B"][use]**(-1/p)
    des = np.c_[np.ones(len(x)),x]
    coeff = np.linalg.lstsq(des,z,rcond=None)[0]
    residual = z-des@coeff
    # Projection onto a quadratic after removing the constant/linear part;
    # a large fraction means smooth bending, rather than unstructured scatter.
    u=(x-x.mean())/(x.max()-x.min())
    v=u*u-np.c_[np.ones(len(u)),u]@np.linalg.lstsq(np.c_[np.ones(len(u)),u],u*u,rcond=None)[0]
    energy = float((v@residual)**2/((v@v)*(residual@residual))) if residual@residual else 0.
    sig=np.sign(residual); sign_changes=int(np.count_nonzero(sig[1:]!=sig[:-1]))
    return dict(intercept=float(coeff[0]),slope=float(coeff[1]),
                trainingPointCount=int(use.sum()),trainingHeightRange=[float(t["height"][use].min()),float(t["height"][use].max())],
                quadraticResidualEnergyFraction=energy,
                transformedTrainingRelativeRMS=rr(residual,z),
                residualSignChanges=sign_changes,
                trainingCoordinate=x.tolist(),transformedResidual=residual.tolist())


def predict_metrics(t, coord, p, fit, lo, hi):
    use=mask(t["height"],lo,hi)
    x = t["height"][use] if coord=="height" else t["m"][use]
    B = t["B"][use]
    z = fit["intercept"]+fit["slope"]*x
    valid=z>0
    e=np.full(len(B),np.nan)
    e[valid]=z[valid]**(-p)/B[valid]-1
    return dict(pointCount=len(B),validCount=int(valid.sum()),
                amplitudeRelativeRMS=rms(e[valid]),amplitudeRelativeMean=float(e[valid].mean()) if valid.any() else None,
                amplitudeRelativeMax=float(np.max(abs(e[valid]))) if valid.any() else None,
                height=t["height"][use].tolist(),relativeResidual=[float(v) if np.isfinite(v) else None for v in e])


def temporal_coordinate_test(frames):
    out=[]
    for track in ["native","inner"]:
        train=next(t for t in frames[0]["tracks"] if t["track"]==track)
        for ref_range in ["middle_only","whole_saved_range"]:
            ref=mask(train["height"],4,64) if ref_range=="middle_only" else train["height"]>=.025
            yh=train["height"][ref];mh=train["m"][ref];bh=train["B"][ref]
            assert np.all(np.diff(mh)>0)
            fy=PchipInterpolator(np.log(yh),np.log(bh),extrapolate=False)
            fm=PchipInterpolator(np.log(mh),np.log(bh),extrapolate=False)
            for frame in frames[1:]:
                t=next(t for t in frame["tracks"] if t["track"]==track)
                for name,lo,hi in DOMAINS:
                    common=mask(t["height"],lo,hi)&(t["height"]>=yh[0])&(t["height"]<=yh[-1])&(t["m"]>=mh[0])&(t["m"]<=mh[-1])
                    if not common.any():
                        continue
                    obs=t["B"][common]
                    py=np.exp(fy(np.log(t["height"][common])))
                    pm=np.exp(fm(np.log(t["m"][common])))
                    out.append(dict(track=track,referenceRange=ref_range,tau=frame["tau"],domain=name,
                                    count=int(common.sum()),heightRange=[float(t["height"][common].min()),float(t["height"][common].max())],
                                    heightCoordinateRMS=rms(py/obs-1),ridgeCoordinateRMS=rms(pm/obs-1),
                                    heightCoordinateMean=float(np.mean(py/obs-1)),ridgeCoordinateMean=float(np.mean(pm/obs-1))))
    return out


def half_height_test(frames):
    """No fitting: one measured half-amplitude anchor per frame."""
    zlog=np.unique(np.r_[np.geomspace(.025,64,161),1.])
    zuniform=np.linspace(0,64,513)
    curves=[];evaluators=[]
    for f in frames:
        t=next(t for t in f["tracks"] if t["track"]=="native")
        fn=PchipInterpolator(t["height"],t["B"],extrapolate=False)
        H=f["halfHeightNative"]
        assert H*64<=t["height"][-1]
        blogs=fn(H*zlog);bunif=fn(H*zuniform)
        curves.append(dict(tau=f["tau"],Hnative=H,Hinner=f["halfHeightInner"],
                           BAtMeasuredAnchor=float(fn(H)),BLogZ=blogs.tolist(),BUniformZ=bunif.tolist()))
        evaluators.append(fn)
    fixedlog=(1+3*zlog)**(-.5);fixeduni=(1+3*zuniform)**(-.5)
    tests=[];native=[]
    for i,f in enumerate(frames):
        c=curves[i];obs=np.array(c["BLogZ"]);uni=np.array(c["BUniformZ"])
        tests.append(dict(tau=f["tau"],tau10RescaledRelativeL2LogZ=rr(np.array(curves[0]["BLogZ"])-obs,obs),
                          tau10RescaledRelativeL2UniformZ=rr(np.array(curves[0]["BUniformZ"])-uni,uni),
                          tau10RescaledLogRMSE=rms(np.log(np.array(curves[0]["BLogZ"])/obs)),
                          fixedProfileRelativeL2LogZ=rr(fixedlog-obs,obs),
                          fixedProfileRelativeL2UniformZ=rr(fixeduni-uni,uni),
                          fixedProfileLogRMSE=rms(np.log(fixedlog/obs)),
                          fixedProfileMeanLogRatio=float(np.log(fixedlog/obs).mean())))
        t=next(t for t in f["tracks"] if t["track"]=="native")
        for name,lo,hi in DOMAINS:
            use=mask(t["height"],lo,hi)
            y=t["height"][use];obs=t["B"][use]
            # Compare fixed-native and scaled predictions at exactly the same
            # queried native heights, respecting reference support for both.
            yref=y*frames[0]["halfHeightNative"]/f["halfHeightNative"]
            valid=(yref>=0)&(yref<=frames[0]["tracks"][0]["height"][-1])
            y,obs,yref=y[valid],obs[valid],yref[valid]
            original=evaluators[0](y);rescaled=evaluators[0](yref)
            pred=(1+3*y/f["halfHeightNative"])**(-.5)
            native.append(dict(tau=f["tau"],domain=name,count=len(y),
                               heightRange=[float(y.min()),float(y.max())],
                               fixedNativeReferenceRelativeRMS=rms(original/obs-1),
                               rescaledReferenceRelativeRMS=rms(rescaled/obs-1),
                               fixedNativeReferenceLogRMSE=rms(np.log(original/obs)),
                               rescaledReferenceLogRMSE=rms(np.log(rescaled/obs)),
                               fixedProfileRelativeRMS=rms(pred/obs-1),fixedProfileLogRMSE=rms(np.log(pred/obs))))
    return dict(definition="H_tau is the directly measured native height with A/A0=1/2; z=Y/H_tau. No fitted parameter. Inner H=H_native/ell gives identical rescaled curve.",
                fixedProfile="B=(1+3z)^(-1/2); B(0)=1,B(1)=1/2,L_native=H_native/3.",
                weighting="uniform-log-z .025..64 emphasizes multiplicative profile shape; uniform-z 0..64 L2 reduces dominance by the dense near-wall part. All finite sampled domains, not infinite-profile norms.",
                zLog=zlog.tolist(),zUniform=zuniform.tolist(),curves=curves,tests=tests,nativeDomainTests=native)


def derivative_relation(frames):
    result=[]
    for track in ["native","inner"]:
        tr=next(t for t in frames[0]["tracks"] if t["track"]==track)
        use=mask(tr["height"],4,64)
        order=np.argsort(tr["B"][use]); bb=tr["B"][use][order]; vv=tr["mPrime"][use][order]
        master=PchipInterpolator(bb,vv,extrapolate=False)
        linear=np.linalg.lstsq(np.c_[np.ones(len(bb)),bb],vv,rcond=None)[0]
        for frame in frames:
            t=next(t for t in frame["tracks"] if t["track"]==track)
            for name,lo,hi in DOMAINS:
                region=mask(t["height"],lo,hi)
                b=t["B"][region];v=t["mPrime"][region]
                common=region&(t["B"]>=bb[0])&(t["B"]<=bb[-1])
                obs=t["mPrime"][common]
                pred=master(t["B"][common])
                normal_obs=obs/t["wallMPrime"]
                normal_pred=pred/tr["wallMPrime"]
                result.append(dict(track=track,tau=frame["tau"],domain=name,
                                   regionCount=int(region.sum()),spearmanMPrimeVsB=float(spearmanr(b,v).statistic),
                                   mPrimeRange=[float(v.min()),float(v.max())],
                                   derivativeCurveRelativeRMS=rr(t["mPrimeCurve"][region]-v,v),
                                   derivativeDecimationRelativeRMS=rr(t["mPrimeDecimated"][region]-v,v),
                                   pmCurveRMSDifference=rms(t["pmCurve"][region]-t["pm"][region]),
                                   pmDecimationRMSDifference=rms(t["pmDecimated"][region]-t["pm"][region]),
                                   amplitudeDecimationRMS=rms(t["amplitudeDecimatedRelative"][region]),
                                   amplitudeCurveConsistencyRMS=rms(t["amplitudeCurveConsistencyRelative"][region]),
                                   pmMedian=float(np.median(t["pm"][region])),pmRange=[float(t["pm"][region].min()),float(t["pm"][region].max())],
                                   pYMedian=float(np.median(t["pY"][region])),pYRange=[float(t["pY"][region].min()),float(t["pY"][region].max())],
                                   commonAmplitudeCount=int(common.sum()),
                                   commonAmplitudeRange=[float(bb[0]),float(bb[-1])],
                                   temporalMasterRelativeRMS=rr(pred-obs,obs),
                                   wallSlopeNormalizedTemporalRelativeRMS=rr(normal_pred-normal_obs,normal_obs),
                                   linearMasterRelativeRMS=rr(linear[0]+linear[1]*t["B"][common]-obs,obs),
                                   maxRidgeStationarity=float(t["ridgeStationarity"][region].max())))
    return result


def serialize(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,dict):return {k:serialize(v) for k,v in x.items()}
    if isinstance(x,list):return [serialize(v) for v in x]
    return x


def run():
    frames=load_curves();affine=[]
    for track in ["native","inner"]:
        for coord in ["height","ridge"]:
            for pname,p in PS:
                for fitframe in frames:
                    train=next(t for t in fitframe["tracks"] if t["track"]==track)
                    fit=fit_affine(train,coord,p)
                    tests=[]
                    # All-current-frame refits are diagnostic; only a tau10
                    # fit followed by tau15/20 is a strict time holdout.
                    targets=frames if fitframe["tau"]==10 else [fitframe]
                    for target in targets:
                        t=next(t for t in target["tracks"] if t["track"]==track)
                        for name,lo,hi in DOMAINS:
                            tests.append(dict(tau=target["tau"],domain=name,**predict_metrics(t,coord,p,fit,lo,hi)))
                    affine.append(dict(track=track,coordinate=coord,pLabel=pname,p=p,fitTau=fitframe["tau"],fit=fit,tests=tests))
    result=dict(title="Joint ridge position and maximum-gradient amplitude relations",
                source=str(SOURCE.relative_to(ROOT)),
                definition="At each height Y, A(Y)=max_X rho_X(X,Y), m_absolute=argmax_X rho_X; this is not max rho or max velocity.",
                coordinateDefinitions={"native":"height=Y, ridge=m_absolute-a=ell*m_inner, amplitude B=A_native/A_native(0)",
                                       "inner":"height=y=Y/ell, ridge=m_inner, B unchanged; evaluate existing measured curves at native Y=ell*y with PCHIP.",
                                       "wallCentering":"Native centered m subtracts the saved three-node peak a. The spline wall maximum has a small residual m(0); no silent re-centering."},
                hypotheses={"fixedP":[x[0] for x in PS],"affine":"B^(-1/p)=b0+b1*height or b0+b1*ridge; two parameters, ordinary least squares in the transformed amplitude.",
                            "timeCoordinateComparison":"PCHIP log B at matched positive Y or m; both errors use exactly the same shared query heights, no extrapolation.",
                            "derivativeGraph":"m'(Y) against B; tau10 middle curve and its linear approximation transfer only on their common B range. A within-frame graph is not evidence for cross-time invariance."},
                splitting="Middle 4..64 trains affine and derivative reference. Near .025..2 and far 128..1024 are range holdouts; tau15 and tau20 are time holdouts for tau10 models. Whole-range PCHIP is used only for cross-time coordinate comparison, never claimed as a height holdout.",
                derivativeMeasurement="At the supplied spline peak, m'=-q_xy/q_xx; A'_native=P*q_y/ell. These are derivatives of the measured maximum, not evolution/PDE equations. p_m=-(m/A)*(A'/m'), p_Y=-(Y/A)*A'. Reused height-curve derivatives and alternate-node q interpolation quantify sensitivity.",
                interpretationLimits="Finite-time data, with unconverged far field. A transformed-affine candidate that fits one range does not establish a power tail or select an asymptotic exponent.",
                frames=serialize(frames),affineTests=affine,
                halfHeightTests=half_height_test(frames),
                temporalCoordinateTests=temporal_coordinate_test(frames),derivativeRelations=derivative_relation(frames))
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"relations.json").write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    write_report(result)
    print(json.dumps({"temporalCoordinateTests":result["temporalCoordinateTests"],"derivativeRelations":result["derivativeRelations"]},ensure_ascii=False,indent=2))


def write_report(d):
    pct=lambda x:"—" if x is None else f"{100*x:.4g}%"
    lines=["# 峰位 m 与每层最大梯度 A 的联合关系", "",
           "这里 A(Y)=max_X ρ_X(X,Y)，m 是该最大值的横向位置。仅分析三帧保存数据，复用 Luna 的高度曲线，不重新提取整剖面；没有方程闭合或动力学推断。", "",
           "当前最强的简化是单个纵向长度锚点：实测半幅高度 H 从 τ10 的 2.14444 降至 τ20 的 1.30160；native 中区的跨时幅值差由固定 Y 的 30.47% 降至 H 缩放后的 2.695%。共同 z=.025..64 上，缩放曲线的等z L2差为2.158%，logRMSE=.018959。纵向收缩解释主要变化，仍留有可辨的形态漂移。", "",
           "峰位提供第二层信息：改用 m 使同一中区跨时差由30.47%降至13.59%。τ20 中区训练的 B^(−5/3) 对 native m 仿射关系，在中区/远域的幅值误差为0.628%/0.723%；但near误差约27.96%，τ10关系迁移到τ20中区仍有13.57%差异，因此只能作为当前中远区的低复杂度候选。", "",
           "## 坐标、留出与比较标准", "",
           "native 用原生 Y、m_centered=m_absolute−a；inner 用 y=Y/ℓ、m_inner=m_centered/ℓ。幅值 B=A(Y)/A(0) 在两坐标中相同。inner 曲线由原生已测 m、B 曲线作保形插值，固定 inner y 跨时对应不同 native Y。所有插值限定保存域，不外推。", "",
           "预指定四个 p：1/2、3/5、5/8、2/3。对每个候选分别拟合 B^(−1/p)=b0+b1Y 与 b0+b1m，中区高度4..64训练，near .025..2、far128..1024留出；τ10拟合还迁移到τ15/20。拟合在变换幅值上最小二乘，评价统一还原为幅值相对 RMS，避免仅凭变换图看起来直就选指数。原生 A^(−1/p) 与 B^(−1/p) 只相差常数 A(0)^(−1/p)。", "",
           "除明示L2的表外，幅值相对RMS指 sqrt(mean((A预测/A实测−1)^2))。近中远区按121层log高度样本等权。", "",
           "## 首要检验：半幅高度归一化是否解释跨时变化", "",
           r"直接采用实测 $H_\tau=Y_{A/A_0=1/2}$，比较 $B_\tau(H_\tau z)$，没有拟合长度。native 与 inner 的 H 分别相差 ℓ，因此归一后的 B(z) 完全相同。固定候选 $B=(1+3z)^{-1/2}$ 自动经过 (0,1)、(1,1/2)，对应 L=H/3；只检验，不调参数。", "",
           "| τ | H native | H inner | 插值 B(H) |", "|---:|---:|---:|---:|"]
    for c in d["halfHeightTests"]["curves"]:
        lines.append(f"| {c['tau']} | {c['Hnative']:.7f} | {c['Hinner']:.7f} | {c['BAtMeasuredAnchor']:.8f} |")
    lines += ["", "共同z=.025..64全部在原保存数据内，使用等log z与等z两种权重；后者还包含z=0。等log z的logRMSE是 log(Bpred/Bobs) 的 RMS，保留远端相对误差，不能与近壁主导的L2混为一谈。", "",
              "| τ | τ10缩放母线 L2 等log z | τ10缩放母线 L2 等z | τ10缩放母线 logRMSE | 固定sqrt候选 L2 等log z | 固定sqrt候选 L2 等z | 固定sqrt候选 logRMSE |", "|---:|---:|---:|---:|---:|---:|---:|"]
    for r in d["halfHeightTests"]["tests"]:
        lines.append(f"| {r['tau']} | {pct(r['tau10RescaledRelativeL2LogZ'])} | {pct(r['tau10RescaledRelativeL2UniformZ'])} | {r['tau10RescaledLogRMSE']:.6f} | {pct(r['fixedProfileRelativeL2LogZ'])} | {pct(r['fixedProfileRelativeL2UniformZ'])} | {r['fixedProfileLogRMSE']:.6f} |")
    lines += ["", "以下在相同native高度逐点比较τ10固定native参照与H缩放参照。样本只取两者都不需要外插的共同范围，JSON给实际高度端点。", "",
              "| τ | native域 | 点数 | τ10固定native A 相对RMS | τ10 H缩放 A 相对RMS | 固定native logRMSE | H缩放 logRMSE | 固定sqrt候选相对RMS |", "|---:|---|---:|---:|---:|---:|---:|---:|"]
    for r in d["halfHeightTests"]["nativeDomainTests"]:
        lines.append(f"| {r['tau']} | {r['domain']} | {r['count']} | {pct(r['fixedNativeReferenceRelativeRMS'])} | {pct(r['rescaledReferenceRelativeRMS'])} | {r['fixedNativeReferenceLogRMSE']:.6f} | {r['rescaledReferenceLogRMSE']:.6f} | {pct(r['fixedProfileRelativeRMS'])} |")
    lines += ["", "如果H缩放显著改善跨时曲线，说明纵向长度变化解释了主要幅值差；残余logRMSE仍检验形态演化。固定候选经过两个锚点本身不算拟合成功，必须看锚点之外的全范围误差。", "",
              "## 换用 m 是否改善跨时幅值重合", "",
           "第一项比较不挑 p：以 τ10 的 log B 对 log Y 或 log m 的保形插值为参照，两者在完全相同且共同可达的目标高度上比较。whole_saved_range 是纯时间迁移，不是高度留出；middle_only 使用训练中区，禁止向外外推。", "",
           "| 轨 | τ | τ10参照范围 | 目标域 | 点数/实际高度范围 | A(Y) 跨时误差 | A(m) 跨时误差 |", "|---|---:|---|---|---|---:|---:|"]
    for r in d["temporalCoordinateTests"]:
        hr=r["heightRange"]
        lines.append(f"| {r['track']} | {r['tau']} | {r['referenceRange']} | {r['domain']} | {r['count']}; {hr[0]:.4g}..{hr[1]:.4g} | {pct(r['heightCoordinateRMS'])} | {pct(r['ridgeCoordinateRMS'])} |")
    lines += ["", "单帧中 B 随高度下降而 m 上升，改用 m 总能形成单值曲线；上表要求跨时误差也下降，才算额外信息。表中不同域的误差不能以不同支持点数互相替代。", "",
              "## 预指定分数幂：严格 τ10 训练的迁移", "",
              "各行两个坐标模型都用相同 τ10 中区高度拟合。表列 τ10中区训练、同帧远域、τ20中区及远域；τ15、near以及逐高度残差全部在JSON。", "",
              "| 轨 | 坐标 | p | τ10中区 | τ10远域留出 | τ20中区时间留出 | τ20远域双重留出 | 训练残差二次弯曲投影 |", "|---|---|---|---:|---:|---:|---:|---:|"]
    for r in d["affineTests"]:
        if r["fitTau"]!=10:continue
        get=lambda tau,dom:next(t for t in r["tests"] if t["tau"]==tau and t["domain"]==dom)
        vals=[get(10,"middle"),get(10,"far"),get(20,"middle"),get(20,"far")]
        lines.append(f"| {r['track']} | {r['coordinate']} | {r['pLabel']} | "+" | ".join(pct(v['amplitudeRelativeRMS']) for v in vals)+f" | {pct(r['fit']['quadraticResidualEnergyFraction'])} |")
    lines += ["", "二次弯曲投影是变换幅值残差在去掉常数/线性部分的二次函数上的能量比例。高比例提示系统弯曲；它不是新的拟合模型，也没有被用于修正预测。仿射在某些 near 点变为非正时，无法还原实数幅值，JSON明确记录validCount，而非裁剪后把误差藏掉。", "",
              "## τ20 重新拟合：高度内结构与误差形状", "",
              "这是同帧中区拟合诊断，不是时间预测。系数是数据拟合值，未认作简单分数。signed mean 为预测A/实测A−1；训练低误差而far有同号系统偏差说明留出仍弯曲。", "",
              "| 轨 | 坐标 | p | 截距 b0 | 斜率 b1 | 中区RMS | far RMS | far signed mean | 弯曲投影 |", "|---|---|---|---:|---:|---:|---:|---:|---:|"]
    for r in d["affineTests"]:
        if r["fitTau"]!=20:continue
        mid=next(t for t in r["tests"] if t["domain"]=="middle")
        far=next(t for t in r["tests"] if t["domain"]=="far")
        lines.append(f"| {r['track']} | {r['coordinate']} | {r['pLabel']} | {r['fit']['intercept']:.6g} | {r['fit']['slope']:.6g} | {pct(mid['amplitudeRelativeRMS'])} | {pct(far['amplitudeRelativeRMS'])} | {pct(far['amplitudeRelativeMean'])} | {pct(r['fit']['quadraticResidualEnergyFraction'])} |")
    lines += ["", "## 局部斜率与 m′ 对幅值的关系", "",
              "p_Y=−dlogA/dlogY，p_m=−dlogA/dlogm。这里 m 用保存壁峰 a 的中心化位置；没有强行减去 spline 重测 m(0)。仅在正高度/正m使用log意义，近壁值不当幂尾指数。直接从已取峰点的保存 q 插值求 m′=−q_xy/q_xx、A′=Pq_y/ℓ，只是最大值曲线的微分公式。", "",
              "| 轨 | τ | 域 | median p_Y | median p_m | p_m 最小..最大 | Spearman(m′,B) | m′最小..最大 |", "|---|---:|---|---:|---:|---|---:|---|"]
    for r in d["derivativeRelations"]:
        pm=r["pmRange"];vv=r["mPrimeRange"]
        lines.append(f"| {r['track']} | {r['tau']} | {r['domain']} | {r['pYMedian']:.5f} | {r['pmMedian']:.5f} | {pm[0]:.5f}..{pm[1]:.5f} | {r['spearmanMPrimeVsB']:.5f} | {vv[0]:.5f}..{vv[1]:.5f} |")
    lines += ["", "只要B在每帧单调，m′(B)单值本身不提供压缩证据。下一表用τ10中区m′(B)经验参照检验同一B处是否跨时重合；同时提供仅以壁面斜率归一化 m′/m′(0) 的参照。只在τ10训练B共同值域内比较，不补远域。", "",
              "| 轨 | τ | 域 | 共同幅值点数 | m′(B)跨时相对误差 | m′/m′(0)误差 | τ10线性参照误差 |", "|---|---:|---|---:|---:|---:|---:|"]
    for r in d["derivativeRelations"]:
        if not r["commonAmplitudeCount"]:continue
        lines.append(f"| {r['track']} | {r['tau']} | {r['domain']} | {r['commonAmplitudeCount']} | {pct(r['temporalMasterRelativeRMS'])} | {pct(r['wallSlopeNormalizedTemporalRelativeRMS'])} | {pct(r['linearMasterRelativeRMS'])} |")
    lines += ["", "## 微分灵敏度与可解释范围", "",
              "曲线导数来自121层既有m、A的三次插值；direct导数来自原2D q插值。隔点抽稀q场后仍在原有峰点比较，作为插值灵敏度，不是独立网格收敛。inner峰位由已测曲线插值，因此另报非零ridge stationarity，避免把插值位置误认为精确峰。", "",
              "| 轨 | τ | 域 | 曲线/direct m′ RMS | 抽稀/direct m′ RMS | p_m曲线差 RMS | p_m抽稀差 RMS | 幅值抽稀RMS | max abs(q_x)·max(1,y)/q |", "|---|---:|---|---:|---:|---:|---:|---:|---:|"]
    for r in d["derivativeRelations"]:
        lines.append(f"| {r['track']} | {r['tau']} | {r['domain']} | {pct(r['derivativeCurveRelativeRMS'])} | {pct(r['derivativeDecimationRelativeRMS'])} | {r['pmCurveRMSDifference']:.5g} | {r['pmDecimationRMSDifference']:.5g} | {pct(r['amplitudeDecimationRMS'])} | {r['maxRidgeStationarity']:.3e} |")
    lines += ["", "τ20 native far 的 p_m median 为0.59670，但q场隔点抽稀造成的p_m RMS差为0.02531，已与3/5和5/8的间隔相当；不得把接近0.6宣称为精确分数。幅值拟合本身的抽稀灵敏度较小，平滑残差的二次投影仍占很大能量，提示未消除的系统弯曲。", "",
              "报告的曲率、斜率漂移与插值灵敏度应一起看。远场尚未收敛，当前拟合失败不能否定极限几何关系；当前拟合成功也不能据单区、单帧选极限指数。高度网格和权重沿用原始121层log采样，本报告RMS对这些样本等权，而非对原生Y积分等权。", "",
              "复现：`OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 PYTHONPATH=/private/tmp/ipm_steady_deps /Users/cbjin/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 solver/steady_analysis/ridge_amplitude_relations.py`。"]
    (OUT/"RELATIONS.md").write_text("\n".join(lines)+"\n")


if __name__=="__main__":run()
