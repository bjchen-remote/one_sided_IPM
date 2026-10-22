"""Held-radius derivative audit of the existing direct-psi polar SVD.

No new mother functions, centres, or rank search. Fits radial coefficients
only on the same tau20 training radii used by streamfunction_fresh_polar.py.
"""
from pathlib import Path
import json
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS","2")
os.environ.setdefault("OMP_NUM_THREADS","2")

import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline
from streamfunction_sources import load_native

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"runs/streamfunction_fresh_20260928"


def error_summary(pred,obs):
    # A row is one held-out radius; combine any component/angle axes.
    axes=tuple(range(1,obs.ndim))
    numerator=np.sum((pred-obs)**2,axis=axes)
    denominator=np.sum(obs**2,axis=axes)
    rel=np.sqrt(numerator/denominator)
    return dict(radiusNormalizedRelativeRms=float(np.sqrt(np.mean(rel**2))),
                perRadiusRelativeL2=rel.tolist(),maxRadiusRelativeL2=float(max(rel)))


def source_values(p,R,theta):
    RR,TT=np.meshgrid(R,theta,indexing="ij")
    x=RR*np.cos(TT)/p["ell"];y=RR*np.sin(TT)/p["ell"]
    if x.min()<p["x"][0] or x.max()>p["x"][-1] or y.min()<0 or y.max()>p["y"][-1]:
        raise ValueError("No extrapolation beyond saved source")
    spl={k:RectBivariateSpline(p["y"],p["x"],p[k]) for k in ("psi","u1","u2","q")}
    psi=p["P"]*p["ell"]**2*spl["psi"].ev(y,x)
    u1=p["P"]*p["ell"]*spl["u1"].ev(y,x)
    u2=p["P"]*p["ell"]*spl["u2"].ev(y,x)
    ur=u1*np.cos(TT)+u2*np.sin(TT)
    ut=-u1*np.sin(TT)+u2*np.cos(TT)
    q=p["P"]*spl["q"].ev(y,x)
    sx=p["P"]*p["ell"]*spl["psi"].ev(y,x,dy=1)
    sy=p["P"]*p["ell"]*spl["psi"].ev(y,x,dx=1)
    qbase=-p["P"]*(spl["psi"].ev(y,x,dy=2)+spl["psi"].ev(y,x,dx=2))
    urbase=-sy*np.cos(TT)+sx*np.sin(TT)
    utbase=sy*np.sin(TT)+sx*np.cos(TT)
    return dict(psi=psi,ur=ur,ut=ut,velocity=np.stack([ur,ut],axis=1),q=q,
                baseline=dict(ur=urbase,ut=utbase,velocity=np.stack([urbase,utbase],axis=1),q=qbase))


def main():
    original=json.loads((OUT/"polar_geometry.json").read_text())
    p=load_native(20)
    records=[]
    with np.load(OUT/"polar_geometry.npz") as z:
        R=z["nativeRadii"];theta=z["thetaRadians"];train=z["trainingRadiusMask"]
        allpsi=z["wall_saddle_tau20_psi"]
        # Omit eight angular midpoint samples on each side (about 3.6 deg).
        interior_theta=theta[8:-8]
        for band in ("middle","outer"):
            meta=original["centers"]["wall_saddle"]["bands"][band]
            lo,hi=meta["nativeRadiusRange"]
            band_idx=np.flatnonzero((R>=lo)&(R<=hi))
            train_idx=band_idx[train[band_idx]]
            held_idx=band_idx[2:-2]
            held_idx=held_idx[~train[held_idx]]
            heldR=R[held_idx]
            obs=source_values(p,heldR,interior_theta)
            degree=meta["amplitudePower"]["exponent"]
            basis=z[f"wall_saddle_{band}_tau20AngularBasis"]
            angular=CubicSpline(theta,basis,axis=1,extrapolate=False)
            v,v1,v2=(angular(interior_theta,nu=i) for i in (0,1,2))
            coefficients=allpsi[train_idx]@basis.T
            scaled=coefficients/R[train_idx,None]**degree
            radial=CubicSpline(np.log(R[train_idx]),scaled,axis=0,extrapolate=False)
            c,c1,c2=(radial(np.log(heldR),nu=i) for i in (0,1,2))
            models=[]
            for rank in (1,2):
                cc,dd,ee=c[:,:rank],c1[:,:rank],c2[:,:rank]
                vv,dv,ev=v[:rank],v1[:rank],v2[:rank]
                psi=heldR[:,None]**degree*(cc@vv)
                ur=-heldR[:,None]**(degree-1)*(cc@dv)
                ut=heldR[:,None]**(degree-1)*((degree*cc+dd)@vv)
                q=-heldR[:,None]**(degree-2)*((degree**2*cc+2*degree*dd+ee)@vv+cc@ev)
                model=dict(rank=rank,psi=error_summary(psi,obs["psi"]),
                    radialVelocity=error_summary(ur,obs["ur"]),
                    angularVelocity=error_summary(ut,obs["ut"]),
                    velocityVector=error_summary(np.stack([ur,ut],axis=1),obs["velocity"]),
                    q=error_summary(q,obs["q"]),
                    qSign=dict(predictedMinimum=float(q.min()),savedMinimum=float(obs["q"].min()),
                        predictedNegativeSamples=int(np.count_nonzero(q<0)),savedNegativeSamples=int(np.count_nonzero(obs["q"]<0)),
                        sampleCount=q.size,perRadiusPredictedMinimum=q.min(axis=1).tolist(),perRadiusSavedMinimum=obs["q"].min(axis=1).tolist()),
                    qCancellation=dict(laplacianRadialTermRms=float(np.sqrt(np.mean((heldR[:,None]**(degree-2)*((degree**2*cc+2*degree*dd+ee)@vv))**2))),
                                       laplacianAngularTermRms=float(np.sqrt(np.mean((heldR[:,None]**(degree-2)*(cc@ev))**2)))))
                models.append(model)
            baseline=dict(velocityVector=error_summary(obs["baseline"]["velocity"],obs["velocity"]),
                          q=error_summary(obs["baseline"]["q"],obs["q"]))
            records.append(dict(band=band,nativeRadiusRange=[lo,hi],trainingRadii=R[train_idx].tolist(),
                heldoutRadii=heldR.tolist(),thetaDegreesRange=(interior_theta[[0,-1]]*180/np.pi).tolist(),
                amplitudePowerUsedOnlyForConditioning=degree,sourcePsiSplineDerivativeBaseline=baseline,models=models))
    result=dict(title="Derivative audit of direct-psi polar rank 1/2 on held-out radii",tau=20,
        source=p["source"],protocol=dict(center="Native centred wall saddle X=a,Y=0; no new centre or mother fit.",
            angularBasis="Use the existing angular SVD modes trained on alternating tau20 radii; cubic-spline theta derivatives, no angular refit.",
            radialFit="Project native centred Psi on the frozen modes at training radii only. Interpolate a_j(R)/R^p by not-a-knot cubic spline in log R; p is the already-fitted amplitude power and does not constrain a_j to powers.",
            holdout="Use the same held-out radii, excluding the first and last two radial samples of each band. Omit eight theta samples at either wall endpoint. No extrapolation.",
            nativeScaling="Psi_native=P ell^2 psi_inner; centred native velocities=P ell u_inner; q_native=P q_inner.",
            equations="u_r=-Psi_theta/R; u_theta=Psi_R; q=-(Psi_RR+Psi_R/R+Psi_thetatheta/R^2).",
            errors="Within each held-out half-circle, take uniform-theta L2 relative to saved field; then equal-radius RMS of these relative errors. The vector velocity combines both components before normalization. No pointwise division by Psi or q.",
            caveat="A radial spline is a descriptive interpolant, not a closed radial evolution equation. Small value error alone does not establish useful Poisson closure."),bands=records)
    (OUT/"polar_derivative.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    lines=["# 直接流函数极坐标骨架的导数审核","",
        "只审核已有 τ20 中区、外区的 1/2 个角向 SVD 模态，不增新模型。角模态与径向系数都只使用原有训练半径；径向 `a_j(R)/R^p` 在 `log R` 上作三次样条。p 仅用于改善插值尺度，并不强迫径向系数为幂律。", "",
        "留出半径去掉各区间首尾两径点；角向去掉两端各八点。误差先取每个半圆的均匀角向相对 L2，再对半径等权合并。真值来自保存的中心化速度和 q；原始 ψ 样条求导给出独立数值基线。", "",
        "| 区间 | 秩 | ψ 误差 | 速度向量误差 | q 误差 | q 负值点 / 总点数 |", "|---|---:|---:|---:|---:|---:|"]
    for row in records:
        for m in row["models"]:
            neg=m["qSign"]
            lines.append(f"| {row['band']} | {m['rank']} | {100*m['psi']['radiusNormalizedRelativeRms']:.3f}% | {100*m['velocityVector']['radiusNormalizedRelativeRms']:.3f}% | {100*m['q']['radiusNormalizedRelativeRms']:.3f}% | {neg['predictedNegativeSamples']} / {neg['sampleCount']} |")
        b=row["sourcePsiSplineDerivativeBaseline"]
        lines.append(f"| {row['band']} 源 ψ 导数基线 | — | — | {100*b['velocityVector']['radiusNormalizedRelativeRms']:.5f}% | {100*b['q']['radiusNormalizedRelativeRms']:.3f}% | — |")
    lines += ["", "公式为 `u_r=−ψ_θ/R`、`u_θ=ψ_R`、`q=−(ψ_RR+ψ_R/R+ψ_θθ/R²)`。归一化严格恢复原生单位：`ψ_native=Pℓ²ψ_inner`、`u_native=Pℓu_inner`、`q_native=Pq_inner`。", "",
        "该检查区分 ψ 的低秩值逼近、速度的一阶导数逼近和 q 的二阶相容性。低秩与径向插值本身不构成 PDE 降维闭合；完整逐半径误差、q 最小值及训练留出点见 `polar_derivative.json`。"]
    (OUT/"POLAR_DERIVATIVE.md").write_text("\n".join(lines)+"\n")
    for row in records:
        print(row["band"],"baseline",row["sourcePsiSplineDerivativeBaseline"])
        for m in row["models"]:
            print("rank",m["rank"],"psi",m["psi"]["radiusNormalizedRelativeRms"],"velocity",m["velocityVector"]["radiusNormalizedRelativeRms"],"q",m["q"]["radiusNormalizedRelativeRms"],"qsign",m["qSign"])


if __name__=="__main__":
    main()
