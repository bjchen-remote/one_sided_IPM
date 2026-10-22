"""Discrete beta sensitivity of the sinh half-width profile alignment."""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline
from scipy.io import loadmat

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/"runs/geometry_loop_20260928"
OUT=ROOT/"runs/profile_geometry_20260928"
BETAS=(1.,1.25,1.5,1.75,2.)
HEIGHTS=np.array([4.,8.,16.,32.,64.,128.,256.])
XI=np.linspace(-1.5,1.5,601)


def frame(tau):
    d=loadmat(DATA/f"frame_tau{tau}.mat",simplify_cells=True)
    return dict(tau=float(d["tau"]),x=np.asarray(d["x"],float),y=np.asarray(d["y"],float),
                q=np.asarray(d["q"],float),Q=np.asarray(d["Q"],float),
                a=float(d["a"]),P=float(d["P"]),ell=float(d["ell"]))


def geometry(qfield,Qfield,x,y):
    qrow=qfield(float(y),x)[0]
    cs=CubicSpline(x,qrow)
    roots=cs.derivative().roots(extrapolate=False)
    roots=roots[(roots>=x[0])&(roots<=x[-1])]
    pts=np.r_[x[0],roots,x[-1]]
    m=float(pts[np.argmax(cs(pts))]); A=float(cs(m))
    hr=CubicSpline(x,qrow-A/2).roots(extrapolate=False)
    left,right=hr[hr<m],hr[hr>m]
    if not len(left) or not len(right):
        raise ValueError(f"half-height crossings unavailable at inner y={y:g}")
    xL,xR=float(left[-1]),float(right[0])
    L=m-xL; R=xR-m; H=(L+R)/2; epsilon=(R-L)/(R+L)
    return dict(m=m,A=A,L=L,R=R,H=H,epsilon=epsilon,Qpeak=float(Qfield(float(y),m)[0,0]))


def run_source(tau):
    p=frame(tau)
    qfield=RectBivariateSpline(p["y"],p["x"],p["q"],kx=3,ky=3)
    Qfield=RectBivariateSpline(p["y"],p["x"],p["Q"],kx=3,ky=3)
    base=[]
    for Y in HEIGHTS:
        y=Y/p["ell"]
        if y<p["y"][0] or y>p["y"][-1]:
            base.append(None); continue
        base.append(geometry(qfield,Qfield,p["x"],y))
    result={"tau":p["tau"],"a":p["a"],"P":p["P"],"ell":p["ell"],"rows":{}}
    for beta in BETAS:
        profiles=[]; rowdata=[]; valid_rows=[]
        for Y,g in zip(HEIGHTS,base):
            if g is None:
                rowdata.append(dict(Ynative=float(Y),available=False,reason="height outside saved y interval; no extrapolation")); continue
            k=g["epsilon"]/np.tanh(beta/2)
            W=g["H"]*beta/np.sinh(beta)
            G=(np.sinh(beta*XI)+k*(np.cosh(beta*XI)-1))/beta
            xmap=g["m"]+W*G
            valid=(xmap>=p["x"][0])&(xmap<=p["x"][-1])
            qshape=np.full_like(XI,np.nan)
            yy=np.full(int(valid.sum()),float(Y/p["ell"]))
            qshape[valid]=qfield(yy,xmap[valid],grid=False)/g["A"]
            valid_rows.append(valid)
            profiles.append(qshape)
            rowdata.append(dict(Ynative=float(Y),yInner=float(Y/p["ell"]),mInner=g["m"],
                AInner=g["A"],HInner=g["H"],epsilon=g["epsilon"],WInner=float(W),
                Wnative=float(W*p["ell"]),k=float(k),kAbsLess1=bool(abs(k)<1),
                mappedXInnerBounds=[float(xmap.min()),float(xmap.max())],
                measuredXDomainInner=[float(p["x"][0]),float(p["x"][-1])],
                xiCoverage=[float(XI[valid][0]),float(XI[valid][-1])] if valid.any() else None))
        common=np.logical_and.reduce(valid_rows) if valid_rows else np.zeros_like(XI,dtype=bool)
        if common.sum()<10:
            shape_spread=None
        else:
            mat=np.asarray(profiles)[:,common]
            master=mat.mean(axis=0)
            shape_spread=dict(relativeL2Spread=float(np.linalg.norm(mat-master[None,:])/np.linalg.norm(mat)),
                rmsAcrossHeight=float(np.sqrt(np.mean((mat-master[None,:])**2))),
                maxPointwiseStd=float(mat.std(axis=0).max()),
                xiCommonSupport=[float(XI[common][0]),float(XI[common][-1])],
                xiSampleCount=int(common.sum()),heightCount=int(mat.shape[0]),
                normalizedProfileMean=master.tolist())
        ks=[r["k"] for r in rowdata if r.get("available",True) and "k" in r]
        result["rows"][str(beta)]=dict(beta=float(beta),
            allRequestedHeightsAvailable=all(r.get("available",True) for r in rowdata),
            allAbsKLess1=(len(ks)==len(HEIGHTS) and all(abs(k)<1 for k in ks)),
            absKLess1Count=int(sum(abs(k)<1 for k in ks)),totalHeightCount=int(len(HEIGHTS)),
            kRange=[float(min(ks)),float(max(ks))] if ks else None,
            rows=rowdata,shapeSpreadAcrossHeights=shape_spread)
    return result


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    sources={str(t):run_source(t) for t in (10,15,20)}
    tau20=sources["20"]
    temporal=[]
    for beta in BETAS:
        bkey=str(beta)
        temporal.append(dict(beta=beta,byTau={t:dict(allAbsKLess1=s["rows"][bkey]["allAbsKLess1"],
            validCount=s["rows"][bkey]["absKLess1Count"],kRange=s["rows"][bkey]["kRange"])
            for t,s in sources.items()}))
    report=dict(title="Discrete sinh-map beta sensitivity on saved near/intermediate-height profiles",
        scope="Data-only coordinate sensitivity; no PDE, continuous beta optimization, or comparison to an independently fixed f0.",
        convention="At each height, peak m and q amplitude A are measured first. H=(L+R)/2 and epsilon=(R-L)/(R+L) come from measured half-widths; for each listed beta, W=H*beta/sinh(beta), k=epsilon/tanh(beta/2), and G=(sinh(beta*xi)+k(cosh(beta*xi)-1))/beta.",
        shapeMetric="q/A is evaluated only at mapped points inside each saved x interval. Cross-height relative L2 spread is computed against the equal-height mean on the intersection of all measured supports within |xi|<=1.5.",
        nativeHeights=HEIGHTS.tolist(),betas=list(BETAS),tau20HeightResults=tau20["rows"],
        timeValidityAtFixedNativeHeights=temporal,
        sourceDetails={t:{"tau":s["tau"],"a":s["a"],"P":s["P"],"ell":s["ell"],"betaResults":s["rows"]} for t,s in sources.items()},
        interpretation="A lower cross-height spread indicates a more consistent coordinate alignment over the chosen heights. It does not establish an exact profile or make beta identifiable independently of the alignment; validity only checks the stated |k|<1 monotonicity condition.")
    (OUT/"beta_sensitivity.json").write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    lines=["# Discrete beta sensitivity","",
        "The same measured peak and half-height widths are remapped for beta in {1, 1.25, 1.5, 1.75, 2}. For each beta, `W=H beta/sinh(beta)` and `k=epsilon/tanh(beta/2)`. Profiles are compared only where every requested height has measured x support, with no extrapolation.","",
        "## tau=20: cross-height alignment at native Y=4 through 256","",
        "| beta | all |k|<1 | k range | common xi interval | q/A relative L2 spread |","|---:|:---:|---:|---:|---:|"]
    for beta in BETAS:
        r=tau20["rows"][str(beta)]; s=r["shapeSpreadAcrossHeights"]
        lines.append(f"| {beta:g} | {r['allAbsKLess1']} | {r['kRange'][0]:.6g} to {r['kRange'][1]:.6g} | {s['xiCommonSupport'][0]:.4g} to {s['xiCommonSupport'][1]:.4g} | {s['relativeL2Spread']:.6g} |")
    lines += ["","## Monotonicity validity at common native heights","","| beta | tau=10 | tau=15 | tau=20 |","|---:|:---:|:---:|:---:|"]
    for r in temporal:
        lines.append(f"| {r['beta']:g} | {r['byTau']['10']['validCount']}/7 | {r['byTau']['15']['validCount']}/7 | {r['byTau']['20']['validCount']}/7 |")
    lines += ["","This is a discrete coordinate-sensitivity check. A smaller spread favors an alignment over this finite set of heights, but does not prove exact shape collapse or select a continuously optimal beta. No fixed profile such as f0 is used in this comparison.",""]
    (OUT/"BETA_SENSITIVITY.md").write_text("\n".join(lines))
    print(json.dumps({"out":str(OUT),"tau20":{b:{"allAbsKLess1":tau20["rows"][str(b)]["allAbsKLess1"],
        "kRange":tau20["rows"][str(b)]["kRange"],"spread":tau20["rows"][str(b)]["shapeSpreadAcrossHeights"]["relativeL2Spread"],
        "xi":tau20["rows"][str(b)]["shapeSpreadAcrossHeights"]["xiCommonSupport"]} for b in BETAS},
        "timeValidity":temporal},indent=2))


if __name__=="__main__": main()
