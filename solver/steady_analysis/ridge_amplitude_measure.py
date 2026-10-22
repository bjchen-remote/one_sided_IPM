"""Native-coordinate ridge position and peak-amplitude curves for saved frames."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator, RectBivariateSpline
from scipy.io import loadmat
from scipy.optimize import brentq, minimize_scalar

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/"runs/geometry_loop_20260928"
OUT=ROOT/"runs/ridge_amplitude_20260928"
YGRID=np.r_[0.,np.geomspace(.025,1024.,121)]
REPRESENTATIVE=np.array([0,.025,.1,.25,.5,1,2,4,8,16,32,64,128,256,512,1024.])
THRESHOLDS=np.array([.75,.5,.25,.125,.0625,.03125])
EXPONENTS=(.5,.6,.625,2/3)


def load_frame(tau):
    d=loadmat(DATA/f"frame_tau{tau}.mat",simplify_cells=True)
    return dict(tau=float(d["tau"]),a=float(d["a"]),P=float(d["P"]),ell=float(d["ell"]),
        x=np.asarray(d["x"],float),y=np.asarray(d["y"],float),
        q=np.asarray(d["q"],float),Q=np.asarray(d["Q"],float))


def peak_row(frame,qfield,y_inner):
    x=frame["x"]
    qrow=qfield(float(y_inner),x)[0]
    cs=CubicSpline(x,qrow)
    roots=cs.derivative().roots(extrapolate=False)
    roots=roots[np.isfinite(roots)&(roots>=x[0])&(roots<=x[-1])]
    points=np.r_[x[0],roots,x[-1]]
    m=float(points[np.argmax(cs(points))]); Ainner=float(cs(m))
    return m,Ainner,float(cs(m,2))


def local_log_slopes(Y,A,window):
    out=[None]*len(Y)
    pos=np.flatnonzero((Y>0)&(A>0))
    half=window//2
    for ii,i in enumerate(pos):
        lo=max(0,ii-half); hi=min(len(pos),ii+half+1)
        if hi-lo<window:
            if lo==0: hi=min(len(pos),window)
            else: lo=max(0,len(pos)-window)
        ids=pos[lo:hi]
        if len(ids)<2: continue
        slope=np.polyfit(np.log(Y[ids]),np.log(A[ids]),1)[0]
        out[i]=float(-slope)
    return out


def wall_decay(Y,A,window):
    # Fit a low-degree polynomial in Y through the wall and the first positive samples.
    n=min(window,len(Y)); yy=Y[:n]; aa=A[:n]
    scale=max(float(yy[-1]),1e-12)
    deg=min(3,n-1)
    coef=np.polynomial.polynomial.polyfit(yy/scale,aa,deg)
    deriv=float(coef[1]/scale)
    return dict(points=n,degree=deg,YRange=[float(yy[0]),float(yy[-1])],
        APrime0=deriv,negativeLogDerivative0=float(-deriv/A[0]),
        effectiveLength0=float(-A[0]/deriv) if deriv<0 else None)


def threshold_crossings(Y,ratio,mabs,a):
    ap=PchipInterpolator(Y,ratio,extrapolate=False)
    mp=PchipInterpolator(Y,mabs,extrapolate=False)
    rows=[]
    for level in THRESHOLDS:
        found=None
        for i in range(len(Y)-1):
            if (ratio[i]-level)==0:
                yroot=float(Y[i]); found=yroot; break
            if (ratio[i]-level)*(ratio[i+1]-level)<0:
                found=float(brentq(lambda z:float(ap(z)-level),Y[i],Y[i+1])); break
        rows.append(dict(fraction=float(level),Ynative=found,
                         mAbsoluteNative=float(mp(found)) if found is not None else None,
                         mCenteredNative=float(mp(found)-a) if found is not None else None,
                         status="crossed_in_saved_interval" if found is not None else "not_reached_before_saved_boundary"))
    return rows


def model_metrics(Y,ratio,p,L,trainmask,holdmask):
    pred=(1+Y/L)**(-p)
    def score(mask):
        obs=ratio[mask]; fit=pred[mask]
        return dict(count=int(mask.sum()),relativeL2=float(np.linalg.norm(fit-obs)/np.linalg.norm(obs)),
                    logRMSE=float(np.sqrt(np.mean((np.log(fit)-np.log(obs))**2))))
    return pred,dict(train=score(trainmask),holdout=score(holdmask))


def fit_amplitude_models(Y,ratio,wall_decay_slope):
    train=(Y>0)&(Y<=16); hold=(Y>16)&(Y<=256)
    results=[]
    for p in EXPONENTS:
        Lwall=float(p/wall_decay_slope) if wall_decay_slope>0 else None
        wall=None
        if Lwall is not None:
            pred,sc=model_metrics(Y,ratio,p,Lwall,train,hold)
            wall=dict(L=Lwall,**sc)
        def objective(logL):
            L=np.exp(logL); pred=(1+Y[train]/L)**(-p)
            return float(np.mean((np.log(pred)-np.log(ratio[train]))**2))
        fit=minimize_scalar(objective,bounds=(np.log(1e-4),np.log(1e8)),method="bounded")
        Ltrain=float(np.exp(fit.x))
        pred,sc=model_metrics(Y,ratio,p,Ltrain,train,hold)
        results.append(dict(p=float(p),wallSlopeFixedL=wall,
                            oneParameterTrainFit=dict(L=Ltrain,**sc)))
    return dict(trainNativeY=[float(Y[train].min()),float(Y[train].max())],
                holdoutNativeY=[float(Y[hold].min()),float(Y[hold].max())],
                exponents=results)


def analyze_frame(tau):
    f=load_frame(tau)
    qfield=RectBivariateSpline(f["y"],f["x"],f["q"],kx=3,ky=3)
    records=[]
    for Y in YGRID:
        y=float(Y/f["ell"])
        if y<f["y"][0] or y>f["y"][-1]:
            records.append(dict(Ynative=float(Y),available=False,reason="outside saved y domain; no extrapolation")); continue
        m,Ainner,qxx=peak_row(f,qfield,y)
        A=f["P"]*Ainner
        curv=f["ell"]*np.sqrt(-2*Ainner/qxx) if qxx<0 and Ainner>0 else None
        records.append(dict(Ynative=float(Y),yInner=y,mInner=m,
            mCenteredNative=float(f["ell"]*m),mAbsoluteNative=float(f["a"]+f["ell"]*m),
            AInner=Ainner,A_native=float(A),qxxInner=qxx,
            curvatureWidthNative=float(curv) if curv is not None else None,available=True))
    good=[r for r in records if r["available"]]
    Y=np.array([r["Ynative"] for r in good]); A=np.array([r["A_native"] for r in good])
    mabs=np.array([r["mAbsoluteNative"] for r in good]); mcenter=np.array([r["mCenteredNative"] for r in good])
    A0=float(A[0]); ratio=A/A0
    for r,a,rat in zip(good,A,ratio):
        r["AOverA0Native"]=float(rat)
        r["AInnerOverA0Inner"]=float(r["AInner"]/good[0]["AInner"])
        r["AOverA0NativeCheckEqualInnerRatio"]=float(r["AInnerOverA0Inner"]-rat)
    # Local log slope of the amplitude and corresponding local effective length Y/p_A.
    slope_sets={str(n):local_log_slopes(Y,A,n) for n in (3,5,9)}
    for i,r in enumerate(good):
        r["pA"]={key:values[i] for key,values in slope_sets.items()}
        r["effectiveLengthFromLocalSlope"]={key:(float(Y[i]/values[i]) if Y[i]>0 and values[i] and values[i]>0 else None)
                                               for key,values in slope_sets.items()}
    wall={str(n):wall_decay(Y,A,n) for n in (3,5,9)}
    primary_wall=wall["5"]["negativeLogDerivative0"]
    crossings=threshold_crossings(Y,ratio,mabs,f["a"])
    models=fit_amplitude_models(Y,ratio,primary_wall)
    representative=[]
    base_lookup={r["Ynative"]:r for r in records if r["available"]}
    for Yrep in REPRESENTATIVE:
        if float(Yrep) in base_lookup:
            representative.append(base_lookup[float(Yrep)])
            continue
        y=float(Yrep/f["ell"])
        if y<f["y"][0] or y>f["y"][-1]:
            representative.append(dict(Ynative=float(Yrep),available=False,reason="outside saved y domain; no extrapolation")); continue
        m,Ainner,qxx=peak_row(f,qfield,y); Arow=f["P"]*Ainner
        positive_ids=np.flatnonzero(Y>0)
        near=positive_ids[np.argsort(abs(np.log(Y[positive_ids]/Yrep)))[:4]]
        fitY=np.r_[Yrep,Y[near]]; fitA=np.r_[Arow,A[near]]
        pa5=float(-np.polyfit(np.log(fitY),np.log(fitA),1)[0])
        representative.append(dict(Ynative=float(Yrep),yInner=y,mInner=m,
            mCenteredNative=float(f["ell"]*m),mAbsoluteNative=float(f["a"]+f["ell"]*m),
            AInner=Ainner,A_native=float(Arow),AOverA0Native=float(Arow/A0),
            AInnerOverA0Inner=float(Ainner/good[0]["AInner"]),qxxInner=qxx,
            curvatureWidthNative=float(f["ell"]*np.sqrt(-2*Ainner/qxx)) if qxx<0 and Ainner>0 else None,
            pA={"5":pa5},effectiveLengthFromLocalSlope={"5":float(Yrep/pa5) if pa5>0 else None},
            available=True,extraRepresentativeSample=True))
    # Preserve compact m/A curves and smooth local slopes in NPZ.
    return dict(tau=f["tau"],a=f["a"],P=f["P"],ell=f["ell"],
        nativeBounds=dict(X=[float(f["a"]+f["ell"]*f["x"][0]),float(f["a"]+f["ell"]*f["x"][-1])],
                          Y=[float(f["ell"]*f["y"][0]),float(f["ell"]*f["y"][-1])]),
        A0Native=A0,A0Inner=float(good[0]["AInner"]),
        wallDerivativeSensitivity=wall,primaryWallDerivativeWindow=5,
        profileRows=records,representativeRows=representative,
        thresholdCrossings=crossings,localAmplitudeModels=models,
        _arrays=dict(Ynative=Y,A_native=A,AOverA0=ratio,mAbsoluteNative=mabs,mCenteredNative=mcenter,
            qxxInner=np.array([r["qxxInner"] for r in good]),
            pA3=np.array([np.nan if r["pA"]["3"] is None else r["pA"]["3"] for r in good]),
            pA5=np.array([np.nan if r["pA"]["5"] is None else r["pA"]["5"] for r in good]),
            pA9=np.array([np.nan if r["pA"]["9"] is None else r["pA"]["9"] for r in good])))


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    analyzed=[analyze_frame(t) for t in (10,15,20)]
    arrays={"requestedNativeY":YGRID,"representativeNativeY":REPRESENTATIVE}
    frames=[]
    for f in analyzed:
        arrays.update({f"tau{int(round(f['tau']))}_{k}":v for k,v in f.pop("_arrays").items()})
        frames.append(f)
    # Join all three time frames on the exact same native-height sampling.
    time=[]
    for i,Y in enumerate(YGRID):
        time.append(dict(Ynative=float(Y),byTau={str(int(round(f["tau"]))):
            {k:f["profileRows"][i].get(k) for k in ("mAbsoluteNative","mCenteredNative","A_native","AOverA0Native","pA")}
            for f in frames}))
    representative_time=[]
    for Y in REPRESENTATIVE:
        by_tau={str(int(round(f["tau"]))):next(r for r in f["representativeRows"] if r["Ynative"]==float(Y)) for f in frames}
        entry=dict(Ynative=float(Y),byTau={t:{k:r.get(k) for k in
            ("mAbsoluteNative","mCenteredNative","A_native","AOverA0Native")} for t,r in by_tau.items()})
        if Y>0:
            entry["AChangeTau20VsTau10Percent"]=float(100*(by_tau["20"]["A_native"]/by_tau["10"]["A_native"]-1))
        representative_time.append(entry)
    result=dict(title="Ridge position and peak amplitude evolution in native coordinates",
        definitions=["X=a+ell*x; Y=ell*y; rho_X=P*q", "A(Y)=max_X rho_X(X,Y)",
            "mAbsoluteNative=a+ell*mInner; mCenteredNative=ell*mInner",
            "A/A0 is reported in native and inner amplitudes; P cancels in the ratio.",
            "Threshold crossings solve A(Y)/A(0)=fraction with a shape-preserving PCHIP interpolant of the 122 saved curve points; they are interpolated curve roots, not direct continuous-field roots. No root is extrapolated beyond the sampled interval."],
        sampling=dict(nativeY=[0.0]+[float(x) for x in np.geomspace(.025,1024,121)],
            fixedRepresentativeY=REPRESENTATIVE.tolist(),noExtrapolation=True),
        timeFrames=frames,commonNativeYTimeEvolution=time,
        representativeNativeYTimeEvolution=representative_time,
        interpretation="The fitted power profiles are diagnostic candidates on prescribed exponents and finite height ranges; no asymptotic exponent is selected.")
    (OUT/"amplitude_detail.json").write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    np.savez_compressed(OUT/"amplitude_detail.npz",**arrays)
    write_markdown(result)
    print(json.dumps({"out":str(OUT),"frames":[dict(tau=f["tau"],A0=f["A0Native"],wallSlope=f["wallDerivativeSensitivity"]["5"]["negativeLogDerivative0"],
        thresholdY={str(r["fraction"]):r["Ynative"] for r in f["thresholdCrossings"]}) for f in frames]},indent=2))


def write_markdown(result):
    lines=["# Ridge position and peak amplitude","",
        "This report isolates the two ridge curves: `m(Y)` and `A(Y)=max_X rho_X(X,Y)`. Each frame uses its own native restoration `X=a+ell*x`, `Y=ell*y`, and `rho_X=P*q`. Absolute `m` and centered `m-a` are both retained.","",
        "## Representative heights","","| tau | Y | m absolute X | m-a | A native | A/A0 | pA (5 point) | effective length Y/pA |","|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for f in result["timeFrames"]:
        for r in f["representativeRows"]:
            if r["Ynative"]==0: continue
            pa=r["pA"]["5"]
            ell=r["effectiveLengthFromLocalSlope"]["5"]
            lines.append(f"| {f['tau']:.6g} | {r['Ynative']:.6g} | {r['mAbsoluteNative']:.6g} | {r['mCenteredNative']:.6g} | {r['A_native']:.6g} | {r['AOverA0Native']:.6g} | {pa:.6g} | {'—' if ell is None else f'{ell:.6g}'} |")
    lines += ["","## Wall amplitude loss and level crossings","","Wall derivative uses polynomial fits through Y=0 and the first 3, 5, or 9 saved height samples. The 5-point cubic estimate is used to set wall-slope model lengths; all window estimates are included in JSON.","",
        "| tau | A0 | -A'(0)/A0, 3 pts | 5 pts | 9 pts | Y at 3/4 | 1/2 | 1/4 | 1/8 | 1/16 | 1/32 |","|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for f in result["timeFrames"]:
        ycross={r["fraction"]:r["Ynative"] for r in f["thresholdCrossings"]}
        vals=[f["wallDerivativeSensitivity"][str(n)]["negativeLogDerivative0"] for n in (3,5,9)]
        lines.append("| "+" | ".join([f"{f['tau']:.6g}",f"{f['A0Native']:.6g}",*[f"{v:.6g}" for v in vals],*["—" if ycross[p] is None else f"{ycross[p]:.6g}" for p in THRESHOLDS]])+" |")
    lines += ["","## Peak-amplitude drift at fixed native heights","","The peak amplitude falls substantially between the saved frames even when each profile is normalized by its own A0. Relative to tau=10, tau=20 A is lower by the following amounts:","",
        "| Y | A(tau10) | A(tau15) | A(tau20) | tau20 change from tau10 | m(tau10) | m(tau20) |","|---:|---:|---:|---:|---:|---:|---:|"]
    for Y in (4.,16.,64.,256.):
        item=next(r for r in result["representativeNativeYTimeEvolution"] if r["Ynative"]==Y)
        rr=item["byTau"]; change=item["AChangeTau20VsTau10Percent"]
        lines.append(f"| {Y:g} | {rr['10']['A_native']:.6g} | {rr['15']['A_native']:.6g} | {rr['20']['A_native']:.6g} | {change:.2f}% | {rr['10']['mAbsoluteNative']:.6g} | {rr['20']['mAbsoluteNative']:.6g} |")
    lines += ["","## Prescribed power profiles","","The family `A/A0=(1+Y/L)^(-p)` is evaluated at p=1/2, 3/5, 5/8, and 2/3. For each p, JSON reports (a) L fixed by the 5-point wall slope and (b) a one-parameter L fit using Y<=16. Both are scored on the training range Y<=16 and held out on 16<Y<=256. Farther heights are shown as data and do not choose an exponent.","",
        "The JSON provides local `pA=-d log(A)/d log(Y)` from centered 3-, 5-, and 9-point log-height windows and `-A/A'` via the corresponding effective length Y/pA. Window variation is explicit; these are finite-time diagnostics, not selected asymptotic laws.",""]
    (OUT/"AMPLITUDE_DETAIL.md").write_text("\n".join(lines))


if __name__=="__main__": main()
