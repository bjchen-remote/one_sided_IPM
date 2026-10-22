"""Detailed native/inner profile geometry measurements for tau10/15/20 exports."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "runs/geometry_loop_20260928"
OUT = ROOT / "runs/profile_geometry_20260928"
NATIVE_Y = np.array([0,.25,.5,1,2,4,8,16,32,64,128,256,512,1024.], float)
COMMON_INNER_Y = np.array([0,.25,.5,1,2,4,8,16,32,64,128,256.], float)
BETA = 1.5
LEVELS = (0.25,0.5,0.75)


def load_frame(tau):
    d=loadmat(SOURCE/f"frame_tau{tau}.mat",simplify_cells=True)
    x=np.asarray(d["x"],float); y=np.asarray(d["y"],float)
    a,P,ell=map(float,(d["a"],d["P"],d["ell"]))
    return dict(tau=float(d["tau"]),a=a,P=P,ell=ell,x=x,y=y,
                X=a+ell*x,Y=ell*y,Q=np.asarray(d["Q"],float),q=np.asarray(d["q"],float),
                source=f"frame_tau{tau}.mat",coordinateNote=str(d.get("coordinateNote","")))


def measure_profile(frame,y_inner,qfield,Qfield):
    x=frame["x"]
    qrow=qfield(float(y_inner),x)[0]
    cs=CubicSpline(x,qrow)
    roots=cs.derivative().roots(extrapolate=False)
    roots=roots[np.isfinite(roots)&(roots>=x[0])&(roots<=x[-1])]
    points=np.r_[x[0],roots,x[-1]]
    m=float(points[np.argmax(cs(points))]); A=float(cs(m)); qxx=float(cs(m,2))
    wcurv=float(np.sqrt(-2*A/qxx)) if qxx<0 and A>0 else None
    Qpeak=float(Qfield(float(y_inner),m)[0,0])
    level_info={}
    half_geometry=None
    for fraction in LEVELS:
        lev=CubicSpline(x,qrow-fraction*A).roots(extrapolate=False)
        left=lev[lev<m]; right=lev[lev>m]
        if len(left) and len(right):
            xl=float(left[-1]); xr=float(right[0]); wl=m-xl; wr=xr-m
            midpoint=.5*(xl+xr)
            skew=(wr-wl)/(wr+wl) if wr+wl else None
            level_info[str(fraction)]=dict(xLeftInner=xl,xRightInner=xr,
                widthLeftInner=wl,widthRightInner=wr,midlineShiftInner=midpoint-m,
                midpointInner=midpoint,skew=skew,
                xLeftNative=frame["a"]+frame["ell"]*xl,
                xRightNative=frame["a"]+frame["ell"]*xr,
                widthLeftNative=frame["ell"]*wl,widthRightNative=frame["ell"]*wr,
                midlineShiftNative=frame["ell"]*(midpoint-m))
        else:
            level_info[str(fraction)]=None
        if fraction==.5 and level_info[str(fraction)] is not None:
            half_geometry=level_info[str(fraction)]
    if half_geometry is not None:
        wl=half_geometry["widthLeftInner"]; wr=half_geometry["widthRightInner"]
        asym=(wr-wl)/(wr+wl)
        k=asym/np.tanh(BETA/2)
        W=(wl+wr)*BETA/(2*np.sinh(BETA))
    else:
        asym=k=W=None
    return dict(yInner=float(y_inner),mInner=m,AInner=A,Qpeak=Qpeak,qxxInner=qxx,
                curvatureWidthInner=wcurv,WInner=W,k=k,halfAsymmetry=asym,levels=level_info)


def window_slopes(rows,metric):
    vals=[]
    for row in rows:
        Y=row["Ynative"]; value=row.get(metric)
        vals.append((Y,value))
    vals=[(Y,v) for Y,v in vals if Y>0 and v is not None and np.isfinite(v) and v>0]
    out=[]
    for size in (3,5,7):
        for i in range(max(0,len(vals)-size+1)):
            group=vals[i:i+size]
            x=np.log([z[0] for z in group]); y=np.log([z[1] for z in group])
            slope=float(np.polyfit(x,y,1)[0])
            out.append(dict(windowPoints=size,YRange=[group[0][0],group[-1][0]],
                            geometricCenter=float(np.exp(np.mean(x))),slope=slope))
    return out


def extract_table(frame, coordinate, heights, qfield,Qfield,wall_offset_inner):
    rows=[]
    for h in heights:
        if coordinate=="nativeY":
            Y=float(h); y=float(Y/frame["ell"])
        else:
            y=float(h); Y=frame["ell"]*y
        if y<frame["y"][0] or y>frame["y"][-1]:
            rows.append(dict(requestedHeight=float(h),Ynative=Y,yInner=y,available=False,
                             reason="outside saved Y domain; no extrapolation")); continue
        g=measure_profile(frame,y,qfield,Qfield)
        m,A,Qpeak=g["mInner"],g["AInner"],g["Qpeak"]
        a,P,ell=frame["a"],frame["P"],frame["ell"]
        cpeak=P*ell*Qpeak
        corigin=P*ell*(Qpeak+wall_offset_inner)
        W=g["WInner"]
        row=dict(requestedHeight=float(h),coordinateRequest=coordinate,Ynative=Y,yInner=y,
            available=True,mInner=m,mCenteredNative=ell*m,mAbsoluteNative=a+ell*m,
            mOverYAbsolute=(a+ell*m)/Y if Y>0 else None,
            mOverYCentered=(ell*m)/Y if Y>0 else None,
            AInner=A,A_native=P*A,QpeakRelativeWallPeak=Qpeak,
            c_peak_rhoMinusRho_a0=cpeak,c_origin_rhoMinusRho00=corigin,
            wallOffset_rho_a0MinusRho00=P*ell*wall_offset_inner,
            curvatureWidthInner=g["curvatureWidthInner"],
            curvatureWidthNative=ell*g["curvatureWidthInner"] if g["curvatureWidthInner"] is not None else None,
            WInner=W,Wnative=ell*W if W is not None else None,k=g["k"],
            halfAsymmetry=g["halfAsymmetry"],qxxInner=g["qxxInner"],
            AW_inner=A*W if W is not None else None,
            AW_native=P*ell*A*W if W is not None else None,
            cOriginOverYA=corigin/(Y*P*A) if Y>0 and A else None,
            fourCOriginOverYA=4*corigin/(Y*P*A) if Y>0 and A else None,
            cPeakOverYA=cpeak/(Y*P*A) if Y>0 and A else None,
            levels=g["levels"])
        rows.append(row)
    return rows


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    frames=[]; arrays={"nativeY":NATIVE_Y,"commonInnerY":COMMON_INNER_Y}
    for tau in (10,15,20):
        frame=load_frame(tau)
        qfield=RectBivariateSpline(frame["y"],frame["x"],frame["q"],kx=3,ky=3)
        Qfield=RectBivariateSpline(frame["y"],frame["x"],frame["Q"],kx=3,ky=3)
        wall_q=CubicSpline(frame["x"],qfield(0.,frame["x"])[0])
        xzero=-frame["a"]/frame["ell"]
        if xzero<frame["x"][0] or xzero>frame["x"][-1]:
            raise ValueError(f"tau{tau}: wall origin x=-a/ell is outside saved x interval")
        wall_offset=float(wall_q.integrate(xzero,0.))
        native_rows=extract_table(frame,"nativeY",NATIVE_Y,qfield,Qfield,wall_offset)
        inner_rows=extract_table(frame,"innerY",COMMON_INNER_Y,qfield,Qfield,wall_offset)
        # The requested W/Y slope uses native width. Attach this derived field before slope calculation.
        for row in native_rows:
            row["WnativeOverY"]=(row["Wnative"]/row["Ynative"]
                                  if row.get("Wnative") is not None and row["Ynative"]>0 else None)
        slopes={metric:window_slopes(native_rows,metric) for metric in
                ("mOverYAbsolute","WnativeOverY","k","cOriginOverYA","AW_native")}
        fdict=dict(tau=frame["tau"],source=frame["source"],coordinateNote=frame["coordinateNote"],
            a=frame["a"],P=frame["P"],ell=frame["ell"],
            innerBounds={"x":[float(frame["x"][0]),float(frame["x"][-1])],"y":[float(frame["y"][0]),float(frame["y"][-1])]},
            nativeBounds={"X":[float(frame["X"][0]),float(frame["X"][-1])],"Y":[float(frame["Y"][0]),float(frame["Y"][-1])]},
            wallOrigin=dict(xAtNativeX0=xzero,normalizedIntegralQ=wall_offset,
            rho_a0MinusRho00=frame["P"]*frame["ell"]*wall_offset),
            nativeYProfiles=native_rows,commonInnerYProfiles=inner_rows,
            localLogSlopesNativeY=slopes)
        frames.append(fdict)
        arrays[f"tau{tau}_native_m_abs"]=np.array([r.get("mAbsoluteNative",np.nan) for r in native_rows])
        arrays[f"tau{tau}_native_A"]=np.array([r.get("A_native",np.nan) for r in native_rows])
        arrays[f"tau{tau}_native_W"]=np.array([r.get("Wnative",np.nan) for r in native_rows])
        arrays[f"tau{tau}_native_k"]=np.array([r.get("k",np.nan) for r in native_rows])
        arrays[f"tau{tau}_native_cOrigin"]=np.array([r.get("c_origin_rhoMinusRho00",np.nan) for r in native_rows])
        arrays[f"tau{tau}_native_cPeak"]=np.array([r.get("c_peak_rhoMinusRho_a0",np.nan) for r in native_rows])
        arrays[f"tau{tau}_native_ratio"]=np.array([r.get("cOriginOverYA",np.nan) for r in native_rows])
        arrays[f"tau{tau}_native_fourCOverYA"]=np.array([r.get("fourCOriginOverYA",np.nan) for r in native_rows])
        arrays[f"tau{tau}_native_curvatureWidth"]=np.array([r.get("curvatureWidthNative",np.nan) for r in native_rows])
        arrays[f"tau{tau}_native_AW"]=np.array([r.get("AW_native",np.nan) for r in native_rows])
        for fraction in LEVELS:
            levelkey=str(fraction)
            for outkey,sourcekey in (("wLeft","widthLeftNative"),("wRight","widthRightNative"),
                                     ("midlineShift","midlineShiftNative"),("skew","skew")):
                arrays[f"tau{tau}_native_level{fraction:g}_{outkey}"]=np.array([
                    (r["levels"][levelkey][sourcekey] if r.get("available") and r["levels"][levelkey] else np.nan)
                    for r in native_rows])
        arrays[f"tau{tau}_innerY_m"]=np.array([r.get("mInner",np.nan) for r in inner_rows])
        inner_fields={"AInner":"AInner","WInner":"WInner","k":"k",
            "curvatureWidthInner":"curvatureWidthInner","Qpeak":"QpeakRelativeWallPeak",
            "cOrigin":"c_origin_rhoMinusRho00","cOriginOverYA":"cOriginOverYA"}
        for outkey,sourcekey in inner_fields.items():
            arrays[f"tau{tau}_innerY_{outkey}"]=np.array([
                (r.get(sourcekey,np.nan) if r.get("available") else np.nan) for r in inner_rows])

    # Direct temporal comparison at identical native Y; same native coordinate means each
    # frame is independently restored with its own (a,P,ell).
    compare=[]
    for i,Y in enumerate(NATIVE_Y):
        by_tau={str(int(round(f["tau"]))):f["nativeYProfiles"][i] for f in frames}
        compare.append(dict(Ynative=float(Y),byTau={t:{key:r.get(key) for key in
            ("mOverYAbsolute","WnativeOverY","k","cOriginOverYA","fourCOriginOverYA","AW_native","mAbsoluteNative","A_native","c_origin_rhoMinusRho00")}
            for t,r in by_tau.items()}))
    tau20=next(f for f in frames if round(f["tau"])==20)
    old_path=SOURCE/"quarter_relation_checks.json"
    old_check=json.loads(old_path.read_text()) if old_path.exists() else None
    consistency=[]
    if old_check is not None:
        old20=next(f for f in old_check["frames"] if round(f["tau"])==20)
        for row in tau20["nativeYProfiles"]:
            old=next((z for z in old20["rows"] if z["Y"]==row["Ynative"]),None)
            if old is not None:
                consistency.append(dict(Ynative=row["Ynative"],
                    peakXAbsDifference=row["mAbsoluteNative"]-old["peakX"],
                    amplitudeDifference=row["A_native"]-old["amplitudeOriginal"],
                    cOriginDifference=row["c_origin_rhoMinusRho00"]-old["cRhoPeakMinusRho00"],
                    fourCOverYADifference=row["fourCOriginOverYA"]-old["fourCOverYA"]))
    result=dict(title="Detailed restored-native profile geometry for tau10/15/20",
        coordinateDefinitions={"inner":"saved frame fields use x,y; q=rho_X/P and Q=(rho-rho(a,0))/(P*ell)",
            "native":"X=a+ell*x, Y=ell*y, rho_X=P*q, rho-rho(a,0)=P*ell*Q",
            "peakPositions":"Both mCenteredNative=ell*mInner and mAbsoluteNative=a+ell*mInner are reported.",
            "densityReferences":"c_peak=P*ell*Q(m,y)=rho(m,Y)-rho(a,0); c_origin=c_peak+P*ell*integral(q(x,0), x=-a/ell..0)=rho(m,Y)-rho(0,0). cOriginOverYA uses native A=P*AInner.",
            "noExtrapolation":"Every requested Y is converted to y=Y/ell and checked against the saved y interval; all level crossings are cubic-spline roots restricted to the saved x interval. Missing crossings are null."},
        nativeYRequested=NATIVE_Y.tolist(),commonInnerYRequested=COMMON_INNER_Y.tolist(),frames=frames,
        fixedNativeYTimeEvolution=compare,
        tau20ConsistencyWithQuarterRelationJSON=dict(source="runs/geometry_loop_20260928/quarter_relation_checks.json",
            matchingNativeHeights=len(consistency),maxAbsoluteDifferences={
                key:max((abs(r[key]) for r in consistency),default=0.) for key in
                ("peakXAbsDifference","amplitudeDifference","cOriginDifference","fourCOverYADifference")},
            rows=consistency),
        interpretation="The outputs describe finite-time geometric drift by height. The listed local log-slope ranges depend on the chosen contiguous logY window and are not asymptotic exponents.")
    (OUT/"geometry_detail.json").write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    np.savez_compressed(OUT/"geometry_detail.npz",**arrays)
    make_markdown(result)
    print(json.dumps({"out":str(OUT),"frames":[{"tau":f["tau"],"ell":f["ell"],"nativeYBounds":f["nativeBounds"]["Y"],
        "fourCOverYAAtY256":next((r["fourCOriginOverYA"] for r in f["nativeYProfiles"] if r["Ynative"]==256),None)} for f in frames],
        "tau20Crosscheck":result["tau20ConsistencyWithQuarterRelationJSON"]["maxAbsoluteDifferences"]},indent=2))


def make_markdown(result):
    tau20=next(f for f in result["frames"] if round(f["tau"])==20)
    lines=["# Detailed profile geometry from saved tau frames","",
        "Coordinates are restored independently for each frame: `X=a+ell*x`, `Y=ell*y`, `rho_X=P*q`, and `rho-rho(a,0)=P*ell*Q`. Both absolute native peak `m=a+ell*mInner` and centered peak `m-a` are retained.","",
        "The table's `c_origin` is `rho(m,Y)-rho(0,0)`, matching the quarter relation. `c_peak` is `rho(m,Y)-rho(a,0)` and differs by the measured wall integral `P*ell*integral(q(x,0), -a/ell..0)`. They are not interchangeable.","",
        "All profiles use cubic splines over saved samples. Requested native heights map to `y=Y/ell`; no height or x extrapolation is performed. Missing crossings are reported as null. `W` and `k` come from the measured left/right half-height widths with beta=1.5; level midlines and skews use each listed fraction of peak height.","",
        "## tau=20 at fixed native heights","",
        "| Y | m absolute X | (m-a)/Y | A native | W native | k | curvature width | c_origin | c_peak | c/(YA) | 4c/(YA) | AW native |","|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in tau20["nativeYProfiles"]:
        if not r["available"]: continue
        def fmt(v): return "—" if v is None else f"{v:.6g}"
        lines.append("| "+" | ".join(fmt(r[k]) for k in ("Ynative","mAbsoluteNative","mOverYCentered","A_native","Wnative","k","curvatureWidthNative","c_origin_rhoMinusRho00","c_peak_rhoMinusRho_a0","cOriginOverYA","fourCOriginOverYA","AW_native"))+" |")
    lines += ["","## Temporal drift at selected common native heights","",
        "Each row below evaluates the same native `Y` independently in all three frames; the frame-specific `a,P,ell` are used for restoration.","",
        "| Y | tau | m/Y absolute | W/Y | k | c_origin/(YA) | 4c_origin/(YA) | AW native |","|---:|---:|---:|---:|---:|---:|---:|---:|"]
    chosen={.25,1.,4.,16.,64.,256.,1024.}
    for item in result["fixedNativeYTimeEvolution"]:
        if item["Ynative"] not in chosen: continue
        for tau,v in item["byTau"].items():
            lines.append(f"| {item['Ynative']:g} | {tau} | {v['mOverYAbsolute']:.6g} | {v['WnativeOverY']:.6g} | {fmt(v['k'])} | {v['cOriginOverYA']:.6g} | {v['fourCOriginOverYA']:.6g} | {v['AW_native']:.6g} |")
    lines += ["","## tau=20 level-set geometry","","For each fixed native Y, `geometry_detail.json` contains left and right widths, both crossing coordinates, equal-height midpoint shift relative to the peak, and normalized left-right skew at q/A=1/4, 1/2, and 3/4. Representative native widths are tabulated below.","",
        "| Y | q/A | left width | right width | midline shift from m | skew |","|---:|---:|---:|---:|---:|---:|"]
    for row in tau20["nativeYProfiles"]:
        if row["Ynative"] not in {0.5,2.,8.,64.,256.}: continue
        for fraction in (0.25,0.5,0.75):
            level=row["levels"][str(fraction)]
            if level:
                lines.append(f"| {row['Ynative']:g} | {fraction:g} | {level['widthLeftNative']:.6g} | {level['widthRightNative']:.6g} | {level['midlineShiftNative']:.6g} | {level['skew']:.6g} |")
    lines += ["","The NPZ keeps compact native-Y arrays for A, W, k, c, m and c/(YA), plus common-inner-height peak positions. The tau20 peak, amplitude, density origin, and 4c/(YA) values at shared native heights match the previous `quarter_relation_checks.json` at stored precision.","",
        "## Local log slopes","","The JSON includes contiguous 3-, 5-, and 7-point windows in log(native Y) for m/Y, W/Y, k, c_origin/(YA), and AW. These slopes vary with the window and height region; they are descriptive finite-time measurements, not asymptotic exponent claims. Fixed common-inner-y tables are also included to separate fixed physical-height comparisons from fixed inner-coordinate comparisons.","",
        "An independent grid-stride audit found small A/W changes near the wall and larger k/curvature sensitivity at far heights. Geometry values are rounded for display; printed precision is not a continuum error bound.",""]
    (OUT/"GEOMETRY_DETAIL.md").write_text("\n".join(lines))


if __name__=="__main__":
    main()
