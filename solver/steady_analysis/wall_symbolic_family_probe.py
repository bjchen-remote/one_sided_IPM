"""Near-wall data test of the one-parameter Laurent-rational wall family."""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.io import loadmat
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/skeleton_closure_20260928/wall_symbolic_family.json"
L4 = ROOT / "runs/steady_20260928_tau20/L4_257x129/state.npz"
L8 = ROOT / "runs/steady_20260928_tau20/L8_257x129_continued/state.npz"
MAT = ROOT / "runs/steady_20260928_tau20/run_analysis.mat"


def get_wall_data():
    result = {}
    for name, path in (("L4_257x129", L4), ("L8_257x129_continued", L8)):
        with np.load(path) as d:
            result[name] = (d["x"].copy(), d["qx"][0].copy(), d["Q"][0].copy())
    d = loadmat(MAT, simplify_cells=True)
    meta = d["metadata"]
    x = (d["xFull"] - float(meta["a"])) / float(meta["ell"])
    Qraw = d["rho"][0]
    origin = np.interp(float(meta["a"]), d["xFull"], Qraw)
    Q = (Qraw-origin)/(float(meta["P"])*float(meta["ell"]))
    q = d["rhoX"][0]/float(meta["P"])
    result["tau20"] = (x, q, Q)
    return result


def model(x, alpha, k):
    p, q = alpha, 1.0
    C = np.sqrt(q*q-p*p)/(p*np.sqrt(2))
    D = np.sqrt(q*q-p*p)/(q*np.sqrt(2))
    ell = q*k/p
    s = np.asarray(x)/D+k
    disc = np.sqrt(s*s+(1-k*k))
    t = (s+disc)/(1+k)
    eta = np.log(t)/q
    Q = C*(np.sinh(p*eta)+ell*(np.cosh(p*eta)-1))
    den = np.cosh(q*eta)+k*np.sinh(q*eta)
    qx = (np.cosh(p*eta)+ell*np.sinh(p*eta))/den
    return qx, Q


def half_width_skew(x, q):
    cs = CubicSpline(x, q)
    peaks = cs.derivative().roots(extrapolate=False)
    peaks = peaks[(peaks >= x[0]) & (peaks <= x[-1])]
    m = float(peaks[np.argmax(cs(peaks))])
    A = float(cs(m))
    roots = CubicSpline(x, q-A/2).roots(extrapolate=False)
    left, right = roots[roots<m], roots[roots>m]
    if not len(left) or not len(right):
        return None
    wl, wr = m-float(left[-1]), float(right[0])-m
    return dict(m=m, A=A, wLeft=wl, wRight=wr, epsilon=(wr-wl)/(wr+wl))


def errors(x, q, Q, alpha, k, indices):
    qp, Qp = model(x[indices], alpha, k)
    q0, Q0 = q[indices], Q[indices]
    return dict(qRelativeL2=float(np.linalg.norm(qp-q0)/np.linalg.norm(q0)),
                QRelativeL2=float(np.linalg.norm(Qp-Q0)/max(np.linalg.norm(Q0),1e-14)),
                qRMSE=float(np.sqrt(np.mean((qp-q0)**2))),
                QRMSE=float(np.sqrt(np.mean((Qp-Q0)**2))))


def simple_fraction(k, alpha):
    candidates = {Fraction(n,d) for d in range(1,13) for n in range(-d,d+1)}
    candidates = [r for r in candidates if -alpha < float(r) < alpha]
    return min(candidates, key=lambda r: abs(float(r)-k))


def jet_diagnostics(x, q, alpha, kfit):
    """Estimate q'''(0), q''''(0), and I4 across local windows/resolutions."""
    estimates=[]
    for radius in (0.30,0.40,0.50,0.65,0.80):
        base=np.flatnonzero(abs(x)<=radius)
        for stride in (1,2,3):
            ids=base[::stride]
            if len(ids)<10 or not (np.any(x[ids]<0) and np.any(x[ids]>0)):
                continue
            z=x[ids]/radius
            V=np.vander(z, N=9, increasing=True)
            coeff=np.linalg.lstsq(V,q[ids],rcond=None)[0]
            q3=float(6*coeff[3]/radius**3)
            q4=float(24*coeff[4]/radius**4)
            I4=float(q4+0.7*q3*q3)
            estimates.append(dict(windowRadius=radius,stride=stride,sampleCount=len(ids),
                                  q3=q3,q4=q4,I4=I4,
                                  kInferredFromQ3=float(q3*np.sqrt(1-alpha*alpha)/(10*np.sqrt(2)))))
    vals=np.asarray([r["I4"] for r in estimates])
    q3vals=np.asarray([r["q3"] for r in estimates])
    predicted=4+32/(1-alpha*alpha)
    return dict(theoreticalI4=float(predicted), theoreticalQ3FromFittedK=float(10*np.sqrt(2)*kfit/np.sqrt(1-alpha*alpha)),
                I4EstimateRange=[float(vals.min()),float(vals.max())] if len(vals) else None,
                I4EstimateMedian=float(np.median(vals)) if len(vals) else None,
                I4SpreadOverWindowsAndResolutions=float(np.ptp(vals)) if len(vals) else None,
                q3EstimateRange=[float(q3vals.min()),float(q3vals.max())] if len(q3vals) else None,
                estimates=estimates)


def fit_domain(x, q, Q, alpha, mask):
    ids = np.flatnonzero(mask)
    train, hold = ids[::2], ids[1::2]
    if len(train)<8 or len(hold)<5:
        return dict(available=False, sampleCount=int(len(ids)))
    def objective(k):
        qh, Qh = model(x[train], alpha, float(k))
        return float(np.mean((qh-q[train])**2)/np.mean(q[train]**2) +
                     np.mean((Qh-Q[train])**2)/max(np.mean(Q[train]**2),1e-14))
    fit = minimize_scalar(objective, bounds=(-alpha+1e-8, alpha-1e-8), method="bounded",
                          options={"xatol":1e-12})
    kfit = float(fit.x)
    frac = simple_fraction(kfit, alpha)
    qfit, Qfit = model(x, alpha, kfit)
    qfrac, Qfrac = model(x, alpha, float(frac))
    measured_shape = half_width_skew(x, q)
    fitted_shape = half_width_skew(x, qfit)
    fraction_shape = half_width_skew(x, qfrac)
    return dict(available=True, sampleCount=int(len(ids)), xBounds=[float(x[ids].min()),float(x[ids].max())],
        trainCount=int(len(train)), holdoutCount=int(len(hold)), kFit=kfit,
        train=errors(x,q,Q,alpha,kfit,train), holdout=errors(x,q,Q,alpha,kfit,hold),
        nearestSimpleFraction=str(frac), nearestSimpleValue=float(frac),
        simpleFractionTrain=errors(x,q,Q,alpha,float(frac),train),
        simpleFractionHoldout=errors(x,q,Q,alpha,float(frac),hold),
        measuredHalfWidth=measured_shape, fittedHalfWidth=fitted_shape,
        fixedFractionHalfWidth=fraction_shape,
        jetDiagnostics=jet_diagnostics(x,q,alpha,kfit),
        wallInvariants=dict(modelCenterQ=float(model(np.array([0.]),alpha,kfit)[1][0]),
                            modelCenterQx=float(model(np.array([0.]),alpha,kfit)[0][0]),
                            modelCenterQxx=-2.0))


def main():
    data = get_wall_data()
    alphas = {"2/5":2/5, "3/8":3/8, "1/3_control":1/3}
    output = {"description":"Near-wall fit of the one-parameter symbolic family; far tails excluded by the two stated masks.",
              "coordinateConvention":"All x,Q,q are in each source's saved or restored inner coordinates, with Q(0)=0.",
              "fitProtocol":"One scalar k is minimized jointly against q and Q on alternating samples; the other parity is held out. A nearest rational k with denominator <=12 is also fixed and scored.",
              "I4Interpretation":"I4EstimateRange is a sensitivity envelope across the listed polynomial windows and downsampling strides, not a confidence interval or rigorous error bound. For L4 the envelope 35.45-40.95 includes the 1/3 prediction 40.0 but excludes the 3/8 and 2/5 predictions 41.24 and 42.10. Window, sampling, and source variation do not provide a robust exponent choice; L4's current window family leans toward 1/3, but this is not convergence evidence.",
              "sources":{}}
    for name,(x,q,Q) in data.items():
        order=np.argsort(x); x,q,Q=x[order],q[order],Q[order]
        masks={"abs_x_le_1p5":abs(x)<=1.5, "q_ge_half":q>=0.5}
        source={"gridPoints":int(len(x)),"fullXBounds":[float(x[0]),float(x[-1])],"wallQAtX0":float(np.interp(0,x,Q)),"wallqAtX0":float(np.interp(0,x,q)),"domains":{}}
        for dname,mask in masks.items():
            source["domains"][dname]={}
            for aname,alpha in alphas.items():
                source["domains"][dname][aname]=fit_domain(x,q,Q,alpha,mask)
        output["sources"][name]=source
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(output,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print(json.dumps({n:{d:{a:{"k":r.get("kFit"),"qhold":r.get("holdout",{}).get("qRelativeL2"),"Qhold":r.get("holdout",{}).get("QRelativeL2"),"frac":r.get("nearestSimpleFraction"),"qfrac":r.get("simpleFractionHoldout",{}).get("qRelativeL2"),"Qfrac":r.get("simpleFractionHoldout",{}).get("QRelativeL2")} for a,r in aa.items()} for d,aa in s["domains"].items()} for n,s in output["sources"].items()},indent=2))


if __name__ == "__main__":
    main()
