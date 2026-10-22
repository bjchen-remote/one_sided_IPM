"""Necessary wall tangency test for a five-function streamfunction skeleton.

Only reads saved source fields. Empirical wall mothers are calibrated at y=0,
then frozen when taking the geometric tangent space; no PDE is evolved.
"""
from pathlib import Path
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import CubicSpline, RectBivariateSpline
from scipy.optimize import least_squares

from fixed_point import fd_matrix
from streamfunction_sources import load_native, load_conditional
from streamfunction_shape_fit import mother, BETA

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/streamfunction_reduction_20260928"
GN, GW = leggauss(128)
PARAMETER_ORDER = ("c", "B", "n", "W", "k")


def G(xi, k):
    return (np.sinh(BETA*xi) + k*(np.cosh(BETA*xi)-1))/BETA


def Gxi(xi, k):
    return np.cosh(BETA*xi) + k*np.sinh(BETA*xi)


def Gk(xi):
    return (np.cosh(BETA*xi)-1)/BETA


def geometry(p, derivative_method="cubic_H"):
    x = p["x"]
    H = CubicSpline(x, -p["u1"][0], extrapolate=False)
    if derivative_method == "cubic_H":
        D = H.derivative()
    else:
        stencil = int(derivative_method.removeprefix("fd"))
        D = CubicSpline(x, fd_matrix(x, 1, stencil) @ (-p["u1"][0]), extrapolate=False)
    roots = D.derivative().roots(extrapolate=False)
    roots = roots[np.isfinite(roots) & (roots > x[0]) & (roots < x[-1])]
    n = float(roots[np.argmax(D(roots))])
    B = float(D(n))
    half = D.solve(B/2, extrapolate=False)
    left, right = half[half<n], half[half>n]
    if not len(left) or not len(right):
        return dict(available=False, n=n, B=B,
                    reason="Both nearest half-height crossings not available inside saved x domain"), H, D
    L, R = n-float(left[-1]), float(right[0])-n
    k = (R-L)/(R+L)/np.tanh(BETA/2)
    W = (R+L)*BETA/(2*np.sinh(BETA))
    return dict(available=True, n=n, B=B, W=W, k=k, c=float(H(n)), L=L, R=R,
                halfCrossings=half.tolist(), globalMapMonotone=bool(abs(k)<1)), H, D


def integrate(g, xi, k, k_derivative=False):
    z = np.asarray(xi)[...,None]*(GN+1)/2
    factor = np.sinh(BETA*z) if k_derivative else Gxi(z,k)
    return np.asarray(xi)/2*np.sum(GW*g(z)*factor, axis=-1)


def relative(error, ref):
    return float(np.linalg.norm(error)/np.linalg.norm(ref))


def scores(pred, obs, train, normalization=None):
    den = float(np.max(abs(obs))) if normalization is None else normalization
    return dict(trainRelativeL2=relative((pred-obs)[train], obs[train]),
                heldoutRelativeL2=relative((pred-obs)[~train], obs[~train]),
                allRelativeL2=relative(pred-obs, obs),
                heldoutLinfOverGlobalPeak=float(np.max(abs(pred-obs)[~train])/den))


def tangent(geo, xi, g):
    B, W, k = (geo[key] for key in ("B", "W", "k"))
    gg = g(xi)
    j, jk = integrate(g,xi,k), integrate(g,xi,k,True)
    return np.column_stack([np.ones_like(xi), W*j, -B*gg,
                            B*(j-gg*G(xi,k)), B*W*(jk-gg*Gk(xi))])


def project(T, target, train):
    scale = np.linalg.norm(T[train],axis=0)
    matrix = T/scale
    coefficient, _, rank, singular = np.linalg.lstsq(matrix[train], target[train], rcond=None)
    parameters = coefficient/scale
    fitted = T@parameters
    allcoef = np.linalg.lstsq(matrix,target,rcond=None)[0]/scale
    return dict(parameterDerivatives=dict(zip(PARAMETER_ORDER,map(float,parameters))),
                allPointParameterDerivatives=dict(zip(PARAMETER_ORDER,map(float,allcoef))),
                scaledMatrixRank=int(rank), scaledSingularValues=singular.tolist(),
                scaledConditionNumber=float(singular[0]/singular[-1]),
                residual=scores(fitted,target,train),
                allPointProjectionResidual=scores(T@allcoef,target,train),
                observedNormalDerivativeMin=float(target.min()),
                observedNormalDerivativeMax=float(target.max()),
                fittedNormalDerivativeMin=float(fitted.min()),
                fittedNormalDerivativeMax=float(fitted.max()),
                remainderRms=float(np.sqrt(np.mean((target-fitted)**2))),
                remainderLinf=float(np.max(abs(target-fitted)))), fitted


def fixed_profile_refit(p, H, D, geo, xi, train, pars):
    """Allow the five geometric parameters to fit H and D on fixed x points.

    It distinguishes inadequate observed half-width calibration from inability
    of the fixed mother family to approximate the actual wall shape.
    """
    x = geo["n"]+geo["W"]*G(xi,geo["k"])
    dobs, hobs = D(x), H(x)
    B0,W0,n0,c0,k0=(geo[key] for key in ("B","W","n","c","k"))
    def evaluate(z):
        n, B, W, k, c = n0+W0*z[0], B0*np.exp(z[1]), W0*np.exp(z[2]), z[3], c0+B0*W0*z[4]
        t=(np.arcsinh((BETA*(x-n)/W+k)/np.sqrt(1-k*k))-np.arctanh(k))/BETA
        f=lambda zz:mother(zz,pars)
        return B*f(t), c+B*W*integrate(f,t,k), dict(n=float(n),B=float(B),W=float(W),k=float(k),c=float(c))
    def residual(z):
        dp,hp,_=evaluate(z)
        return np.r_[(dp-dobs)[train]/B0,(hp-hobs)[train]/(B0*W0)]
    fit=least_squares(residual,[0,0,0,k0,0],bounds=([-1,-1,-1,-.98,-2],[1,1,1,.98,2]),
                      ftol=1e-11,xtol=1e-11,gtol=1e-11)
    dp,hp,gg=evaluate(fit.x)
    t=(np.arcsinh((BETA*(x-gg["n"])/gg["W"]+gg["k"])/np.sqrt(1-gg["k"]**2))-np.arctanh(gg["k"]))/BETA
    T=tangent(gg,t,lambda zz:mother(zz,pars))
    q=CubicSpline(p["x"],p["q"][0],extrapolate=False)(x)
    gate,_=project(T,-q/2,train)
    return dict(geometry=gg,D=scores(dp,dobs,train,B0),H=scores(hp,hobs,train,B0*W0),
                normalProjection=gate,
                note="Five geometry values jointly fit saved H and D on training fixed x only; frozen mother parameters unchanged. Validation keeps original physical x samples.")


def normal_difference_check(p,x,train):
    spline=RectBivariateSpline(p["y"],p["x"],p["psi"])
    saved_target=-CubicSpline(p["x"],p["q"][0],extrapolate=False)(x)/2
    h0=spline.ev(np.zeros_like(x),x,dx=1)
    wall_velocity=-CubicSpline(p["x"],p["u1"][0],extrapolate=False)(x)
    rows=[]
    for h in (.005,.01,.02):
        if 2*h>p["y"][-1]:
            continue
        h1=spline.ev(np.full_like(x,h),x)/h
        h2=spline.ev(np.full_like(x,2*h),x)/(2*h)
        estimate=(-3*h0+4*h1-h2)/(2*h)
        from_saved_wall=(-3*wall_velocity+4*h1-h2)/(2*h)
        rows.append(dict(innerYStep=h,consistentSplineWall=scores(estimate,saved_target,train),
                         savedVelocityWall=scores(from_saved_wall,saved_target,train)))
    analytic=spline.ev(np.zeros_like(x),x,dx=2)/2
    return dict(method="Second-order one-sided difference of H=psi/y. The consistent wall uses psi_y(x,0); savedVelocityWall uses -u1(x,0) and reveals interpolation mismatch divided by h.",
        wallTraceDifferenceLinf=float(max(abs(h0-wall_velocity))),
        directSplineHNormalVsMinusQHalf=scores(analytic,saved_target,train),steps=rows)


def one_source(p,frozen):
    name=p["name"]
    native=p["kind"]=="finite_time_native"
    extent,trainlimit=(1.5,1.) if native else (1.,.7)
    xi=np.linspace(-extent,extent,801)
    train=abs(xi)<=trainlimit
    rows=[]
    for method in ("cubic_H","fd5","fd7"):
        geo,H,D=geometry(p,method)
        if not geo["available"]:
            rows.append(dict(derivativeMethod=method,geometry=geo));continue
        x=geo["n"]+geo["W"]*G(xi,geo["k"])
        if x.min()<p["x"][0] or x.max()>p["x"][-1]:
            raise ValueError(f"No extrapolation: {name}/{method}")
        q=CubicSpline(p["x"],p["q"][0],extrapolate=False)(x)
        gemp=lambda zz:D(geo["n"]+geo["W"]*G(zz,geo["k"]))/geo["B"]
        models=[]
        mothers={**{family:(lambda zz, pars=pars:mother(zz,pars)) for family,pars in frozen.items()},
                 "empirical_wall":gemp}
        for family,g in mothers.items():
            dp=geo["B"]*g(xi)
            hp=geo["c"]+geo["B"]*geo["W"]*integrate(g,xi,geo["k"])
            T=tangent(geo,xi,g)
            gate,_=project(T,-q/2,train)
            rec=dict(family=family,geometrySource="measured wall peak and half widths; c=H(n)",
                     wallD=scores(dp,D(x),train,geo["B"]),wallH=scores(hp,H(x),train,geo["B"]*geo["W"]),
                     normalProjection=gate)
            if family!="empirical_wall" and method=="cubic_H":
                rec["refittedWallGeometry"]=fixed_profile_refit(p,H,D,geo,xi,train,frozen[family])
            models.append(rec)
        rows.append(dict(derivativeMethod=method,geometry=geo,models=models,
                    physicalXInterval=[float(x.min()),float(x.max())],
                    normalFiniteDifference=(normal_difference_check(p,x,train) if method=="cubic_H" else None)))
    return dict(name=name,source=p["source"],kind=p["kind"],
                savedXInterval=[float(p["x"][0]),float(p["x"][-1])],
                normProtocol=dict(coordinate="uniform measured-wall xi, not physical x quadrature",xiRange=[-extent,extent],
                    trainingXiRange=[-trainlimit,trainlimit],holdout="outer xi bands, both sides",samples=len(xi),
                    derivativeConvention="All parameter derivatives and tangent vectors take y derivatives at fixed physical/inner x."),
                derivativeVariants=rows)


def main():
    shape=json.loads((ROOT/"runs/streamfunction_geometry_20260928/shape_fit.json").read_text())
    frozen={family:np.asarray(next(r["parameters"] for r in shape["fits"] if r["coordinate"]=="sinh" and r["family"]==family))
            for family in ("student","student_skew_even")}
    result=dict(title="Necessary wall tangency of a five-geometry streamfunction skeleton", 
        definition="H=c+B W J(xi,k), x=n+W G(xi,k), J=int_0^xi g(t) G_t(t,k)dt; g is held fixed.",
        parameterOrder=list(PARAMETER_ORDER),
        tangentBasis=["1","W J","-B g","B (J-g G)","B W (J_k-g G_k)"],
        necessaryEquation="At y=0, H_y=-q/2. Project -q/2 onto the five columns. The residual must be retained as a normal shape mode unless it vanishes at source accuracy.",
        caveats=["A small one-step tangent residual does not establish invariant finite-dimensional closure, elliptic globality, or transport compatibility.",
                 "Finite-height mothers can fail to represent the wall; empirical_wall freezes the observed wall profile and separates mother mismatch from tangency failure.",
                 "All projection errors are held-out x-band tests of the same wall, not time or height prediction.",
                 "Column-normalized condition numbers diagnose whether numerical parameter derivatives are identifiable.",
                 "Conditional finite-box L4 results depend on frozen exterior; nearest left half-height crossing lies close to the boundary."],
        sources=[one_source(load_native(20),frozen),one_source(load_conditional("L4_257x129"),frozen)])
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"wall_gate.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    for s in result["sources"]:
        print(s["name"])
        for v in s["derivativeVariants"]:
            print(v["derivativeMethod"],v["geometry"])
            for m in v.get("models",[]):
                a=m["normalProjection"]
                print(m["family"],"H",m["wallH"]["heldoutRelativeL2"],"D",m["wallD"]["heldoutRelativeL2"],
                      "target",a["residual"],"cond",a["scaledConditionNumber"],"derivatives",a["parameterDerivatives"])
                if "refittedWallGeometry" in m:
                    r=m["refittedWallGeometry"]
                    print("refit","H",r["H"],"D",r["D"],"target",r["normalProjection"]["residual"])


if __name__=="__main__":
    main()
