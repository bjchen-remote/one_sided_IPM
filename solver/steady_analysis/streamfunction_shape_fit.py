"""Finite-family streamfunction shape probe on the saved tau=20 field.

Uses measured geometry at every height, including held-out heights. Only the
common transverse mother shape is trained. This is not height-law prediction,
a closed solution, or a PDE evolution. Every query stays in the saved domain.
"""
from pathlib import Path
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import CubicSpline, RectBivariateSpline
from scipy.optimize import least_squares

from streamfunction_sources import load_native

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/streamfunction_geometry_20260928"
SOURCE = ROOT / "runs/steady_20260928_tau20/run_analysis.mat"
HEIGHTS = np.array([4., 8., 16., 32., 64., 128., 256.])
TRAIN = np.array([True, False, True, False, True, False, True])
XI = np.linspace(-1.5, 1.5, 601)
BETA = 1.5
GN, GW = leggauss(64)
MAPS = ("affine_halfwidth", "piecewise_halfwidth", "sinh")
FAMILIES = ("student", "student_skew", "student_skew_even")


def load_source(tau=20):
    raw = load_native(tau)
    x, y = raw["x"], raw["y"]
    return dict(x=x, y=y, a=raw["a"], P=raw["P"], ell=raw["ell"],
                tau=raw["tau"], source=raw["source"],
                fields={k: RectBivariateSpline(y, x, raw[k]) for k in ("psi", "u1", "u2", "q")})


def field(source, name, y, x):
    x = np.asarray(x)
    if y < source["y"][0] or y > source["y"][-1] or x.min() < source["x"][0] or x.max() > source["x"][-1]:
        raise ValueError("No extrapolation beyond saved field")
    return source["fields"][name](y, x, grid=False)


def measure(source, Y):
    y = float(Y / source["ell"])
    x = source["x"]
    d = field(source, "u2", y, x) / y
    cs = CubicSpline(x, d, extrapolate=False)
    extrema = cs.derivative().roots(extrapolate=False)
    extrema = extrema[np.isfinite(extrema)]
    n = float(extrema[np.argmax(cs(extrema))])
    B = float(cs(n))
    cross = CubicSpline(x, d - B / 2, extrapolate=False).roots(extrapolate=False)
    L, R = n - float(cross[cross < n][-1]), float(cross[cross > n][0]) - n
    zeros = CubicSpline(x, field(source, "psi", y, x) / y, extrapolate=False).roots(extrapolate=False)
    # The saved finite box has an additional distant zero; choose the core one.
    s = float(zeros[np.argmin(abs(zeros - n))])
    k = (R - L) / (R + L) / np.tanh(BETA / 2)
    if abs(k) >= 1:
        raise ValueError("Sinh map ceases to be globally monotone")
    return dict(Ynative=float(Y), yInner=y, n=n, B=B, L=L, R=R,
                W=(L + R) * BETA / (2 * np.sinh(BETA)), k=k, s=s,
                nNative=float(source["a"] + source["ell"] * n),
                sNative=float(source["a"] + source["ell"] * s),
                allSavedHZeroes=zeros.tolist())


def coordinate(xi, geo, name):
    if name == "affine_halfwidth":
        return geo["n"] + (geo["L"] + geo["R"]) / 2 * xi
    if name == "piecewise_halfwidth":
        return geo["n"] + np.where(xi < 0, geo["L"], geo["R"]) * xi
    return geo["n"] + geo["W"] * (np.sinh(BETA * xi) + geo["k"] * (np.cosh(BETA * xi) - 1)) / BETA


def inverse(x, geo, name):
    z = x - geo["n"]
    if name == "affine_halfwidth":
        return 2 * z / (geo["L"] + geo["R"])
    if name == "piecewise_halfwidth":
        return z / np.where(z < 0, geo["L"], geo["R"])
    k = geo["k"]
    return (np.arcsinh((BETA * z / geo["W"] + k) / np.sqrt(1 - k*k)) - np.arctanh(k)) / BETA


def jacobian(xi, geo, name):
    if name == "affine_halfwidth":
        return np.full_like(xi, (geo["L"] + geo["R"]) / 2)
    if name == "piecewise_halfwidth":
        return np.where(xi < 0, geo["L"], geo["R"])
    return geo["W"] * (np.cosh(BETA * xi) + geo["k"] * np.sinh(BETA * xi))


def mother(xi, parameters):
    p = parameters[0]
    c = np.expm1(np.log(2) / p)
    logshape = -p * np.log1p(c * xi*xi)
    if len(parameters) >= 2:
        logshape += parameters[1] * xi**3 * (xi*xi - 1) / (1 + xi*xi)**2.5
    if len(parameters) >= 3:
        logshape += parameters[2] * xi*xi * (xi*xi - 1) / (1 + xi*xi)**2
    return np.exp(logshape)


def primitive(xi, geo, name, parameters):
    xi = np.asarray(xi)
    z = xi[..., None] * (GN + 1) / 2
    return xi / 2 * np.sum(GW * mother(z, parameters) * jacobian(z, geo, name), axis=-1)


def predict_psi(x, geo, name, parameters):
    xi = inverse(np.asarray(x), geo, name)
    xizero = inverse(geo["s"], geo, name)
    return geo["yInner"] * geo["B"] * (primitive(xi, geo, name, parameters) - primitive(xizero, geo, name, parameters))


def relative(prediction, observed):
    return float(np.linalg.norm(prediction - observed) / np.linalg.norm(observed))


def metrics(prediction, observed, scale=1.):
    return dict(relativeL2=relative(prediction, observed),
                linfOverScale=float(np.max(abs(prediction - observed)) / scale))


def derivative_audit(source, parameters):
    """Conditional measured-geometry derivative test, not geometry prediction."""
    result = []
    for Y in HEIGHTS[~TRAIN]:
        center = measure(source, Y)
        # Common core in the selected smooth coordinates, safely away from box.
        xi = np.linspace(-1.5, 1.5, 301)
        x = coordinate(xi, center, "sinh")
        y = center["yInner"]
        u1 = field(source, "u1", y, x)
        q = field(source, "q", y, x)
        row = dict(Ynative=float(Y), stepChecks=[])
        for fraction in (.005, .01, .02):
            hY = Y * fraction
            h = hY / source["ell"]
            psi = np.array([predict_psi(x, measure(source, Y + j * hY), "sinh", parameters) for j in (-2, -1, 0, 1, 2)])
            u1fit = -(psi[0] - 8 * psi[1] + 8 * psi[3] - psi[4]) / (12 * h)
            psiyy = (-psi[4] + 16 * psi[3] - 30 * psi[2] + 16 * psi[1] - psi[0]) / (12 * h*h)
            # Analytic psi_x=y B F; stable central derivative of known F.
            eps = 1e-4
            fp = (mother(xi + eps, parameters) - mother(xi - eps, parameters)) / (2 * eps)
            psixx = y * center["B"] * fp / jacobian(xi, center, "sinh")
            qfit = -psixx - psiyy
            row["stepChecks"].append(dict(nativeYStepFraction=fraction,
                u1=metrics(u1fit, u1, max(abs(u1))), q=metrics(qfit, q, max(abs(q))),
                qfitMinimum=float(qfit.min()), savedQMinimum=float(q.min()),
                qfitNegativeSamples=int(np.count_nonzero(qfit < 0)),
                savedQNegativeSamples=int(np.count_nonzero(q < 0)),
                qfitMinimumAtXi=float(xi[np.argmin(qfit)]),
                sampleCount=len(xi)))
        sf = source["fields"]["psi"]
        direct_u1 = -sf(y, x, dx=1, dy=0, grid=False)
        direct_q = -sf(y, x, dx=2, dy=0, grid=False) - sf(y, x, dx=0, dy=2, grid=False)
        row["savedPsiSplineDerivativeBaseline"] = dict(u1=metrics(direct_u1, u1, max(abs(u1))), q=metrics(direct_q, q, max(abs(q))))
        result.append(row)
    return result


def temporal_transfer(fits, arrays):
    """Freeze tau20 mothers; remeasure only each source's geometry."""
    frozen = {family: next(f for f in fits if f["coordinate"] == "sinh" and f["family"] == family)
              for family in ("student", "student_skew_even")}
    result = []
    for tau in (10, 15, 20):
        source = load_source(tau)
        geometry = [measure(source, Y) for Y in HEIGHTS]
        mapped = np.array([coordinate(XI, g, "sinh") for g in geometry])
        observed = np.array([field(source, "u2", g["yInner"], x) / (g["yInner"] * g["B"])
                             for x, g in zip(mapped, geometry)])
        psi = np.array([field(source, "psi", g["yInner"], x) for x, g in zip(mapped, geometry)])
        psiscale = np.array([g["yInner"] * g["B"] * (g["L"] + g["R"]) / 2 for g in geometry])[:, None]
        models = []
        for family, frozen_model in frozen.items():
            pars = np.asarray(frozen_model["parameters"])
            prediction = np.broadcast_to(mother(XI, pars), observed.shape)
            psifit = np.array([predict_psi(x, g, "sinh", pars) for x, g in zip(mapped, geometry)])
            models.append(dict(family=family, frozenTau20Parameters=pars.tolist(),
                shapeAllHeights=metrics(prediction, observed),
                psiAllHeightsHeightNormalized=metrics(psifit / psiscale, psi / psiscale),
                shapeTau20HeldoutHeights=metrics(prediction[~TRAIN], observed[~TRAIN]),
                psiTau20HeldoutHeightsHeightNormalized=metrics((psifit / psiscale)[~TRAIN], (psi / psiscale)[~TRAIN]),
                perHeight=[dict(Ynative=g["Ynative"], shapeAndU2=metrics(prediction[i], observed[i]),
                                psi=metrics(psifit[i], psi[i], psiscale[i, 0])) for i, g in enumerate(geometry)]))
            arrays[f"transfer_tau{tau}_{family}_prediction_psi"] = psifit
        arrays[f"transfer_tau{tau}_x"] = mapped
        arrays[f"transfer_tau{tau}_observed_D_over_B"] = observed
        arrays[f"transfer_tau{tau}_observed_psi"] = psi
        result.append(dict(tau=tau, source=source["source"], geometry=geometry, models=models))
    return dict(protocol="Tau20 student and student_skew_even mother parameters are frozen; no refit at tau10/15. Each frame has separately measured n,B,L,R,s at the same native Y=4..256. This tests transfer of transverse shape, not prediction of height geometry or time evolution.",
                norm="Equal height and uniform xi sampling, with psi divided by y B (L+R)/2 before pooling. Within each height, shape relative L2 equals u2 relative L2.",
                frames=result)


def main():
    source = load_source()
    geometry = [measure(source, Y) for Y in HEIGHTS]
    arrays = dict(nativeY=HEIGHTS, xi=XI, trainHeightMask=TRAIN)
    alignments, fits = {}, []
    for name in MAPS:
        mapped = np.array([coordinate(XI, g, name) for g in geometry])
        observed = np.array([field(source, "u2", g["yInner"], x) / (g["yInner"] * g["B"]) for x, g in zip(mapped, geometry)])
        mean = observed.mean(axis=0)
        trainmean = observed[TRAIN].mean(axis=0)
        alignments[name] = dict(meanCurveSpreadRelativeL2=relative(np.broadcast_to(mean, observed.shape), observed),
            trainMeanToHoldoutRelativeL2=relative(np.broadcast_to(trainmean, observed[~TRAIN].shape), observed[~TRAIN]),
            maxPointwiseStd=float(observed.std(axis=0).max()),
            mappedXBounds=[float(mapped.min()), float(mapped.max())],
            note="Equal height/xi weights; physical comparison bands differ between maps.")
        arrays[f"{name}_x"] = mapped
        arrays[f"{name}_observed_D_over_B"] = observed
        for npar, family in enumerate(FAMILIES, start=1):
            init = [1.] + [0.] * (npar - 1)
            fit = least_squares(lambda par: (mother(XI, par)[None, :] - observed[TRAIN]).ravel(), init,
                                bounds=([.08] + [-4.] * (npar - 1), [30.] + [4.] * (npar - 1)),
                                ftol=1e-12, xtol=1e-12, gtol=1e-12)
            curve = mother(XI, fit.x)
            pred = np.broadcast_to(curve, observed.shape)
            psifit = np.array([predict_psi(x, g, name, fit.x) for x, g in zip(mapped, geometry)])
            psi = np.array([field(source, "psi", g["yInner"], x) for x, g in zip(mapped, geometry)])
            # Scale each height before pooling to avoid large-height dominance.
            psiscale = np.array([g["yInner"] * g["B"] * (g["L"] + g["R"]) / 2 for g in geometry])[:, None]
            record = dict(coordinate=name, family=family, parameters=fit.x.tolist(),
                shapeTrain=metrics(pred[TRAIN], observed[TRAIN]), shapeHoldout=metrics(pred[~TRAIN], observed[~TRAIN]),
                psiTrainHeightNormalized=metrics((psifit/psiscale)[TRAIN], (psi/psiscale)[TRAIN]),
                psiHoldoutHeightNormalized=metrics((psifit/psiscale)[~TRAIN], (psi/psiscale)[~TRAIN]),
                positivityOnFitInterval=bool(np.all(curve > 0)),
                monotoneFromPeakOnFitInterval=bool(np.all(np.diff(curve[:301]) >= 0) and np.all(np.diff(curve[300:]) <= 0)),
                perHeight=[dict(Ynative=g["Ynative"], training=bool(TRAIN[i]),
                    u2=metrics(pred[i], observed[i]), psi=metrics(psifit[i], psi[i], psiscale[i, 0])) for i, g in enumerate(geometry)])
            fits.append(record)
            arrays[f"{name}_{family}_prediction_D_over_B"] = curve
            arrays[f"{name}_{family}_prediction_psi"] = psifit
        arrays[f"{name}_observed_psi"] = psi
    # Diagnostic explicitly chosen a priori: smooth sinh map, richest preset
    # family (three parameters). All nine fit results remain reported.
    chosen = next(f for f in fits if f["coordinate"] == "sinh" and f["family"] == "student_skew_even")
    single = next(f for f in fits if f["coordinate"] == "sinh" and f["family"] == "student")
    output = dict(source=str(SOURCE), tau=20, coordinateDefinitions=dict(
        inner="x=(X-a)/ell; y=Y/ell; psi_c=(psi_native+uAtPeak*Y)/(P*ell^2)",
        fields="H=psi_c/y, D=H_x=u2_inner/y, n=argmax D, B=max D, s=near-core H zero", 
        gauge="Centered velocity gauge fixed by saved uAtPeak; wall psi=0. Far-box H zero is not selected.",
        nativeParameters={k:source[k] for k in ("a", "P", "ell")}),
        protocol=dict(nativeY=HEIGHTS.tolist(), trainingY=HEIGHTS[TRAIN].tolist(), holdoutY=HEIGHTS[~TRAIN].tolist(),
            xiRange=[-1.5,1.5], xiSamples=len(XI), beta=BETA,
            conditionalGeometry="n,B,L,R,s are measured independently at every height, including held-out heights. Only shared mother-shape parameters are trained.",
            families="F=(1+(2^(1/p)-1)xi^2)^(-p) * exp(a xi^3(xi^2-1)/(1+xi^2)^(5/2) + b xi^2(xi^2-1)/(1+xi^2)^2); successive presets use p, (p,a), (p,a,b). Positive on all real xi; no global monotonicity claim.",
            identities="F(0)=1, F'(0)=0, F(+/-1)=1/2. Shapes were not borrowed from q. For asymmetric affine scaling the last identity is deliberately restrictive.",
            errorNorms="Relative L2 uses equally spaced xi samples, not physical x quadrature. Pooled D/B and psi/[y B (L+R)/2] scores give heights comparable scale; u2 relative error equals D/B error at each fixed height.",
            reconstruction="psi_fit=y B [integral(0..xi)F(t) x_t dt - integral(0..xi_s)F(t) x_t dt], xi_s inverse map of measured s.",
            caveat="This is finite-time descriptive fitting, not a global elliptic or transport closure. Coordinate maps cover different physical bands; scores compare collapse, not a common-x model contest.",
            noExtrapolation=True), geometry=geometry, alignments=alignments, fits=fits,
        derivativeAudit=dict(model="sinh/student_skew_even, chosen a priori for smoothness and three-parameter preset", 
            note="Geometry and zero level remeasured at Y+-h,Y+-2h, then five-point y derivatives of fitted psi at fixed physical/inner x, not fixed xi. Metrics and sign checks sample uniform xi in [-1.5,1.5] at the central height, not uniform x quadrature. Conditional on observed height dependence, not independent height prediction. q is -Laplacian psi, tested against saved rho_X/P; positivity is asserted only on these samples if observed.",
            rows=derivative_audit(source, np.asarray(chosen["parameters"]))),
        singleParameterDerivativeAudit=dict(model="sinh/student, frozen tau20 p from the same four training heights", 
            parameters=single["parameters"],
            note="Same fixed-x derivative and uniform-xi scoring protocol as derivativeAudit. No extra models or refits. Geometry remains measured at every neighboring height; this is conditional derivative consistency.",
            rows=derivative_audit(source, np.asarray(single["parameters"]))),
        temporalTransfer=temporal_transfer(fits, arrays))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "shape_fit.json").write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    np.savez_compressed(OUT / "shape_fit.npz", **arrays)
    print(json.dumps(dict(alignments=alignments,
        fits=[{k:r[k] for k in ("coordinate", "family", "parameters", "shapeTrain", "shapeHoldout", "psiHoldoutHeightNormalized")} for r in fits],
        derivativeAudit=output["derivativeAudit"],
        temporalTransfer=[dict(tau=f["tau"], models=[{k:m[k] for k in ("family", "shapeAllHeights", "psiAllHeightsHeightNormalized")} for m in f["models"]]) for f in output["temporalTransfer"]["frames"]]), indent=2))


if __name__ == "__main__":
    main()
