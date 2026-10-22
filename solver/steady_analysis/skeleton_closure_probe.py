"""Prepare near-wall skeleton geometry and actual velocity samples.

The probe only reads the saved L4 fine-grid root and the tau20 trajectory.
It does not fit a closure or evolve the PDE.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline
from scipy.integrate import quad
from scipy.integrate import solve_ivp
from scipy.io import loadmat

from fixed_point import fd_matrix


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/skeleton_closure_20260928"
STEADY = ROOT / "runs/steady_20260928_tau20/L4_257x129/state.npz"
ANALYSIS = ROOT / "runs/steady_20260928_tau20/run_analysis.mat"
BETA = 1.5
Y_SAMPLES = np.linspace(0.05, 0.4, 15)
XI_FULL = np.linspace(-1.5, 1.5, 121)


def shape_at_height(x, qrow, Qfield, y, height, xi_grid=XI_FULL):
    cs = CubicSpline(x, qrow)
    roots = cs.derivative().roots(extrapolate=False)
    roots = roots[np.isfinite(roots) & (roots >= x[0]) & (roots <= x[-1])]
    points = np.r_[x[0], roots, x[-1]]
    m = float(points[np.argmax(cs(points))])
    A = float(cs(m))
    half = CubicSpline(x, qrow - A / 2).roots(extrapolate=False)
    left, right = half[half < m], half[half > m]
    if not len(left) or not len(right):
        raise ValueError(f"Half-height crossings missing at y={height:g}")
    wL, wR = m - float(left[-1]), float(right[0]) - m
    epsilon = (wR - wL) / (wR + wL)
    k = epsilon / np.tanh(BETA / 2)
    if abs(k) >= 1:
        raise ValueError(f"Sinh geometry is not globally monotone at y={height:g}: k={k:g}")
    W = (wL + wR) * BETA / (2 * np.sinh(BETA))
    G = (np.sinh(BETA * xi_grid) + k * (np.cosh(BETA * xi_grid) - 1)) / BETA
    xmap = m + W * G
    valid = (xmap >= x[0]) & (xmap <= x[-1])
    c = float(Qfield(height, m)[0, 0])
    return dict(m=m, A=A, W=float(W), k=float(k), epsilon=float(epsilon),
                c=c, wL=float(wL), wR=float(wR), G=G, x=xmap, valid=valid)


def candidate_data(name="L4_257x129"):
    state_path = ROOT / f"runs/steady_20260928_tau20/{name}/state.npz"
    with np.load(state_path) as d:
        x, y = d["x"].copy(), d["y"].copy()
        Q, q = d["Q"].copy(), d["qx"].copy()
        psi, rates = d["psi"].copy(), d["rates"].copy()
    Dx = fd_matrix(x, 1, 5)
    Dy = fd_matrix(y, 1, 5)
    u1 = -(Dy @ psi)
    u2 = psi @ Dx.T
    qy = Dy @ Q
    u1_center = float(RectBivariateSpline(y, x, u1)(0, 0)[0, 0])
    u1 -= u1_center
    return dict(name=name, kind="conditional_frozen_exterior_root",
                x=x, y=y, Q=Q, q=q, qy=qy, psi=psi, u1=u1, u2=u2,
                rates=rates, u1CenterBeforeSubtract=u1_center,
                velocityConvention="psi is normalized streamfunction; u1=-D_y psi, u2=D_x psi using the solver's five-point matrices; subtract center u1.")


def trajectory_data():
    d = loadmat(ANALYSIS, simplify_cells=True)
    meta = d["metadata"]
    a, P, ell = float(meta["a"]), float(meta["P"]), float(meta["ell"])
    x = (d["xFull"] - a) / ell
    y = d["yFull"] / ell
    u1center = float(meta["uAtPeak"])
    u1 = (d["u1"] - u1center) / (P * ell)
    u2 = d["u2"] / (P * ell)
    return dict(name="tau20", kind="native_finite_time", x=x, y=y,
                Q=(d["rho"] - RectBivariateSpline(d["yFull"], d["xFull"], d["rho"])(0, a)[0, 0]) / (P * ell),
                q=d["rhoX"] / P, u1=u1, u2=u2, rates=np.asarray(meta["innerRates"], float),
                a=a, P=P, ell=ell, rawCenterVelocity=u1center,
                velocityConvention="u1,u2 are native physical velocities converted by P*ell; subtract u1(a,0) before normalizing.")


def extract(p, xi_grid=XI_FULL):
    qfield = RectBivariateSpline(p["y"], p["x"], p["q"])
    Qfield = RectBivariateSpline(p["y"], p["x"], p["Q"])
    Qyfield = RectBivariateSpline(p["y"], p["x"], p.get("qy", np.gradient(p["Q"], p["y"], axis=0, edge_order=2)))
    u1field = RectBivariateSpline(p["y"], p["x"], p["u1"])
    u2field = RectBivariateSpline(p["y"], p["x"], p["u2"])
    lam, mu, v = map(float, p["rates"][:3])
    records = []
    arrays = {k: [] for k in ("x", "q", "Q", "Qy", "u1relative", "u2relative", "V1", "V2")}
    masks = []
    for h in Y_SAMPLES:
        shape = shape_at_height(p["x"], qfield(float(h), p["x"])[0], Qfield, float(h), float(h), xi_grid)
        xx = shape["x"]; valid = shape["valid"]
        qvalues = np.full_like(xi_grid, np.nan); Qvalues = qvalues.copy()
        Qyvalues = qvalues.copy(); ux = qvalues.copy(); uy = qvalues.copy()
        qvalues[valid] = qfield(float(h), xx[valid])[0]
        Qvalues[valid] = Qfield(float(h), xx[valid])[0]
        Qyvalues[valid] = Qyfield(float(h), xx[valid])[0]
        ux[valid] = u1field(float(h), xx[valid])[0]
        uy[valid] = u2field(float(h), xx[valid])[0]
        V1 = ux + lam * xx + v
        V2 = uy + lam * float(h)
        records.append(dict(y=float(h), **{k: shape[k] for k in ("m", "A", "W", "k", "epsilon", "c", "wL", "wR")},
                            maxAbsXMapped=float(np.max(abs(xx))),
                            mappedXBounds=[float(xx.min()), float(xx.max())],
                            validXiBounds=[float(xi_grid[valid][0]), float(xi_grid[valid][-1])] if valid.any() else None,
                            maxAbsRelativeU1=float(np.max(abs(ux[valid]))),
                            maxAbsRelativeU2=float(np.max(abs(uy[valid]))),
                            maxAbsV1=float(np.max(abs(V1[valid]))), maxAbsV2=float(np.max(abs(V2[valid])))))
        arrays["x"].append(np.where(valid, xx, np.nan))
        arrays["q"].append(qvalues)
        arrays["Q"].append(Qvalues)
        arrays["Qy"].append(Qyvalues)
        arrays["u1relative"].append(ux)
        arrays["u2relative"].append(uy)
        arrays["V1"].append(V1)
        arrays["V2"].append(V2)
        masks.append(valid)
    common = np.logical_and.reduce(masks)
    if common.sum() < 20 or not np.any(common & (abs(xi_grid) < 1e-12)):
        raise ValueError(f"Too little common measured xi support for {p['name']}: {int(common.sum())}")
    arrays = {k: np.asarray(v)[:, common] for k, v in arrays.items()}
    return dict(name=p["name"], kind=p["kind"],
                grid=[len(p["x"]), len(p["y"])], xBounds=[float(p["x"][0]), float(p["x"][-1])],
                yBounds=[float(p["y"][0]), float(p["y"][-1])], rates=[float(z) for z in p["rates"]],
                velocityConvention=p["velocityConvention"],
                residualVelocityConvention="V1=u1_rel+lambda*x+v; V2=u2_rel+lambda*y, matching the centered fixed-point operator.",
                xi=xi_grid[common], xiSupport=[float(xi_grid[common][0]), float(xi_grid[common][-1])],
                records=records, arrays=arrays)


def primitive(values, xi):
    antiderivative = CubicSpline(xi, values).antiderivative()
    return antiderivative(xi) - antiderivative(0.0)


def wall_learned_shape():
    p = candidate_data("L4_257x129")
    qfield = RectBivariateSpline(p["y"], p["x"], p["q"])
    Qfield = RectBivariateSpline(p["y"], p["x"], p["Q"])
    shape = shape_at_height(p["x"], qfield(0, p["x"])[0], Qfield, 0.0, 0.0)
    return qfield(0, shape["x"])[0] / shape["A"], shape


def f0_rational(xi):
    return 32 / (xi**6 - 3*xi**5 + 3*xi**4 + 3*xi**3 + 28*xi**2 + 32)


def dae_projection(profile, f, train_phase=0):
    xi, recs, arrays = profile["xi"], profile["records"], profile["arrays"]
    fs = np.asarray(f(xi) if callable(f) else f, float)
    n = len(xi); train = np.arange(n) % 2 == train_phase; test = ~train
    lam, mu, v, _ = profile["rates"]
    per_height = []
    M_rows, C_rows, r0_rows, E_rows, zprime_rows = [], [], [], [], []
    qhat_rows, Qhat_rows, baseline_rows = [], [], []
    for j, rec in enumerate(recs):
        k, A, W, c = rec["k"], rec["A"], rec["W"], rec["c"]
        G = (np.sinh(BETA*xi) + k*(np.cosh(BETA*xi)-1))/BETA
        Gk = (np.cosh(BETA*xi)-1)/BETA
        h = primitive(fs*(np.cosh(BETA*xi)+k*np.sinh(BETA*xi)), xi)
        hk = primitive(fs*np.sinh(BETA*xi), xi)
        V1, V2 = arrays["V1"][j], arrays["V2"][j]
        D = A*W
        cols = np.column_stack([-A*fs, D*h, D*(h-fs*G), D*(hk-fs*Gk), np.ones(n)])
        M = V2[:, None]*cols
        r0 = A*fs*V1 - mu*(c+D*h) - v
        scales = np.maximum(np.linalg.norm(M[train], axis=0), np.finfo(float).tiny)
        Ms = M[train]/scales[None, :]
        coeff, _, rank, _ = np.linalg.lstsq(Ms, -r0[train], rcond=None)
        zprime = coeff/scales
        residual = M@zprime + r0
        hold = residual[test]
        base = V1*arrays["q"][j] + V2*arrays["Qy"][j] - mu*arrays["Q"][j] - v
        qhat = A*fs
        Qhat = c+D*h
        measured_scale = np.sqrt(np.linalg.norm(V1*arrays["q"][j])**2+
                                 np.linalg.norm(V2*arrays["Qy"][j])**2+
                                 np.linalg.norm(mu*arrays["Q"][j])**2+(abs(v)*np.sqrt(n))**2)
        Qyhat = cols@zprime
        skeleton_scale = np.sqrt(np.linalg.norm(V1*qhat)**2+
                                 np.linalg.norm(V2*Qyhat)**2+
                                 np.linalg.norm(mu*Qhat)**2+(abs(v)*np.sqrt(n))**2)
        hold_scale = np.sqrt(np.linalg.norm((V1*qhat)[test])**2+
                             np.linalg.norm((V2*Qyhat)[test])**2+
                             np.linalg.norm((mu*Qhat)[test])**2+(abs(v)*np.sqrt(test.sum()))**2)
        singular = np.linalg.svd(Ms, compute_uv=False)
        per_height.append(dict(y=float(rec["y"]), rank=int(rank),
            condition=float(singular[0]/singular[-1]) if singular[-1] else None,
            singularValues=singular.tolist(), scaledColumnNorms=scales.tolist(), zprime=zprime.tolist(),
            trainResidualRelativeL2=float(np.linalg.norm(residual[train])/np.linalg.norm(r0[train])),
            holdoutResidualRelativeL2=float(np.linalg.norm(hold)/np.linalg.norm(r0[test])),
            holdoutResidualLinf=float(np.max(abs(hold))),
            fullResidualRelativeL2=float(np.linalg.norm(residual)/np.linalg.norm(r0)),
            holdoutTermScaleRelativeL2=float(np.linalg.norm(hold)/max(hold_scale,np.finfo(float).tiny)),
            fullTermScaleRelativeL2=float(np.linalg.norm(residual)/max(skeleton_scale,np.finfo(float).tiny)),
            measuredTransportBaselineRelativeL2=float(np.linalg.norm(base)/max(measured_scale,np.finfo(float).tiny)),
            measuredTransportBaselineLinf=float(np.max(abs(base))),
            measuredTransportRms=float(np.sqrt(np.mean(base**2))),
            densityReconstructionRelativeL2=float(np.linalg.norm(Qhat-arrays["Q"][j])/np.linalg.norm(arrays["Q"][j])),
            densityReconstructionRelativeToAW=float(np.linalg.norm(Qhat-arrays["Q"][j])/max(D*np.sqrt(n),np.finfo(float).tiny)),
            wallCenterIdentity=float(Qhat[int(np.argmin(abs(xi)))]-c)))
        M_rows.append(M); C_rows.append(cols); r0_rows.append(r0); E_rows.append(residual); zprime_rows.append(zprime)
        qhat_rows.append(qhat); Qhat_rows.append(Qhat); baseline_rows.append(base)
    return dict(metrics=per_height,
                arrays=dict(M=np.asarray(M_rows), fiveBasisColumns=np.asarray(C_rows),
                            r0=np.asarray(r0_rows), Eperp=np.asarray(E_rows),
                            zprime=np.asarray(zprime_rows), qhat=np.asarray(qhat_rows),
                            Qhat=np.asarray(Qhat_rows), measuredTransportResidual=np.asarray(baseline_rows),
                            trainMask=train, xi=xi))


def integrate_fixed_velocity_l4(profile, f, label):
    """Integrate the five-geometry DAE while holding saved L4 velocity fixed."""
    p = candidate_data("L4_257x129")
    xi = profile["xi"]
    dense_xi = np.linspace(-1.5, 1.5, 1201)
    fdense = np.asarray(f(dense_xi) if callable(f) else f, float)
    f_spline = CubicSpline(dense_xi, fdense, extrapolate=False)
    # h(k)=int f(cosh+k*sinh); hk=int f*sinh. High-resolution cubic
    # antiderivatives keep primitive error below the xi-grid error.
    hc_spline = CubicSpline(dense_xi, fdense*np.cosh(BETA*dense_xi)).antiderivative()
    hs_spline = CubicSpline(dense_xi, fdense*np.sinh(BETA*dense_xi)).antiderivative()
    u1field = RectBivariateSpline(p["y"], p["x"], p["u1"])
    u2field = RectBivariateSpline(p["y"], p["x"], p["u2"])
    lam, mu, v = map(float, p["rates"][:3])
    xlo, xhi = float(p["x"][0]), float(p["x"][-1])
    ylo, yhi = float(p["y"][0]), float(p["y"][-1])
    train = np.arange(len(xi)) % 2 == 0
    test = ~train
    first = profile["records"][0]
    z0 = np.array([first["m"], np.log(first["A"]), np.log(first["W"]), first["k"], first["c"]])
    domain_margins = []
    def geometry(z):
        m, logA, logW, k, c = z
        A, W = np.exp(logA), np.exp(logW)
        G = (np.sinh(BETA*xi)+k*(np.cosh(BETA*xi)-1))/BETA
        Gk = (np.cosh(BETA*xi)-1)/BETA
        xx = m+W*G
        return m, A, W, k, c, G, Gk, xx
    def rhs(y, z):
        m, A, W, k, c, G, Gk, xx = geometry(z)
        if abs(k) >= 1 or xx.min() < xlo or xx.max() > xhi or y < ylo or y > yhi:
            raise ValueError("ODE stage left the measured coordinate domain")
        fxi = f_spline(xi)
        hc = hc_spline(xi)-hc_spline(0.0)
        hs = hs_spline(xi)-hs_spline(0.0)
        h = hc+k*hs
        hk = hs
        V1 = u1field(float(y), xx)[0]+lam*xx+v
        V2 = u2field(float(y), xx)[0]+lam*float(y)
        D = A*W
        cols = np.column_stack([-A*fxi, D*h, D*(h-fxi*G), D*(hk-fxi*Gk), np.ones(len(xi))])
        M = V2[:, None]*cols
        r0 = A*fxi*V1-mu*(c+D*h)-v
        scales = np.maximum(np.linalg.norm(M[train], axis=0), np.finfo(float).tiny)
        coeff, _, rank, _ = np.linalg.lstsq(M[train]/scales, -r0[train], rcond=None)
        if rank < 5:
            raise ValueError(f"ODE DAE rank dropped to {rank}")
        domain_margins.append([float(y), float(xx.min()-xlo), float(xhi-xx.max()), float(1-abs(k))])
        return coeff/scales
    def boundary_event(y, z):
        try:
            _, _, _, k, _, _, _, xx = geometry(z)
            return min(float(xx.min()-xlo), float(xhi-xx.max()), float(1-abs(k)), float(y-ylo), float(yhi-y))
        except Exception:
            return -1.0
    boundary_event.terminal = True
    boundary_event.direction = -1
    ode_error = None
    try:
        sol = solve_ivp(rhs, (float(Y_SAMPLES[0]), float(Y_SAMPLES[-1])), z0,
                        method="RK45", t_eval=Y_SAMPLES, events=boundary_event,
                        max_step=0.003, rtol=2e-7, atol=1e-9)
    except (ValueError, FloatingPointError) as exc:
        sol = None
        ode_error = str(exc)
    if sol is None or not sol.success or len(sol.t) != len(Y_SAMPLES):
        return dict(label=label, status="terminated_or_failed", error=ode_error or (sol.message if sol else None),
                    completedHeights=[] if sol is None else sol.t.tolist(),
                    fixedVelocityNecessaryTest=True)
    metrics = []
    for i, (y,z) in enumerate(zip(sol.t, sol.y.T)):
        m, A, W, k, c, G, Gk, xx = geometry(z)
        rec = profile["records"][i]
        h = (hc_spline(xi)-hc_spline(0.0))+k*(hs_spline(xi)-hs_spline(0.0))
        qhat = A*f_spline(xi); Qhat = c+A*W*h
        qfield = RectBivariateSpline(p["y"], p["x"], p["q"])
        Qfield = RectBivariateSpline(p["y"], p["x"], p["Q"])
        qmeas = qfield(float(y), xx)[0]; Qmeas = Qfield(float(y), xx)[0]
        m_err = (m-rec["m"])/rec["W"]
        metrics.append(dict(y=float(y), mErrorInW=float(m_err), ARelativeError=float(A/rec["A"]-1),
                            WRelativeError=float(W/rec["W"]-1), kError=float(k-rec["k"]),
                            cErrorInAW=float((c-rec["c"])/(rec["A"]*rec["W"])),
                            qReconstructionRelativeL2=float(np.linalg.norm(qhat-qmeas)/np.linalg.norm(qmeas)),
                            QReconstructionRelativeL2=float(np.linalg.norm(Qhat-Qmeas)/np.linalg.norm(Qmeas)),
                            densityCenterError=float(Qhat[np.argmin(abs(xi))]-c),
                            xMappedBounds=[float(xx.min()), float(xx.max())],
                            validA=bool(A>0), validW=bool(W>0), validK=bool(abs(k)<1)))
    return dict(label=label, status="success", fixedVelocityNecessaryTest=True,
                interpretation="Geometry ODE driven by saved L4 velocities held fixed; not a full Biot-Savart closure.",
                metrics=metrics,
                minimumObservedMargins={"xLeft":float(min(m[1] for m in domain_margins)),
                                        "xRight":float(min(m[2] for m in domain_margins)),
                                        "monotonicity":float(min(m[3] for m in domain_margins))})


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    profiles = [candidate_data("L4_257x129"), candidate_data("L8_257x129_continued"), trajectory_data()]
    results = [extract(p) for p in profiles]
    wall_f, wall_geometry = wall_learned_shape()
    wall_spline = CubicSpline(XI_FULL, wall_f, extrapolate=False)
    shape_families = {
        "fixed_integer_rational_f0": f0_rational,
        "learned_L4_wall_f": lambda x: wall_spline(x),
    }
    projections = {}
    projection_arrays = {}
    for profile in results:
        projections[profile["name"]] = {}
        for label, f in shape_families.items():
            projected = dae_projection(profile, f, train_phase=0)
            projections[profile["name"]][label] = projected["metrics"]
            projection_arrays[f"{profile['name']}_{label}"] = projected["arrays"]
    l4_ode = [integrate_fixed_velocity_l4(results[0], f, label) for label, f in shape_families.items()]
    refinement = {}
    l4raw = profiles[0]
    l4_dense = extract(l4raw, np.linspace(-1.5, 1.5, 241))
    l4_index = int(np.argmin(abs(Y_SAMPLES-0.2)))
    for label, f in shape_families.items():
        detail = {}
        for density, prof in ((121, results[0]), (241, l4_dense)):
            for phase in (0, 1):
                projected = dae_projection(prof, f, train_phase=phase)
                row = projected["metrics"][l4_index]
                detail[f"grid{density}_trainParity{phase}"] = {k:row[k] for k in (
                    "rank", "condition", "holdoutResidualRelativeL2", "holdoutTermScaleRelativeL2",
                    "fullTermScaleRelativeL2", "measuredTransportBaselineRelativeL2",
                    "densityReconstructionRelativeL2", "densityReconstructionRelativeToAW")}
        refinement[label] = detail

    # Independently verify the analytic f0 primitives against adaptive quadrature.
    quad_check = {}
    k_mid = l4_dense["records"][l4_index]["k"]
    for density, prof in ((121, results[0]), (241, l4_dense)):
        xi_grid = prof["xi"]
        f_grid = f0_rational(xi_grid)
        h_grid = primitive(f_grid*(np.cosh(BETA*xi_grid)+k_mid*np.sinh(BETA*xi_grid)), xi_grid)
        hk_grid = primitive(f_grid*np.sinh(BETA*xi_grid), xi_grid)
        points = xi_grid[np.linspace(0, len(xi_grid)-1, 17).round().astype(int)]
        h_error, hk_error = [], []
        for z in points:
            hq = quad(lambda s:f0_rational(s)*(np.cosh(BETA*s)+k_mid*np.sinh(BETA*s)), 0, z, epsabs=1e-12, epsrel=1e-12)[0]
            hkq = quad(lambda s:f0_rational(s)*np.sinh(BETA*s), 0, z, epsabs=1e-12, epsrel=1e-12)[0]
            j = int(np.argmin(abs(xi_grid-z)))
            h_error.append(float(h_grid[j]-hq)); hk_error.append(float(hk_grid[j]-hkq))
        quad_check[str(density)] = dict(sampleCount=len(points), hMaxAbsError=float(max(abs(np.asarray(h_error)))),
                                         hkMaxAbsError=float(max(abs(np.asarray(hk_error)))))
    quad_check["y"] = float(Y_SAMPLES[l4_index]); quad_check["k"] = float(k_mid)
    tau20_direct = loadmat(ANALYSIS, simplify_cells=True)
    keys = {"state_npz": sorted(np.load(STEADY).files),
            "run_analysis_mat": sorted(k for k in tau20_direct if not k.startswith("__"))}
    summaries = []
    for r in results:
        summary = {k:v for k,v in r.items() if k not in ("arrays", "xi")}
        summary["xi"] = r["xi"].tolist()
        summaries.append(summary)
    result = dict(title="Near-wall geometry skeleton samples and relative velocity data",
                  inspectedFields=keys,
                  interpretation={
                      "candidate": "state.npz saves psi but no u1/u2; velocities were reconstructed from psi using the same five-point derivative matrices as the solver.",
                      "trajectory": "run_analysis.mat saves native u1/u2; they were converted to inner units with P*ell and centered by u1(a,0).",
                      "sampleDomain": "0.05<=y<=0.4; per-source common measured xi support is used, with no extrapolation. L8 continued and tau20 cover the full requested xi interval; L4 fine is clipped by its x box.",
                      "rateConvention": "rates=(lambda,mu,v,b), with b=v. The candidate rates are the converged finite-box root; tau20 rates are the final-frame instantaneous estimates.",
                      "velocityConvention": "u_rel is the center-subtracted incompressible velocity; V adds the similarity frame terms used in the discrete steady residual.",
                      "dae":"fixed_velocity_necessary_test: fit five geometry derivatives on alternating xi nodes and score on held-out nodes; this is a necessary kinematic test with saved velocities, not a full Biot-Savart closure."},
                  data=summaries, daeProjection=projections,
                  fixedVelocityL4ODE=l4_ode,
                  xiRefinementAndParityExchange=refinement,
                  f0PrimitiveQuadratureCheck=quad_check,
                  learnedWall=dict(m=wall_geometry["m"], A=wall_geometry["A"], W=wall_geometry["W"],
                                   k=wall_geometry["k"], xi=XI_FULL.tolist(), f=wall_f.tolist()))
    (OUT / "skeleton_samples.json").write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    payload = {"y": Y_SAMPLES, "xi_full": XI_FULL}
    for r in results:
        payload[f"{r['name']}_xi"] = r["xi"]
        for k,v in r["arrays"].items():
            payload[f"{r['name']}_{k}"] = v
    np.savez_compressed(OUT / "skeleton_samples.npz", **payload)
    flat_projection = {}
    for name, families in projection_arrays.items():
        for key, value in families.items():
            flat_projection[f"{name}_{key}"] = value
    np.savez_compressed(OUT / "skeleton_projection_samples.npz", **flat_projection)
    lines = ["# Skeleton closure sample inventory", "",
             "Candidate `state.npz` includes `psi`, `Q`, `qx`, rates, and residual, but no velocity fields. Native `run_analysis.mat` includes `u1`, `u2`, and `psi`.", "",
             "The candidate velocity was reconstructed as `u1=-D_y psi`, `u2=D_x psi` with the solver's five-point matrices, then centered at the wall peak. The trajectory velocities were converted by `P*ell` and centered by the saved `uAtPeak`.", "",
             "The DAE projection tests `q=A f(xi)`, `Q=c+AW h` with a fixed rational master and the measured L4 wall profile. L4 fine uses its common measured xi interval; L8 continued and tau20 cover xi=[-1.5,1.5]. Projection results include scaled rank/condition, alternating-node holdout residuals, the measured-field transport residual baseline, and direct density reconstruction errors.", "",
             "`skeleton_projection_samples.npz` stores the five basis columns, `M`, `r0`, the full orthogonal residual `Eperp`, fitted parameter derivatives `zprime`, and predicted/measured Q fields for every height and both fixed shape functions. A 121/241-point xi refinement and train/holdout parity exchange are recorded at y≈0.2; the f0 primitive is checked against adaptive quadrature.", "",
             "The L4 five-geometry ODE is integrated with saved velocities held fixed. Both shape choices complete y=.05 to .4 while remaining in-domain; geometry drift and direct q/Q reconstruction errors are recorded. This is a fixed-velocity necessary test only.", "",
             "`wall_symbolic_family.json` is produced by `solver/steady_analysis/wall_symbolic_family_probe.py` and reports q,Q train/holdout errors on |x|≤1.5 and q≥1/2 masks for L4, L8, and tau20.", "",
             "`fixed_velocity_necessary_test` uses saved velocities as fixed spatial fields. It is not a full Biot-Savart closure and the finite-time tau20 field is only a control.", "",
             "`V1=u1_rel+lambda*x+v` and `V2=u2_rel+lambda*y` are saved for comparison with the centered fixed-point transport operator.", "",
             "The JSON records contain measured `m,A,W,k,c`, `lambda,mu,v,b`, coordinate bounds, and velocity magnitudes. The NPZ stores point samples of `q,Q,u1_rel,u2_rel,V1,V2`."]
    (OUT / "README.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"out": str(OUT), "inventory": keys,
                      "available": [{"name": r["name"], "grid": r["grid"],
                                     "rates": r["rates"], "maxSampleHeight": max(d["y"] for d in r["records"]),
                                     "xiSupport": r["xiSupport"]}
                                    for r in results],
                      "projection": {name:{f:{"condition":max(r["condition"] for r in rows),
                                            "holdout":float(np.mean([r["holdoutResidualRelativeL2"] for r in rows])),
                                            "baseline":float(np.mean([r["measuredTransportBaselineRelativeL2"] for r in rows]))}
                                        for f,rows in fields.items()} for name,fields in projections.items()},
                      "quadCheck":quad_check},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
