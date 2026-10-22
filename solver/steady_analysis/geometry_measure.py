"""Inventory saved steady candidates and measure peak geometry consistently.

This is a lightweight analysis of saved fields; it does not evolve the PDE.
All candidate states use the saved inner coordinates x,y and q=Q_x.  The
native tau=20 field is converted with the centered curvature scale stored in
run_analysis.mat.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[2]
STEADY = ROOT / "runs/steady_20260928_tau20"
RUN = ROOT / "runs/run_20260927_133548_651"
OUT = ROOT / "runs/geometry_loop_20260928"
TARGETS = ["L4_129x65", "L4_257x129", "L6_193x97", "L8_257x129_continued"]
HEIGHTS = np.array([0.0, 0.25, 0.5, 1.0, 1.5, 2.0])
PROFILE_HEIGHTS = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
XI = np.linspace(-1.5, 1.5, 401)
BETA = 1.5


def trajectory_frame(tau):
    d = loadmat(OUT / f"frame_tau{tau}.mat", simplify_cells=True)
    rates = np.full(4, np.nan)
    if tau == 20:
        rates = np.array([0.60735924, 0.2219324, 0.17473202, 0.17473202])
    return {
        "name": f"tau{tau}", "kind": "native_finite_time", "source": f"frame_tau{tau}.mat",
        "x": np.asarray(d["x"]), "y": np.asarray(d["y"]),
        "Q": np.asarray(d["Q"]), "q": np.asarray(d["q"]), "rates": rates,
        "tau": float(d["tau"]),
        "coordinate_note": str(d["coordinateNote"]),
    }


def candidate(name):
    with np.load(STEADY / name / "state.npz") as d:
        return {
            "name": name, "kind": "conditional_frozen_exterior_root",
            "source": f"{name}/state.npz", "x": d["x"].copy(), "y": d["y"].copy(),
            "Q": d["Q"].copy(), "q": d["qx"].copy(),
            "rates": d["rates"].copy(),
            "coordinate_note": "saved inner x,y,Q with q=Q_x; finite-box conditional root",
        }


def interp_fields(p):
    return (RectBivariateSpline(p["y"], p["x"], p["q"], kx=3, ky=3),
            RectBivariateSpline(p["y"], p["x"], p["Q"], kx=3, ky=3))


def peak_measure(p, height, qfield, Qfield):
    x = p["x"]
    row = qfield(height, x)[0]
    cs = CubicSpline(x, row)
    roots = cs.derivative().roots(extrapolate=False)
    roots = roots[np.isfinite(roots) & (roots >= x[0]) & (roots <= x[-1])]
    points = np.r_[x[0], roots, x[-1]]
    m = float(points[np.argmax(cs(points))])
    A = float(cs(m))
    qxx = float(cs(m, 2))
    width = float(np.sqrt(-2 * A / qxx)) if qxx < 0 else float("nan")

    # A five-node local quadratic gives an interpolation sensitivity check.
    i = int(np.searchsorted(x, m))
    lo = max(0, min(i - 2, len(x) - 5))
    xx = x[lo:lo + 5]
    yy = row[lo:lo + 5]
    coef = np.polyfit(xx - m, yy, 2)
    qxx5 = float(2 * coef[0])
    w5 = float(np.sqrt(-2 * A / qxx5)) if qxx5 < 0 else float("nan")

    half_spline = CubicSpline(x, row - A / 2)
    half_roots = half_spline.roots(extrapolate=False)
    left = half_roots[half_roots < m]
    right = half_roots[half_roots > m]
    if len(left) and len(right):
        wl, wr = m - float(left[-1]), float(right[0]) - m
        delta = (wr - wl) / (wr + wl)
        skew = delta / np.tanh(BETA / 2)
        sinh_width = (wl + wr) * BETA / (2 * np.sinh(BETA))
    else:
        wl = wr = skew = sinh_width = float("nan")
    js = lambda z: float(z) if np.isfinite(z) else None
    return dict(y=float(height), m=m, A=A, qxx=qxx, w=width,
                c=float(Qfield(height, m)[0, 0]),
                halfWidthGeometry=dict(left=js(wl), right=js(wr), skew=js(skew), W=js(sinh_width)),
                interp=dict(fiveNodeQuadratic=dict(m=float(m), qxx=qxx5, w=w5)))


def sinh_geometry_profile(p, qfield, Qfield, heights=PROFILE_HEIGHTS):
    """Evaluate q/A and (Q-c)/(A W) without extrapolating beyond any field."""
    x, out = p["x"], []
    for h in heights:
        row = qfield(float(h), x)[0]
        cs = CubicSpline(x, row)
        roots = cs.derivative().roots(extrapolate=False)
        roots = roots[np.isfinite(roots) & (roots >= x[0]) & (roots <= x[-1])]
        points = np.r_[x[0], roots, x[-1]]
        m = float(points[np.argmax(cs(points))]); A = float(cs(m))
        half = CubicSpline(x, row - A / 2).roots(extrapolate=False)
        left, right = half[half < m], half[half > m]
        if not len(left) or not len(right):
            out.append(dict(y=float(h), valid=np.zeros(len(XI), dtype=bool),
                            qOverA=np.full_like(XI, np.nan),
                            QCenteredOverAW=np.full_like(XI, np.nan)))
            continue
        wl, wr = m - float(left[-1]), float(right[0]) - m
        skew = ((wr - wl) / (wr + wl)) / np.tanh(BETA / 2)
        if abs(skew) >= 1:
            out.append(dict(y=float(h), valid=np.zeros(len(XI), dtype=bool),
                            qOverA=np.full_like(XI, np.nan),
                            QCenteredOverAW=np.full_like(XI, np.nan)))
            continue
        W = (wl + wr) * BETA / (2 * np.sinh(BETA))
        G = (np.sinh(BETA * XI) + skew * (np.cosh(BETA * XI) - 1)) / BETA
        xx = m + W * G
        valid = (xx >= x[0]) & (xx <= x[-1])
        qv = np.full_like(XI, np.nan); Qv = np.full_like(XI, np.nan)
        qv[valid] = qfield(float(h), xx[valid])[0] / A
        c = float(Qfield(float(h), m)[0, 0])
        Qv[valid] = (Qfield(float(h), xx[valid])[0] - c) / (A * W)
        Gxi = np.cosh(BETA * XI) + skew * np.sinh(BETA * XI)
        out.append(dict(y=float(h), m=m, A=A, c=c, W=W, skew=float(skew),
                        wl=wl, wr=wr, cOverAW=float(c / (A * W)),
                        vOverMu=(float(p["rates"][2] / p["rates"][1])
                                  if np.isfinite(p["rates"][1]) and p["rates"][1] != 0 else None),
                        x=xx, G=G, Gxi=Gxi, qOverA=qv, QCenteredOverAW=Qv,
                        valid=valid))
    return out


def collapse_audit(all_shapes):
    train = [r for r in all_shapes["L4_257x129"] if r["valid"].any()]
    train_mask = np.logical_and.reduce([r["valid"] for r in train])
    if not train_mask.any():
        return dict(available=False, reason="training geometry has no shared physical xi support")
    master = np.full(len(XI), np.nan)
    master[train_mask] = np.mean(np.array([r["qOverA"][train_mask] for r in train]), axis=0)
    per_dataset = []
    for name, dataset in all_shapes.items():
        q_num = q_den = Q_num = Q_den = 0.0
        rows_out = []
        for row in dataset:
            mask = train_mask & row["valid"]
            if not mask.any() or not (mask & (abs(XI) < 1e-14)).any():
                rows_out.append(dict(y=row["y"], available=False,
                                     reason="no shared measured xi support with training master"))
                continue
            xi = XI[mask]
            qobs = row["qOverA"][mask]
            qhat = master[mask]
            q_num += float(np.sum((qobs - qhat)**2)); q_den += float(np.sum(qobs**2))
            integrand = qhat * row["Gxi"][mask]
            ih = np.zeros_like(integrand)
            iz = int(np.argmin(abs(xi)))
            for j in range(iz + 1, len(xi)):
                ih[j] = ih[j-1] + (integrand[j] + integrand[j-1])*(xi[j]-xi[j-1])/2
            for j in range(iz - 1, -1, -1):
                ih[j] = ih[j+1] - (integrand[j] + integrand[j+1])*(xi[j+1]-xi[j])/2
            Qobs = row["QCenteredOverAW"][mask]
            Q_num += float(np.sum((Qobs - ih)**2)); Q_den += float(np.sum(Qobs**2))
            rows_out.append(dict(y=row["y"], available=True,
                                 xiSupport=[float(xi[0]), float(xi[-1])],
                                 gradientRelativeL2=float(np.linalg.norm(qobs-qhat)/np.linalg.norm(qobs)),
                                 densityRelativeL2WithIndependentC=float(np.linalg.norm(Qobs-ih)/np.linalg.norm(Qobs)),
                                 cOverAW=float(row["cOverAW"])))
        per_dataset.append(dict(name=name,
                                aggregateGradientRelativeL2=float(np.sqrt(q_num/q_den)) if q_den else None,
                                aggregateDensityRelativeL2WithIndependentC=float(np.sqrt(Q_num/Q_den)) if Q_den else None,
                                rows=rows_out))
    return dict(available=True, training="L4_257x129; all valid measured y<=1 rows; no fitted shape parameters",
                trainingHeights=[float(r["y"]) for r in train],
                trainingXiSupport=[float(XI[train_mask][0]), float(XI[train_mask][-1])],
                masterProfile=master[train_mask].tolist(), xi=XI[train_mask].tolist(),
                note="Each test row is scored only on intersection of its measured support with the training support; no curve is extrapolated.",
                datasets=per_dataset,
                vOverMu={name:next((r["vOverMu"] for r in rows if r["vOverMu"] is not None), None)
                         for name, rows in all_shapes.items()})


def inverse_rank(p, qfield, Qfield, common_heights=HEIGHTS):
    """SVD of X(Q,y)-X(Q,0) on a common, monotone core density interval."""
    x, y = p["x"], p["y"]
    rows = []
    for h in common_heights:
        qrow = qfield(float(h), x)[0]
        cs = CubicSpline(x, qrow)
        roots = cs.derivative().roots(extrapolate=False)
        roots = roots[np.isfinite(roots) & (roots >= x[0]) & (roots <= x[-1])]
        points = np.r_[x[0], roots, x[-1]]
        m = float(points[np.argmax(cs(points))]); A = float(cs(m)); qxx = float(cs(m, 2))
        w = float(np.sqrt(-2*A/qxx)) if qxx < 0 else 0.0
        mask = (x >= m - 2*w) & (x <= m + 2*w)
        rr = Qfield(float(h), x)[0]
        xs, rs = x[mask], rr[mask]
        keep = np.r_[True, np.diff(rs) > 1e-12]
        xs, rs = xs[keep], rs[keep]
        rows.append((xs, rs))
    lo = max(row[1][0] for row in rows)
    hi = min(row[1][-1] for row in rows)
    if not hi > lo:
        return {"available": False, "reason": "no shared increasing density interval"}
    levels = np.linspace(lo, hi, 241)
    XQ = np.array([np.interp(levels, r, xx) for xx, r in rows])
    displacement = XQ - XQ[0]
    sv = np.linalg.svd(displacement, compute_uv=False)
    norm = np.linalg.norm(sv)
    return dict(available=True, densityInterval=[float(lo), float(hi)],
                relativeTailAfterRank1=float(np.linalg.norm(sv[1:]) / norm),
                relativeTailAfterRank2=float(np.linalg.norm(sv[2:]) / norm),
                relativeTailAfterRank3=float(np.linalg.norm(sv[3:]) / norm),
                singularValues=sv.tolist(), heights=common_heights.tolist())


def rel_l2(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def tail_master_audit():
    path = OUT / "tail_master_arrays.npz"
    if not path.exists():
        return {"available": False, "reason": "tail_master_arrays.npz not found"}
    with np.load(path) as d:
        xi = d["outer_xi"].copy(); heights = d["heights"].copy()
        rows = d["rows"].copy(); valid = d["valid"].copy()
        geom = d["geometry"].copy()
    near = (heights <= 1.0)
    train = heights <= 8.0
    holdout = heights >= 64.0
    master_mask = np.logical_and.reduce(valid[near])
    master = np.mean(rows[near][:, master_mask], axis=0)
    cs = CubicSpline(xi[master_mask], master)
    fprime = cs.derivative()(xi[master_mask])
    xic = xi[master_mask]
    delta = rows[:, master_mask] - master[None, :]
    raw_train = delta[train]
    raw_u, raw_s, raw_vt = np.linalg.svd(raw_train, full_matrices=False)
    raw_relative = {}
    for subset, mask in (("nearTrain", train), ("intermediate", (heights > 8) & (heights < 64)),
                         ("farHoldout", holdout)):
        mat = delta[mask]
        norm = np.linalg.norm(mat)
        raw_relative[subset] = {
            "rows": int(mask.sum()),
            "rank1RemainingRelativeL2": float(np.linalg.norm(mat - (mat @ raw_vt[:1].T) @ raw_vt[:1]) / norm),
            "rank2RemainingRelativeL2": float(np.linalg.norm(mat - (mat @ raw_vt[:2].T) @ raw_vt[:2]) / norm),
        }

    projected = np.zeros_like(delta)
    gauge_coeff = np.zeros((len(heights), 4))
    core = np.abs(xic) <= 1.0
    for j, k in enumerate(geom[:, 3]):
        G = (np.sinh(BETA*xic) + k*(np.cosh(BETA*xic)-1))/BETA
        Gxi = np.cosh(BETA*xic) + k*np.sinh(BETA*xic)
        Gk = (np.cosh(BETA*xic)-1)/BETA
        basis = np.column_stack([master, fprime/Gxi, G*fprime/Gxi, Gk*fprime/Gxi])
        coeff, *_ = np.linalg.lstsq(basis[core], delta[j, core], rcond=None)
        gauge_coeff[j] = coeff
        projected[j] = delta[j] - basis @ coeff
    proj_train = projected[train]
    _, proj_s, proj_vt = np.linalg.svd(proj_train, full_matrices=False)
    proj_relative = {}
    proj_coeffs = projected @ proj_vt[:2].T
    for subset, mask in (("nearTrain", train), ("intermediate", (heights > 8) & (heights < 64)),
                         ("farHoldout", holdout)):
        mat = projected[mask]
        norm = np.linalg.norm(mat)
        proj_relative[subset] = {
            "rows": int(mask.sum()),
            "rank1RemainingRelativeL2": float(np.linalg.norm(mat-(mat@proj_vt[:1].T)@proj_vt[:1])/norm),
            "rank2RemainingRelativeL2": float(np.linalg.norm(mat-(mat@proj_vt[:2].T)@proj_vt[:2])/norm),
        }

    # Fit only train coefficients; test the held-out far rows without refitting.
    coef_models = []
    for mode in range(2):
        y = heights[train]; z = proj_coeffs[train, mode]
        log_design = np.column_stack([np.ones_like(y), np.log1p(y)])
        log_coef, *_ = np.linalg.lstsq(log_design, z, rcond=None)
        log_pred = np.column_stack([np.ones_like(heights), np.log1p(heights)]) @ log_coef
        # Saturating law c0+c1*y/(y+s); search positive s, refit c0,c1 linearly.
        candidates = np.geomspace(1e-3, 1e3, 300)
        fits = []
        for scale in candidates:
            v = y / (y + scale)
            design = np.column_stack([np.ones_like(v), v])
            c, *_ = np.linalg.lstsq(design, z, rcond=None)
            fits.append((np.linalg.norm(design @ c - z), scale, c))
        _, scale, sat_coef = min(fits, key=lambda x: x[0])
        v_all = heights / (heights + scale)
        sat_pred = sat_coef[0] + sat_coef[1]*v_all
        coef_models.append(dict(mode=mode+1,
            logTrainRms=float(np.sqrt(np.mean((log_design @ log_coef-z)**2))),
            logHeldoutRms=float(np.sqrt(np.mean((log_pred[holdout]-proj_coeffs[holdout,mode])**2))),
            logHeldoutRelativeL2=float(np.linalg.norm(log_pred[holdout]-proj_coeffs[holdout,mode])/np.linalg.norm(proj_coeffs[holdout,mode])),
            saturationScale=float(scale),
            saturationTrainRms=float(np.sqrt(np.mean((sat_coef[0]+sat_coef[1]*(y/(y+scale))-z)**2))),
            saturationHeldoutRms=float(np.sqrt(np.mean((sat_pred[holdout]-proj_coeffs[holdout,mode])**2))),
            saturationHeldoutRelativeL2=float(np.linalg.norm(sat_pred[holdout]-proj_coeffs[holdout,mode])/np.linalg.norm(proj_coeffs[holdout,mode]))))
    return dict(available=True, source=path.name,
        caveat="All rows are one finite-time tau20 field. Far-height residuals describe this saved field and do not establish a limiting tail mode.",
        nearReference=dict(heightMax=1.0,xiSupport=[float(xic[0]),float(xic[-1])],values=master.tolist()),
        gauge=dict(fitOnlyOn="|xi|<=1",basis=["f", "f'/G_xi", "G f'/G_xi", "G_k f'/G_xi"],
                   residualSvdTrainingHeightMax=8.0, heldoutHeightMin=64.0,
                   rawSvdRelativeRemainder=raw_relative,
                   gaugeProjectedSvdRelativeRemainder=proj_relative,
                   gaugeCoefficients=gauge_coeff.tolist(),
                   projectedModeHeightCoefficients=proj_coeffs.tolist(),
                   coefficientHeightModels=coef_models),
        heights=heights.tolist(), xi=xic.tolist(), singularValuesRaw=raw_s.tolist(),
        singularValuesGaugeProjected=proj_s.tolist())


def time_tail_comparison(profiles):
    targets = [0.5, 2.0, 8.0, 64.0, 256.0]
    xit = np.linspace(-3.0, 3.0, 601)
    frames = {p["name"]: p for p in profiles if p["name"].startswith("tau")}
    rows, curves = [], {}
    for name, p in frames.items():
        qfield, Qfield = interp_fields(p)
        sampled = {}
        for h in targets:
            if h > p["y"][-1]:
                sampled[h] = None
                continue
            x = p["x"]
            qrow = qfield(h, x)[0]
            cs = CubicSpline(x, qrow)
            roots = cs.derivative().roots(extrapolate=False)
            roots = roots[np.isfinite(roots) & (roots >= x[0]) & (roots <= x[-1])]
            points = np.r_[x[0], roots, x[-1]]
            m = float(points[np.argmax(cs(points))]); A = float(cs(m))
            half = CubicSpline(x, qrow-A/2).roots(extrapolate=False)
            left, right = half[half < m], half[half > m]
            if not len(left) or not len(right):
                sampled[h] = None
                continue
            wl, wr = m-float(left[-1]), float(right[0])-m
            k = ((wr-wl)/(wr+wl))/np.tanh(BETA/2)
            W = (wl+wr)*BETA/(2*np.sinh(BETA))
            G = (np.sinh(BETA*xit)+k*(np.cosh(BETA*xit)-1))/BETA
            xx = m + W*G
            valid = (xx >= x[0]) & (xx <= x[-1])
            qnorm = np.full_like(xit, np.nan)
            qnorm[valid] = qfield(h, xx[valid])[0]/A
            sampled[h] = dict(m=m,A=A,W=W,k=float(k),qOverA=qnorm,valid=valid)
        curves[name] = sampled
    for h in targets:
        present = [name for name in ("tau10","tau15","tau20") if curves[name][h] is not None]
        pairwise = []
        for a, b in (("tau10","tau15"),("tau15","tau20"),("tau10","tau20")):
            if curves[a][h] is None or curves[b][h] is None:
                pairwise.append(dict(pair=[a,b],available=False)); continue
            A, B = curves[a][h], curves[b][h]
            mask = A["valid"] & B["valid"]
            d = A["qOverA"][mask]-B["qOverA"][mask]
            ref = A["qOverA"][mask]
            region_data = {}
            for region, use in (("core",abs(xit)<=1),("leftTail",(xit>=-3)&(xit<=-1.5)),
                                ("rightTail",(xit>=1.5)&(xit<=3)),("all",np.ones_like(xit,dtype=bool))):
                rmask=mask&use
                region_data[region]=float(np.linalg.norm(A["qOverA"][rmask]-B["qOverA"][rmask])/np.linalg.norm(A["qOverA"][rmask])) if rmask.any() else None
            pairwise.append(dict(pair=[a,b],available=True,xiSupport=[float(xit[mask][0]),float(xit[mask][-1])],
                                 relativeShapeL2ByRegion=region_data,
                                 geometryChange=dict(deltaM=B["m"]-A["m"],relativeA=B["A"]/A["A"]-1,
                                                     relativeW=B["W"]/A["W"]-1,deltaK=B["k"]-A["k"])))
        rows.append(dict(y=h,availableFrames=present,pairwise=pairwise))
    return dict(xi=[float(x) for x in xit],rows=rows,
                caveat="Each time frame uses its own saved peak, amplitude, half widths and inner scale; changes are measured after this frame-specific geometry normalization.")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    profiles = [trajectory_frame(t) for t in (10, 15, 20)] + [candidate(n) for n in TARGETS]
    output = []
    grids = {}
    geometries = {}
    for p in profiles:
        qfield, Qfield = interp_fields(p)
        heights = HEIGHTS[HEIGHTS <= min(2.0, p["y"][-1])]
        geo = [peak_measure(p, float(h), qfield, Qfield) for h in heights]
        profile_rows = sinh_geometry_profile(p, qfield, Qfield)
        rank = inverse_rank(p, qfield, Qfield, heights)
        output.append(dict(name=p["name"], kind=p["kind"], source=p["source"],
                           coordinates=p["coordinate_note"],
                           tau=p.get("tau"), rates=[float(v) if np.isfinite(v) else None for v in p["rates"]],
                           xBounds=[float(p["x"][0]), float(p["x"][-1])],
                           yBounds=[float(p["y"][0]), float(p["y"][-1])],
                           grid=[int(len(p["x"])), int(len(p["y"]))],
                           heights=geo, inverseLevelSetRank=rank))
        grids[p["name"]] = profile_rows
        geometries[p["name"]] = geo

    collapse = collapse_audit(grids)

    manifest_path = RUN / "checkpoints/checkpoint_manifest.jsonl"
    records = [json.loads(line) for line in manifest_path.read_text().splitlines() if line.strip()]
    checkpoint_frames = {}
    for target in (10.0, 15.0, 20.0):
        rec = min(records, key=lambda r: abs(float(r["canonicalTime"]) - target))
        checkpoint_frames[str(int(target))] = dict(canonicalTime=float(rec["canonicalTime"]),
                                                   step=int(rec["step"]),
                                                   file=Path(rec["checkpointFile"]).name,
                                                   exportedToSmallMat=True)

    result = dict(title="Geometry data inventory and first-pass extraction",
                  temporalFrames=checkpoint_frames,
                  limitations=[
                      "tau=10,15,20 full fields were exported from existing checkpoints into small MAT files in unified inner coordinates; tau20 was pointwise cross-checked with run_analysis.mat.",
                      "Candidate solves are conditional finite-box roots, not trajectory time frames.",
                      "Cubic-spline peak geometry is primary; local five-node quadratic is an interpolation sensitivity check.",
                      "The rank statistic uses the common increasing Q interval at the listed core heights; it is geometric, not a differential residual."],
                  profiles=output, normalizedSineWidthCollapse=collapse,
                  sameTrajectoryTailDrift=time_tail_comparison(profiles),
                  singleFrameTailResidualModes=tail_master_audit())

    # Correctly classify fixed-step mesh/box pairs from actual saved grids.
    byname = {p["name"]: p for p in profiles}
    def dx(name): return float(np.median(np.diff(byname[name]["x"])))
    result["decomposedComparisons"] = {
                  "time": "same original run and 401x321 mesh; the three exported frames use each frame's own peak-curvature scale.",
        "mesh": "L4_129x65 vs L4_257x129: same L=4 box, dx changes from %.6g to %.6g." % (dx("L4_129x65"), dx("L4_257x129")),
        "region": "L4_129x65 vs L6_193x97: same dx=%.6g, box L changes 4 to 6; L6_193x97 vs L8_257x129_continued: same dx=%.6g, box L changes 6 to 8." % (dx("L4_129x65"), dx("L6_193x97")),
        "continuedNewton": "L8_257x129 and L8_257x129_continued are same mesh/box with additional Newton iterations; included in existing inventory but only continued solution was requested for initial extraction.",
        "interpolation": "At each peak compare cubic-spline derivatives with local five-node quadratic; differences remain separately attached to each geometry row.",
    }
    json_path = OUT / "geometry_inventory.json"
    json_path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")

    arrays = {}
    for name, rows in grids.items():
        arrays[f"{name}_qOverA"] = np.array([r["qOverA"] for r in rows])
        arrays[f"{name}_QCenteredOverAW"] = np.array([r["QCenteredOverAW"] for r in rows])
        arrays[f"{name}_valid"] = np.array([r["valid"] for r in rows])
    np.savez_compressed(OUT / "geometry_profiles.npz", xi=XI, heights=PROFILE_HEIGHTS,
                        **arrays)
    lines = ["# Geometry data inventory and first-pass extraction", "",
             "All fields are represented in the shared inner coordinates. The native trajectory uses the centered peak/curvature scaling; candidate arrays already use inner x,y,Q and q=Q_x.", "",
             "| Data | source type | grid | x box | y max |", "|---|---|---:|---:|---:|"]
    for p in result["profiles"]:
        lines.append(f"| {p['name']} | {p['kind']} | {p['grid'][0]}×{p['grid'][1]} | [{p['xBounds'][0]:g},{p['xBounds'][1]:g}] | {p['yBounds'][1]:g} |")
    lines += ["", "Checkpoint frames: " + "; ".join(f"τ={k}: t={v['canonicalTime']:.8f}, step={v['step']}" for k, v in checkpoint_frames.items()), "",
              "Primary geometry and interpolation alternatives are in `geometry_inventory.json`; sampled q/A profiles are in `geometry_profiles.npz`.", "",
              "The sinh half-width collapse compares q/A on strictly common measured xi intervals. Density residuals use (Q-c)/(A W) with each measured c retained independently.", "",
              "The first level-set SVD compares X(Q,y)-X(Q,0) on a core interval of ±2 curvature widths. It summarizes geometry only.", "",
              "The time-drift table below reports q/A after each frame's measured sinh half-width remapping; regional errors are core |xi|<=1 and tails 1.5<=|xi|<=3.", "",
              "Temporal evolution, fixed-box mesh change, fixed-spacing box change, and peak interpolation sensitivity are kept as separate comparisons."]
    rows_by_y = {r["y"]: r for r in result["sameTrajectoryTailDrift"]["rows"]}
    lines += ["", "| y | tau10→15 core / left / right | tau15→20 core / left / right |", "|---:|---:|---:|"]
    for h in (0.5, 2.0, 8.0, 64.0, 256.0):
        pairs = rows_by_y[h]["pairwise"]
        vals = []
        for item in pairs[:2]:
            reg = item.get("relativeShapeL2ByRegion", {})
            vals.append(" / ".join(f"{100*reg.get(k,0):.2f}%" for k in ("core","leftTail","rightTail")) if reg else "unavailable")
        lines.append(f"| {h:g} | {vals[0]} | {vals[1]} |")
    lines += ["", "The tail-residual SVD is exploratory and uses one finite-time tau20 field; high-y rows have not converged in time, so this does not establish a persistent extra mode."]
    (OUT / "README.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"output": str(OUT), "frames": checkpoint_frames,
                      "profileSummary": [{"name": p["name"], "grid": p["grid"],
                                          "first": p["heights"][0], "rank": p["inverseLevelSetRank"]}
                                         for p in output]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
