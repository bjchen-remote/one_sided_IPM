"""Check whether 4*c(Y)/(Y*A(Y)) approaches one in native coordinates."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline
from scipy.io import loadmat

from geometry_measure import OUT, STEADY, interp_fields, trajectory_frame


HEIGHTS = np.array([0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0])


def peak_at(p, y_inner, qfield, Qfield):
    x = p["x"]
    qrow = qfield(float(y_inner), x)[0]
    cs = CubicSpline(x, qrow)
    roots = cs.derivative().roots(extrapolate=False)
    roots = roots[np.isfinite(roots) & (roots >= x[0]) & (roots <= x[-1])]
    points = np.r_[x[0], roots, x[-1]]
    m = float(points[np.argmax(cs(points))])
    return m, float(cs(m)), float(Qfield(float(y_inner), m)[0, 0])


def run():
    native = loadmat(STEADY / "run_analysis.mat", simplify_cells=True)
    raw = RectBivariateSpline(native["yFull"], native["xFull"], native["rho"])
    direct_center = float(raw(0.0, 0.0)[0, 0])
    direct_peak = float(raw(0.0, float(native["metadata"]["a"]))[0, 0])
    direct_wall_offset = direct_peak - direct_center
    frames = []
    for tau in (10, 15, 20):
        p = trajectory_frame(tau)
        qfield, Qfield = interp_fields(p)
        a = float(loadmat(OUT / f"frame_tau{tau}.mat", simplify_cells=True)["a"])
        P = float(loadmat(OUT / f"frame_tau{tau}.mat", simplify_cells=True)["P"])
        ell = float(loadmat(OUT / f"frame_tau{tau}.mat", simplify_cells=True)["ell"])

        # Recover rho(a,0)-rho(0,0) from the saved wall derivative in the same frame.
        wall_q = CubicSpline(p["x"], qfield(0.0, p["x"])[0])
        xzero = -a / ell
        wall_offset_normalized = float(wall_q.integrate(xzero, 0.0))
        wall_offset_density = P * ell * wall_offset_normalized
        rows = []
        for Y in HEIGHTS:
            y_inner = float(Y / ell)
            m, A, Qpeak = peak_at(p, y_inner, qfield, Qfield)
            c_reconstructed = P * ell * (Qpeak + wall_offset_normalized)
            ratio = 4 * c_reconstructed / (float(Y) * P * A)
            row = dict(Y=float(Y), yInner=y_inner, peakX=float(a + ell*m),
                       mInner=m, amplitudeNormalized=A, amplitudeOriginal=P*A,
                       QatPeak=Qpeak, cFromWallPeakNormalized=Qpeak,
                       rhoA0MinusRho00=float(wall_offset_density),
                       rhoA0MinusRho00Normalized=wall_offset_normalized,
                       cRhoPeakMinusRho00=float(c_reconstructed),
                       fourCOverYA=float(ratio))
            if tau == 20:
                x_original = a + ell*m
                rho_peak = float(raw(float(Y), x_original)[0, 0])
                c_direct = rho_peak - direct_center
                row.update(dict(cRhoPeakMinusRho00Direct=float(c_direct),
                                fourCOverYADirect=float(4*c_direct/(float(Y)*P*A)),
                                directVsReconstructedC=float(c_direct-c_reconstructed)))
            rows.append(row)
        frames.append(dict(tau=float(p["tau"]), source=p["source"],
                           a=a, P=P, ell=ell,
                           wallOrigin=dict(rho00DirectTau20=direct_center if tau == 20 else None,
                                           rhoA0MinusRho00DirectTau20=direct_wall_offset if tau == 20 else None,
                                           rhoA0MinusRho00FromWallIntegral=wall_offset_density,
                                           normalizedIntegral=wall_offset_normalized),
                           rows=rows))

    result = dict(title="Quarter-density relation in restored native coordinates",
                  definition="c(Y)=rho(m(Y),Y)-rho(0,0), A(Y)=max_X rho_X(X,Y), ratio=4*c/(Y*A)",
                  coordinateRestoration="X=a+ell*x, Y=ell*y, rho-rho(a,0)=P*ell*Q, rho_X=P*q",
                  originHandling="For each frame recover rho(a,0)-rho(0,0)=P*ell*integral(q(x,0), x=-a/ell..0). At tau20 this is checked against the native run_analysis rho field; the requested ratio uses the measured rho(0,0) origin.",
                  tau20NativeOrigin=dict(rho00=direct_center,
                                         rhoA0=float(direct_peak),
                                         rhoA0MinusRho00=float(direct_wall_offset)),
                  heights=HEIGHTS.tolist(), frames=frames,
                  interpretation="This is a finite-time geometry check. Values trending toward one at large native Y support the quarter relation over the sampled range; near-wall values diagnose its range of validity.")
    (OUT / "quarter_relation_checks.json").write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n")

    lines = ["# Quarter-density relation check", "",
             r"Measured quantity: $4c(Y)/(Y A(Y))$, with $c(Y)=\rho(m(Y),Y)-\rho(0,0)$ and $A(Y)=\max_X\rho_X(X,Y)$.", "",
             r"The saved inner fields were restored with each frame's own $(a,P,\ell)$. The offset from the wall-peak reference to the measured origin was recovered by integrating the saved wall gradient; for τ=20 it was also checked against the original native density field.", "",
             "| τ | Y | m(Y) in native X | A(Y) | c(Y) | 4c/(YA) |", "|---:|---:|---:|---:|---:|---:|"]
    for frame in frames:
        for row in frame["rows"]:
            lines.append(f"| {frame['tau']:.3f} | {row['Y']:g} | {row['peakX']:.7g} | {row['amplitudeOriginal']:.7g} | {row['cRhoPeakMinusRho00']:.7g} | {row['fourCOverYA']:.6f} |")
    lines += ["", "At τ=20, `fourCOverYADirect` and `cRhoPeakMinusRho00Direct` in the JSON give the direct-native cross-check. The near-wall rows are range-of-validity diagnostics; the larger-Y rows show whether the ratio moves toward one. These saved frames remain finite-time data."]
    (OUT / "QUARTER_RELATION.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"output": str(OUT), "tau20OriginOffset": direct_wall_offset,
                      "ratios": {str(f["tau"]): [round(r["fourCOverYA"], 6) for r in f["rows"]]
                                for f in frames}}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    run()
