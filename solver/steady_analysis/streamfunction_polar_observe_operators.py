"""Observe the streamfunction and the velocity it generates on native circles.

No angular basis or radial law is fitted. This is finite-time observation,
not a steady closure test. Reuses the saved streamfunction at tau=20.
"""
from pathlib import Path
import json
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ipm_observation_mpl")

import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import brentq, minimize_scalar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from streamfunction_sources import load_native


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/streamfunction_polar_observation_20260929"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    p = load_native(20)
    f = RectBivariateSpline(p["y"], p["x"], p["psi"])
    radii = np.array([.25, 1., 4., 16., 64., 256., 1024., 4096.])
    theta = np.linspace(0, np.pi, 801)
    rr, tt = np.meshgrid(radii, theta, indexing="ij")
    xx, yy = rr*np.cos(tt)/p["ell"], rr*np.sin(tt)/p["ell"]
    assert xx.min() >= p["x"][0] and xx.max() <= p["x"][-1]
    assert yy.min() >= p["y"][0] and yy.max() <= p["y"][-1]
    psi = p["P"]*p["ell"]**2*f.ev(yy, xx)
    psi_x = p["P"]*p["ell"]*f.ev(yy, xx, dy=1)
    psi_y = p["P"]*p["ell"]*f.ev(yy, xx, dx=1)
    psi_r = np.cos(tt)*psi_x+np.sin(tt)*psi_y
    psi_theta = rr*(-np.sin(tt)*psi_x+np.cos(tt)*psi_y)
    angular_mean = lambda value: np.trapezoid(value, theta, axis=1)/np.pi
    amp = np.sqrt(angular_mean(psi**2))
    phi = psi/amp[:, None]
    degree = angular_mean(psi*rr*psi_r)/amp**2
    euler = rr*psi_r/amp[:, None]
    drift = euler-degree[:, None]*phi
    ur_scaled = -psi_theta/amp[:, None]
    ut_shape = euler/degree[:, None]

    np.savez_compressed(OUT/"polar_velocity_observation.npz",
        native_radii=radii, theta=theta, amplitude=amp, degree=degree,
        phi=phi, radial_velocity_scaled=ur_scaled,
        azimuthal_velocity_shape=ut_shape, angular_drift=drift)
    result = dict(source=p["source"], tau=20, ell=p["ell"], P=p["P"],
        gauge=p["gauge"], noPdeEvolution=True, noFittedBasis=True,
        definitions=dict(phi="Psi/A, A^2=mean_theta(Psi^2)",
            p="d log A / d log R = mean(Psi R Psi_R)/A^2",
            radial="R u_R/A = -phi_theta",
            azimuthal="R u_theta/(p A) = phi + phi_logR/p",
            drift="phi_logR = R Psi_R/A - p phi"),
        rows=[dict(nativeR=float(r),innerR=float(r/p["ell"]),
            amplitude=float(amp[i]), degree=float(degree[i]),
            driftOrthogonality=float(angular_mean(phi*drift)[i]),
            radialVelocitySignChanges=int(np.count_nonzero(
                ur_scaled[i, :-1]*ur_scaled[i, 1:] < 0)))
            for i,r in enumerate(radii)],
        scope="Spline derivatives of a saved finite-time field. Endpoints are sampled, not prescribed. No extrapolation or exact closure claim.")
    (OUT/"polar_velocity_observation.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False)+"\n")

    # An observation-selected comparison: the same measured zero angle occurs
    # twice. Query the original Cartesian spline independently of polar-shell
    # interpolation, and compare shape quantities without fitting a family.
    zero = np.deg2rad(36.)
    def at(radius, angle):
        return f.ev(radius*np.sin(angle)/p["ell"],
                    radius*np.cos(angle)/p["ell"])
    comparison = []
    for bracket in ((.25,.5), (64.,128.)):
        radius = brentq(lambda r: float(at(r, zero)), *bracket, xtol=1e-12)
        positive = minimize_scalar(lambda t: -float(at(radius,t)),
            bounds=(0,zero),method="bounded",options={"xatol":1e-13}).x
        negative = minimize_scalar(lambda t: float(at(radius,t)),
            bounds=(zero,np.pi),method="bounded",options={"xatol":1e-13}).x
        def kernel(angle):
            return at(radius,angle)/(np.sin(angle)*np.sin(zero-angle))
        comparison.append(dict(nativeR=float(radius),innerR=radius/p["ell"],
            zeroAngleDegrees=36.,positivePeakAngle=float(np.degrees(positive)),
            negativePeakAngle=float(np.degrees(negative)),
            positiveOverNegativePeak=float(-at(radius,positive)/at(radius,negative)),
            G15OverG90=float(kernel(np.deg2rad(15))/kernel(np.pi/2)),
            G165OverG90=float(kernel(np.deg2rad(165))/kernel(np.pi/2))))
    (OUT/"same_zero_direct_check.json").write_text(json.dumps(dict(
        method="Direct evaluation of the saved tau20 Cartesian streamfunction spline; root of Psi(R,36 deg)=0. Independent of interpolation between stored polar shells.",
        source=p["source"],gauge=p["gauge"],rows=comparison),indent=2)+"\n")

    plt.rcParams.update({"font.size":10,"axes.spines.top":False,
                         "axes.spines.right":False})
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.5), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(.04, .93, len(radii)))
    quantities = [phi, ur_scaled, ut_shape, drift]
    titles = ["Streamfunction: one small positive lobe, one large negative lobe",
              "Radial velocity: angular derivative of the observed shape",
              "Azimuthal velocity: shape plus its radial evolution",
              "Angular shape evolution with log radius"]
    labels = [r"$\phi=\Psi/A$",r"$R u_R/A=-\phi_\theta$",
              r"$R u_\theta/(pA)$",r"$\partial_{\log R}\phi$"]
    for ax, values, title, label in zip(axes.flat, quantities, titles, labels):
        for i,r in enumerate(radii):
            ax.plot(np.degrees(theta),values[i],color=colors[i],
                    label=f"R={r:g}",lw=1.6)
        ax.axhline(0,color=".65",lw=.7)
        ax.set(xlabel="Angle (degrees)",ylabel=label,title=title,xlim=(0,180))
        ax.grid(alpha=.15)
    axes[0,0].legend(ncol=2,fontsize=8)
    fig.suptitle("Direct polar observation | tau=20 | native R | no fitted angular functions",fontsize=14)
    fig.savefig(OUT/"polar_velocity_observation.png",dpi=180)
    plt.close(fig)
    print("Saved direct polar streamfunction/velocity observations for",len(radii),"radii.")
    print("Maximum normalized amplitude/drift inner product:",
          max(abs(angular_mean(phi*drift))))


if __name__ == "__main__":
    main()
