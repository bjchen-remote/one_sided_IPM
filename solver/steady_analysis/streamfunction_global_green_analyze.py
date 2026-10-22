"""Analyse the centered-Green branch and its critical Mellin closure.

This script is deliberately downstream of the nonlinear solves.  It does not
alter or re-solve a state.  It checks three claims that are needed for the
infinite-domain interpretation:

1. the finite-box eigenvalue approaches alpha=1/2 with the logarithmic
   detuning expected at the N=2, m=1 Mellin resonance;
2. rays from the physical wall core admit a bounded core-to-far-field
   Mellin--Pade representation, including the resonant r log r term;
3. shifted-centre far-field rays share one r log r sin(theta) coefficient,
   rather than requiring an independent logarithmic coefficient on every ray.

The output is a machine-readable certificate plus publication-style figures.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import numpy as np
from scipy.interpolate import RectBivariateSpline, UnivariateSpline
from scipy.optimize import minimize_scalar


ROOT = Path(__file__).resolve().parents[2]
RUN_ROOT = ROOT / "runs" / "streamfunction_global_green_20260929"


def read_json(path: Path):
    with path.open() as stream:
        return json.load(stream)


def root_path(length: int, nx: int, ny: int, stencil: int, version: int = 1):
    return RUN_ROOT / (
        f"greenroot_L{length}_{nx}x{ny}_stencil{stencil}_q4_v{version}")


def load_sequence(entries):
    rows = []
    for length, path in entries:
        report = read_json(path / "report.json")
        metrics = report["finalMetrics"]
        rows.append(dict(length=float(length), path=str(path), **metrics))
    return rows


def linear_fit(x, y, intercept=True):
    x, y = np.asarray(x, float), np.asarray(y, float)
    matrix = (np.column_stack([np.ones_like(x), x])
              if intercept else x[:, None])
    coefficients = np.linalg.lstsq(matrix, y, rcond=None)[0]
    prediction = matrix @ coefficients
    return coefficients, prediction, float(np.sqrt(np.mean((prediction - y)**2)))


def finite_size_analysis(sequences):
    finest = sequences["257x129_s7"]
    length = np.array([row["length"] for row in finest])
    alpha = np.array([row["alpha"] for row in finest])
    log_length = np.log(length)
    coordinate = log_length**-2
    constrained, constrained_prediction, constrained_rms = linear_fit(
        coordinate, alpha - 0.5, intercept=False)
    unconstrained, unconstrained_prediction, unconstrained_rms = linear_fit(
        coordinate, alpha, intercept=True)
    compensated = (alpha - 0.5) * log_length**2

    rate_fits = {}
    for name, key in [("lambda", "lambda_"), ("mu", "mu"), ("v", "v")]:
        values = np.array([row[key] for row in finest])
        all_fit, _, all_rms = linear_fit(coordinate, values)
        last_fit, _, last_rms = linear_fit(coordinate[-3:], values[-3:])
        rate_fits[name] = dict(
            allLimit=float(all_fit[0]), allSlope=float(all_fit[1]),
            allRms=all_rms, lastThreeLimit=float(last_fit[0]),
            lastThreeSlope=float(last_fit[1]), lastThreeRms=last_rms)

    leave_one_out = []
    for omitted in range(len(length)):
        keep = np.arange(len(length)) != omitted
        fit, _, _ = linear_fit(coordinate[keep], alpha[keep])
        leave_one_out.append(float(fit[0]))

    return dict(
        length=length, alpha=alpha, logLength=log_length,
        coordinate=coordinate, compensated=compensated,
        constrainedConstant=float(constrained[0]),
        constrainedPrediction=0.5 + constrained_prediction,
        constrainedRms=constrained_rms,
        unconstrainedLimit=float(unconstrained[0]),
        unconstrainedSlope=float(unconstrained[1]),
        unconstrainedPrediction=unconstrained_prediction,
        unconstrainedRms=unconstrained_rms,
        leaveOneOutLimits=leave_one_out,
        rateFits=rate_fits)


def load_state(path: Path):
    with np.load(path) as data:
        return dict(
            x=np.asarray(data["x"], float),
            y=np.asarray(data["y"], float),
            Q=np.asarray(data["Q"], float),
            psi=np.asarray(data["psi"], float),
            parameters=np.asarray(data["parameters"], float))


def sample_rays(state, angles, radii, shifted):
    x, y = state["x"], state["y"]
    q_spline = RectBivariateSpline(y, x, state["Q"])
    psi_spline = RectBivariateSpline(y, x, state["psi"])
    lam, mu, drift, offset = state["parameters"]
    shift = drift / lam if shifted else 0.0
    theta = np.asarray(angles, float)[:, None]
    radius = np.asarray(radii, float)[None, :]
    query_x = radius * np.cos(theta) - shift
    query_y = radius * np.sin(theta)
    q = q_spline.ev(query_y.ravel(), query_x.ravel()).reshape(query_x.shape)
    psi = psi_spline.ev(
        query_y.ravel(), query_x.ravel()).reshape(query_x.shape)
    p = q + offset / mu
    return dict(Q=q, P=p, psi=psi, x=query_x, y=query_y,
                shift=float(shift))


def fit_core_rays(state, angles, radii):
    """Fit a bounded resonant Mellin--Pade model on physical-core rays."""
    sampled = sample_rays(state, angles, radii, shifted=False)
    psi = sampled["psi"]
    radii = np.asarray(radii, float)
    angles = np.asarray(angles, float)
    lam, mu, _, _ = state["parameters"]
    sigma = lam - mu
    core_coefficient = (
        sigma * np.sin(angles) * np.cos(angles)
        - 0.5 * np.sin(angles)**2)
    observed = psi / radii[None, :]**2
    train = np.arange(radii.size) % 2 == 0
    validate = ~train
    scales = np.sqrt(np.mean(observed[:, train]**2, axis=1))
    scales = np.maximum(scales, 1e-10 * np.max(scales))

    def basis(length, resonant):
        s = radii / length
        leading = (1.0 + s * s)**(-0.25)
        regular = (1.0 + s * s)**(-0.5)
        if not resonant:
            return leading, regular, None
        logarithmic = np.ones_like(s)
        nonzero = s != 0.0
        logarithmic[nonzero] = np.arcsinh(s[nonzero]) / s[nonzero]
        return leading, regular, logarithmic

    def solve_at(length, resonant, mask):
        leading, regular, logarithmic = basis(length, resonant)
        coefficients = []
        fitted = np.empty_like(observed)
        for index, coefficient in enumerate(core_coefficient):
            if resonant:
                design = np.column_stack(
                    [leading - logarithmic, regular - logarithmic])
                target = observed[index] - coefficient * logarithmic
                local = np.linalg.lstsq(
                    design[mask] / scales[index],
                    target[mask] / scales[index], rcond=None)[0]
                model = (coefficient * logarithmic
                         + design @ local)
                coefficients.append([local[0], local[1],
                                     coefficient - local.sum()])
            else:
                design = (leading - regular)[:, None]
                target = observed[index] - coefficient * regular
                local = np.linalg.lstsq(
                    design[mask] / scales[index],
                    target[mask] / scales[index], rcond=None)[0]
                model = coefficient * regular + design[:, 0] * local[0]
                coefficients.append([local[0], coefficient - local[0]])
            fitted[index] = model
        return np.asarray(coefficients), fitted

    def objective(log_length, resonant):
        length = np.exp(log_length)
        _, fitted = solve_at(length, resonant, train)
        residual = (fitted[:, train] - observed[:, train]) / scales[:, None]
        return float(np.sqrt(np.mean(residual**2)))

    fits = {}
    for name, resonant in [("nonresonant", False), ("resonant", True)]:
        optimum = minimize_scalar(
            lambda value: objective(value, resonant),
            bounds=(np.log(0.05), np.log(256.0)), method="bounded",
            options={"xatol": 1e-8})
        length = float(np.exp(optimum.x))
        coefficients, fitted = solve_at(length, resonant, train)
        train_residual = ((fitted[:, train] - observed[:, train])
                          / scales[:, None])
        validation_residual = ((fitted[:, validate] - observed[:, validate])
                               / scales[:, None])
        physical_relative = np.linalg.norm(
            fitted * radii[None, :]**2 - psi) / np.linalg.norm(psi)
        fits[name] = dict(
            crossoverLength=length, coefficients=coefficients,
            fittedNormalized=fitted,
            fittedPsi=fitted * radii[None, :]**2,
            trainNormalizedRms=float(np.sqrt(np.mean(train_residual**2))),
            validationNormalizedRms=float(
                np.sqrt(np.mean(validation_residual**2))),
            physicalRelativeL2=float(physical_relative))
    return dict(
        angles=angles, radii=radii, sampled=sampled,
        observedNormalized=observed, coreCoefficient=core_coefficient,
        fits=fits)


def block_ray_fit(values, radii, angles, kind, train, alpha):
    """Fit shared-coefficient Mellin models to a matrix of ray values."""
    number_angles, number_radii = values.shape
    rows = number_angles * number_radii
    angle_index = np.repeat(np.arange(number_angles), number_radii)
    radius = np.tile(radii, number_angles)
    theta = np.repeat(angles, number_radii)
    decay = 1.0 - alpha
    if kind == "psi":
        target = (values / radii[None, :]**(1.0 + alpha)).ravel()
        matrix = np.zeros((rows, 3 * number_angles + 1))
        matrix[np.arange(rows), angle_index] = 1.0
        matrix[np.arange(rows), number_angles + angle_index] = radius**(-decay)
        matrix[np.arange(rows), 2 * number_angles + angle_index] = radius**(-2.0 * decay)
        # The only growing Dirichlet harmonic compatible with the selected
        # class is c*y.  At alpha close to 1/2 it is almost collinear with the
        # m=1 part of the N=2 Mellin term; keeping it global avoids assigning a
        # separate fake logarithmic coefficient to every ray.
        matrix[:, -1] = np.sin(theta) * radius**(-alpha)
        no_log_columns = 3 * number_angles
    elif kind == "density":
        target = (values / radii[None, :]**alpha).ravel()
        matrix = np.zeros((rows, 3 * number_angles))
        matrix[np.arange(rows), angle_index] = 1.0
        matrix[np.arange(rows), number_angles + angle_index] = radius**(-decay)
        matrix[np.arange(rows), 2 * number_angles + angle_index] = radius**(-2.0 * decay)
        no_log_columns = matrix.shape[1]
    else:
        raise ValueError(kind)
    train_rows = np.tile(train, number_angles)
    validation_rows = ~train_rows
    coefficient = np.linalg.lstsq(
        matrix[train_rows], target[train_rows], rcond=None)[0]
    fitted = matrix @ coefficient
    scale = max(float(np.sqrt(np.mean(target**2))), 1e-300)
    result = dict(
        coefficients=coefficient,
        fitted=fitted.reshape(number_angles, number_radii),
        trainRelativeRms=float(np.sqrt(np.mean(
            (fitted[train_rows] - target[train_rows])**2)) / scale),
        validationRelativeRms=float(np.sqrt(np.mean(
            (fitted[validation_rows] - target[validation_rows])**2)) / scale))
    if kind == "psi":
        no_log = np.linalg.lstsq(
            matrix[train_rows, :no_log_columns], target[train_rows],
            rcond=None)[0]
        no_log_fitted = matrix[:, :no_log_columns] @ no_log
        result.update(
            noHarmonicValidationRelativeRms=float(np.sqrt(np.mean(
                (no_log_fitted[validation_rows]
                 - target[validation_rows])**2)) / scale),
            harmonicCoefficient=float(coefficient[-1]),
            B1=coefficient[:number_angles],
            B2=coefficient[number_angles:2 * number_angles],
            B3=coefficient[2 * number_angles:3 * number_angles])
    else:
        result.update(
            A1=coefficient[:number_angles],
            A2=coefficient[number_angles:2 * number_angles])
    return result


def angular_poisson_prediction(theta, angular_density, beta):
    """Solve one angular Poisson equation for density exponent beta."""
    dense = np.linspace(0.0, np.pi, 4001)
    # A very small smoothing level suppresses radial-fit noise while retaining
    # the narrow quadratic minimum near theta=pi/2.
    scale = max(float(np.linalg.norm(angular_density)), 1e-12)
    spline = UnivariateSpline(
        theta, angular_density, k=3, s=2e-5 * scale**2)
    dense_a = spline(dense)
    dense_ap = spline.derivative()(dense)
    forcing = beta * np.cos(dense) * dense_a - np.sin(dense) * dense_ap
    modes = np.arange(1, 81)
    sine = np.sin(np.outer(dense, modes))
    forcing_modes = (2.0 / np.pi) * np.trapezoid(
        forcing[:, None] * sine, dense, axis=0)
    kappa = beta + 1.0
    b_modes = forcing_modes / (modes * modes - kappa * kappa)
    dense_b = sine @ b_modes
    dense_bp = (np.cos(np.outer(dense, modes))
                @ (modes * b_modes))
    return dict(
        denseTheta=dense, A1=dense_a, A1Prime=dense_ap,
        beta=float(beta), kappa=float(kappa), forcing=forcing,
        B=dense_b, BPrime=dense_bp,
        forcingModes=forcing_modes, BModes=b_modes)


def fit_far_rays(state, lambda_limit, radius_minimum=24.0,
                 radius_maximum=512.0):
    lam, mu, _, _ = state["parameters"]
    alpha = mu / lam
    decay = 1.0 - alpha
    angles = np.linspace(0.0, np.pi, 73)
    radii = np.geomspace(radius_minimum, radius_maximum, 121)
    sampled = sample_rays(state, angles, radii, shifted=True)
    train = np.arange(radii.size) % 2 == 0
    density = block_ray_fit(
        sampled["P"], radii, angles, "density", train, alpha)
    stream = block_ray_fit(
        sampled["psi"], radii, angles, "psi", train, alpha)

    poisson = angular_poisson_prediction(angles, density["A1"], alpha)
    b1_predicted = np.interp(
        angles, poisson["denseTheta"], poisson["B"])
    b1_data = stream["B1"]
    b1_relative = float(
        np.linalg.norm(b1_predicted - b1_data)
        / max(np.linalg.norm(b1_data), 1e-300))

    a1 = poisson["A1"]
    a1p = poisson["A1Prime"]
    b1 = poisson["B"]
    b1p = poisson["BPrime"]
    a2 = (-alpha * b1p * a1 + (1.0 + alpha) * b1 * a1p) / (
        lam * decay)
    a2_data = density["A2"]
    a2_predicted_on_rays = np.interp(
        angles, poisson["denseTheta"], a2)
    a2_relative = float(
        np.linalg.norm(a2_predicted_on_rays - a2_data)
        / max(np.linalg.norm(a2_data), 1e-300))

    beta2 = 2.0 * alpha - 1.0
    poisson2 = angular_poisson_prediction(
        poisson["denseTheta"], a2, beta2)
    forced_projection = float(poisson2["forcingModes"][0])
    b2_sine = float(poisson2["BModes"][0])
    predicted_log = -0.5 * forced_projection

    theta_matrix = angles[:, None]
    radius_matrix = radii[None, :]
    b1_on_rays = np.interp(
        angles, poisson["denseTheta"], poisson["B"])
    b2_on_rays = np.interp(
        angles, poisson2["denseTheta"], poisson2["B"])
    leading = (radius_matrix**(1.0 + alpha)
               * b1_on_rays[:, None])
    second = radius_matrix**(2.0 * alpha) * b2_on_rays[:, None]
    harmonic = radius_matrix * np.sin(theta_matrix)
    train_matrix = np.tile(train, (angles.size, 1))

    def fit_harmonic(base):
        coefficient = float(np.dot(
            harmonic[train_matrix],
            (sampled["psi"] - base)[train_matrix])
            / np.dot(harmonic[train_matrix], harmonic[train_matrix]))
        model = base + coefficient * harmonic
        validation = ~train_matrix
        relative = float(
            np.linalg.norm((model - sampled["psi"])[validation])
            / max(np.linalg.norm(sampled["psi"][validation]), 1e-300))
        return coefficient, model, relative

    c_leading, leading_model, leading_error = fit_harmonic(leading)
    c_second, second_model, second_error = fit_harmonic(leading + second)
    return dict(
        alpha=float(alpha), decay=float(decay), beta2=float(beta2),
        angles=angles, radii=radii, sampled=sampled,
        density=density, stream=stream, poisson=poisson, poisson2=poisson2,
        B1Predicted=b1_predicted, B1PoissonRelativeL2=b1_relative,
        A2Predicted=a2, A2PredictedOnRays=a2_predicted_on_rays,
        A2RecurrenceRelativeL2=a2_relative,
        forcedProjection=forced_projection,
        B2SineCoefficient=b2_sine,
        predictedLogCoefficient=float(predicted_log),
        leadingHarmonicCoefficient=c_leading,
        leadingValidationRelativeL2=leading_error,
        secondHarmonicCoefficient=c_second,
        secondValidationRelativeL2=second_error,
        leadingModel=leading_model, secondModel=second_model)


def plot_finite_size(output, sequences, analysis, arclength):
    figure, axes = plt.subplots(2, 2, figsize=(13.5, 9.5), constrained_layout=True)
    colors = dict(s7="#0068b5", s5="#d04a35", coarse="#888888")
    labels = {
        "257x129_s7": "257x129, 7-point",
        "321x161_s7": "321x161, 7-point (check)",
        "193x97_s7": "193x97, 7-point",
        "193x97_s5": "193x97, 5-point",
    }
    for key, rows in sequences.items():
        length = np.array([row["length"] for row in rows])
        alpha = np.array([row["alpha"] for row in rows])
        style = ("#148f77" if key.startswith("321") else
                 colors["s7"] if "s7" in key else colors["s5"])
        marker = ("*" if key.startswith("321") else
                  "o" if key.startswith("257") else "s")
        axes[0, 0].plot(np.log(length)**-2, alpha, marker=marker,
                        color=style, label=labels[key], lw=1.3)
        axes[0, 1].plot(length, (alpha - 0.5) * np.log(length)**2,
                        marker=marker, color=style, label=labels[key], lw=1.3)
    xline = np.linspace(0.0, 1.04 * analysis["coordinate"].max(), 200)
    axes[0, 0].plot(
        xline, 0.5 + analysis["constrainedConstant"] * xline,
        "k--", lw=1.6,
        label=r"$1/2+C/(\log L)^2$")
    axes[0, 0].axhline(0.5, color="0.5", lw=1.0, ls=":")
    axes[0, 0].set(xlabel=r"$1/(\log L)^2$", ylabel=r"$\alpha(L)=\mu/\lambda$",
                   title="Critical exponent selection")
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].grid(alpha=0.25)
    axes[0, 1].axhline(
        analysis["constrainedConstant"], color="k", ls="--", lw=1.4,
        label=f"C={analysis['constrainedConstant']:.4f}")
    axes[0, 1].set_xscale("log", base=2)
    axes[0, 1].set(xlabel="box half-size L",
                   ylabel=r"$(\alpha-1/2)(\log L)^2$",
                   title="Compensated logarithmic drift")
    axes[0, 1].legend(fontsize=8)
    axes[0, 1].grid(alpha=0.25)

    finest = sequences["257x129_s7"]
    length = np.array([row["length"] for row in finest])
    coordinate = np.log(length)**-2
    lam = np.array([row["lambda_"] for row in finest])
    mu = np.array([row["mu"] for row in finest])
    drift = np.array([row["v"] for row in finest])
    axes[1, 0].plot(coordinate, lam, "o-", label=r"$\lambda$")
    axes[1, 0].plot(coordinate, 2.0 * mu, "s-", label=r"$2\mu$")
    axes[1, 0].plot(coordinate, 5.0 * drift, "^-", label=r"$5v$")
    axes[1, 0].set(xlabel=r"$1/(\log L)^2$", ylabel="scaled rate",
                   title="Rates (extrapolation is less rigid than alpha)")
    axes[1, 0].grid(alpha=0.25)
    axes[1, 0].legend()

    branch = arclength["branch"]
    axes[1, 1].plot(
        [row["lambda_"] for row in branch],
        [row["tau"] for row in branch], "o-", ms=3.0, lw=1.2,
        color="#6a3d9a")
    maximum = max(branch, key=lambda row: row["tau"])
    axes[1, 1].plot(maximum["lambda_"], maximum["tau"], "r*", ms=13,
                    label=f"fold: tau={maximum['tau']:.6f}")
    axes[1, 1].axhline(1.0, color="0.4", ls=":")
    axes[1, 1].set(xlabel=r"$\lambda$", ylabel=r"homotopy $\tau$",
                   title="Native t=20 branch folds before the Green root")
    axes[1, 1].grid(alpha=0.25)
    axes[1, 1].legend(fontsize=8)
    figure.savefig(output, dpi=190)
    plt.close(figure)


def plot_global_rays(output, state, core_fit):
    x, y, psi = state["x"], state["y"], state["psi"]
    lam, mu, drift, _ = state["parameters"]
    figure = plt.figure(figsize=(14.5, 10.0), constrained_layout=True)
    grid = figure.add_gridspec(2, 2)
    axis = figure.add_subplot(grid[:, 0])
    radius = np.hypot(x[None, :], y[:, None])
    normalized = psi / (1.0 + radius**1.5)
    limit = np.nanpercentile(np.abs(normalized), 99.0)
    mesh = axis.pcolormesh(
        x, y, normalized, shading="auto", cmap="coolwarm",
        norm=SymLogNorm(max(limit * 2e-3, 1e-8), vmin=-limit, vmax=limit))
    axis.set_xscale("symlog", linthresh=1.0)
    axis.set_yscale("symlog", linthresh=1.0)
    axis.set_xlim(x[0], x[-1])
    axis.set_ylim(0.0, y[-1])
    for angle in core_fit["angles"]:
        end = 0.72 * min(x[-1], y[-1])
        axis.plot([0.0, end * np.cos(angle)], [0.0, end * np.sin(angle)],
                  color="k", lw=0.35, alpha=0.35)
    axis.plot(0.0, 0.0, "ko", ms=4, label="physical core")
    axis.plot(-drift / lam, 0.0, "k*", ms=8, label="similarity centre")
    axis.set(xlabel="x (symlog)", ylabel="y (symlog)",
             title=r"Global Green solution: $\psi/(1+r^{3/2})$")
    axis.legend(loc="upper right", fontsize=8)
    figure.colorbar(mesh, ax=axis, shrink=0.78)

    chosen = [0, 2, 4, 6, 8, 10]
    colors = plt.cm.turbo(np.linspace(0.05, 0.95, len(chosen)))
    signed_axis = figure.add_subplot(grid[0, 1])
    error_axis = figure.add_subplot(grid[1, 1])
    for panel, quantity in [(signed_axis, "signed"),
                            (error_axis, "error")]:
        for color, index in zip(colors, chosen):
            angle = np.rad2deg(core_fit["angles"][index])
            radius_values = core_fit["radii"]
            observed = core_fit["observedNormalized"][index]
            fitted = core_fit["fits"]["resonant"]["fittedNormalized"][index]
            if quantity == "signed":
                panel.plot(radius_values, observed, color=color, lw=1.35,
                           label=f"{angle:.0f} deg")
                panel.plot(radius_values, fitted, color=color, ls="--", lw=1.1)
            else:
                scale = max(np.sqrt(np.mean(observed**2)), 1e-12)
                panel.plot(radius_values, (fitted - observed) / scale,
                           color=color, lw=1.1, label=f"{angle:.0f} deg")
        panel.set_xscale("log")
        panel.grid(alpha=0.25)
    signed_axis.set(
        ylabel=r"$\psi/r^2$",
        title="Core rays: data (solid), bounded resonant fit (dashed)")
    signed_axis.legend(ncol=2, fontsize=8)
    error_axis.axhline(0.0, color="k", lw=0.7)
    error_axis.set(
        xlabel="radius from physical core", ylabel="normalized fit defect",
        title="Same fit, signed defect across the full radial window")
    figure.savefig(output, dpi=190)
    plt.close(figure)


def plot_angular_closure(output, far_fit):
    figure, axes = plt.subplots(2, 2, figsize=(13.5, 9.5), constrained_layout=True)
    degrees = np.rad2deg(far_fit["angles"])
    axes[0, 0].plot(degrees, far_fit["density"]["A1"], "o-", ms=3,
                    label=r"radial fit $A_1$")
    axes[0, 0].plot(
        np.rad2deg(far_fit["poisson"]["denseTheta"]),
        far_fit["poisson"]["A1"], "k--", lw=1.1,
        label="smoothed angular profile")
    axes[0, 0].set(title=r"$P=r^\alpha A_1+r^{2\alpha-1}A_2+\cdots$",
                   xlabel="theta (deg)", ylabel=r"$A_1$")
    axes[0, 0].legend(fontsize=8)

    axes[0, 1].plot(degrees, far_fit["stream"]["B1"], "o", ms=3,
                    label=r"radial fit $B_1$")
    axes[0, 1].plot(degrees, far_fit["B1Predicted"], "-", lw=1.4,
                    label="Poisson prediction from A1")
    axes[0, 1].set(title="Leading angular Poisson check (finite window)",
                   xlabel="theta (deg)", ylabel=r"$B_1$")
    axes[0, 1].legend(fontsize=8)

    chosen = [8, 16, 24, 32, 40, 48, 56, 64]
    colors = plt.cm.viridis(np.linspace(0.05, 0.95, len(chosen)))
    radii = far_fit["radii"]
    for color, index in zip(colors, chosen):
        observed = (far_fit["sampled"]["psi"][index]
                    / radii**(1.0 + far_fit["alpha"]))
        fitted = far_fit["stream"]["fitted"][index]
        axes[1, 0].plot(radii, observed, color=color, lw=1.1,
                        label=f"{degrees[index]:.0f} deg")
        axes[1, 0].plot(radii, fitted, color=color, ls="--", lw=1.0)
    axes[1, 0].set_xscale("log")
    axes[1, 0].set(title="Shifted-centre rays: finite-alpha Mellin fit",
                   xlabel="r", ylabel=r"$\psi/r^{1+\alpha}$")
    axes[1, 0].legend(ncol=2, fontsize=8)

    residual_no_log = far_fit["leadingValidationRelativeL2"]
    residual_log = far_fit["secondValidationRelativeL2"]
    axes[1, 1].bar(
        ["N=1 + c y", "N=1,2 + c y"], [residual_no_log, residual_log],
        color=["#999999", "#2b8cbe"])
    axes[1, 1].set_yscale("log")
    axes[1, 1].set(ylabel="held-out relative RMS",
                   title=("Mellin PDE reconstruction, held-out radii\n"
                          f"B2,m=1={far_fit['B2SineCoefficient']:.3g}, "
                          f"c={far_fit['secondHarmonicCoefficient']:.3g}"))
    for axis in axes.ravel():
        axis.grid(alpha=0.23)
    figure.savefig(output, dpi=190)
    plt.close(figure)


def serializable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {key: serializable(item) for key, item in value.items()
                if key not in {"sampled", "poisson", "fittedNormalized",
                               "fittedPsi", "observedNormalized"}}
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out", type=Path,
        default=RUN_ROOT / "final_closure")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    sequences = {
        "257x129_s7": load_sequence([
            (length, root_path(length, 257, 129, 7))
            for length in [256, 512, 1024, 2048, 4096]]),
        "321x161_s7": load_sequence([
            (4096, root_path(4096, 321, 161, 7))]),
        "193x97_s7": load_sequence([
            (length, root_path(length, 193, 97, 7, 2 if length == 4096 else 1))
            for length in [64, 128, 256, 512, 1024, 2048, 4096]]),
        "193x97_s5": load_sequence([
            (length, root_path(
                length, 193, 97, 5, 2 if length == 512 else 1))
            for length in [64, 128, 256, 512, 1024, 2048, 4096]]),
    }
    finite = finite_size_analysis(sequences)
    high_resolution = sequences["321x161_s7"][0]
    finite["highResolutionL4096"] = dict(
        alpha=float(high_resolution["alpha"]),
        compensated=float(
            (high_resolution["alpha"] - 0.5)
            * np.log(high_resolution["length"])**2),
        residualInf=float(high_resolution["residualInf"]))
    final_path = root_path(4096, 321, 161, 7)
    state = load_state(final_path / "state.npz")
    core_angles = np.deg2rad([10, 25, 40, 55, 70, 85, 100, 120, 140, 160, 175])
    core_radii = np.geomspace(0.125, 512.0, 401)
    core_fit = fit_core_rays(state, core_angles, core_radii)
    lambda_limit = finite["rateFits"]["lambda"]["lastThreeLimit"]
    # At larger fractions of the finite box the omitted exterior source is
    # already visible in psi.  The 4--64 annulus is the widest tested window
    # on which the first two PDE-generated Mellin orders improve the held-out
    # reconstruction; it is therefore reported explicitly, not hidden by an
    # automatically selected favourable interval.
    far_fit = fit_far_rays(state, lambda_limit, 4.0, 64.0)
    arclength = read_json(
        RUN_ROOT / "native_L64_129x65_stencil7_q2_arclength_v6"
        / "report.json")

    plot_finite_size(
        args.out / "critical_scaling.png", sequences, finite, arclength)
    plot_global_rays(
        args.out / "global_solution_core_rays.png", state, core_fit)
    plot_angular_closure(
        args.out / "mellin_resonant_closure.png", far_fit)

    with (args.out / "ray_global_fit.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "angle_deg", "core_coefficient", "leading_amplitude",
            "regular_r_amplitude", "log_amplitude", "crossover_length"])
        coefficients = core_fit["fits"]["resonant"]["coefficients"]
        length = core_fit["fits"]["resonant"]["crossoverLength"]
        for angle, core, coefficient in zip(
                core_angles, core_fit["coreCoefficient"], coefficients):
            writer.writerow([
                np.rad2deg(angle), core, coefficient[0], coefficient[1],
                coefficient[2], length])

    maximum_fold = max(arclength["branch"], key=lambda row: row["tau"])
    certificate = dict(
        formulation=(
            "full transport + Poisson + centered half-plane Green closure"),
        finalRoot=dict(
            path=str(final_path),
            parameters=state["parameters"].tolist(),
            report=read_json(final_path / "report.json")),
        finiteSize=serializable(finite),
        coreRayFit=dict(
            anglesDeg=np.rad2deg(core_angles).tolist(),
            radiusWindow=[float(core_radii[0]), float(core_radii[-1])],
            nonresonant=serializable(core_fit["fits"]["nonresonant"]),
            resonant=serializable(core_fit["fits"]["resonant"])),
        farRayFit=dict(
            shiftedCentre=-float(state["parameters"][2]
                                 / state["parameters"][0]),
            radiusWindow=[float(far_fit["radii"][0]),
                          float(far_fit["radii"][-1])],
            densityTrainRelativeRms=far_fit["density"]["trainRelativeRms"],
            densityValidationRelativeRms=(
                far_fit["density"]["validationRelativeRms"]),
            directStreamFitTrainRelativeRms=(
                far_fit["stream"]["trainRelativeRms"]),
            directStreamFitValidationRelativeRms=(
                far_fit["stream"]["validationRelativeRms"]),
            fittedHarmonicCoefficient=(
                far_fit["stream"]["harmonicCoefficient"]),
            predictedLogCoefficient=far_fit["predictedLogCoefficient"],
            forcedProjection=far_fit["forcedProjection"],
            B1PoissonRelativeL2=far_fit["B1PoissonRelativeL2"],
            A2RecurrenceRelativeL2=far_fit["A2RecurrenceRelativeL2"],
            B2SineCoefficient=far_fit["B2SineCoefficient"],
            leadingHarmonicCoefficient=(
                far_fit["leadingHarmonicCoefficient"]),
            leadingValidationRelativeL2=(
                far_fit["leadingValidationRelativeL2"]),
            secondHarmonicCoefficient=(
                far_fit["secondHarmonicCoefficient"]),
            secondValidationRelativeL2=(
                far_fit["secondValidationRelativeL2"])),
        nativeBranchFold=dict(
            tau=maximum_fold["tau"], lambda_=maximum_fold["lambda_"],
            alpha=maximum_fold["alpha"], reachesPhysicalTauOne=False),
        figures=[
            "critical_scaling.png", "global_solution_core_rays.png",
            "mellin_resonant_closure.png"],
    )
    with (args.out / "certificate.json").open("w") as stream:
        json.dump(certificate, stream, indent=2)
        stream.write("\n")
    print(json.dumps(certificate, indent=2))


if __name__ == "__main__":
    main()
