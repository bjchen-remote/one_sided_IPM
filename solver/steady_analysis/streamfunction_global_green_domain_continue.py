"""Continue a centered-Green self-similar IPM root to another box/grid.

The source is interpolated in its resolved interior.  Outside that region it
is extended along rays from the shifted similarity centre with the leading
homogeneities Q+b/mu ~ r**alpha and psi ~ r**(1+alpha).  This extension is
only a Newton seed; the target always solves the full transport, Poisson, and
centered half-plane Green boundary equations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
from scipy.interpolate import RectBivariateSpline

from fixed_point import fd_matrix
from streamfunction_global_green_solve import CenteredGreenBVP, solve
from streamfunction_global_selfsimilar_solve import (
    NewtonOptions,
    dump,
    stretched_grid,
)


def load_state(path: Path):
    with np.load(path) as data:
        return (
            np.asarray(data["x"], float),
            np.asarray(data["y"], float),
            np.asarray(data["Q"], float),
            np.asarray(data["psi"], float),
            np.asarray(data["parameters"], float),
        )


def source_center_derivative(x, y, psi, stencil):
    origin = int(np.argmin(np.abs(x)))
    row = fd_matrix(y, 1, stencil).getrow(0)
    return float((row @ psi[:, origin]).item())


def enforce_core_gauges(problem, q, psi, parameters):
    envelope = np.exp(-(problem.X / 1.5)**8 - (problem.Y / 1.5)**4)
    basis = np.column_stack([
        (envelope * problem.X**degree).ravel() for degree in range(4)
    ])
    qflat = q.ravel().copy()
    qflat += basis @ np.linalg.solve(
        problem.C @ basis, problem.target - problem.C @ qflat)
    psiflat = psi.ravel().copy()
    center = float((problem.Dy.getrow(problem.origin) @ psiflat).item())
    psiflat -= center * problem.yy
    adjusted = np.asarray(parameters, float).copy()
    adjusted[2] -= center
    return qflat, psiflat, adjusted, center


def transferred_seed(problem, source_path, blend_start, blend_end, stencil):
    sx, sy, sq, spsi, parameters = load_state(source_path)
    if parameters.size != 4:
        raise ValueError("centered-Green source must have four parameters")
    lam, mu, drift, offset = map(float, parameters)
    if not (lam > mu > 0):
        raise ValueError("source rates must satisfy 0 < mu < lambda")

    # Remove the finite-difference remnant of the forbidden linear harmonic
    # mode while preserving u_x + drift exactly in the seed.
    source_center = source_center_derivative(sx, sy, spsi, stencil)
    spsi = spsi - source_center * sy[:, None]
    parameters = parameters.copy()
    parameters[2] -= source_center
    lam, mu, drift, offset = map(float, parameters)

    alpha = mu / lam
    power = 1.0 + alpha
    shift = drift / lam
    source_radius = min(
        -sx[0] - abs(shift), sx[-1] - abs(shift), sy[-1])
    if source_radius <= 0:
        raise ValueError("source box does not surround the shifted centre")
    r0 = blend_start * source_radius
    r1 = blend_end * source_radius
    if not (0.0 < r0 < r1 < source_radius):
        raise ValueError("need 0 < blend-start < blend-end < 1")

    qspline = RectBivariateSpline(sy, sx, sq)
    psispline = RectBivariateSpline(sy, sx, spsi)
    z = problem.X + shift
    radius = np.hypot(z, problem.Y)
    angle = np.arctan2(problem.Y, z)
    x0 = r0 * np.cos(angle) - shift
    y0 = r0 * np.sin(angle)
    q0 = qspline.ev(y0.ravel(), x0.ravel()).reshape(problem.shape)
    psi0 = psispline.ev(y0.ravel(), x0.ravel()).reshape(problem.shape)
    scale = np.maximum(radius / r0, 1e-300)
    q_outer = scale**alpha * (q0 + offset / mu) - offset / mu
    psi_outer = scale**power * psi0

    q_direct = qspline.ev(
        problem.Y.ravel(), problem.X.ravel()).reshape(problem.shape)
    psi_direct = psispline.ev(
        problem.Y.ravel(), problem.X.ravel()).reshape(problem.shape)
    target_inside = bool(
        problem.x[0] >= sx[0] and problem.x[-1] <= sx[-1]
        and problem.y[-1] <= sy[-1])
    if target_inside:
        q, psi = q_direct, psi_direct
        mode = "interpolation"
    else:
        blend = np.clip((radius - r0) / (r1 - r0), 0.0, 1.0)
        blend = blend * blend * (3.0 - 2.0 * blend)
        q = (1.0 - blend) * q_direct + blend * q_outer
        psi = (1.0 - blend) * psi_direct + blend * psi_outer
        mode = "interpolation_plus_homogeneous_ray_extension"
    qflat, psiflat, parameters, target_center = enforce_core_gauges(
        problem, q, psi, parameters)
    state = np.r_[qflat, psiflat, parameters]
    metadata = dict(
        source=str(source_path.resolve()),
        sourceBox=[float(sx[0]), float(sx[-1]), float(sy[-1])],
        sourceCenterDerivative=source_center,
        targetInterpolationCenterCorrection=target_center,
        shiftedCenter=float(-shift), usableSourceRadius=float(source_radius),
        blendRadius=[float(r0), float(r1)], transferMode=mode)
    return state, metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--length", type=float, required=True)
    parser.add_argument("--nx", type=int, default=129)
    parser.add_argument("--ny", type=int, default=65)
    parser.add_argument("--core-step", type=float, default=0.125)
    parser.add_argument("--stencil", type=int, default=5)
    parser.add_argument("--quadrature-order", type=int, default=4)
    parser.add_argument("--blend-start", type=float, default=0.65)
    parser.add_argument("--blend-end", type=float, default=0.85)
    parser.add_argument("--iterations", type=int, default=40)
    parser.add_argument("--tolerance", type=float, default=2e-8)
    parser.add_argument("--taus", default="")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError(f"refusing to overwrite {args.out}")
    args.out.mkdir(parents=True)

    sx, sy, _, _, source_parameters = load_state(args.source)
    x, y, ax, ay = stretched_grid(
        args.length, args.nx, args.ny, args.core_step)
    problem = CenteredGreenBVP(
        x, y, source_parameters[:4], args.stencil,
        args.quadrature_order)
    state, transfer = transferred_seed(
        problem, args.source, args.blend_start, args.blend_end,
        args.stencil)
    initial = problem.metrics(state)
    config = dict(
        transfer=transfer, length=args.length, nx=args.nx, ny=args.ny,
        requestedCoreStep=args.core_step,
        actualCoreStepX=float(x[problem.origin + 1]),
        actualCoreStepY=float(y[1]), stretchX=ax, stretchY=ay,
        stencil=args.stencil, quadratureOrder=args.quadrature_order,
        sourceGrid=[len(sx), len(sy)], initialMetrics=initial,
        scriptSha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    dump(args.out / "config.json", config)
    print(json.dumps(config, indent=2), flush=True)

    initial_residual = problem.evaluate(state)[0]
    taus = ([float(value) for value in args.taus.split(",")]
            if args.taus else [1.0])
    if args.taus and (
            taus[0] != 0.0 or taus[-1] != 1.0
            or any(b <= a for a, b in zip(taus, taus[1:]))):
        raise ValueError("taus must increase strictly from 0 to 1")
    rows = []
    status = "not_started"
    final_tau = None
    start = time.time()
    with (args.out / "iterations.jsonl").open("w") as stream:
        for tau in taus:
            rhs = ((1.0 - tau) * initial_residual
                   if args.taus else np.zeros(problem.size))

            def callback(row):
                tagged = dict(homotopyTau=tau, **row)
                rows.append(tagged)
                stream.write(json.dumps(tagged) + "\n")
                stream.flush()
                print(
                    f"tau={tau:.5g} iter={row['iteration']:02d} "
                    f"res={row['residualInf']:.3e} "
                    f"rates=({row['lambda_']:.9f},{row['mu']:.9f},"
                    f"{row['v']:.9f}) alpha={row['alpha']:.9f}",
                    flush=True)

            state, status, _ = solve(
                problem, state,
                NewtonOptions(args.iterations, args.tolerance),
                right_hand_side=rhs, callback=callback)
            final_tau = tau
            np.savez_compressed(
                args.out / f"stage_tau_{tau:.5f}.npz",
                x=x, y=y,
                Q=state[:problem.N].reshape(problem.shape),
                psi=state[problem.N:2 * problem.N].reshape(problem.shape),
                parameters=state[2 * problem.N:], homotopyTau=tau)
            if status != "converged":
                break
    final = problem.metrics(state)
    report = dict(
        status=status,
        converged=bool(status == "converged" and final_tau == 1.0),
        finalHomotopyTau=final_tau, elapsedSeconds=time.time() - start,
        iterations=len(rows) - len(taus), finalMetrics=final)
    dump(args.out / "report.json", report)
    q = state[:problem.N].reshape(problem.shape)
    np.savez_compressed(
        args.out / "state.npz", x=x, y=y, Q=q,
        psi=state[problem.N:2 * problem.N].reshape(problem.shape),
        parameters=state[2 * problem.N:],
        qx=(problem.Dx @ q.ravel()).reshape(problem.shape),
        residual=problem.evaluate(state)[0])
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
