"""Continue a solved self-similar IPM BVP to a larger box or finer grid.

The interior of the old solution is interpolated.  When the new box extends
beyond it, the seed is continued along rays from the shifted similarity
centre using the exact leading homogeneities

    Q + b/mu ~ r**alpha A(theta),
    psi - h*y ~ r**(1+alpha) B(theta).

The extrapolation is only an initial guess; the full transport--Poisson BVP
is then polished by Newton on the new grid.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.interpolate import RectBivariateSpline

from streamfunction_global_selfsimilar_solve import (
    NewtonOptions,
    dump,
    native_seed,
    solve,
)


ROOT = Path(__file__).resolve().parents[2]


def load_state(path: Path):
    with np.load(path) as data:
        return (
            np.asarray(data["x"], float),
            np.asarray(data["y"], float),
            np.asarray(data["Q"], float),
            np.asarray(data["psi"], float),
            np.asarray(data["parameters"], float),
        )


def enforce_core_gauges(problem, q: np.ndarray, psi: np.ndarray):
    """Remove interpolation defects without changing the far field."""
    envelope = np.exp(-(problem.X / 1.5)**8 - (problem.Y / 1.5)**4)
    basis = np.column_stack([
        (envelope * problem.X**degree).ravel() for degree in range(4)
    ])
    qflat = q.ravel().copy()
    qflat += basis @ np.linalg.solve(
        problem.C @ basis, problem.target - problem.C @ qflat)
    psiflat = psi.ravel().copy()
    center = float((problem.center_row @ psiflat).item())
    psiflat -= center * problem.yy
    return qflat, psiflat, center


def transferred_seed(problem, source_path: Path, blend_start: float,
                     blend_end: float, source_resonant_log: bool = False):
    sx, sy, sq, spsi, parameters = load_state(source_path)
    if parameters.size != 5:
        raise ValueError("source state must contain exactly five parameters")
    if source_resonant_log:
        lam, drift, offset, harmonic, log_coefficient = map(float, parameters)
        mu = 0.5 * lam
    else:
        lam, mu, drift, offset, harmonic = map(float, parameters)
        log_coefficient = 0.0
    if problem.resonant_log:
        target_parameters = np.array(
            [lam, drift, offset, harmonic, log_coefficient])
    elif not source_resonant_log:
        target_parameters = parameters
    else:
        raise ValueError("cannot transfer a resonant-log state to a nonresonant problem")
    if not (lam > mu > 0):
        raise ValueError("source rates must satisfy 0 < mu < lambda")
    alpha = mu / lam
    power = 1.0 + alpha
    shift = drift / lam
    source_radius = min(-sx[0] - abs(shift), sx[-1] - abs(shift), sy[-1])
    if source_radius <= 0:
        raise ValueError("source box does not surround the shifted centre")
    r0 = blend_start * source_radius
    r1 = blend_end * source_radius
    if not (0 < r0 < r1 < source_radius):
        raise ValueError("need 0 < blend-start < blend-end < 1")

    qspline = RectBivariateSpline(sy, sx, sq)
    psispline = RectBivariateSpline(sy, sx, spsi)
    z = problem.X + shift
    radius = np.hypot(z, problem.Y)
    angle = np.arctan2(problem.Y, z)

    # The circle r=r0 lies strictly inside the old rectangle.
    x0 = r0 * np.cos(angle) - shift
    y0 = r0 * np.sin(angle)
    q0 = qspline.ev(y0.ravel(), x0.ravel()).reshape(problem.shape)
    psi0 = psispline.ev(y0.ravel(), x0.ravel()).reshape(problem.shape)
    scale = np.maximum(radius / r0, 1e-300)
    q_outer = scale**alpha * (q0 + offset / mu) - offset / mu
    resonant0 = harmonic * y0 + log_coefficient * y0 * np.log(r0)
    resonant = (harmonic * problem.Y
                + log_coefficient * problem.Y
                * np.log(np.maximum(radius, 1e-300)))
    psi_outer = resonant + scale**power * (psi0 - resonant0)

    # Direct interpolation is used only where it is guaranteed to remain
    # inside the source box.  The smooth blend removes a derivative kink.
    q_direct = qspline.ev(problem.Y.ravel(), problem.X.ravel()).reshape(
        problem.shape)
    psi_direct = psispline.ev(
        problem.Y.ravel(), problem.X.ravel()).reshape(problem.shape)
    target_inside = bool(
        problem.x[0] >= sx[0] and problem.x[-1] <= sx[-1]
        and problem.y[-1] <= sy[-1])
    if target_inside:
        q, psi = q_direct, psi_direct
        transfer_mode = "interpolation"
    else:
        blend = np.clip((radius - r0) / (r1 - r0), 0.0, 1.0)
        blend = blend * blend * (3.0 - 2.0 * blend)
        q = (1.0 - blend) * q_direct + blend * q_outer
        psi = (1.0 - blend) * psi_direct + blend * psi_outer
        transfer_mode = "interpolation_plus_homogeneous_ray_extension"
    qflat, psiflat, center = enforce_core_gauges(problem, q, psi)
    state = np.r_[qflat, psiflat, target_parameters]
    transfer = dict(
        source=str(source_path), sourceLengthX=[float(sx[0]), float(sx[-1])],
        sourceLengthY=float(sy[-1]), shiftedCenter=float(-shift),
        usableSourceRadius=float(source_radius), blendRadius=[float(r0), float(r1)],
        transferMode=transfer_mode,
        sourceResonantLog=bool(source_resonant_log),
        targetResonantLog=bool(problem.resonant_log),
        interpolationCenterCorrection=center,
    )
    return state, transfer


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--length", type=float, required=True)
    parser.add_argument("--nx", type=int, default=129)
    parser.add_argument("--ny", type=int, default=65)
    parser.add_argument("--core-step", type=float, default=0.125)
    parser.add_argument("--stencil", type=int, default=5)
    parser.add_argument("--resonant-log", action="store_true")
    parser.add_argument("--source-resonant-log", action="store_true")
    parser.add_argument("--blend-start", type=float, default=0.65)
    parser.add_argument("--blend-end", type=float, default=0.85)
    parser.add_argument("--iterations", type=int, default=30)
    parser.add_argument("--tolerance", type=float, default=2e-8)
    parser.add_argument(
        "--taus", default="",
        help="optional residual-withdrawal homotopy, e.g. 0,.02,.05,.1,...,1")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError(f"refusing to overwrite {args.out}")
    args.out.mkdir(parents=True)

    problem, _, grid_metadata = native_seed(
        args.length, args.nx, args.ny, args.core_step, args.stencil,
        resonant_log=args.resonant_log)
    state, transfer = transferred_seed(
        problem, args.source, args.blend_start, args.blend_end,
        source_resonant_log=args.source_resonant_log)
    initial = problem.metrics(state)
    config = dict(
        grid=grid_metadata, transfer=transfer, initialMetrics=initial,
        stencil=args.stencil,
        scriptSha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    dump(args.out / "config.json", config)
    print(json.dumps(config, indent=2), flush=True)

    history_rows = []
    with (args.out / "iterations.jsonl").open("w") as stream:
        start = time.time()
        options = NewtonOptions(args.iterations, args.tolerance)
        initial_residual = problem.evaluate(state)[0]
        taus = ([float(value) for value in args.taus.split(",")]
                if args.taus else [1.0])
        if args.taus and (
                taus[0] != 0.0 or taus[-1] != 1.0
                or any(b <= a for a, b in zip(taus, taus[1:]))):
            raise ValueError("taus must increase strictly from 0 to 1")
        status = "not_started"
        final_tau = None
        for tau in taus:
            rhs = ((1.0 - tau) * initial_residual
                   if args.taus else np.zeros(problem.size))

            def callback(row):
                tagged = dict(homotopyTau=tau, **row)
                history_rows.append(tagged)
                stream.write(json.dumps(tagged) + "\n")
                stream.flush()
                print(
                    f"tau={tau:.5g} iter={row['iteration']:02d} "
                    f"res={row['residualInf']:.3e} "
                    f"rates=({row['lambda_']:.9f},{row['mu']:.9f},"
                    f"{row['v']:.9f}) alpha={row['alpha']:.9f} "
                    f"out={row['minimumOutwardVelocity']:.6g}", flush=True)

            state, status, history = solve(
                problem, state, options, right_hand_side=rhs,
                callback=callback)
            final_tau = tau
            np.savez_compressed(
                args.out / f"stage_tau_{tau:.5f}.npz",
                x=problem.x, y=problem.y,
                Q=state[:problem.N].reshape(problem.shape),
                psi=state[problem.N:2 * problem.N].reshape(problem.shape),
                parameters=state[2 * problem.N:], homotopyTau=tau)
            if status != "converged":
                break
        elapsed = time.time() - start
    final = problem.metrics(state)
    report = dict(
        status=status,
        converged=bool(status == "converged" and final_tau == 1.0),
        finalHomotopyTau=final_tau,
        elapsedSeconds=elapsed, iterations=len(history_rows) - len(taus),
        finalMetrics=final,
    )
    dump(args.out / "report.json", report)
    q = state[:problem.N].reshape(problem.shape)
    psi = state[problem.N:2 * problem.N].reshape(problem.shape)
    np.savez_compressed(
        args.out / "state.npz", x=problem.x, y=problem.y, Q=q, psi=psi,
        parameters=state[2 * problem.N:],
        qx=(problem.Dx @ q.ravel()).reshape(problem.shape),
        residual=problem.evaluate(state)[0])
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
