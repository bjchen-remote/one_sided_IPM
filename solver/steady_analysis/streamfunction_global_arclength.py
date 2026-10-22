"""Pseudo-arclength continuation for the global self-similar IPM BVP.

The continuation equation is

    F(Z) - (1-chi) F(Z_native) = 0,

with chi treated as an unknown.  A scaled secant phase condition crosses folds
where fixed-chi Newton loses invertibility.  Once the branch brackets chi=1,
the script interpolates the state and polishes the physical equation F(Z)=0.
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
from scipy import sparse
from scipy.sparse.linalg import splu

from streamfunction_global_selfsimilar_solve import (
    NewtonOptions,
    dump,
    native_seed,
    solve,
)


ROOT = Path(__file__).resolve().parents[2]


def load_stage(path: Path) -> tuple[np.ndarray, float]:
    with np.load(path) as data:
        state = np.r_[data["Q"].ravel(), data["psi"].ravel(),
                      data["parameters"]]
        tau = float(data["homotopyTau"])
    return state, tau


def continuation_scale(problem, native_state: np.ndarray) -> np.ndarray:
    (q, psi, lam, mu, drift, offset, harmonic,
     log_coefficient) = problem.unpack(native_state)
    alpha = mu / lam
    power = 1.0 + alpha
    radius = np.hypot(problem.xx + drift / lam, problem.yy)
    qscale = (1.0 + radius)**(-alpha)
    psiscale = (1.0 + radius)**(-power)
    parameter_scale = 1.0 / np.maximum(
        abs(native_state[2 * problem.N:]), 0.1)
    return np.r_[qscale, psiscale, parameter_scale]


def correct(problem, initial_state: np.ndarray, initial_tau: float,
            native_residual: np.ndarray, tangent: np.ndarray,
            state_scale: np.ndarray, predicted_state: np.ndarray,
            predicted_tau: float, tolerance: float, maximum_iterations: int):
    state = initial_state.copy()
    tau = float(initial_tau)
    tangent_state = tangent[:-1]
    tangent_tau = float(tangent[-1])
    phase_row = sparse.csr_matrix((tangent_state * state_scale)[None, :])
    history = []
    status = "maximum_iterations"
    for iteration in range(maximum_iterations + 1):
        physical, _, jacobian = problem.evaluate(state, jacobian=True)
        equation = physical - (1.0 - tau) * native_residual
        phase = float(
            tangent_state @ (state_scale * (state - predicted_state))
            + tangent_tau * (tau - predicted_tau))
        residual = np.r_[equation, phase]
        metrics = problem.metrics(state)
        physical_inf = metrics.pop("residualInf")
        physical_rms = metrics.pop("residualRms")
        row = dict(
            iteration=iteration, tau=tau,
            residualInf=float(np.max(abs(residual))),
            residualRms=float(np.linalg.norm(residual) / np.sqrt(residual.size)),
            physicalResidualInf=physical_inf,
            physicalResidualRms=physical_rms,
            phase=phase, **metrics)
        history.append(row)
        if row["residualInf"] < tolerance:
            status = "converged"
            break
        if iteration == maximum_iterations:
            break
        augmented = sparse.bmat([
            [jacobian, sparse.csc_matrix(native_residual[:, None])],
            [phase_row, sparse.csc_matrix([[tangent_tau]])],
        ], format="csc")
        try:
            step = splu(augmented, permc_spec="COLAMD").solve(-residual)
        except RuntimeError as exc:
            row["linearError"] = str(exc)
            status = "factorization_failed"
            break
        row["linearResidualRelative"] = float(
            np.linalg.norm(augmented @ step + residual)
            / max(np.linalg.norm(residual), 1e-300))
        length = 1.0
        old_norm = np.linalg.norm(residual)
        accepted = False
        while length >= 2.0**-16:
            trial_state = state + length * step[:-1]
            trial_tau = tau + length * step[-1]
            try:
                trial_physical, _ = problem.evaluate(trial_state)
            except FloatingPointError:
                length *= 0.5
                continue
            trial_equation = (
                trial_physical - (1.0 - trial_tau) * native_residual)
            trial_phase = float(
                tangent_state @ (state_scale * (trial_state - predicted_state))
                + tangent_tau * (trial_tau - predicted_tau))
            trial_residual = np.r_[trial_equation, trial_phase]
            if (np.isfinite(trial_residual).all()
                    and np.linalg.norm(trial_residual)
                    < (1.0 - 1e-4 * length) * old_norm):
                accepted = True
                break
            length *= 0.5
        row["stepLength"] = float(length)
        row["stepInf"] = float(np.max(abs(step)))
        row["stepTau"] = float(step[-1])
        if not accepted:
            status = "line_search_failed"
            break
        state, tau = trial_state, trial_tau
    return state, tau, status, history


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length", type=float, default=64.0)
    parser.add_argument("--nx", type=int, default=129)
    parser.add_argument("--ny", type=int, default=65)
    parser.add_argument("--core-step", type=float, default=0.125)
    parser.add_argument("--stencil", type=int, default=5)
    parser.add_argument(
        "--base", type=Path,
        help="optional tau=0 state defining the withdrawn residual")
    parser.add_argument("--density-robin", action="store_true")
    parser.add_argument("--resonant-log", action="store_true")
    parser.add_argument("--first", type=Path, required=True)
    parser.add_argument("--second", type=Path, required=True)
    parser.add_argument("--points", type=int, default=40)
    parser.add_argument("--step-factor", type=float, default=0.7)
    parser.add_argument("--tolerance", type=float, default=2e-8)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError(f"refusing to overwrite {args.out}")
    args.out.mkdir(parents=True)

    problem, native_state, metadata = native_seed(
        args.length, args.nx, args.ny, args.core_step, args.stencil,
        density_robin=args.density_robin,
        resonant_log=args.resonant_log)
    if args.base is not None:
        native_state, base_tau = load_stage(args.base)
        if abs(base_tau) > 1e-13:
            raise ValueError("base state must have homotopyTau=0")
        if native_state.size != problem.size:
            raise ValueError("base state grid does not match target problem")
    native_residual = problem.evaluate(native_state)[0]
    first_state, first_tau = load_stage(args.first)
    second_state, second_tau = load_stage(args.second)
    for label, state, tau in [
            ("first", first_state, first_tau),
            ("second", second_state, second_tau)]:
        defect = problem.evaluate(state)[0] - (1.0 - tau) * native_residual
        if np.max(abs(defect)) > 5e-7:
            raise ValueError(f"{label} start is not converged: {np.max(abs(defect)):g}")
    state_scale = continuation_scale(problem, native_state)
    config = dict(
        metadata=metadata, base=(str(args.base) if args.base else "native_seed"),
        densityRobin=bool(args.density_robin),
        resonantLog=bool(args.resonant_log),
        first=str(args.first), second=str(args.second),
        firstTau=first_tau, secondTau=second_tau,
        points=args.points, stepFactor=args.step_factor,
        tolerance=args.tolerance,
        scriptSha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        stateScaleRange=[float(state_scale.min()), float(state_scale.max())],
    )
    dump(args.out / "config.json", config)
    branch = [(first_state, first_tau), (second_state, second_tau)]
    initial_secant = np.r_[
        state_scale * (second_state - first_state), second_tau - first_tau]
    base_arclength = args.step_factor * np.linalg.norm(initial_secant)
    records = []
    start = time.time()
    bracket = None
    with (args.out / "arclength.jsonl").open("w") as stream:
        for index in range(args.points):
            z0, t0 = branch[-2]
            z1, t1 = branch[-1]
            secant = np.r_[state_scale * (z1 - z0), t1 - t0]
            norm = np.linalg.norm(secant)
            tangent = secant / norm
            trial_arclength = base_arclength
            solved = None
            while trial_arclength >= base_arclength / 32.0:
                factor = trial_arclength / norm
                predicted_state = z1 + factor * (z1 - z0)
                predicted_tau = t1 + factor * (t1 - t0)
                candidate = correct(
                    problem, predicted_state, predicted_tau, native_residual,
                    tangent, state_scale, predicted_state, predicted_tau,
                    args.tolerance, 14)
                if candidate[2] == "converged":
                    solved = candidate
                    break
                trial_arclength *= 0.5
            if solved is None:
                records.append(dict(index=index, status="step_failed"))
                break
            state, tau, status, history = solved
            metrics = problem.metrics(state)
            record = dict(
                index=index + 2, status=status, tau=float(tau),
                predictorFactor=float(factor),
                predictorArcLength=float(trial_arclength),
                correctorIterations=len(history) - 1,
                lambda_=metrics["lambda_"], mu=metrics["mu"],
                v=metrics["v"], alpha=metrics["alpha"],
                harmonic=metrics["harmonic"],
                homotopyResidualInf=history[-1]["residualInf"],
                physicalResidualInf=metrics["residualInf"],
            )
            records.append(record)
            stream.write(json.dumps(record) + "\n")
            stream.flush()
            print(json.dumps(record), flush=True)
            np.savez_compressed(
                args.out / f"point_{index + 2:03d}.npz",
                x=problem.x, y=problem.y,
                Q=state[:problem.N].reshape(problem.shape),
                psi=state[problem.N:2 * problem.N].reshape(problem.shape),
                parameters=state[2 * problem.N:], homotopyTau=tau)
            branch.append((state, tau))
            if (t1 - 1.0) * (tau - 1.0) <= 0 and max(t1, tau) >= 1.0:
                bracket = ((z1, t1), (state, tau))
                break

    polished = None
    if bracket is not None:
        (za, ta), (zb, tb) = bracket
        weight = (1.0 - ta) / (tb - ta)
        guess = za + weight * (zb - za)
        def callback(row):
            print(json.dumps(dict(polish=True, **row)), flush=True)
        state, status, history = solve(
            problem, guess, NewtonOptions(24, args.tolerance), callback=callback)
        polished = dict(status=status, metrics=problem.metrics(state),
                        iterations=len(history) - 1)
        np.savez_compressed(
            args.out / "state.npz", x=problem.x, y=problem.y,
            Q=state[:problem.N].reshape(problem.shape),
            psi=state[problem.N:2 * problem.N].reshape(problem.shape),
            parameters=state[2 * problem.N:])
    report = dict(
        elapsedSeconds=time.time() - start, branch=records,
        bracketedPhysicalTarget=bool(bracket is not None), polished=polished)
    dump(args.out / "report.json", report)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
