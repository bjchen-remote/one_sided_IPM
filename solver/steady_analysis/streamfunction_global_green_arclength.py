"""Pseudo-arclength continuation for the centered Green-closed IPM BVP.

The homotopy equation is

    F(Z) - (1-tau) F(Z_native) = 0.

Treating ``tau`` as an unknown and adding a scaled secant phase condition
allows the branch to pass folds at which a fixed-tau Newton solve fails.  If
the branch crosses tau=1, the interpolated state is polished against the
physical equation F(Z)=0.
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
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import LinearOperator, gmres, splu

from streamfunction_global_green_solve import native_seed, solve
from streamfunction_global_selfsimilar_solve import NewtonOptions, dump


def load_stage(path: Path) -> tuple[np.ndarray, float]:
    with np.load(path) as data:
        state = np.r_[data["Q"].ravel(), data["psi"].ravel(),
                      data["parameters"]]
        tau = float(data["homotopyTau"])
    return state, tau


def continuation_scale(problem, native_state: np.ndarray) -> np.ndarray:
    _, _, lam, mu, drift, _ = problem.unpack(native_state)
    alpha = mu / lam
    radius = np.hypot(problem.xx + drift / lam, problem.yy)
    qscale = (1.0 + radius)**(-alpha)
    psiscale = (1.0 + radius)**(-(1.0 + alpha))
    parameter_scale = 1.0 / np.maximum(
        np.abs(native_state[2 * problem.N:]), 0.1)
    return np.r_[qscale, psiscale, parameter_scale]


def bordered_linear_solve(problem, state, jacobian,
                          homotopy_column, phase_state, phase_tau,
                          equation, phase):
    """Direct Woodbury solve of the Green-updated augmented system."""
    size = state.size
    local = problem.local_jacobian(state)
    factorization = splu(local, permc_spec="COLAMD")

    def solve_local(right_hand_side):
        solution = factorization.solve(right_hand_side)
        for _ in range(2):
            defect = local @ solution - right_hand_side
            solution += factorization.solve(-defect)
        return solution

    local_homotopy = solve_local(homotopy_column)
    denominator = float(phase_tau - phase_state @ local_homotopy)
    if abs(denominator) < 1e-13:
        raise RuntimeError("singular local bordered preconditioner")

    # Full augmented matrix = local bordered matrix + U V.  U selects the
    # artificial-boundary streamfunction equations and V contains their
    # weighted centered-Green dependence on Q.
    boundary = problem.boundary_indices
    rank = len(boundary)
    local_rhs = np.zeros((size, rank + 1), order="F")
    rhs = -np.r_[equation, phase]
    local_rhs[:, 0] = rhs[:size]
    local_rhs[problem.N + boundary, 1 + np.arange(rank)] = 1.0
    direct = solve_local(local_rhs)
    phase_rhs = np.r_[rhs[size], np.zeros(rank)]
    tau_parts = (
        phase_rhs - phase_state @ direct) / denominator
    state_parts = direct - local_homotopy[:, None] * tau_parts[None, :]
    local_bordered_solutions = np.vstack([state_parts, tau_parts])
    local_solution = local_bordered_solutions[:, 0]
    woodbury_basis = local_bordered_solutions[:, 1:]

    green_rows = problem.green_q[boundary, :].toarray()
    green_update_q = (
        -problem.psi_weight[boundary, None] * green_rows)
    projected_solution = green_update_q @ local_solution[:problem.N]
    projected_basis = green_update_q @ woodbury_basis[:problem.N, :]
    woodbury_matrix = np.eye(rank) + projected_basis
    try:
        woodbury_factor = lu_factor(woodbury_matrix)
        coefficients = lu_solve(woodbury_factor, projected_solution)
    except (ValueError, np.linalg.LinAlgError) as exc:
        raise RuntimeError(f"singular Green Woodbury matrix: {exc}") from exc
    step = local_solution - woodbury_basis @ coefficients

    def apply_local_bordered(vector):
        direct_vector = solve_local(vector[:size])
        tau_vector = float(
            (vector[size] - phase_state @ direct_vector) / denominator)
        return np.r_[
            direct_vector - local_homotopy * tau_vector, tau_vector]

    def apply_full_inverse(vector):
        local_vector = apply_local_bordered(vector)
        projected = green_update_q @ local_vector[:problem.N]
        coefficient = lu_solve(woodbury_factor, projected)
        return local_vector - woodbury_basis @ coefficient

    def augmented_matvec(vector):
        state_vector, tau_vector = vector[:size], vector[size]
        return np.r_[
            jacobian @ state_vector + homotopy_column * tau_vector,
            phase_state @ state_vector + phase_tau * tau_vector,
        ]

    augmented = LinearOperator(
        (size + 1, size + 1), matvec=augmented_matvec)
    woodbury_preconditioner = LinearOperator(
        (size + 1, size + 1), matvec=apply_full_inverse)
    krylov_history = []
    step, info = gmres(
        augmented, rhs, x0=step, M=woodbury_preconditioner,
        rtol=2e-11, atol=0.0, restart=20, maxiter=10,
        callback=krylov_history.append, callback_type="pr_norm")
    defect = augmented @ step - rhs
    relative_residual = float(
        np.linalg.norm(defect) / max(np.linalg.norm(rhs), 1e-300))
    if info != 0 and relative_residual > 2e-8:
        raise RuntimeError(
            f"Woodbury-GMRES failed: info={info}, "
            f"relative residual={relative_residual:g}")
    diagnostics = dict(
        woodburyRank=rank,
        woodburyCondition=float(np.linalg.cond(woodbury_matrix)),
        krylovIterations=len(krylov_history),
        relativeResidual=relative_residual)
    return step, diagnostics


def correct(problem, initial_state: np.ndarray, initial_tau: float,
            native_residual: np.ndarray, tangent: np.ndarray,
            state_scale: np.ndarray, predicted_state: np.ndarray,
            predicted_tau: float, tolerance: float,
            maximum_iterations: int):
    state = initial_state.copy()
    tau = float(initial_tau)
    tangent_state = tangent[:-1]
    tangent_tau = float(tangent[-1])
    phase_state = tangent_state * state_scale
    history = []
    status = "maximum_iterations"
    for iteration in range(maximum_iterations + 1):
        try:
            physical, _, jacobian = problem.evaluate(state, jacobian=True)
        except FloatingPointError:
            return state, tau, "invalid_predictor", history
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
            residualInf=float(np.max(np.abs(residual))),
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
        try:
            step, linear_diagnostics = bordered_linear_solve(
                problem, state, jacobian, native_residual,
                phase_state, tangent_tau, equation, phase)
        except RuntimeError as exc:
            row["linearError"] = str(exc)
            status = "factorization_failed"
            break
        row["linearResidualRelative"] = linear_diagnostics["relativeResidual"]
        row["woodburyRank"] = linear_diagnostics["woodburyRank"]
        row["woodburyCondition"] = linear_diagnostics["woodburyCondition"]
        row["linearKrylovIterations"] = linear_diagnostics["krylovIterations"]
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
                tangent_state @ (
                    state_scale * (trial_state - predicted_state))
                + tangent_tau * (trial_tau - predicted_tau))
            trial_residual = np.r_[trial_equation, trial_phase]
            if (np.isfinite(trial_residual).all()
                    and np.linalg.norm(trial_residual)
                    < (1.0 - 1e-4 * length) * old_norm):
                accepted = True
                break
            length *= 0.5
        row["stepLength"] = float(length)
        row["stepInf"] = float(np.max(np.abs(step)))
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
    parser.add_argument("--stencil", type=int, default=7)
    parser.add_argument("--quadrature-order", type=int, default=2)
    parser.add_argument("--base", type=Path)
    parser.add_argument("--first", type=Path, required=True)
    parser.add_argument("--second", type=Path, required=True)
    parser.add_argument("--points", type=int, default=60)
    parser.add_argument("--step-factor", type=float, default=0.7)
    parser.add_argument("--tolerance", type=float, default=2e-8)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError(f"refusing to overwrite {args.out}")
    args.out.mkdir(parents=True)

    problem, native_state, metadata = native_seed(
        args.length, args.nx, args.ny, args.core_step,
        args.stencil, args.quadrature_order)
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
        if np.max(np.abs(defect)) > 5e-7:
            raise ValueError(
                f"{label} start is not converged: "
                f"{np.max(np.abs(defect)):g}")

    state_scale = continuation_scale(problem, native_state)
    initial_secant = np.r_[
        state_scale * (second_state - first_state), second_tau - first_tau]
    base_arclength = args.step_factor * np.linalg.norm(initial_secant)
    config = dict(
        metadata=metadata,
        base=(str(args.base) if args.base else "native_seed"),
        first=str(args.first), second=str(args.second),
        firstTau=first_tau, secondTau=second_tau,
        points=args.points, stepFactor=args.step_factor,
        baseArcLength=float(base_arclength), tolerance=args.tolerance,
        stencil=args.stencil, quadratureOrder=args.quadrature_order,
        stateScaleRange=[float(state_scale.min()), float(state_scale.max())],
        scriptSha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    dump(args.out / "config.json", config)

    branch = [(first_state, first_tau), (second_state, second_tau)]
    records = []
    start = time.time()
    bracket = None
    with (args.out / "arclength.jsonl").open("w") as stream:
        for index in range(args.points):
            z0, t0 = branch[-2]
            z1, t1 = branch[-1]
            secant = np.r_[state_scale * (z1 - z0), t1 - t0]
            secant_norm = np.linalg.norm(secant)
            tangent = secant / secant_norm
            trial_arclength = base_arclength
            solved = None
            attempts = []
            while trial_arclength >= base_arclength / 64.0:
                factor = trial_arclength / secant_norm
                predicted_state = z1 + factor * (z1 - z0)
                predicted_tau = t1 + factor * (t1 - t0)
                candidate = correct(
                    problem, predicted_state, predicted_tau, native_residual,
                    tangent, state_scale, predicted_state, predicted_tau,
                    args.tolerance, 16)
                attempts.append(dict(
                    arcLength=float(trial_arclength), status=candidate[2],
                    iterations=max(len(candidate[3]) - 1, 0),
                    linearError=(candidate[3][-1].get("linearError")
                                 if candidate[3] else None)))
                if candidate[2] == "converged":
                    solved = candidate
                    break
                trial_arclength *= 0.5
            if solved is None:
                record = dict(index=index + 2, status="step_failed",
                              attempts=attempts)
                records.append(record)
                stream.write(json.dumps(record) + "\n")
                stream.flush()
                print(json.dumps(record), flush=True)
                break
            state, tau, status, history = solved
            metrics = problem.metrics(state)
            record = dict(
                index=index + 2, status=status, tau=float(tau),
                predictorFactor=float(factor),
                predictorArcLength=float(trial_arclength),
                correctorIterations=len(history) - 1,
                lambda_=metrics["lambda_"], mu=metrics["mu"],
                v=metrics["v"], b=metrics["b"], alpha=metrics["alpha"],
                centerDerivative=metrics["centerDerivative"],
                minimumOutwardVelocity=metrics["minimumOutwardVelocity"],
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
            if (t1 - 1.0) * (tau - 1.0) <= 0.0 and max(t1, tau) >= 1.0:
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
            problem, guess, NewtonOptions(30, args.tolerance),
            callback=callback)
        polished = dict(
            status=status, metrics=problem.metrics(state),
            iterations=len(history) - 1)
        np.savez_compressed(
            args.out / "state.npz", x=problem.x, y=problem.y,
            Q=state[:problem.N].reshape(problem.shape),
            psi=state[problem.N:2 * problem.N].reshape(problem.shape),
            parameters=state[2 * problem.N:])
    report = dict(
        elapsedSeconds=time.time() - start, branch=records,
        bracketedPhysicalTarget=bool(bracket is not None),
        polished=polished)
    dump(args.out / "report.json", report)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
