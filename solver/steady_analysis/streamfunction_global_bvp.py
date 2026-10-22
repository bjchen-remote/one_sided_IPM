"""Solve the centered steady IPM profile with self-similar outer conditions.

This is a coupled Q/psi boundary-value problem.  Unlike the earlier box
solver, no exterior velocity is frozen: Poisson is solved together with the
transport equation, and the artificial boundary is closed by the leading
self-similar radial relations.  A homotopy starts at a saved finite-box root
and replaces its outer Dirichlet/Cauchy treatment by the asymptotic boundary
conditions.

Unknowns are Q, an asymptotically gauged psi, and (lambda, mu, v, b).  The
velocity entering transport is explicitly centered by subtracting u_1(0,0).
The four Q gauges make the system square, and the origin transport equation
then enforces b=v.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path
import time

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from fixed_point import fd_matrix


ROOT = Path(__file__).resolve().parents[2]


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


@dataclass
class NewtonOptions:
    max_iterations: int = 16
    tolerance: float = 2e-9
    minimum_step: float = 2.0**-14


class GlobalProfileBVP:
    """Finite-box discretization of the one-patch global BVP."""

    def __init__(self, x: np.ndarray, y: np.ndarray, psi_reference: np.ndarray,
                 stencil: int = 5):
        self.x = np.asarray(x, float)
        self.y = np.asarray(y, float)
        self.nx = len(self.x)
        self.ny = len(self.y)
        self.shape = (self.ny, self.nx)
        self.N = self.nx * self.ny
        self.X, self.Y = np.meshgrid(self.x, self.y)
        self.xx = self.X.ravel()
        self.yy = self.Y.ravel()
        if abs(self.y[0]) > 1e-13 or min(abs(self.x)) > 1e-13:
            raise ValueError("grid must contain y=0 and x=0")
        self.origin = int(np.argmin(abs(self.x)))
        ix = sparse.eye(self.nx, format="csr")
        iy = sparse.eye(self.ny, format="csr")
        self.Dx = sparse.kron(iy, fd_matrix(self.x, 1, stencil), format="csr")
        self.Dy = sparse.kron(fd_matrix(self.y, 1, stencil), ix, format="csr")
        self.Lap = (
            sparse.kron(iy, fd_matrix(self.x, 2, stencil), format="csr")
            + sparse.kron(fd_matrix(self.y, 2, stencil), ix, format="csr")
        )
        wall = np.zeros(self.shape, bool)
        wall[0, :] = True
        artificial = np.zeros(self.shape, bool)
        artificial[-1, :] = True
        artificial[:, [0, -1]] = True
        self.wall = wall.ravel()
        self.artificial_q = artificial.ravel()
        self.artificial_psi = (artificial & ~wall).ravel()
        self.poisson = ~(self.wall | self.artificial_psi)
        self.transport = ~self.artificial_q
        self.I = sparse.eye(self.N, format="csr")
        self.Stransport = sparse.diags(self.transport.astype(float), format="csr")
        self.Stailq = sparse.diags(self.artificial_q.astype(float), format="csr")
        self.Spoisson = sparse.diags(self.poisson.astype(float), format="csr")
        self.Swall = sparse.diags(self.wall.astype(float), format="csr")
        self.Stailpsi = sparse.diags(self.artificial_psi.astype(float), format="csr")

        # The same four centered normalization conditions used by fixed_point.py.
        x2 = fd_matrix(self.x, 2, max(stencil, 7)).getrow(self.origin)
        x3 = fd_matrix(self.x, 3, max(stencil, 7)).getrow(self.origin)
        wall_pick = sparse.csr_matrix(
            ([1.0], ([0], [0])), shape=(1, self.ny))
        self.C = sparse.vstack([
            sparse.csr_matrix(
                ([1.0], ([0], [self.origin])), shape=(1, self.N)),
            self.Dx.getrow(self.origin),
            sparse.kron(wall_pick, x2, format="csr"),
            sparse.kron(wall_pick, x3, format="csr"),
        ], format="csr")
        self.target = np.array([0.0, 1.0, 0.0, -2.0])
        cnorm = np.sqrt(np.asarray(self.C.multiply(self.C).sum(axis=1))).ravel()
        self.constraint_scale = 1.0 / np.maximum(cnorm, 1.0)
        self.center_row = self.Dy.getrow(self.origin)
        self.center_velocity_operator = (
            -self.Dy
            + sparse.csr_matrix(np.ones((self.N, 1))) @ self.center_row
        )
        self.psi_reference = np.asarray(psi_reference, float).reshape(-1)

    @property
    def size(self) -> int:
        return 2 * self.N + 4

    def unpack(self, state: np.ndarray):
        q = state[:self.N]
        psi = state[self.N:2 * self.N]
        lam, mu, drift, offset = state[2 * self.N:]
        return q, psi, lam, mu, drift, offset

    def evaluate(self, state: np.ndarray, tau: float, jacobian: bool = False):
        q, psi, lam, mu, drift, offset = self.unpack(state)
        if not (lam > 0 and mu > 0):
            raise FloatingPointError("rates left the positive branch")
        alpha = mu / lam
        power = 1.0 + alpha
        zcoord = self.xx + drift / lam
        shifted_q = q + offset / mu
        qx = self.Dx @ q
        qy = self.Dy @ q
        psix = self.Dx @ psi
        psiy = self.Dy @ psi
        raw_ux = -psiy
        center_ux = raw_ux[self.origin]
        ux = raw_ux - center_ux
        uy = psix
        vx = ux + lam * self.xx + drift
        vy = uy + lam * self.yy
        transport_residual = vx * qx + vy * qy - mu * q - offset
        q_tail = zcoord * qx + self.yy * qy - alpha * shifted_q
        q_residual = transport_residual.copy()
        q_residual[self.artificial_q] = (
            (1.0 - tau) * transport_residual[self.artificial_q]
            + tau * q_tail[self.artificial_q]
        )

        poisson_residual = -self.Lap @ psi - qx
        psi_tail = zcoord * psix + self.yy * psiy - power * psi
        psi_residual = poisson_residual.copy()
        psi_residual[self.wall] = psi[self.wall]
        psi_residual[self.artificial_psi] = (
            (1.0 - tau)
            * (psi[self.artificial_psi]
               - self.psi_reference[self.artificial_psi])
            + tau * psi_tail[self.artificial_psi]
        )
        gauges = (self.C @ q - self.target) * self.constraint_scale
        residual = np.r_[q_residual, psi_residual, gauges]

        cache = dict(
            q=q, psi=psi, qx=qx, qy=qy, psix=psix, psiy=psiy,
            vx=vx, vy=vy, alpha=alpha, power=power, zcoord=zcoord,
            shifted_q=shifted_q, transport=transport_residual,
            q_tail=q_tail, poisson=poisson_residual, psi_tail=psi_tail,
            lam=lam, mu=mu, drift=drift, offset=offset,
            rawCenterVelocity=center_ux, gauges=gauges,
        )
        if not jacobian:
            return residual, cache
        return residual, cache, self.jacobian(cache, tau)

    def jacobian(self, c: dict, tau: float) -> sparse.csc_matrix:
        lam = c["lam"]
        mu = c["mu"]
        drift = c["drift"]
        alpha = c["alpha"]
        q = c["q"]
        shifted_q = c["shifted_q"]
        qx, qy = c["qx"], c["qy"]
        psi, psix = c["psi"], c["psix"]
        zcoord = c["zcoord"]

        transport_qq = (
            sparse.diags(c["vx"]) @ self.Dx
            + sparse.diags(c["vy"]) @ self.Dy
            - mu * self.I
        )
        transport_qpsi = (
            sparse.diags(qx) @ self.center_velocity_operator
            + sparse.diags(qy) @ self.Dx
        )
        tail_qq = (
            sparse.diags(zcoord) @ self.Dx
            + sparse.diags(self.yy) @ self.Dy
            - alpha * self.I
        )
        blend_q = self.transport.astype(float) + (
            self.artificial_q.astype(float) * (1.0 - tau))
        Jqq = sparse.diags(blend_q) @ transport_qq + tau * self.Stailq @ tail_qq
        Jqpsi = sparse.diags(blend_q) @ transport_qpsi

        q_params = np.column_stack([
            self.xx * qx + self.yy * qy,
            -q,
            qx,
            -np.ones(self.N),
        ])
        outer = self.artificial_q
        tail_params = np.column_stack([
            -(drift / lam**2) * qx + (alpha / lam) * shifted_q,
            -q / lam,
            qx / lam,
            -np.ones(self.N) / lam,
        ])
        q_params[outer] = (
            (1.0 - tau) * q_params[outer] + tau * tail_params[outer]
        )

        Jpsi_q = -self.Spoisson @ self.Dx
        tail_psipsi = (
            sparse.diags(zcoord) @ self.Dx
            + sparse.diags(self.yy) @ self.Dy
            - c["power"] * self.I
        )
        Jpsi_psi = (
            -self.Spoisson @ self.Lap + self.Swall
            + (1.0 - tau) * self.Stailpsi
            + tau * self.Stailpsi @ tail_psipsi
        )
        psi_params = np.zeros((self.N, 4))
        outer_psi = self.artificial_psi
        tail_psi_params = np.column_stack([
            -(drift / lam**2) * psix
            + (alpha / lam) * psi,
            -psi / lam,
            psix / lam,
            np.zeros(self.N),
        ])
        psi_params[outer_psi] = tau * tail_psi_params[outer_psi]

        zeros_4n = sparse.csr_matrix((4, self.N))
        zeros_44 = sparse.csr_matrix((4, 4))
        return sparse.bmat([
            [Jqq, Jqpsi, sparse.csr_matrix(q_params)],
            [Jpsi_q, Jpsi_psi, sparse.csr_matrix(psi_params)],
            [sparse.diags(self.constraint_scale) @ self.C,
             zeros_4n, zeros_44],
        ], format="csc")

    def metrics(self, state: np.ndarray, tau: float) -> dict:
        residual, c = self.evaluate(state, tau)
        q_active = ~self.artificial_q
        psi_active = self.poisson
        radial = np.hypot(c["zcoord"], self.yy)
        outer_radius = radial[self.artificial_q]
        core = (abs(self.xx) <= 2.0) & (self.yy <= 2.0)
        return dict(
            tau=float(tau),
            residualInf=float(np.max(abs(residual))),
            residualRms=float(np.linalg.norm(residual) / np.sqrt(residual.size)),
            transportInf=float(np.max(abs(c["transport"][q_active]))),
            coreTransportInf=float(np.max(abs(c["transport"][core]))),
            poissonInf=float(np.max(abs(c["poisson"][psi_active]))),
            densityTailInf=float(np.max(abs(c["q_tail"][self.artificial_q]))),
            streamfunctionTailInf=float(
                np.max(abs(c["psi_tail"][self.artificial_psi]))),
            gaugeInf=float(np.max(abs(c["gauges"]))),
            rawCenterVelocity=float(c["rawCenterVelocity"]),
            lambda_=float(c["lam"]),
            mu=float(c["mu"]),
            v=float(c["drift"]),
            b=float(c["offset"]),
            alpha=float(c["alpha"]),
            offsetDriftDefect=float(c["offset"] - c["drift"]),
            minimumOuterRadius=float(np.min(outer_radius)),
            maxAbsQ=float(np.max(abs(c["q"]))),
            maxAbsPsi=float(np.max(abs(c["psi"]))),
        )


def load_seed(state_path: Path, stencil: int):
    with np.load(state_path) as data:
        x = data["x"].copy()
        y = data["y"].copy()
        q = data["Q"].copy()
        psi = data["psi"].copy()
        rates = data["rates"].copy()
    ix = sparse.eye(len(x), format="csr")
    dy = sparse.kron(fd_matrix(y, 1, stencil), ix, format="csr")
    origin = int(np.argmin(abs(x)))
    horizontal_velocity = -(dy.getrow(origin) @ psi.ravel()).item()
    initial = np.r_[q.ravel(), psi.ravel(), rates]
    return x, y, psi, initial, horizontal_velocity


def directional_check(problem: GlobalProfileBVP, state: np.ndarray,
                      tau: float) -> float:
    residual, _, matrix = problem.evaluate(state, tau, jacobian=True)
    rng = np.random.default_rng(290926)
    direction = rng.standard_normal(state.size)
    direction /= np.linalg.norm(direction)
    epsilon = 2e-6
    plus = problem.evaluate(state + epsilon * direction, tau)[0]
    minus = problem.evaluate(state - epsilon * direction, tau)[0]
    numerical = (plus - minus) / (2 * epsilon)
    exact = matrix @ direction
    return float(np.linalg.norm(numerical - exact)
                 / max(np.linalg.norm(exact), 1e-300))


def newton(problem: GlobalProfileBVP, initial: np.ndarray, tau: float,
           options: NewtonOptions, callback=None):
    state = np.asarray(initial, float).copy()
    history = []
    status = "maximum_iterations"
    for iteration in range(options.max_iterations + 1):
        residual, cache, matrix = problem.evaluate(state, tau, jacobian=True)
        metrics = problem.metrics(state, tau)
        entry = dict(iteration=iteration, **metrics)
        history.append(entry)
        if callback:
            callback(entry)
        if metrics["residualInf"] < options.tolerance:
            status = "converged"
            break
        if iteration == options.max_iterations:
            break
        try:
            factor = splu(matrix, permc_spec="COLAMD")
            step = factor.solve(-residual)
        except RuntimeError as exc:
            status = "factorization_failed"
            entry["linearError"] = str(exc)
            break
        entry["linearResidualRelative"] = float(
            np.linalg.norm(matrix @ step + residual)
            / max(np.linalg.norm(residual), 1e-300))
        length = 1.0
        current_norm = np.linalg.norm(residual)
        current_inf = np.max(abs(residual))
        accepted = False
        while length >= options.minimum_step:
            trial = state + length * step
            _, _, trial_lam, trial_mu, _, _ = problem.unpack(trial)
            if not (trial_lam > trial_mu > 0):
                length *= 0.5
                continue
            try:
                trial_residual, _ = problem.evaluate(trial, tau)
            except FloatingPointError:
                length *= 0.5
                continue
            if (np.isfinite(trial_residual).all()
                    and np.linalg.norm(trial_residual)
                    < (1.0 - 1e-4 * length) * current_norm
                    and np.max(abs(trial_residual))
                    <= max(1.05 * current_inf, 2 * options.tolerance)):
                accepted = True
                break
            length *= 0.5
        entry["stepLength"] = float(length)
        entry["stepInf"] = float(np.max(abs(step)))
        if not accepted:
            status = "line_search_failed"
            break
        state = trial
    return state, status, history


def parse_taus(text: str) -> list[float]:
    values = [float(item) for item in text.split(",") if item.strip()]
    if not values or values[0] != 0.0 or values[-1] != 1.0:
        raise ValueError("homotopy schedule must start at 0 and end at 1")
    if any(b <= a for a, b in zip(values, values[1:])):
        raise ValueError("homotopy schedule must be strictly increasing")
    return values


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--seed", type=Path,
        default=ROOT / "runs/steady_20260928_tau20/L4_129x65/state.npz")
    parser.add_argument(
        "--out", type=Path,
        default=ROOT / "runs/streamfunction_global_bvp_20260929/L4_129x65_v3")
    parser.add_argument("--stencil", type=int, default=5)
    parser.add_argument(
        "--taus", default="0,0.05,0.1,0.2,0.35,0.5,0.7,0.85,1")
    parser.add_argument("--iterations", type=int, default=16)
    parser.add_argument("--tolerance", type=float, default=2e-9)
    parser.add_argument("--audit-only", action="store_true")
    args = parser.parse_args()

    x, y, psi_reference, state, center_shift = load_seed(
        args.seed, args.stencil)
    problem = GlobalProfileBVP(x, y, psi_reference, args.stencil)
    check = directional_check(problem, state, 0.37)
    if check > 2e-6:
        raise AssertionError(f"Jacobian directional error {check:g}")
    initial_metrics = problem.metrics(state, 0.0)
    initial_target_metrics = problem.metrics(state, 1.0)
    print(json.dumps(dict(
        seed=str(args.seed), shape=problem.shape,
        jacobianDirectionalRelativeError=check,
        centeredVelocityShift=center_shift,
        initialHomotopyMetrics=initial_metrics,
        initialTargetMetrics=initial_target_metrics,
    ), indent=2), flush=True)
    if args.audit_only:
        return
    if args.out.exists():
        raise ValueError(f"refusing to overwrite {args.out}")
    args.out.mkdir(parents=True)
    config = dict(
        seed=str(args.seed.resolve()),
        seedSha256=hashlib.sha256(args.seed.read_bytes()).hexdigest(),
        scriptSha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        shape=list(problem.shape), stencil=args.stencil,
        homotopy=parse_taus(args.taus),
        jacobianDirectionalRelativeError=check,
        centeredVelocityShift=center_shift,
        initialHomotopyMetrics=initial_metrics,
        initialTargetMetrics=initial_target_metrics,
        formulation=(
            "Coupled Q/psi solve; transport and Poisson in the box; "
            "leading shifted self-similar Robin conditions on all artificial "
            "boundaries; the Poisson gauge excludes the linear y mode and "
            "transport explicitly subtracts u_1(0,0)."),
    )
    write_json(args.out / "config.json", config)
    all_history = []
    start = time.time()
    status = "not_started"
    options = NewtonOptions(args.iterations, args.tolerance)
    with (args.out / "iterations.jsonl").open("w") as stream:
        for tau in config["homotopy"]:
            def callback(row):
                tagged = dict(stageTau=tau, **row)
                all_history.append(tagged)
                stream.write(json.dumps(tagged) + "\n")
                stream.flush()
                print(
                    f"tau={tau:.3f} iter={row['iteration']:02d} "
                    f"res={row['residualInf']:.3e} "
                    f"rates=({row['lambda_']:.7f},{row['mu']:.7f},"
                    f"{row['v']:.7f},{row['b']:.7f}) "
                    f"alpha={row['alpha']:.7f}", flush=True)
            state, status, _ = newton(problem, state, tau, options, callback)
            np.savez_compressed(
                args.out / f"stage_tau_{tau:.3f}.npz",
                x=x, y=y,
                Q=state[:problem.N].reshape(problem.shape),
                psi=state[problem.N:2 * problem.N].reshape(problem.shape),
                parameters=state[2 * problem.N:], tau=tau)
            if status != "converged":
                break
    final_tau = all_history[-1]["stageTau"]
    final_metrics = problem.metrics(state, final_tau)
    result = dict(
        status=status,
        completedTarget=bool(status == "converged" and final_tau == 1.0),
        finalTau=final_tau,
        seconds=time.time() - start,
        finalMetrics=final_metrics,
        stages=[row for row in all_history if row["iteration"] == 0],
    )
    write_json(args.out / "report.json", result)
    np.savez_compressed(
        args.out / "state.npz", x=x, y=y,
        Q=state[:problem.N].reshape(problem.shape),
        psi=state[problem.N:2 * problem.N].reshape(problem.shape),
        parameters=state[2 * problem.N:], tau=final_tau)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
