"""Global self-similar IPM solve with the centered half-plane Green closure.

The density transport equation is imposed at every node.  The streamfunction
solves Poisson in the interior, vanishes on the physical wall, and on the
artificial faces is supplied by the centered half-plane Green integral

  psi(X) = int [G_H(X,S) - y*eta/(pi*(xi^2+eta^2))] Q_x(S) dS.

The subtraction fixes the velocity at the wall core and removes the otherwise
free linear harmonic mode.  Thus no fitted far-field Robin or harmonic
coefficient is present.  The continuum closure is exact; its numerical
approximations are finite-box truncation, quadrature, and finite differences.
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
from scipy.interpolate import RectBivariateSpline
from scipy.sparse.linalg import LinearOperator, gmres, splu

from fixed_point import fd_matrix
from streamfunction_sources import load_native
from streamfunction_global_selfsimilar_solve import NewtonOptions, dump, stretched_grid


ROOT = Path(__file__).resolve().parents[2]


class CenteredGreenBVP:
    def __init__(self, x: np.ndarray, y: np.ndarray,
                 reference_rates: np.ndarray, stencil: int = 5,
                 quadrature_order: int = 2):
        self.x, self.y = np.asarray(x, float), np.asarray(y, float)
        self.nx, self.ny = len(self.x), len(self.y)
        self.shape = self.ny, self.nx
        self.N = self.nx * self.ny
        self.X, self.Y = np.meshgrid(self.x, self.y)
        self.xx, self.yy = self.X.ravel(), self.Y.ravel()
        self.origin = int(np.argmin(abs(self.x)))
        ix = sparse.eye(self.nx, format="csr")
        iy = sparse.eye(self.ny, format="csr")
        self.Dx = sparse.kron(
            iy, fd_matrix(self.x, 1, stencil), format="csr")
        self.Dy = sparse.kron(
            fd_matrix(self.y, 1, stencil), ix, format="csr")
        self.Lap = (
            sparse.kron(iy, fd_matrix(self.x, 2, stencil), format="csr")
            + sparse.kron(fd_matrix(self.y, 2, stencil), ix, format="csr"))
        self.I = sparse.eye(self.N, format="csr")

        wall = np.zeros(self.shape, bool)
        wall[0, :] = True
        artificial = np.zeros(self.shape, bool)
        artificial[-1, :] = True
        artificial[:, [0, -1]] = True
        self.wall = wall.ravel()
        self.artificial = (artificial & ~wall).ravel()
        self.poisson = ~(self.wall | self.artificial)
        self.boundary_indices = np.flatnonzero(self.artificial)
        self.Spoisson = sparse.diags(self.poisson.astype(float))
        self.Swall = sparse.diags(self.wall.astype(float))
        self.Sartificial = sparse.diags(self.artificial.astype(float))

        wall_pick = sparse.csr_matrix(
            ([1.0], ([0], [0])), shape=(1, self.ny))
        x2 = fd_matrix(self.x, 2, max(7, stencil)).getrow(self.origin)
        x3 = fd_matrix(self.x, 3, max(7, stencil)).getrow(self.origin)
        self.C = sparse.vstack([
            sparse.csr_matrix(
                ([1.0], ([0], [self.origin])), shape=(1, self.N)),
            self.Dx.getrow(self.origin),
            sparse.kron(wall_pick, x2, format="csr"),
            sparse.kron(wall_pick, x3, format="csr"),
        ], format="csr")
        self.target = np.array([0.0, 1.0, 0.0, -2.0])
        cnorm = np.sqrt(np.asarray(
            self.C.multiply(self.C).sum(axis=1))).ravel()
        self.gauge_scale = 1.0 / np.maximum(cnorm, 1.0)

        lam0, mu0, drift0, _ = map(float, reference_rates)
        alpha0 = mu0 / lam0
        power0 = 1.0 + alpha0
        radius = np.hypot(self.xx + drift0 / lam0, self.yy)
        self.q_weight = (1.0 + radius)**(-alpha0)
        self.psi_weight = np.empty(self.N)
        self.psi_weight[self.poisson] = (
            1.0 + radius[self.poisson])**(1.0 - alpha0)
        self.psi_weight[self.wall | self.artificial] = (
            1.0 + radius[self.wall | self.artificial])**(-power0)

        self.quadrature_order = quadrature_order
        self.green_q = self._green_density_matrix(quadrature_order)

    @property
    def size(self):
        return 2 * self.N + 4

    def unpack(self, state):
        q = state[:self.N]
        psi = state[self.N:2 * self.N]
        lam, mu, drift, offset = state[2 * self.N:]
        return q, psi, lam, mu, drift, offset

    def _green_density_matrix(self, order: int) -> sparse.csr_matrix:
        """Map nodal Q directly to centered Green boundary values."""
        if order < 1:
            raise ValueError("quadrature order must be positive")
        xb = self.xx[self.artificial][:, None]
        yb = self.yy[self.artificial][:, None]
        nb = len(self.boundary_indices)
        nodal_source = np.zeros((nb, self.N), dtype=float)
        dx = np.diff(self.x)
        dy = np.diff(self.y)
        area = (dy[:, None] * dx[None, :]).ravel()
        bottom_left = (
            np.arange(self.ny - 1)[:, None] * self.nx
            + np.arange(self.nx - 1)[None, :]).ravel()
        corners = [
            bottom_left,
            bottom_left + 1,
            bottom_left + self.nx,
            bottom_left + self.nx + 1,
        ]
        nodes, weights = np.polynomial.legendre.leggauss(order)
        fractions = 0.5 * (nodes + 1.0)
        mapped_weights = 0.5 * weights
        points = [
            (a, b, wa * wb)
            for b, wb in zip(fractions, mapped_weights)
            for a, wa in zip(fractions, mapped_weights)
        ]

        for tx, ty, quadrature_weight in points:
            xi = (self.x[:-1][None, :]
                  + tx * dx[None, :])
            eta = (self.y[:-1][:, None]
                   + ty * dy[:, None])
            xi = np.broadcast_to(xi, (self.ny - 1, self.nx - 1)).ravel()
            eta = np.broadcast_to(eta, (self.ny - 1, self.nx - 1)).ravel()
            direct = (xb - xi[None, :])**2 + (yb - eta[None, :])**2
            kernel = np.log1p(4.0 * yb * eta[None, :] / direct) / (4.0 * np.pi)
            kernel -= yb * eta[None, :] / (
                np.pi * (xi[None, :]**2 + eta[None, :]**2))
            weighted = kernel * (quadrature_weight * area)[None, :]
            shapes = [
                (1.0 - tx) * (1.0 - ty),
                tx * (1.0 - ty),
                (1.0 - tx) * ty,
                tx * ty,
            ]
            for indices, shape_weight in zip(corners, shapes):
                nodal_source[:, indices] += shape_weight * weighted

        # The quadrature first maps nodal Q_x to boundary psi.  Compose with
        # the same discrete Q_x used in Poisson, preserving exact linearity.
        boundary_q = (self.Dx.T @ nodal_source.T).T
        rows = np.repeat(self.boundary_indices, self.N)
        cols = np.tile(np.arange(self.N), nb)
        return sparse.csr_matrix(
            (boundary_q.ravel(), (rows, cols)), shape=(self.N, self.N))

    def evaluate(self, state, jacobian=False):
        q, psi, lam, mu, drift, offset = self.unpack(state)
        if not (lam > mu > 0):
            raise FloatingPointError("rates left 0 < mu < lambda")
        qx, qy = self.Dx @ q, self.Dy @ q
        psix, psiy = self.Dx @ psi, self.Dy @ psi
        ux, uy = -psiy, psix
        vx = ux + lam * self.xx + drift
        vy = uy + lam * self.yy
        transport = vx * qx + vy * qy - mu * q - offset
        poisson = -self.Lap @ psi - qx
        green_values = self.green_q @ q
        raw_psi = poisson.copy()
        raw_psi[self.wall] = psi[self.wall]
        raw_psi[self.artificial] = (
            psi[self.artificial] - green_values[self.artificial])
        gauges = self.C @ q - self.target
        residual = np.r_[
            self.q_weight * transport,
            self.psi_weight * raw_psi,
            self.gauge_scale * gauges,
        ]
        cache = dict(
            q=q, psi=psi, qx=qx, qy=qy, psix=psix, psiy=psiy,
            vx=vx, vy=vy, transport=transport, poisson=poisson,
            greenValues=green_values, rawPsi=raw_psi, gauges=gauges,
            lam=lam, mu=mu, drift=drift, offset=offset,
            alpha=mu / lam,
        )
        if not jacobian:
            return residual, cache
        return residual, cache, self._assemble_jacobian(cache, True)

    def _assemble_jacobian(self, cache, include_green: bool):
        q = cache["q"]
        qx, qy = cache["qx"], cache["qy"]
        vx, vy = cache["vx"], cache["vy"]
        mu = cache["mu"]
        transport_qq = (
            sparse.diags(vx) @ self.Dx
            + sparse.diags(vy) @ self.Dy - mu * self.I)
        transport_qpsi = (
            -sparse.diags(qx) @ self.Dy + sparse.diags(qy) @ self.Dx)
        q_params = np.column_stack([
            self.xx * qx + self.yy * qy,
            -q,
            qx,
            -np.ones(self.N),
        ])
        psi_q = -self.Spoisson @ self.Dx
        if include_green:
            psi_q = psi_q - self.green_q
        psi_psi = (
            -self.Spoisson @ self.Lap + self.Swall + self.Sartificial)
        matrix = sparse.bmat([
            [sparse.diags(self.q_weight) @ transport_qq,
             sparse.diags(self.q_weight) @ transport_qpsi,
             sparse.diags(self.q_weight) @ sparse.csr_matrix(q_params)],
            [sparse.diags(self.psi_weight) @ psi_q,
             sparse.diags(self.psi_weight) @ psi_psi,
             sparse.csr_matrix((self.N, 4))],
            [sparse.diags(self.gauge_scale) @ self.C,
             sparse.csr_matrix((4, self.N)), sparse.csr_matrix((4, 4))],
        ], format="csc")
        return matrix

    def local_jacobian(self, state):
        """Sparse PDE Jacobian with the nonlocal Green update removed."""
        _, cache = self.evaluate(state)
        return self._assemble_jacobian(cache, False)

    def metrics(self, state):
        residual, c = self.evaluate(state)
        core = (abs(self.xx) <= 2.0) & (self.yy <= 2.0)
        vx = c["vx"].reshape(self.shape)
        vy = c["vy"].reshape(self.shape)
        outward = np.r_[-vx[1:-1, 0], vx[1:-1, -1], vy[-1, :]]
        boundary_error = (
            c["psi"][self.artificial]
            - c["greenValues"][self.artificial])
        return dict(
            residualInf=float(np.max(abs(residual))),
            residualRms=float(np.linalg.norm(residual) / np.sqrt(residual.size)),
            weightedQInf=float(np.max(abs(self.q_weight * c["transport"]))),
            weightedPsiInf=float(np.max(abs(self.psi_weight * c["rawPsi"]))),
            transportInf=float(np.max(abs(c["transport"]))),
            coreTransportInf=float(np.max(abs(c["transport"][core]))),
            poissonInf=float(np.max(abs(c["poisson"][self.poisson]))),
            greenBoundaryInf=float(np.max(abs(boundary_error))),
            minimumOutwardVelocity=float(np.min(outward)),
            maximumOutwardVelocity=float(np.max(outward)),
            inflowBoundaryNodes=int(np.count_nonzero(outward <= 0.0)),
            artificialBoundaryFaceNodes=int(outward.size),
            gaugeInf=float(np.max(abs(c["gauges"]))),
            centerDerivative=float((self.Dy.getrow(self.origin) @ c["psi"]).item()),
            lambda_=float(c["lam"]), mu=float(c["mu"]),
            v=float(c["drift"]), b=float(c["offset"]),
            alpha=float(c["alpha"]), harmonic=0.0, logCoefficient=0.0,
            offsetDriftDefect=float(c["offset"] - c["drift"]),
            maxAbsQ=float(np.max(abs(c["q"]))),
            maxAbsPsi=float(np.max(abs(c["psi"]))),
        )


def native_seed(length, nx, ny, core_step, stencil, quadrature_order):
    profile = load_native(20)
    x, y, ax, ay = stretched_grid(length, nx, ny, core_step)
    q = RectBivariateSpline(
        profile["y"], profile["x"], profile["Q"])(y, x)
    psi = RectBivariateSpline(
        profile["y"], profile["x"], profile["psi"])(y, x)
    problem = CenteredGreenBVP(
        x, y, profile["rates"], stencil, quadrature_order)
    envelope = np.exp(-(problem.X / 1.5)**8 - (problem.Y / 1.5)**4)
    basis = np.column_stack([
        (envelope * problem.X**degree).ravel() for degree in range(4)])
    qflat = q.ravel()
    qflat += basis @ np.linalg.solve(
        problem.C @ basis, problem.target - problem.C @ qflat)
    state = np.r_[qflat, psi.ravel(), profile["rates"]]
    metadata = dict(
        source=profile["source"], length=length, nx=nx, ny=ny,
        requestedCoreStep=core_step,
        actualCoreStepX=float(x[problem.origin + 1]),
        actualCoreStepY=float(y[1]), stretchX=ax, stretchY=ay,
        quadratureOrder=quadrature_order,
        greenBoundaryRows=int(len(problem.boundary_indices)),
        greenMatrixNonzeros=int(problem.green_q.nnz),
    )
    return problem, state, metadata


def directional_check(problem, state):
    _, _, matrix = problem.evaluate(state, jacobian=True)
    rng = np.random.default_rng(290929)
    direction = rng.standard_normal(state.size)
    direction /= np.linalg.norm(direction)
    step = 1e-4
    plus = problem.evaluate(state + step * direction)[0]
    minus = problem.evaluate(state - step * direction)[0]
    numerical = (plus - minus) / (2.0 * step)
    exact = matrix @ direction
    return float(np.linalg.norm(numerical - exact)
                 / max(np.linalg.norm(exact), 1e-300))


def linear_solve(problem, state, matrix, right_hand_sides,
                 relative_tolerance=1e-9):
    """Solve Green-coupled systems with the local PDE Jacobian as a preconditioner."""
    local = problem.local_jacobian(state)
    factorization = splu(local, permc_spec="COLAMD")
    preconditioner = LinearOperator(matrix.shape, matvec=factorization.solve)
    right_hand_sides = np.asarray(right_hand_sides, float)
    vector_input = right_hand_sides.ndim == 1
    if vector_input:
        right_hand_sides = right_hand_sides[:, None]
    solutions = np.empty_like(right_hand_sides)
    diagnostics = []
    for column in range(right_hand_sides.shape[1]):
        history = []
        rhs = right_hand_sides[:, column]
        solution, info = gmres(
            matrix, rhs, M=preconditioner,
            rtol=relative_tolerance, atol=0.0,
            restart=80, maxiter=12,
            callback=history.append, callback_type="pr_norm")
        defect = matrix @ solution - rhs
        diagnostics.append(dict(
            info=int(info), iterations=len(history),
            preconditionedResidual=(float(history[-1])
                                    if history else 0.0),
            relativeResidual=float(
                np.linalg.norm(defect)
                / max(np.linalg.norm(rhs), 1e-300))))
        if (info != 0
                and diagnostics[-1]["relativeResidual"]
                > max(1e-6, 20.0 * relative_tolerance)):
            raise RuntimeError(
                f"GMRES failed for RHS {column}: info={info}, "
                f"relative residual={diagnostics[-1]['relativeResidual']:g}")
        solutions[:, column] = solution
    if vector_input:
        solutions = solutions[:, 0]
    return solutions, diagnostics


def solve(problem, initial, options, right_hand_side=None, callback=None):
    """Damped Newton solve using the Green-aware iterative linear solver."""
    state = np.asarray(initial, float).copy()
    rhs = (np.zeros(problem.size) if right_hand_side is None
           else np.asarray(right_hand_side, float))
    history = []
    status = "maximum_iterations"
    for iteration in range(options.maximum_iterations + 1):
        physical_residual, _, matrix = problem.evaluate(state, jacobian=True)
        residual = physical_residual - rhs
        metrics = problem.metrics(state)
        row = dict(
            iteration=iteration,
            residualInf=float(np.max(np.abs(residual))),
            residualRms=float(np.linalg.norm(residual) / np.sqrt(residual.size)),
            physicalResidualInf=metrics.pop("residualInf"),
            physicalResidualRms=metrics.pop("residualRms"),
            **metrics,
        )
        history.append(row)
        if callback:
            callback(row)
        if row["residualInf"] < options.tolerance:
            status = "converged"
            break
        if iteration == options.maximum_iterations:
            break
        try:
            step, diagnostics = linear_solve(
                problem, state, matrix, -residual)
        except RuntimeError as exc:
            row["linearError"] = str(exc)
            status = "linear_solve_failed"
            break
        row["linearIterations"] = diagnostics[0]["iterations"]
        row["linearResidualRelative"] = diagnostics[0]["relativeResidual"]
        length = 1.0
        old_norm = np.linalg.norm(residual)
        accepted = False
        while length >= options.minimum_step:
            trial = state + length * step
            try:
                trial_physical, _ = problem.evaluate(trial)
            except FloatingPointError:
                length *= 0.5
                continue
            trial_residual = trial_physical - rhs
            if (np.isfinite(trial_residual).all()
                    and np.linalg.norm(trial_residual)
                    < (1.0 - 1e-4 * length) * old_norm):
                accepted = True
                break
            length *= 0.5
        row["stepLength"] = float(length)
        row["stepInf"] = float(np.max(np.abs(step)))
        if not accepted:
            status = "line_search_failed"
            break
        state = trial
    return state, status, history


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length", type=float, default=64.0)
    parser.add_argument("--nx", type=int, default=129)
    parser.add_argument("--ny", type=int, default=65)
    parser.add_argument("--core-step", type=float, default=0.125)
    parser.add_argument("--stencil", type=int, default=5)
    parser.add_argument("--quadrature-order", type=int, default=2)
    parser.add_argument("--iterations", type=int, default=24)
    parser.add_argument("--tolerance", type=float, default=2e-8)
    parser.add_argument(
        "--taus", default="0,0.02,0.05,0.1,0.2,0.35,0.5,0.65,0.8,0.9,1")
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--out", type=Path, required=False,
                        default=ROOT / "runs/streamfunction_global_green_20260929/L64")
    args = parser.parse_args()
    problem, state, metadata = native_seed(
        args.length, args.nx, args.ny, args.core_step,
        args.stencil, args.quadrature_order)
    check = directional_check(problem, state)
    summary = dict(
        metadata=metadata, jacobianDirectionalRelativeError=check,
        initialMetrics=problem.metrics(state))
    print(json.dumps(summary, indent=2), flush=True)
    if check > 3e-6:
        raise AssertionError(f"Jacobian directional error {check:g}")
    if args.audit_only:
        return
    if args.out.exists():
        raise ValueError(f"refusing to overwrite {args.out}")
    args.out.mkdir(parents=True)
    config = dict(
        **summary, stencil=args.stencil,
        scriptSha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        formulation="Full transport plus centered half-plane Green boundary closure")
    dump(args.out / "config.json", config)
    initial_residual = problem.evaluate(state)[0]
    taus = [float(value) for value in args.taus.split(",")]
    if (not taus or taus[0] != 0.0 or taus[-1] != 1.0
            or any(b <= a for a, b in zip(taus, taus[1:]))):
        raise ValueError("taus must increase strictly from 0 to 1")
    options = NewtonOptions(args.iterations, args.tolerance)
    history_rows = []
    start = time.time()
    status = "not_started"
    final_tau = None
    with (args.out / "iterations.jsonl").open("w") as stream:
        for tau in taus:
            rhs = (1.0 - tau) * initial_residual

            def callback(row):
                tagged = dict(homotopyTau=tau, **row)
                history_rows.append(tagged)
                stream.write(json.dumps(tagged) + "\n")
                stream.flush()
                print(
                    f"tau={tau:.3f} iter={row['iteration']:02d} "
                    f"res={row['residualInf']:.3e} "
                    f"rates=({row['lambda_']:.8f},{row['mu']:.8f},"
                    f"{row['v']:.8f}) alpha={row['alpha']:.8f} "
                    f"out={row['minimumOutwardVelocity']:.5g}", flush=True)

            state, status, _ = solve(
                problem, state, options, right_hand_side=rhs,
                callback=callback)
            final_tau = tau
            np.savez_compressed(
                args.out / f"stage_tau_{tau:.3f}.npz",
                x=problem.x, y=problem.y,
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
        iterations=len(history_rows) - len(taus), finalMetrics=final)
    dump(args.out / "report.json", report)
    np.savez_compressed(
        args.out / "state.npz", x=problem.x, y=problem.y,
        Q=state[:problem.N].reshape(problem.shape),
        psi=state[problem.N:2 * problem.N].reshape(problem.shape),
        parameters=state[2 * problem.N:], residual=problem.evaluate(state)[0])
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
