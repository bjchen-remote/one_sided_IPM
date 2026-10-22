"""Direct global self-similar BVP solve on a stretched half-plane box.

The grid is fine near the normalized core and approximately logarithmic in
the far field.  Q and the *centered* streamfunction are solved simultaneously.
The density equation is retained on every grid point, including the
artificial boundary.  This is the correct characteristic closure when the
effective similarity velocity points outward there.  Only the elliptic
streamfunction needs a far-field condition:

    z psi_x + y psi_y = (1+alpha) psi - alpha h y.

The scalar h retains the admissible linear harmonic y mode, while
psi_y(0,0)=0 fixes the centered gauge.  Together with the four Q gauges this
gives a square nonlinear system for Q, psi, lambda, mu, v, b, and h.
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
from scipy.interpolate import RectBivariateSpline
from scipy.sparse.linalg import splu

from fixed_point import fd_matrix
from streamfunction_sources import load_native


ROOT = Path(__file__).resolve().parents[2]


def dump(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def stretch_parameter(length: float, intervals: int, core_step: float) -> float:
    """Choose a in x=L*sinh(a*s)/sinh(a) for the requested first step."""
    delta = 2.0 / intervals
    target = core_step / length
    lo, hi = 1e-8, 40.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        value = np.sinh(mid * delta) / np.sinh(mid)
        if value > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def stretched_grid(length: float, nx: int, ny: int, core_step: float):
    if nx % 2 != 1 or ny < 9:
        raise ValueError("nx must be odd and both dimensions must be nontrivial")
    ax = stretch_parameter(length, nx - 1, core_step)
    # y has half the parameter interval, so its normalized spacing is 1/(ny-1).
    ay = stretch_parameter(length, 2 * (ny - 1), core_step)
    sx = np.linspace(-1.0, 1.0, nx)
    sy = np.linspace(0.0, 1.0, ny)
    x = length * np.sinh(ax * sx) / np.sinh(ax)
    y = length * np.sinh(ay * sy) / np.sinh(ay)
    return x, y, ax, ay


@dataclass
class NewtonOptions:
    maximum_iterations: int = 24
    tolerance: float = 2e-8
    minimum_step: float = 2.0**-18


class CenteredSelfSimilarBVP:
    def __init__(self, x: np.ndarray, y: np.ndarray, reference_rates: np.ndarray,
                 stencil: int = 5, density_robin: bool = False,
                 resonant_log: bool = False):
        self.x = np.asarray(x, float)
        self.y = np.asarray(y, float)
        self.nx, self.ny = len(x), len(y)
        self.shape = (self.ny, self.nx)
        self.N = self.nx * self.ny
        self.X, self.Y = np.meshgrid(self.x, self.y)
        self.xx, self.yy = self.X.ravel(), self.Y.ravel()
        self.origin = int(np.argmin(abs(self.x)))
        if abs(self.x[self.origin]) > 1e-13 or abs(self.y[0]) > 1e-13:
            raise ValueError("grid must contain the wall core")
        ix = sparse.eye(self.nx, format="csr")
        iy = sparse.eye(self.ny, format="csr")
        self.Dx = sparse.kron(iy, fd_matrix(self.x, 1, stencil), format="csr")
        self.Dy = sparse.kron(fd_matrix(self.y, 1, stencil), ix, format="csr")
        self.Lap = (
            sparse.kron(iy, fd_matrix(self.x, 2, stencil), format="csr")
            + sparse.kron(fd_matrix(self.y, 2, stencil), ix, format="csr")
        )
        self.I = sparse.eye(self.N, format="csr")
        wall = np.zeros(self.shape, bool)
        wall[0, :] = True
        artificial = np.zeros(self.shape, bool)
        artificial[-1, :] = True
        artificial[:, [0, -1]] = True
        self.wall = wall.ravel()
        self.artificial_q = artificial.ravel()
        self.artificial_psi = (artificial & ~wall).ravel()
        self.transport = ~self.artificial_q
        self.poisson = ~(self.wall | self.artificial_psi)
        self.Stransport = sparse.diags(self.transport.astype(float))
        self.Stailq = sparse.diags(self.artificial_q.astype(float))
        self.Spoisson = sparse.diags(self.poisson.astype(float))
        self.Swall = sparse.diags(self.wall.astype(float))
        self.Stailpsi = sparse.diags(self.artificial_psi.astype(float))
        self.density_robin = bool(density_robin)
        self.resonant_log = bool(resonant_log)
        # At the m=1 Mellin resonance alpha=1/2 is imposed exactly.  Replacing
        # the dependent rate mu by the logarithmic coefficient keeps the
        # bordered problem square and avoids a singular constraint row.
        self.parameter_count = 5

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
        cnorm = np.sqrt(np.asarray(self.C.multiply(self.C).sum(axis=1))).ravel()
        self.gauge_scale = 1.0 / np.maximum(cnorm, 1.0)
        self.center_row = self.Dy.getrow(self.origin)
        center_norm = np.sqrt(self.center_row.multiply(self.center_row).sum())
        self.center_scale = 1.0 / max(float(center_norm), 1.0)

        lam0, mu0, drift0, _ = map(float, reference_rates)
        alpha0 = mu0 / lam0
        power0 = 1.0 + alpha0
        radius = np.hypot(self.xx + drift0 / lam0, self.yy)
        self.q_weight = (1.0 + radius)**(-alpha0)
        self.psi_weight = np.empty(self.N)
        self.psi_weight[self.poisson] = (1.0 + radius[self.poisson])**(1.0 - alpha0)
        self.psi_weight[self.wall] = (1.0 + radius[self.wall])**(-power0)
        self.psi_weight[self.artificial_psi] = (
            (1.0 + radius[self.artificial_psi])**(-power0))

    @property
    def size(self) -> int:
        return 2 * self.N + self.parameter_count

    def unpack(self, state: np.ndarray):
        q = state[:self.N]
        psi = state[self.N:2 * self.N]
        parameters = state[2 * self.N:]
        if self.resonant_log:
            lam, drift, offset, harmonic, log_coefficient = parameters
            mu = 0.5 * lam
        else:
            lam, mu, drift, offset, harmonic = parameters
            log_coefficient = 0.0
        return (q, psi, lam, mu, drift, offset, harmonic,
                log_coefficient)

    def evaluate(self, state: np.ndarray, jacobian: bool = False):
        (q, psi, lam, mu, drift, offset, harmonic,
         log_coefficient) = self.unpack(state)
        if not (lam > mu > 0):
            raise FloatingPointError("rates left 0 < mu < lambda")
        alpha = 0.5 if self.resonant_log else mu / lam
        power = 1.0 + alpha
        zcoord = self.xx + drift / lam
        shifted_q = q + offset / mu
        qx, qy = self.Dx @ q, self.Dy @ q
        psix, psiy = self.Dx @ psi, self.Dy @ psi
        ux, uy = -psiy, psix
        vx = ux + lam * self.xx + drift
        vy = uy + lam * self.yy
        transport = vx * qx + vy * qy - mu * q - offset
        q_tail = zcoord * qx + self.yy * qy - alpha * shifted_q
        raw_q = transport.copy()
        if self.density_robin:
            raw_q[self.artificial_q] = q_tail[self.artificial_q]
        poisson = -self.Lap @ psi - qx
        radius = np.hypot(zcoord, self.yy)
        psi_tail = (
            zcoord * psix + self.yy * psiy - power * psi
            + alpha * harmonic * self.yy
        )
        if self.resonant_log:
            psi_tail += (log_coefficient * self.yy
                         * (alpha * np.log(radius) - 1.0))
        raw_psi = poisson.copy()
        raw_psi[self.wall] = psi[self.wall]
        raw_psi[self.artificial_psi] = psi_tail[self.artificial_psi]
        gauges = self.C @ q - self.target
        center = (self.center_row @ psi).item()
        scalar_residual = [
            *(self.gauge_scale * gauges), self.center_scale * center,
        ]
        residual = np.r_[
            self.q_weight * raw_q,
            self.psi_weight * raw_psi,
            scalar_residual,
        ]
        cache = dict(
            q=q, psi=psi, qx=qx, qy=qy, psix=psix, psiy=psiy,
            vx=vx, vy=vy, alpha=alpha, power=power, zcoord=zcoord,
            shifted_q=shifted_q, transport=transport, q_tail=q_tail,
            poisson=poisson, psi_tail=psi_tail, rawQ=raw_q,
            rawPsi=raw_psi, gauges=gauges, center=center,
            lam=lam, mu=mu, drift=drift, offset=offset,
            harmonic=harmonic,
            logCoefficient=log_coefficient, radius=radius,
        )
        if not jacobian:
            return residual, cache
        return residual, cache, self.jacobian(cache)

    def jacobian(self, c: dict) -> sparse.csc_matrix:
        lam, mu = c["lam"], c["mu"]
        drift, offset = c["drift"], c["offset"]
        harmonic = c["harmonic"]
        log_coefficient = c["logCoefficient"]
        alpha = c["alpha"]
        q, qx, qy = c["q"], c["qx"], c["qy"]
        psi, psix = c["psi"], c["psix"]
        zcoord, shifted_q = c["zcoord"], c["shifted_q"]
        transport_qq = (
            sparse.diags(c["vx"]) @ self.Dx
            + sparse.diags(c["vy"]) @ self.Dy - mu * self.I)
        transport_qpsi = (
            -sparse.diags(qx) @ self.Dy + sparse.diags(qy) @ self.Dx)
        tail_qq = (
            sparse.diags(zcoord) @ self.Dx
            + sparse.diags(self.yy) @ self.Dy - alpha * self.I)
        if self.density_robin:
            Jqq = self.Stransport @ transport_qq + self.Stailq @ tail_qq
            Jqpsi = self.Stransport @ transport_qpsi
        else:
            Jqq = transport_qq
            Jqpsi = transport_qpsi
        if self.resonant_log:
            q_params = np.column_stack([
                self.xx * qx + self.yy * qy - 0.5 * q,
                qx,
                -np.ones(self.N),
                np.zeros(self.N),
                np.zeros(self.N),
            ])
            q_tail_params = np.column_stack([
                -(drift / lam**2) * qx + offset / lam**2,
                qx / lam,
                -np.ones(self.N) / lam,
                np.zeros(self.N),
                np.zeros(self.N),
            ])
        else:
            q_params = np.column_stack([
                self.xx * qx + self.yy * qy,
                -q,
                qx,
                -np.ones(self.N),
                np.zeros(self.N),
            ])
            q_tail_params = np.column_stack([
                -(drift / lam**2) * qx + (alpha / lam) * shifted_q,
                -q / lam,
                qx / lam,
                -np.ones(self.N) / lam,
                np.zeros(self.N),
            ])
        if self.density_robin:
            q_params[self.artificial_q] = q_tail_params[self.artificial_q]

        Jpsi_q = -self.Spoisson @ self.Dx
        tail_psipsi = (
            sparse.diags(zcoord) @ self.Dx
            + sparse.diags(self.yy) @ self.Dy - c["power"] * self.I)
        Jpsi_psi = (
            -self.Spoisson @ self.Lap + self.Swall
            + self.Stailpsi @ tail_psipsi)
        psi_params = np.zeros((self.N, self.parameter_count))
        if self.resonant_log:
            radius2 = c["radius"]**2
            log_radius = np.log(c["radius"])
            tail_psi_params = np.column_stack([
                -(drift / lam**2) * psix
                - log_coefficient * self.yy * alpha * drift * zcoord
                / (lam**2 * radius2),
                psix / lam + log_coefficient * self.yy * alpha * zcoord
                / (lam * radius2),
                np.zeros(self.N),
                alpha * self.yy,
                self.yy * (alpha * log_radius - 1.0),
            ])
        else:
            tail_psi_params = np.column_stack([
                -(drift / lam**2) * psix
                + (alpha / lam) * (psi - harmonic * self.yy),
                (-psi + harmonic * self.yy) / lam,
                psix / lam,
                np.zeros(self.N),
                alpha * self.yy,
            ])
        psi_params[self.artificial_psi] = (
            tail_psi_params[self.artificial_psi])

        Wq = sparse.diags(self.q_weight)
        Wpsi = sparse.diags(self.psi_weight)
        zeros_4n = sparse.csr_matrix((4, self.N))
        zeros_4p = sparse.csr_matrix((4, self.parameter_count))
        center_params = sparse.csr_matrix((1, self.parameter_count))
        blocks = [
            [Wq @ Jqq, Wq @ Jqpsi, Wq @ sparse.csr_matrix(q_params)],
            [Wpsi @ Jpsi_q, Wpsi @ Jpsi_psi,
             Wpsi @ sparse.csr_matrix(psi_params)],
            [sparse.diags(self.gauge_scale) @ self.C, zeros_4n, zeros_4p],
            [sparse.csr_matrix((1, self.N)),
             self.center_scale * self.center_row, center_params],
        ]
        return sparse.bmat(blocks, format="csc")

    def metrics(self, state: np.ndarray) -> dict:
        residual, c = self.evaluate(state)
        core = (abs(self.xx) <= 2.0) & (self.yy <= 2.0)
        vx = c["vx"].reshape(self.shape)
        vy = c["vy"].reshape(self.shape)
        # Face-normal characteristic speeds on the three artificial faces.
        # Corners are counted on the top face and the physical wall is omitted.
        outward = np.r_[
            -vx[1:-1, 0], vx[1:-1, -1], vy[-1, :],
        ]
        return dict(
            residualInf=float(np.max(abs(residual))),
            residualRms=float(np.linalg.norm(residual) / np.sqrt(residual.size)),
            weightedQInf=float(np.max(abs(self.q_weight * c["rawQ"]))),
            weightedPsiInf=float(np.max(abs(self.psi_weight * c["rawPsi"]))),
            transportInf=float(np.max(abs(c["transport"]))),
            coreTransportInf=float(np.max(abs(c["transport"][core]))),
            poissonInf=float(np.max(abs(c["poisson"][self.poisson]))),
            densityTailInf=float(np.max(abs(c["q_tail"][self.artificial_q]))),
            streamfunctionTailInf=float(
                np.max(abs(c["psi_tail"][self.artificial_psi]))),
            minimumOutwardVelocity=float(np.min(outward)),
            maximumOutwardVelocity=float(np.max(outward)),
            inflowBoundaryNodes=int(np.count_nonzero(outward <= 0.0)),
            artificialBoundaryFaceNodes=int(outward.size),
            gaugeInf=float(np.max(abs(c["gauges"]))),
            centerDerivative=float(c["center"]),
            lambda_=float(c["lam"]), mu=float(c["mu"]),
            v=float(c["drift"]), b=float(c["offset"]),
            alpha=float(c["alpha"]), harmonic=float(c["harmonic"]),
            logCoefficient=float(c["logCoefficient"]),
            alphaConstraint=0.0 if self.resonant_log else None,
            offsetDriftDefect=float(c["offset"] - c["drift"]),
            maxAbsQ=float(np.max(abs(c["q"]))),
            maxAbsPsi=float(np.max(abs(c["psi"]))),
        )


def native_seed(length: float, nx: int, ny: int, core_step: float,
                stencil: int, density_robin: bool = False,
                resonant_log: bool = False):
    profile = load_native(20)
    x, y, ax, ay = stretched_grid(length, nx, ny, core_step)
    q_spline = RectBivariateSpline(profile["y"], profile["x"], profile["Q"])
    psi_spline = RectBivariateSpline(
        profile["y"], profile["x"], profile["psi"])
    q = q_spline(y, x)
    psi = psi_spline(y, x)
    problem = CenteredSelfSimilarBVP(
        x, y, profile["rates"], stencil, density_robin=density_robin,
        resonant_log=resonant_log)

    # Enforce the four Q gauges by a localized correction only.
    envelope = np.exp(-(problem.X / 1.5)**8 - (problem.Y / 1.5)**4)
    basis = np.column_stack(
        [(envelope * problem.X**degree).ravel() for degree in range(4)])
    qflat = q.ravel()
    qflat += basis @ np.linalg.solve(
        problem.C @ basis, problem.target - problem.C @ qflat)

    # Remove the small interpolation/FD centering defect by an exact y mode.
    psiflat = psi.ravel()
    center = (problem.center_row @ psiflat).item()
    psiflat -= center * problem.yy

    lam, mu, drift, offset = map(float, profile["rates"])
    alpha = 0.5 if resonant_log else mu / lam
    power = 1.0 + alpha
    psix, psiy = problem.Dx @ psiflat, problem.Dy @ psiflat
    zcoord = problem.xx + drift / lam
    base = zcoord * psix + problem.yy * psiy - power * psiflat
    direction = alpha * problem.yy
    mask = problem.artificial_psi
    weight = problem.psi_weight[mask]
    harmonic = -float(
        (weight * direction[mask]) @ (weight * base[mask])) / float(
        (weight * direction[mask]) @ (weight * direction[mask]))
    if resonant_log:
        state = np.r_[qflat, psiflat, lam, drift, offset, harmonic, 0.0]
    else:
        state = np.r_[qflat, psiflat, profile["rates"], harmonic]
    metadata = dict(
        source=profile["source"], length=length, nx=nx, ny=ny,
        requestedCoreStep=core_step,
        actualCoreStepX=float(x[problem.origin + 1]),
        actualCoreStepY=float(y[1]), stretchX=ax, stretchY=ay,
        interpolationCenterCorrection=float(center),
        initialHarmonic=float(harmonic),
        densityBoundary=("leading_homogeneous_robin" if density_robin
                         else "full_transport_outflow"),
        resonantLog=bool(resonant_log),
    )
    return problem, state, metadata


def directional_check(problem: CenteredSelfSimilarBVP,
                      state: np.ndarray) -> float:
    _, _, matrix = problem.evaluate(state, jacobian=True)
    rng = np.random.default_rng(290929)
    direction = rng.standard_normal(state.size)
    direction /= np.linalg.norm(direction)
    # The stretched-grid Laplacian has coefficients spanning many decades;
    # a larger directional step avoids roundoff in the far-field rows.
    step = 1e-3
    plus = problem.evaluate(state + step * direction)[0]
    minus = problem.evaluate(state - step * direction)[0]
    numerical = (plus - minus) / (2 * step)
    exact = matrix @ direction
    return float(np.linalg.norm(numerical - exact)
                 / max(np.linalg.norm(exact), 1e-300))


def solve(problem: CenteredSelfSimilarBVP, initial: np.ndarray,
          options: NewtonOptions, right_hand_side: np.ndarray | None = None,
          callback=None):
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
            residualInf=float(np.max(abs(residual))),
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
            factor = splu(matrix, permc_spec="COLAMD")
            step = factor.solve(-residual)
        except RuntimeError as exc:
            row["linearError"] = str(exc)
            status = "factorization_failed"
            break
        row["linearResidualRelative"] = float(
            np.linalg.norm(matrix @ step + residual)
            / max(np.linalg.norm(residual), 1e-300))
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
        row["stepInf"] = float(np.max(abs(step)))
        if not accepted:
            status = "line_search_failed"
            break
        state = trial
    return state, status, history


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length", type=float, default=1024.0)
    parser.add_argument("--nx", type=int, default=129)
    parser.add_argument("--ny", type=int, default=65)
    parser.add_argument("--core-step", type=float, default=0.125)
    parser.add_argument("--stencil", type=int, default=5)
    parser.add_argument("--iterations", type=int, default=24)
    parser.add_argument("--tolerance", type=float, default=2e-8)
    parser.add_argument(
        "--taus",
        default="0,0.02,0.05,0.1,0.2,0.35,0.5,0.65,0.8,0.9,1")
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--density-robin", action="store_true")
    parser.add_argument("--resonant-log", action="store_true")
    parser.add_argument(
        "--out", type=Path,
        default=ROOT / "runs/streamfunction_global_bvp_20260929/native_L1024_129x65")
    args = parser.parse_args()
    problem, state, metadata = native_seed(
        args.length, args.nx, args.ny, args.core_step, args.stencil,
        density_robin=args.density_robin,
        resonant_log=args.resonant_log)
    check = directional_check(problem, state)
    initial = problem.metrics(state)
    summary = dict(metadata=metadata,
                   jacobianDirectionalRelativeError=check,
                   initialMetrics=initial)
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
        formulation=(
            "One-patch coupled Q/centered-psi BVP on a stretched box; "
            "large-radius self-similar Robin closure; explicit linear y "
            "mode; no frozen exterior velocity."),
    )
    dump(args.out / "config.json", config)
    options = NewtonOptions(args.iterations, args.tolerance)
    start = time.time()
    initial_residual = problem.evaluate(state)[0]
    taus = [float(value) for value in args.taus.split(",")]
    if (not taus or taus[0] != 0.0 or taus[-1] != 1.0
            or any(b <= a for a, b in zip(taus, taus[1:]))):
        raise ValueError("taus must increase strictly from 0 to 1")
    config["residualWithdrawalHomotopy"] = taus
    dump(args.out / "config.json", config)
    all_history = []
    status = "not_started"
    with (args.out / "iterations.jsonl").open("w") as stream:
        for tau in taus:
            rhs = (1.0 - tau) * initial_residual
            def callback(row):
                tagged = dict(homotopyTau=tau, **row)
                all_history.append(tagged)
                stream.write(json.dumps(tagged) + "\n")
                stream.flush()
                print(
                    f"tau={tau:.3f} iter={row['iteration']:02d} "
                    f"res={row['residualInf']:.3e} "
                    f"physical={row['physicalResidualInf']:.3e} "
                    f"rates=({row['lambda_']:.8f},{row['mu']:.8f},"
                    f"{row['v']:.8f},{row['b']:.8f}) "
                    f"alpha={row['alpha']:.8f} h={row['harmonic']:.5g}",
                    flush=True)
            state, status, history = solve(
                problem, state, options, right_hand_side=rhs,
                callback=callback)
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
        converged=bool(status == "converged" and all_history[-1]["homotopyTau"] == 1.0),
        finalHomotopyTau=all_history[-1]["homotopyTau"],
        seconds=time.time() - start, finalMetrics=final,
        iterations=len(all_history) - len(taus),
    )
    dump(args.out / "report.json", report)
    q = state[:problem.N].reshape(problem.shape)
    psi = state[problem.N:2 * problem.N].reshape(problem.shape)
    np.savez_compressed(
        args.out / "state.npz", x=problem.x, y=problem.y,
        Q=q, psi=psi, parameters=state[2 * problem.N:],
        qx=(problem.Dx @ q.ravel()).reshape(problem.shape),
        residual=problem.evaluate(state)[0])
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
