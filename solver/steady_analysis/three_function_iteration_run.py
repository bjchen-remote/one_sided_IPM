"""Bounded physical IPM test of the three-function primal/dual iteration.

The sibling agent's finite-box evaluator is imported READ ONLY. This runner
owns a separate output directory. It ties the width exponent to the CURRENT
rates, uses piecewise Gaussian integration, and never treats a small update
or a vanishing projected gradient as a solution of the full PDE.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import solve as dense_solve

ROOT = Path(__file__).resolve().parents[2]
BACKEND = ROOT / 'runs/three_function_iteration_20260928/backend_snapshot.py'
OUT = ROOT / 'runs/three_function_iteration_20260928'
spec = importlib.util.spec_from_file_location('three_function_pilot_backend', BACKEND)
backend = importlib.util.module_from_spec(spec)
spec.loader.exec_module(backend)
BACKEND_HASH = hashlib.sha256(BACKEND.read_bytes()).hexdigest()


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')


class PhysicalModel(backend.ThreeFunctions):
    """Same strict three functions; no y-dependent mother correction modes."""
    def __init__(self, *args, tie_power=True, **kwargs):
        super().__init__(*args, **kwargs)
        self.tie_power = tie_power
        self.gauss = {n: leggauss(n) for n in (12, 24)}

    def select_power(self, z):
        if self.fit_power:
            self.power = float(np.exp(z[-1]))
        elif self.tie_power:
            lam, mu, _, _ = self.selected_rates(z)
            self.power = lam / (lam-mu)

    def evaluate(self, z, cache=False):
        self.select_power(z)
        return super().evaluate(z, cache)

    def admissible(self, z):
        self.select_power(z)
        return super().admissible(z)

    def primitive(self, x, hs, order=12):
        """Integrate each cubic-spline piece separately, anchored at xi=0."""
        gx, gw = self.gauss[order] if order in self.gauss else leggauss(order)
        k = self.xknots
        width = np.diff(k)
        points = k[:-1, None]+width[:, None]*(gx+1)/2
        spans = width*np.sum(self.mother(points, hs)*gw, axis=-1)/2
        prefix = np.r_[0., np.cumsum(spans)]
        index = np.clip(np.searchsorted(k, x, side='right')-1, 0, len(k)-2)
        left = k[index]
        points = left[..., None]+(x-left)[..., None]*(gx+1)/2
        partial = (x-left)*np.sum(self.mother(points, hs)*gw, axis=-1)/2
        return prefix[index]-prefix[self.zero]+partial

    def jacobian(self, z, r=None, epsilon=1e-5):
        r = self.evaluate(z) if r is None else r
        J = np.empty((len(r), len(z)))
        for j in range(len(z)):
            dz = epsilon*max(1., abs(z[j]))
            zp, zm = z.copy(), z.copy()
            zp[j] += dz
            zm[j] -= dz
            J[:, j] = (self.evaluate(zp)-self.evaluate(zm))/(2*dz)
        self.select_power(z)
        return J


def diagnostic(model, z):
    r, state = model.evaluate(z, cache=True)
    J = model.jacobian(z, r)
    direction = np.random.default_rng(280928).normal(size=len(z))
    direction /= np.linalg.norm(direction)
    eps = 3e-6
    direct = (model.evaluate(z+eps*direction)-model.evaluate(z-eps*direction))/(2*eps)
    jac_error = np.linalg.norm(J@direction-direct)/np.linalg.norm(direct)
    hs, aa, mm = model.unpack(z)
    model.select_power(z)
    xi = (model.problem.x[None, :]-mm(model.problem.y)[:, None])*np.exp(model.power*aa(model.problem.y))[:, None]
    quadrature = np.max(abs(model.primitive(xi, hs, 12)-model.primitive(xi, hs, 24)))
    mu = model.selected_rates(z)[1]
    matrix = (state['vy'][:, :, None]*model.D[:, None, :]-mu*model.E[:, None, :]).reshape(model.problem.N, -1)
    c_check = np.linalg.lstsq(matrix, -state['raw0'].ravel(), rcond=None)[0]
    density_error = np.max(abs(c_check-state['ccoef']))
    U, s, _ = np.linalg.svd(J, full_matrices=False)
    rank = int(np.sum(s > s[0]*1e-8))
    normal = r-U[:, :rank]@(U[:, :rank].T@r)
    result = dict(jacobian_direction_relative=float(jac_error),
                  primitive_12_24_max_difference=float(quadrature),
                  density_projection_difference=float(density_error),
                  evaluator_agreement=state['metrics']['evaluator_agreement'],
                  linearized_rank=rank, unknowns=len(z),
                  jacobian_singular_values=s.tolist(),
                  residual_norm=float(np.linalg.norm(r)),
                  normal_to_jacobian_residual_norm=float(np.linalg.norm(normal)),
                  normal_residual_fraction=float(np.linalg.norm(normal)/np.linalg.norm(r)))
    assert jac_error < .003 and density_error < 1e-7 and quadrature < 1e-7
    return result


def iterate(model, z0, count, label, mode, tau, seconds):
    z = z0.copy()
    r0 = model.evaluate(z)
    dual = np.zeros_like(r0)
    J0 = model.jacobian(z, r0)
    # Fixed positive metric, retaining cross-function coupling. The floor
    # removes unobservable outer-spline directions, not residual equations.
    gram = J0.T@J0
    floor = 1e-5*max(float(np.linalg.eigvalsh(gram)[-1]), 1e-8)
    metric = gram+floor*np.eye(len(z))
    start = time.monotonic()
    history = []
    status = 'iteration_limit'
    for it in range(count+1):
        r, state = model.evaluate(z, cache=True)
        row = dict(iteration=it, **state['metrics'],
                   residual_norm=float(np.linalg.norm(r)),
                   dual_norm=float(np.linalg.norm(dual)),
                   dual_increment_if_taken=float(np.linalg.norm(r)),
                   exponent=float(model.power), elapsed=time.monotonic()-start)
        history.append(row)
        if it % 5 == 0:
            print(label, it, 'rms', f"{row['rms']:.9g}", 'inf', f"{row['inf']:.9g}",
                  'dual', f"{row['dual_norm']:.9g}", flush=True)
        if row['inf'] < 1e-7 and max(abs(v) for v in row['discrete_gauges']) < 1e-6:
            status = 'full_residual_root'
            break
        if it == count:
            break
        if time.monotonic()-start > seconds:
            status = 'time_bound'
            break
        J = J0 if it == 0 else model.jacobian(z, r)
        gradient = J.T@r
        row['projected_gradient'] = float(np.linalg.norm(gradient))
        row['metric_gradient'] = float(np.sqrt(max(0., gradient@dense_solve(metric, gradient, assume_a='pos'))))
        accepted = False
        if mode == 'primal':
            for damping in (1e-4, 1e-2, 1., 100.):
                step = dense_solve(J.T@J+damping*metric, -gradient, assume_a='pos')
                for alpha in (1., .5, .25, .125, .03125, .0078125):
                    trial = z+alpha*step
                    if not model.admissible(trial):
                        continue
                    rt = model.evaluate(trial)
                    if np.linalg.norm(rt) < np.linalg.norm(r)*(1-1e-6*alpha):
                        accepted = True
                        delta = alpha*step
                        row.update(alpha=alpha, damping=damping)
                        break
                if accepted:
                    break
        else:
            # gamma=1; a fixed tau is used whenever the declared safeguards
            # permit it. No residual-decrease test: dual dynamics can oscillate.
            for retry in range(9):
                effective_tau = tau*4.**(-retry)
                delta = dense_solve(metric/effective_tau+J.T@J,
                                    -J.T@(dual+r), assume_a='pos')
                trial = z+delta
                if np.linalg.norm(delta) > 2. or not model.admissible(trial):
                    continue
                rt = model.evaluate(trial)
                if np.isfinite(rt).all() and np.linalg.norm(rt) < max(.05, 5*np.linalg.norm(r)):
                    row.update(tau=effective_tau, safeguard_retries=retry)
                    accepted = True
                    break
        if not accepted:
            status = 'safeguard_blocked' if mode == 'dual' else 'primal_stationarity_or_line_search_stall'
            break
        row['step_norm'] = float(np.linalg.norm(delta))
        row['metric_step_norm'] = float(np.sqrt(delta@metric@delta))
        z = trial
        if mode == 'dual':
            dual += rt
        np.savez_compressed(OUT/f'{label}_checkpoint.npz', z=z, dual=dual, metric=metric, iteration=it+1)
        dump(OUT/f'{label}_history.json', history)
    r, state = model.evaluate(z, cache=True)
    final = dict(label=label, mode=mode, requested_tau=tau, gamma=1., status=status,
                 iterations=len(history)-1, seconds=time.monotonic()-start,
                 initial=history[0], final=history[-1], history=history,
                 diagnostic=diagnostic(model, z))
    np.savez_compressed(OUT/f'{label}_state.npz', z=z, dual=dual, metric=metric,
                        x=model.problem.x, y=model.problem.y,
                        rates=model.selected_rates(z),
                        **{k:v for k, v in state.items() if isinstance(v, np.ndarray) and k != 'selected_rates'})
    dump(OUT/f'{label}.json', final)
    print('FINISHED', label, status, final['iterations'],
          'rms', final['final']['rms'], 'dual', final['final']['dual_norm'], flush=True)
    return z, final


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nx', type=int, default=65)
    ap.add_argument('--ny', type=int, default=33)
    ap.add_argument('--shape', type=int, default=12)
    ap.add_argument('--height', type=int, default=7)
    ap.add_argument('--primal', type=int, default=25)
    ap.add_argument('--dual', type=int, default=60)
    ap.add_argument('--tau', type=float, default=1.)
    ap.add_argument('--seconds', type=float, default=60.)
    ap.add_argument('--fixed-rates', action='store_true')
    ap.add_argument('--calibrated-power', action='store_true')
    ap.add_argument('--free-power', action='store_true')
    ap.add_argument('--direct', action='store_true')
    ap.add_argument('--seed', type=Path)
    ap.add_argument('--label', required=True)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    p, Qref, rates, qfun, provenance = backend.source_problem(args.nx, args.ny)
    control = backend.solve(p, np.r_[p.normalize_seed(Qref), rates],
                            backend.SolverOptions(max_iterations=12))
    zcontrol = control.pop('z')
    dump(OUT/f'{args.label}_control.json', control)
    if control['box_root']:
        rates = zcontrol[p.N:]
    yy = np.linspace(0, 1, 33)
    geo = backend.geometry(qfun, yy)
    lx = -np.log(geo[:, 1]/geo[0, 1])
    ly = np.log(geo[:, 2]/geo[0, 2])
    calibrated = float(lx@ly/(lx@lx))
    power = calibrated if args.calibrated_power or args.free_power else rates[0]/(rates[0]-rates[1])
    model = PhysicalModel(p, rates, qfun, args.shape, args.height, power,
                          free_rates=not args.fixed_rates, smooth_c=True,
                          fit_power=args.free_power,
                          tie_power=not (args.calibrated_power or args.free_power))
    z = model.z0 if args.seed is None else np.load(args.seed)['z']
    assert len(z) == model.size
    config = dict(provenance=provenance, backend=str(BACKEND), backend_sha256=BACKEND_HASH,
                  runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  seed=None if args.seed is None else str(args.seed.resolve()),
                  arguments={k:str(v) if isinstance(v, Path) else v for k,v in vars(args).items()},
                  fixed_rates=args.fixed_rates, free_power=args.free_power,
                  dynamic_theory_exponent=not (args.calibrated_power or args.free_power),
                  power=power, calibrated_power=calibrated, coefficients=model.size,
                  variables='h(xi), log A(y), m(y), plus finite rate parameters when free',
                  normalization='f(0)=1, f_prime(0)=0, f_second(0)=-2, A(0)=1, m(0)=C(0)=0',
                  density_trace='C(y) globally eliminated in the same cubic height space',
                  residual='full box PDE/sqrt(N) and .5 times the three discrete derivative gauges',
                  exterior='frozen response, retained source response updated, not half-plane closure',
                  control_root=control['box_root'], control_metrics=control['metrics'],
                  initial_checks=diagnostic(model, z))
    dump(OUT/f'{args.label}_config.json', config)
    print('CONFIG', args.label, 'control', config['control_root'],
          'control_inf', config['control_metrics']['active_pde_residual_inf'],
          'unknowns', model.size, 'power', power, flush=True)
    if not args.direct and args.primal:
        z, _ = iterate(model, z, args.primal, args.label+'_primal', 'primal', args.tau, args.seconds)
    if args.dual:
        iterate(model, z, args.dual, args.label+'_dual', 'dual', args.tau, args.seconds)


if __name__ == '__main__':
    main()
