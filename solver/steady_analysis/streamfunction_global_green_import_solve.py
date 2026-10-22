"""Polish an imported same-grid state against the centered Green IPM BVP."""
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

from streamfunction_global_green_solve import (
    CenteredGreenBVP,
    directional_check,
    solve,
)
from streamfunction_global_selfsimilar_solve import NewtonOptions, dump


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=Path, required=True)
    parser.add_argument("--stencil", type=int, default=5)
    parser.add_argument("--quadrature-order", type=int, default=4)
    parser.add_argument("--iterations", type=int, default=50)
    parser.add_argument("--tolerance", type=float, default=2e-8)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError(f"refusing to overwrite {args.out}")
    args.out.mkdir(parents=True)

    with np.load(args.seed) as data:
        x = np.asarray(data["x"], float)
        y = np.asarray(data["y"], float)
        q = np.asarray(data["Q"], float)
        psi = np.asarray(data["psi"], float)
        parameters = np.asarray(data["parameters"], float)
    if q.shape != (len(y), len(x)) or psi.shape != q.shape:
        raise ValueError("seed arrays do not match their coordinate vectors")
    if parameters.size < 4:
        raise ValueError("seed must contain lambda, mu, drift, and offset")

    problem = CenteredGreenBVP(
        x, y, parameters[:4], args.stencil, args.quadrature_order)
    state = np.r_[q.ravel(), psi.ravel(), parameters[:4]]
    audit = directional_check(problem, state)
    initial = problem.metrics(state)
    config = dict(
        seed=str(args.seed.resolve()), stencil=args.stencil,
        quadratureOrder=args.quadrature_order,
        discardedSeedParameters=parameters[4:].tolist(),
        jacobianDirectionalRelativeError=audit,
        initialMetrics=initial,
        scriptSha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    dump(args.out / "config.json", config)
    print(json.dumps(config, indent=2), flush=True)

    rows = []
    start = time.time()
    with (args.out / "iterations.jsonl").open("w") as stream:
        def callback(row):
            rows.append(row)
            stream.write(json.dumps(row) + "\n")
            stream.flush()
            print(
                f"iter={row['iteration']:02d} res={row['residualInf']:.3e} "
                f"rates=({row['lambda_']:.9g},{row['mu']:.9g},"
                f"{row['v']:.9g}) alpha={row['alpha']:.9g} "
                f"out={row['minimumOutwardVelocity']:.6g}", flush=True)

        state, status, _ = solve(
            problem, state, NewtonOptions(args.iterations, args.tolerance),
            callback=callback)
    final = problem.metrics(state)
    report = dict(
        status=status, converged=bool(status == "converged"),
        elapsedSeconds=time.time() - start,
        iterations=max(len(rows) - 1, 0), finalMetrics=final)
    dump(args.out / "report.json", report)
    np.savez_compressed(
        args.out / "state.npz", x=x, y=y,
        Q=state[:problem.N].reshape(problem.shape),
        psi=state[problem.N:2 * problem.N].reshape(problem.shape),
        parameters=state[2 * problem.N:],
        residual=problem.evaluate(state)[0])
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
