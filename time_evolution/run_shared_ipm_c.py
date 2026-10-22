#!/usr/bin/env python3
"""Bounded regression of shared IPM time integration with polar operators."""
from pathlib import Path
import argparse
import fcntl
import hashlib
import importlib.util
import json
import shutil
import signal

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--restart-only', action='store_true')
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location('shared_ipm_executor', HERE / 'run_ipm_c.py')
    executor = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(executor)
    signal.signal(signal.SIGTERM, executor.interrupted)
    root = args.output.resolve()
    q = executor.quote
    with (HERE / '.ipm_c_runner.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if not args.restart_only:
            root.mkdir(parents=True, exist_ok=False)
            files = sorted((HERE / '+ipm').rglob('*.m')) + sorted((HERE / '+ipm_c').rglob('*.m'))
            files += sorted((HERE / 'tests_ipm_c').glob('*.m'))
            files += [HERE / 'tests' / 'check_shared_spatial_sdirk.m',
                      HERE / 'tests' / 'fixtures' / 'stepRk_before_spatial.m',
                      HERE / 'main_c.m', Path(__file__).resolve()]
            records = []
            for path in files:
                relative = path.relative_to(HERE)
                target = root / 'source_snapshot' / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, target)
                records.append(dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
            executor.atomic(root / 'source_inventory.json', records)
            body = ("addpath('" + q(HERE / 'tests') + "','" + q(HERE / 'tests_ipm_c') + "');"
                    "r=run_shared_components();ipm_c.output.atomicJson('" + q(root / 'components.json') + "',r);"
                    "r=test_shared_workflow('" + q(root / 'native') + "');assert(r.passed);"
                    "r=test_shared_default_step();ipm_c.output.atomicJson('" + q(root / 'default_step.json') + "',r);")
            executor.run_matlab(root, 'shared_components_native', body, seconds=300)
        records = json.loads((root / 'source_inventory.json').read_text())
        for record in records:
            if hashlib.sha256(Path(record['path']).read_bytes()).hexdigest() != record['sha256']:
                raise RuntimeError('Source changed during regression: ' + record['path'])
        body = ("addpath('" + q(HERE / 'tests_ipm_c') + "');"
                "r=test_shared_restart('" + q(root / 'native') + "');assert(r.passed);")
        executor.run_matlab(root, 'shared_cold_restart', body, seconds=180)
        for record in records:
            if hashlib.sha256(Path(record['path']).read_bytes()).hexdigest() != record['sha256']:
                raise RuntimeError('Source changed during restart regression: ' + record['path'])
        summary = dict(passed=True, architecture='shared_ipm_sdirk2_with_polar_spatial_operators',
                       components=str(root / 'components.json'), native=str(root / 'native' / 'native_report.json'),
                       defaultStep=str(root / 'default_step.json'),
                       restart=str(root / 'native' / 'restart_report.json'), longTrajectoryQualified=False)
        executor.atomic(root / 'summary.json', summary)
        print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
