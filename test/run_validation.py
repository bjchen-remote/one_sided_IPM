#!/usr/bin/env python3
"""Run the bounded, single-thread ipm_cr elliptic qualification suite."""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matlab', default='/Applications/MATLAB_R2026a.app/bin/matlab')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--timeout', type=int, default=180)
    args = parser.parse_args()
    folder = Path(__file__).resolve().parent
    stamp = dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    output = (args.output or folder / 'results' / stamp).resolve()
    output.mkdir(parents=True, exist_ok=False)
    originals = list((folder.parent / 'time_evolution' / '+ipm').rglob('*.m'))
    originals += list((folder.parent / 'time_evolution' / '+ipm_c').rglob('*.m'))
    before = {str(p): digest(p) for p in originals}
    inputs = list((folder / '+ipm_cr').rglob('*.m'))
    inputs += list((folder / 'tests').rglob('*.m')) + [folder / 'main_cr.m', Path(__file__).resolve()]
    receipt = dict(status='prepared', startedUtc=dt.datetime.now(dt.timezone.utc).isoformat(),
                   threads=1, hardWallSeconds=args.timeout,
                   inputHashes={str(p.relative_to(folder)): digest(p) for p in inputs},
                   productionSourceHashes=before, trajectoryAdvanced=False)

    def persist():
        (output / 'execution.json').write_text(json.dumps(receipt, indent=2) + '\n')

    def mq(path):
        return str(path).replace("'", "''")

    command = [args.matlab, '-batch',
               "restoredefaultpath;maxNumCompThreads(1);addpath('" + mq(folder) + "','" +
               mq(folder / 'tests') + "');r=test_ipm_cr('" + mq(output) +
               "',struct('includeIsolatedRefinement',true));assert(r.passed);[~,~,d]=main_cr;assert(d.accepted);"]
    receipt['command'] = command
    persist()
    start = time.monotonic()
    child = None
    try:
        with (output / 'console.log').open('w') as log:
            child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                     start_new_session=True,
                                     env=dict(os.environ, OMP_NUM_THREADS='1',
                                              OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))
            receipt.update(status='running', pid=child.pid)
            persist()
            try:
                code = child.wait(timeout=args.timeout)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGTERM)
                try:
                    child.wait(timeout=8)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                    child.wait(timeout=8)
                raise RuntimeError('MATLAB exceeded bounded wall time')
        receipt['exitCode'] = code
        receipt['productionSourceUnchanged'] = all(digest(Path(p)) == h for p, h in before.items())
        receipt['testInputsUnchanged'] = all(digest(folder / p) == h
                                             for p, h in receipt['inputHashes'].items())
        if code != 0:
            raise RuntimeError('MATLAB failed; inspect console.log')
        if not receipt['productionSourceUnchanged'] or not receipt['testInputsUnchanged']:
            raise RuntimeError('Source changed while validation ran')
        report = json.loads((output / 'report.json').read_text())
        if not report['passed']:
            raise RuntimeError('Qualification did not pass')
        receipt['status'] = 'completed'
    except BaseException as error:
        receipt.update(status='failed', error=f'{type(error).__name__}: {error}')
        raise
    finally:
        if child is not None and child.poll() is None:
            os.killpg(child.pid, signal.SIGTERM)
            child.wait(timeout=8)
        receipt.update(seconds=time.monotonic() - start,
                       endedUtc=dt.datetime.now(dt.timezone.utc).isoformat())
        persist()
    print(json.dumps(dict(status=receipt['status'], output=str(output), seconds=receipt['seconds'])))


if __name__ == '__main__':
    main()
