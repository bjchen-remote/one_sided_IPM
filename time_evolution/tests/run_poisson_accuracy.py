#!/usr/bin/env python3
"""Bounded, one-thread MATLAB analytic Poisson test; never advances a run."""
import argparse
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = PACKAGE.parent
GIB = 1024**3


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resources():
    vm = subprocess.check_output(['vm_stat'], text=True)
    page = int(re.search(r'page size of (\d+) bytes', vm).group(1))
    names = ('Pages free', 'Pages inactive', 'Pages speculative', 'Pages purgeable')
    available = page * sum(int(re.search(r'^'+n+r':\s+(\d+)', vm, re.M).group(1)) for n in names)
    rows = subprocess.check_output(['ps', '-axo', 'pid=,rss=,comm='], text=True).splitlines()
    native = [r.split(None, 2) for r in rows if '/Contents/MacOS/MATLAB_' in r]
    return {'reclaimableBytes': available,
            'matlab': [{'pid': int(r[0]), 'rssBytes': int(r[1])*1024} for r in native]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    inputs = list((PACKAGE / '+ipm').rglob('*.m')) + [
        PACKAGE / 'main.m', Path(__file__).resolve(), PACKAGE / 'tests/check_poisson_accuracy.m']
    pins = {str(p): sha(p) for p in inputs}
    receipt = {'kind': 'bounded_standard_poisson_accuracy_execution_v1', 'status': 'prepared',
               'startedUtc': dt.datetime.now(dt.timezone.utc).isoformat(), 'inputHashes': pins,
               'threads': 1, 'hardWallSeconds': 180, 'maximumOwnedRssBytes': int(1.5*GIB),
               'minimumInitialReclaimableBytes': int(2.5*GIB),
               'minimumRunningReclaimableBytes': int(.75*GIB), 'maximumTotalMatlab': 2,
               'sharedLock': str(ROOT / 'experiments/fast_slow_high_order/.runner.lock'),
               'trajectoryAdvanced': False}

    def persist():
        (out / 'execution.json').write_text(json.dumps(receipt, indent=2)+'\n')

    persist()
    child = None
    began = time.monotonic()
    with open(receipt['sharedLock'], 'a+') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            before = resources(); receipt['resourcesBefore'] = before
            if len(before['matlab']) > 1 or before['reclaimableBytes'] < 2.5*GIB:
                raise RuntimeError('Insufficient resources for the bounded small-grid test')
            quote = lambda p: str(p).replace("'", "''")
            expression = ("restoredefaultpath;maxNumCompThreads(1);addpath('"+quote(PACKAGE)+
                "','-begin');addpath('"+quote(PACKAGE / 'tests')+"','-begin');"+
                "assert(startsWith(which('ipm.field.poisson'),'"+quote(PACKAGE)+"/'));"+
                "r=check_poisson_accuracy('"+quote(out)+"');assert(r.completed);")
            command = ['/Applications/MATLAB_R2026a.app/bin/matlab', '-batch', expression]
            env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
            with (out / 'console.log').open('w') as log:
                child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                         start_new_session=True, env=env)
                receipt.update(status='running', pid=child.pid, command=command); persist()
                peak = 0
                while child.poll() is None:
                    snapshot = resources()
                    owned = [p for p in snapshot['matlab'] if p['pid'] == child.pid]
                    peak = max([peak] + [p['rssBytes'] for p in owned])
                    reason = None
                    if time.monotonic()-began > 180: reason = 'hard wall timeout'
                    if peak > 1.5*GIB: reason = 'owned MATLAB RSS limit'
                    if snapshot['reclaimableBytes'] < .75*GIB: reason = 'memory reserve'
                    if len(snapshot['matlab']) > 2: reason = 'global MATLAB count'
                    if reason:
                        receipt['guardStop'] = reason
                        raise RuntimeError(reason)
                    time.sleep(1)
            receipt.update(exitCode=child.returncode, peakObservedRssBytes=peak,
                           protectedInputsUnchanged=all(sha(Path(p)) == h for p, h in pins.items()))
            assert child.returncode == 0, 'MATLAB exited with an error; inspect console.log'
            assert receipt['protectedInputsUnchanged'], 'A test input changed'
            report = json.loads((out / 'report.json').read_text())
            assert report['completed']
            receipt['status'] = 'completed'
        except Exception as error:
            receipt.update(status='failed', error=type(error).__name__+': '+str(error))
            raise
        finally:
            if child is not None and child.poll() is None:
                os.killpg(child.pid, signal.SIGTERM)
                try: child.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL); child.wait(timeout=10)
                receipt['exitCode'] = child.returncode
            receipt['seconds'] = time.monotonic()-began
            receipt['endedUtc'] = dt.datetime.now(dt.timezone.utc).isoformat()
            persist()
            fcntl.flock(lock, fcntl.LOCK_UN)
    print(json.dumps({'status': receipt['status'], 'output': str(out),
                      'seconds': receipt['seconds'], 'peakRss': receipt['peakObservedRssBytes']}))


if __name__ == '__main__':
    main()
