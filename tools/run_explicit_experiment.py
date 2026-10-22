#!/usr/bin/env python3
"""Freeze and supervise the public ipm_cr.run trajectory on macOS or Linux."""
from pathlib import Path
import argparse
import datetime as dt
import fcntl
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from common import (atomic, available_memory, freeze, initial_input, matlab_executable,
                    now, process_rss, quote, require_kernel, source_files,
                    stop_child, verify_snapshot)

ROOT = Path(__file__).resolve().parent.parent


def prepare(folder, args):
    options = json.loads(args.options.read_text()) if args.options else {}
    if not isinstance(options, dict):
        raise ValueError('Numerical options must be a JSON object')
    authority = initial_input(ROOT)
    inputs = source_files(ROOT) + [authority, Path('tools/common.py'), Path('tools/run_explicit_experiment.py')]
    if not args.prepare_only:
        require_kernel(ROOT, options)
    folder.mkdir(parents=True, exist_ok=False)
    inventory = freeze(ROOT, folder / 'runtime', inputs)
    initial_input(folder / 'runtime')
    atomic(folder / 'source_inventory.json', inventory)
    atomic(folder / 'explicit_options.json', options)
    # Options are immutable input, independent of mutable mission/progress records.
    from common import sha
    atomic(folder / 'options_identity.json', {'sha256': sha(folder / 'explicit_options.json')})
    deadline = dt.datetime.now(dt.timezone.utc) + dt.timedelta(hours=args.hours)
    atomic(folder / 'mission.json', dict(
        kind='ipm_cr_frozen_run_v2', status='prepared', targetTau=args.target,
        createdUtc=now(), deadlineUtc=deadline.isoformat(), segmentSeconds=args.segment_seconds,
        threads=1, rssLimitBytes=2 * 1024**3, minimumFreeDiskBytes=8 * 1024**3,
        matlab=args.matlab, sourceFrozen=True,
        timeIntegrator=options.get('timeIntegrator', 'rk4'),
        timeControl='all_active_both_patches_all_stages_physical_CFL',
        implicitStages=0, polarOriginPresent=False, trajectoryAccuracyCertified=False,
        options=options, segments=[]))


def dispatch(folder, matlab):
    entry = folder / 'runtime/tools/run_explicit_experiment.py'
    with (folder / 'supervisor.log').open('ab') as log:
        child = subprocess.Popen([sys.executable, str(entry), '_supervise', '--output', str(folder),
                                  '--matlab', matlab], stdout=log, stderr=subprocess.STDOUT,
                                 start_new_session=True)
    print(json.dumps(dict(output=str(folder), supervisorPid=child.pid)), flush=True)


def supervise(folder, matlab):
    def interrupted(_signum, _frame):
        raise InterruptedError('Supervisor interrupted')
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    with (folder / 'owner.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        mission = json.loads((folder / 'mission.json').read_text())
        if mission.get('kind') != 'ipm_cr_frozen_run_v2':
            raise RuntimeError('Resume requires this v2 frozen run protocol')
        from common import sha
        expected_options = json.loads((folder / 'options_identity.json').read_text())['sha256']
        if sha(folder / 'explicit_options.json') != expected_options:
            raise RuntimeError('Frozen numerical options changed')
        if json.loads((folder / 'explicit_options.json').read_text()) != mission['options']:
            raise RuntimeError('Mission options differ from the frozen numerical input')
        require_kernel(folder / 'runtime', mission['options'])
        frozen = json.loads((folder / 'source_inventory.json').read_text())
        verify_snapshot(folder / 'runtime', frozen)
        expected_sources = sorted(Path(path) for path in frozen if path.startswith('src/') or path == 'data/initial.json')
        if source_files(folder / 'runtime') != expected_sources:
            raise RuntimeError('Frozen source file set differs from its inventory')
        previous_wall = mission.get('supervisorWallSeconds', 0)
        native = folder / 'native'
        native.mkdir(exist_ok=True)
        deadline = dt.datetime.fromisoformat(mission['deadlineUtc'])
        child = None
        awake = None
        if sys.platform == 'darwin' and Path('/usr/bin/caffeinate').is_file():
            awake = subprocess.Popen(['/usr/bin/caffeinate', '-i', '-w', str(os.getpid())],
                                     stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        mission['idleSleepInhibitedDuringRun'] = awake is not None
        started = time.monotonic()
        try:
            while dt.datetime.now(dt.timezone.utc) < deadline:
                if (folder / 'STOP').exists() or (native / 'STOP').exists():
                    mission['status'] = 'stopped'
                    break
                assert shutil.disk_usage(folder).free >= mission['minimumFreeDiskBytes'], 'Insufficient free disk'
                assert available_memory() >= 1024**3, 'Insufficient available memory'
                verify_snapshot(folder / 'runtime', frozen)
                verify_snapshot(folder, {'explicit_options.json': expected_options})
                if source_files(folder / 'runtime') != expected_sources:
                    raise RuntimeError('Frozen source file set changed')
                index = len(mission['segments']) + 1
                script = folder / 'segment.m'
                script.write_text(
                    "restoredefaultpath;maxNumCompThreads(1);\n"
                    "addpath('" + quote(folder / 'runtime/src') + "');\n"
                    "opts=struct('outputDir','" + quote(native) + "','tauEnd'," + repr(mission['targetTau']) + ","
                    "'wallSeconds'," + repr(mission['segmentSeconds']) + ",'checkpointEveryTau',5,'checkpointWallSeconds',60,'logEveryTau',5);\n"
                    "opts.explicitOptions=jsondecode(fileread('" + quote(folder / 'explicit_options.json') + "'));\n"
                    "summary=ipm_cr.run(opts);disp(summary.status);\n")
                entry = dict(index=index, startedUtc=now(), maximumRssBytes=0)
                mission['segments'].append(entry)
                mission.update(status='running', supervisorPid=os.getpid(), updatedUtc=now())
                atomic(folder / 'mission.json', mission)
                segment_started = time.monotonic()
                reason = None
                with (folder / ('segment_%05d.log' % index)).open('xb') as log:
                    child = subprocess.Popen([matlab, '-singleCompThread', '-batch', "run('" + quote(script) + "')"],
                                             stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                                             env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))
                    mission['matlabPid'] = child.pid
                    while child.poll() is None:
                        entry['maximumRssBytes'] = max(entry['maximumRssBytes'], process_rss(child.pid))
                        if (folder / 'STOP').exists():
                            reason = 'stopped'
                        elif dt.datetime.now(dt.timezone.utc) >= deadline:
                            reason = 'wall_budget'
                        elif entry['maximumRssBytes'] > mission['rssLimitBytes'] or available_memory() < 512 * 1024**2:
                            reason = 'memory_budget'
                        elif shutil.disk_usage(folder).free < mission['minimumFreeDiskBytes']:
                            reason = 'disk_budget'
                        elif time.monotonic() - segment_started > mission['segmentSeconds'] + 180:
                            reason = 'segment_timeout'
                        if reason:
                            stop_child(child, native)
                            break
                        cursor = native / 'latest.json'
                        if cursor.is_file():
                            cp = json.loads(cursor.read_text())
                            mission.update(tau=cp['tau'], steps=cp['sequence'])
                        mission['updatedUtc'] = now()
                        atomic(folder / 'mission.json', mission)
                        time.sleep(5)
                    entry.update(returncode=child.wait(), elapsedSeconds=time.monotonic() - segment_started, endedUtc=now())
                mission.pop('matlabPid', None)
                verify_snapshot(folder / 'runtime', frozen)
                verify_snapshot(folder, {'explicit_options.json': expected_options})
                if source_files(folder / 'runtime') != expected_sources:
                    raise RuntimeError('Frozen source file set changed')
                if reason:
                    # A final checkpoint can be committed while MATLAB is
                    # still shutting down. Prefer that durable completion.
                    final_report = native / 'summary.json'
                    final_cursor = native / 'latest.json'
                    completed = False
                    if final_report.is_file() and final_cursor.is_file():
                        report = json.loads(final_report.read_text())
                        cursor = json.loads(final_cursor.read_text())
                        completed = (report.get('status') == 'target_reached'
                                     and cursor.get('tau', 0) >= mission['targetTau'])
                        mission.update(tau=cursor['tau'], steps=cursor['sequence'])
                    mission['status'] = 'target_reached' if completed else reason
                    mission['stopReason'] = reason
                    break
                if entry['returncode']:
                    mission['status'] = 'numerical_or_runtime_failure'
                    break
                report = json.loads((native / 'summary.json').read_text())
                mission.update(tau=report['tau'], steps=report['steps'], updatedUtc=now())
                if report['status'] != 'wall_budget':
                    mission['status'] = report['status']
                    break
                atomic(folder / 'mission.json', mission)
            else:
                mission['status'] = 'wall_budget'
        except BaseException as error:
            mission.update(status='supervisor_failure', error=repr(error))
            raise
        finally:
            stop_child(child, native)
            if awake is not None:
                awake.terminate()
                awake.wait(timeout=5)
            mission.update(updatedUtc=now(), invocationSupervisorWallSeconds=time.monotonic() - started,
                           supervisorWallSeconds=previous_wall + time.monotonic() - started)
            atomic(folder / 'mission.json', mission)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('start', 'resume', 'status', '_supervise'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--options', type=Path)
    parser.add_argument('--matlab', help='MATLAB executable; otherwise discover it on this host')
    parser.add_argument('--target', type=float, default=200)
    parser.add_argument('--hours', type=float, default=50/60)
    parser.add_argument('--segment-seconds', type=float, default=600)
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    if args.prepare_only and args.action != 'start':
        parser.error('--prepare-only is only valid with start')
    folder = args.output.resolve()
    if args.action == 'start':
        if not (args.target > 0 and args.hours > 0 and args.segment_seconds > 0):
            parser.error('Target and budgets must be positive')
        if not args.prepare_only:
            args.matlab = matlab_executable(args.matlab)
        prepare(folder, args)
        if args.prepare_only:
            print(json.dumps(dict(status='prepared', output=str(folder), dispatched=False)))
        else:
            dispatch(folder, args.matlab)
    elif args.action == 'resume':
        if (folder / 'STOP').exists() or (folder / 'native/STOP').exists():
            raise RuntimeError('Remove the requested STOP before resuming')
        mission = json.loads((folder / 'mission.json').read_text())
        executable = matlab_executable(args.matlab or mission.get('matlab'))
        dispatch(folder, executable)
    elif args.action == '_supervise':
        supervise(folder, matlab_executable(args.matlab))
    else:
        print(json.dumps(json.loads((folder / 'mission.json').read_text()), ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
