#!/usr/bin/env python3
"""Freeze and supervise the public ipm_cr.run trajectory on macOS or Linux."""
from pathlib import Path
import argparse
import datetime as dt
import fcntl
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from common import (atomic, available_memory, freeze, initial_input, matlab_executable,
                    now, process_rss, quote, relative_file, require_kernel, source_files,
                    sha, stop_child, verify_snapshot)

ROOT = Path(__file__).resolve().parent.parent


def read_options(path):
    def reject(value):
        raise ValueError('Non-finite JSON constant is not allowed: ' + value)
    return json.loads(path.read_text(), parse_constant=reject)


def seed_input(options, case_format):
    """Validate an optional actual-time accepted-state seed, without loading MATLAB."""
    if not case_format:
        return None
    run = options.get('run', {})
    if not isinstance(run, dict):
        raise ValueError('case.run must be an object')
    seed = run.get('seed')
    if seed is None or seed == []:
        return None
    if (not isinstance(seed, dict) or set(seed) != {'file', 'sha256', 'tau'}
            or not isinstance(seed['file'], str) or not seed['file']
            or not isinstance(seed['sha256'], str)
            or re.fullmatch('[0-9a-f]{64}', seed['sha256']) is None
            or not isinstance(seed['tau'], (int, float)) or isinstance(seed['tau'], bool)
            or not math.isfinite(seed['tau']) or seed['tau'] < 0):
        raise ValueError('Seed requires a file, lowercase SHA256 and finite nonnegative actual tau')
    return seed


def freeze_seed(folder, options, case_format):
    seed = seed_input(options, case_format)
    if seed is None:
        return None
    original = Path(seed['file']).expanduser().resolve()
    if not original.is_file() or sha(original) != seed['sha256']:
        raise RuntimeError('Seed file differs from its declared SHA256')
    destination = folder / 'input_seed.mat'
    temporary = destination.with_name(destination.name + '.writing')
    shutil.copy2(original, temporary)
    if sha(original) != seed['sha256'] or sha(temporary) != seed['sha256']:
        raise RuntimeError('Seed changed while preparing its immutable copy')
    os.replace(temporary, destination)
    if sha(destination) != seed['sha256']:
        raise RuntimeError('Frozen seed readback differs')
    # The immutable option uses a run-relative name; each segment binds it
    # to this folder, so moving a complete run does not lose the seed.
    record = dict(kind='ipm_cr_frozen_seed_v1', originalFile=str(original),
                  frozenFile=destination.name, sha256=seed['sha256'], tau=seed['tau'])
    seed['file'] = destination.name
    return record


def verify_frozen(folder, mission, frozen, expected_options):
    """Check all frozen inputs before/after each segment and at cold resume."""
    runtime = folder / 'runtime'
    verify_snapshot(runtime, frozen)
    verify_snapshot(folder, {'explicit_options.json': expected_options})
    if read_options(folder / 'explicit_options.json') != mission['options']:
        raise RuntimeError('Mission options differ from the frozen numerical input')
    # source_files is the only source-set definition, including mesh.json.
    # Missing/changed original files fail above; newly added sources fail here.
    if not set(source_files(runtime)).issubset(Path(name) for name in frozen):
        raise RuntimeError('Frozen source file set differs from its inventory')
    seed = seed_input(mission['options'], mission.get('caseFormat'))
    record = mission.get('seedInput')
    if seed is None:
        if record is not None:
            raise RuntimeError('Frozen seed provenance differs from numerical input')
        return
    if (not isinstance(record, dict) or record.get('kind') != 'ipm_cr_frozen_seed_v1'
            or record.get('frozenFile') != 'input_seed.mat' or seed['file'] != 'input_seed.mat'
            or record.get('sha256') != seed['sha256'] or record.get('tau') != seed['tau']):
        raise RuntimeError('Frozen seed identity differs from numerical input')
    if not (folder / 'input_seed.mat').is_file() or sha(folder / 'input_seed.mat') != seed['sha256']:
        raise RuntimeError('Frozen seed file differs from its declared SHA256')


def prepare(folder, args):
    option_file = args.case or args.options
    options = read_options(option_file) if option_file else {}
    if not isinstance(options, dict):
        raise ValueError('Numerical options must be a JSON object')
    method = options.get('method', {}) if args.case else options
    if not isinstance(method, dict):
        raise ValueError('Numerical method must be a JSON object')
    initial = options.get('initial', method.get('initial', 'data/initial.json'))
    authority = initial_input(ROOT, initial)
    inputs = source_files(ROOT) + [initial_input(ROOT), authority, relative_file(ROOT, initial),
                                   Path('tools/common.py'), Path('tools/run_explicit_experiment.py')]
    if not args.prepare_only:
        require_kernel(ROOT, method)
    folder.mkdir(parents=True, exist_ok=False)
    inventory = freeze(ROOT, folder / 'runtime', inputs)
    initial_input(folder / 'runtime', initial)
    seed = freeze_seed(folder, options, bool(args.case))
    atomic(folder / 'source_inventory.json', inventory)
    atomic(folder / 'explicit_options.json', options)
    # Options are immutable input, independent of mutable mission/progress records.
    atomic(folder / 'options_identity.json', {'sha256': sha(folder / 'explicit_options.json')})
    deadline = dt.datetime.now(dt.timezone.utc) + dt.timedelta(hours=args.hours)
    atomic(folder / 'mission.json', dict(
        kind='ipm_cr_frozen_run_v2', status='prepared', targetTau=args.target,
        createdUtc=now(), deadlineUtc=deadline.isoformat(), segmentSeconds=args.segment_seconds,
        threads=1, rssLimitBytes=2 * 1024**3, minimumFreeDiskBytes=8 * 1024**3,
        matlab=args.matlab, sourceFrozen=True,
        caseFormat=bool(args.case), initial=initial, seedInput=seed,
        timeIntegrator=method.get('timeIntegrator', 'rk4'),
        timeControl='initial_rate_CFL_formula_stage_rates_observed',
        timePartition='ordinary_original_IPC_cubic_spline_knots_fixed_four_secant',
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
        expected_options = json.loads((folder / 'options_identity.json').read_text())['sha256']
        frozen = json.loads((folder / 'source_inventory.json').read_text())
        verify_frozen(folder, mission, frozen, expected_options)
        numerical = mission['options'].get('method', {}) if mission.get('caseFormat') else mission['options']
        require_kernel(folder / 'runtime', numerical)
        initial_input(folder / 'runtime', mission.get('initial', 'data/initial.json'))
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
                verify_frozen(folder, mission, frozen, expected_options)
                index = len(mission['segments']) + 1
                script = folder / 'segment.m'
                if mission.get('caseFormat'):
                    invocation = (
                        "opts=jsondecode(fileread('" + quote(folder / 'explicit_options.json') + "'));\n"
                        "if ~isfield(opts,'run'),opts.run=struct();end\n"
                        + ("opts.run.seed.file='" + quote(folder / 'input_seed.mat') + "';\n"
                           if mission.get('seedInput') is not None else "") +
                        "opts.run.outputDir='" + quote(native) + "';opts.run.tauEnd=" + repr(mission['targetTau']) + ";\n"
                        "opts.run.wallSeconds=" + repr(min(mission['segmentSeconds'], max(0.001, (deadline-dt.datetime.now(dt.timezone.utc)).total_seconds()))) + ";\n")
                else:
                    invocation = (
                        "opts=struct('outputDir','" + quote(native) + "','tauEnd'," + repr(mission['targetTau']) + ","
                        "'wallSeconds'," + repr(mission['segmentSeconds']) + ",'checkpointEveryTau',5,'checkpointWallSeconds',60,'logEveryTau',5);\n"
                        "opts.explicitOptions=jsondecode(fileread('" + quote(folder / 'explicit_options.json') + "'));\n")
                script.write_text(
                    "restoredefaultpath;maxNumCompThreads(1);\n"
                    "addpath('" + quote(folder / 'runtime/src') + "');\n"
                    + invocation +
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
                verify_frozen(folder, mission, frozen, expected_options)
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
    parser.add_argument('--output', type=Path, help='New run folder (default: case run.outputDir or result)')
    options = parser.add_mutually_exclusive_group()
    options.add_argument('--options', type=Path, help='Legacy numerical-method JSON')
    options.add_argument('--case', type=Path, help='Portable config(\'case\') JSON')
    parser.add_argument('--matlab', help='MATLAB executable; otherwise discover it on this host')
    parser.add_argument('--target', type=float)
    parser.add_argument('--hours', type=float)
    parser.add_argument('--segment-seconds', type=float, default=600)
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    if args.prepare_only and args.action != 'start':
        parser.error('--prepare-only is only valid with start')
    case_run = {}
    if args.case:
        case = read_options(args.case)
        if not isinstance(case, dict) or case.get('format') != 'ipm_cr_case_v1':
            parser.error('--case requires an ipm_cr_case_v1 JSON object')
        case_run = case.get('run', {})
        if not isinstance(case_run, dict):
            parser.error('case.run must be an object')
    seconds = case_run.get('wallSeconds', 600 if args.case else 3000)
    if not isinstance(seconds, (int,float)) or isinstance(seconds,bool) or not math.isfinite(seconds) or seconds<=0:
        parser.error('case.run.wallSeconds must be a finite positive number')
    args.target = args.target if args.target is not None else case_run.get('tauEnd', 200)
    args.hours = args.hours if args.hours is not None else seconds / 3600
    folder = (args.output or Path(case_run.get('outputDir', ROOT / 'result'))).resolve()
    if args.action == 'start':
        if not all(isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v) and v>0
                   for v in (args.target,args.hours,args.segment_seconds)):
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
