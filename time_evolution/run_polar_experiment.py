#!/usr/bin/env python3
"""Frozen-source, checkpointed polar experiment: start, status, or resume."""
from pathlib import Path
import argparse
import datetime as dt
import fcntl
import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
MATLAB = Path('/Applications/MATLAB_R2026a.app/bin/matlab')


def atomic(path, value):
    temporary = path.with_suffix('.writing')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')
    os.replace(temporary, path)


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def quote(value):
    return str(value).replace("'", "''")


def interrupted(signum, _frame):
    raise InterruptedError('Supervisor received signal %d' % signum)


def stop_owned(child, native):
    if child is None or child.poll() is not None:
        return
    (native / 'STOP').touch()
    try:
        child.wait(timeout=30)
    except subprocess.TimeoutExpired:
        os.killpg(child.pid, signal.SIGTERM)
        try:
            child.wait(timeout=15)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGKILL)
            child.wait()


def memory():
    result = subprocess.run(['vm_stat'], capture_output=True, text=True, check=True).stdout
    page = int(re.search(r'page size of (\d+)', result).group(1))
    pages = {k: int(v) for k, v in re.findall(r'^(Pages [^:]+):\s+(\d+)\.', result, re.M)}
    return page * sum(pages.get(k, 0) for k in ('Pages free', 'Pages inactive', 'Pages speculative'))


def prepare(root, args):
    root.mkdir(parents=True, exist_ok=False)
    runtime = root / 'runtime'
    runtime.mkdir()
    for name in ('+ipm', '+ipm_c'):
        shutil.copytree(HERE / name, runtime / name)
    shutil.copy2(Path(__file__), runtime / Path(__file__).name)
    # An immutable copy keeps the original wall authority available independently.
    authority = HERE.parent / 'runs/local_validation_20261013/gaussian26_fresh_warm_01/data_frames/frame_step000000000000_tau0000000000000000.mat'
    shutil.copy2(authority, root / 'authority.mat')
    inventory = [dict(path=str(p.relative_to(root)), sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                 for p in sorted(runtime.rglob('*')) if p.is_file()]
    inventory.append(dict(path='authority.mat', sha256=hashlib.sha256((root / 'authority.mat').read_bytes()).hexdigest()))
    atomic(root / 'source_inventory.json', inventory)
    mission = dict(kind='ipm_c_frozen_direct_cfl_experiment_v1', status='prepared',
                   createdUtc=now(), nr=args.nr, ntheta=args.ntheta, targetTau=args.target,
                   fieldEveryTau=args.field_every, firstSpacing=args.first_spacing,
                   maxSteps=2000000, segmentSeconds=300,
                   wallBudgetHours=args.hours, rssLimitBytes=int(1.5 * 1024**3),
                   minimumFreeDiskBytes=10 * 1024**3, sourceFrozen=True,
                   spatialTrajectoryQualified=False, longTimeQualified=False,
                   timeIntegrator='ssprk3', timeControl='direct CFL', stepRetryUsed=False,
                   automaticRemesh=True, poissonDecompositionCache=True, segments=[])
    mission['deadlineUtc'] = (dt.datetime.now(dt.timezone.utc) + dt.timedelta(hours=args.hours)).isoformat()
    atomic(root / 'mission.json', mission)


def dispatch(root):
    with (root / 'supervisor.log').open('ab') as log:
        child = subprocess.Popen([sys.executable, str(root / 'runtime/run_polar_experiment.py'),
                                  '_supervise', '--output', str(root)],
                                 stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
    print(json.dumps(dict(supervisorPid=child.pid, output=str(root)), ensure_ascii=False), flush=True)


def supervise(root):
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    mission = json.loads((root / 'mission.json').read_text())
    with (root / 'supervisor.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for item in json.loads((root / 'source_inventory.json').read_text()):
            assert hashlib.sha256((root / item['path']).read_bytes()).hexdigest() == item['sha256'], item['path']
        native = root / 'native'
        native.mkdir(exist_ok=True)
        mission.update(status='starting', supervisorPid=os.getpid(), lastStartedUtc=now())
        atomic(root / 'mission.json', mission)
        started = time.monotonic()
        # A resume keeps the original wall-time deadline rather than granting
        # the same budget again after every diagnostic pause.
        deadline = dt.datetime.fromisoformat(mission.get('deadlineUtc',
            (dt.datetime.now(dt.timezone.utc) + dt.timedelta(hours=mission['wallBudgetHours'])).isoformat()))
        child = None
        try:
            while dt.datetime.now(dt.timezone.utc) < deadline:
                if (root / 'STOP').exists() or (native / 'STOP').exists():
                    mission['status'] = 'stopped'
                    break
                if shutil.disk_usage(root).free < mission['minimumFreeDiskBytes']:
                    mission['status'] = 'disk_budget_reached'
                    break
                if memory() < 2 * 1024**3:
                    mission.update(status='waiting_for_memory', updatedUtc=now())
                    atomic(root / 'mission.json', mission)
                    time.sleep(30)
                    continue
                index = max((entry['index'] for entry in mission['segments']), default=0) + 1
                script = root / 'segment.m'
                script.write_text(
                    "addpath('" + quote(root / 'runtime') + "');\n"
                    "folder='" + quote(native) + "';\n"
                    "if isfile(fullfile(folder,'latest.json')),cp=ipm_c.output.loadCheckpoint(folder);cfg=cp.cfg;\n"
                    "else,cfg=ipm_c.config();cfg.grid.nr=" + str(mission['nr']) + ";cfg.grid.ntheta=" + str(mission['ntheta']) + ";\n"
                    "cfg.authorityFile='" + quote(root / 'authority.mat') + "';cfg.runtime.threads=1;\n"
                    "cfg.grid.firstSpacing=" + repr(mission.get('firstSpacing',.03)) + ";\n"
                    "cfg.remesh.minimumFirstSpacing=cfg.grid.firstSpacing;\n"
                    "cfg.time.outputEvery=" + repr(mission['fieldEveryTau']) + ";end\n"
                    "cfg.targetTau=" + repr(mission['targetTau']) + ";cfg.maxSteps=" + str(mission['maxSteps']) + ";\n"
                    "cfg.wallSeconds=" + str(mission['segmentSeconds']) + ";\n"
                    "report=ipm_c.solve(cfg,folder);disp(report);\n")
                logpath = root / ('segment_%05d.log' % index)
                segment = dict(index=index, status='reserved', startedUtc=now(), log=logpath.name, maximumRssBytes=0)
                mission['segments'].append(segment)
                mission.update(status='running', updatedUtc=now(), activeSegment=index)
                atomic(root / 'mission.json', mission)
                with logpath.open('xb') as log:
                    child = subprocess.Popen([str(MATLAB), '-singleCompThread', '-batch',
                                              "run('" + quote(script) + "')"],
                                             stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                    segment['status'] = 'running'
                    mission['matlabPid'] = child.pid
                    atomic(root / 'mission.json', mission)
                    begin = time.monotonic()
                    reason = None
                    while child.poll() is None:
                        sample = subprocess.run(['ps', '-p', str(child.pid), '-o', 'rss='], capture_output=True, text=True)
                        rss = int(sample.stdout.strip() or 0) * 1024
                        segment['maximumRssBytes'] = max(segment['maximumRssBytes'], rss)
                        if (root / 'STOP').exists():
                            (native / 'STOP').touch()
                        if rss > mission['rssLimitBytes'] or memory() < 512 * 1024**2:
                            reason = 'memory_budget_reached'
                        elif dt.datetime.now(dt.timezone.utc) >= deadline:
                            reason = 'wall_budget_reached'
                        elif shutil.disk_usage(root).free < mission['minimumFreeDiskBytes']:
                            reason = 'disk_budget_reached'
                        elif time.monotonic() - begin > mission['segmentSeconds'] + 180:
                            reason = 'segment_timeout'
                        if reason:
                            stop_owned(child, native)
                            break
                        time.sleep(5)
                    segment.update(status='finished', returncode=child.wait(), elapsedSeconds=time.monotonic() - begin, endedUtc=now())
                mission.pop('matlabPid', None)
                if reason:
                    mission['status'] = reason
                    break
                if segment['returncode']:
                    mission['status'] = 'numerical_or_runtime_failure'
                    break
                report = json.loads((native / 'report.json').read_text())
                mission.update(tau=report['tau'], physicalTime=report['physicalTime'], steps=report['steps'],
                               remeshCount=report['remeshCount'], updatedUtc=now())
                if report['status'] == 'target_reached':
                    mission['status'] = 'target_reached'
                    break
                if report['status'] != 'segment_budget' or report['steps'] >= mission['maxSteps']:
                    mission['status'] = report['status']
                    break
                atomic(root / 'mission.json', mission)
                # A segment resumes only its last accepted checkpoint, never a failed step.
            else:
                mission['status'] = 'wall_budget_reached'
        except BaseException as error:
            mission.update(status='interrupted' if isinstance(error, (InterruptedError, KeyboardInterrupt))
                           else 'supervisor_failure', error=repr(error))
            raise
        finally:
            stop_owned(child, native)
            mission.update(updatedUtc=now(), supervisorElapsedSeconds=time.monotonic() - started)
            atomic(root / 'mission.json', mission)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('start', 'status', 'resume', '_supervise'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--nr', type=int, default=385)
    parser.add_argument('--ntheta', type=int, default=128)
    parser.add_argument('--target', type=float, default=200)
    parser.add_argument('--field-every', type=float, default=.5)
    parser.add_argument('--first-spacing', type=float, default=.03)
    parser.add_argument('--hours', type=float, default=48)
    args = parser.parse_args()
    root = args.output.resolve()
    if args.action == 'start':
        assert args.nr >= 65 and args.ntheta >= 32 and args.target > 0 and args.field_every > 0 and args.hours > 0
        assert args.first_spacing > 0
        assert shutil.disk_usage(root.parent).free >= 10 * 1024**3
        prepare(root, args)
        dispatch(root)
    elif args.action == 'resume':
        assert not (root / 'STOP').exists() and not (root / 'native/STOP').exists(), 'Remove STOP explicitly before resuming.'
        dispatch(root)
    elif args.action == '_supervise':
        supervise(root)
    else:
        result = dict(mission=json.loads((root / 'mission.json').read_text()))
        if (root / 'native/report.json').exists():
            result['latestReport'] = json.loads((root / 'native/report.json').read_text())
        print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
