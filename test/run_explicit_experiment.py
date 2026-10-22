#!/usr/bin/env python3
"""Frozen-source supervisor for the actual mixed-grid explicit IPM trajectory."""
from pathlib import Path, PurePosixPath
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
MATLAB = '/Applications/MATLAB_R2026a.app/bin/matlab'
AUTHORITY = Path('runs/local_validation_20261013/gaussian26_fresh_warm_01/data_frames/frame_step000000000000_tau0000000000000000.mat')
SHARED_MANIFEST = Path('test/SHARED_DEPENDENCIES.json')


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def atomic(path, value):
    temporary = path.with_name(path.name + '.writing')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')
    os.replace(temporary, path)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def shared_dependencies(project):
    """Resolve the declared shared inputs without scanning legacy packages."""
    project = project.resolve()
    manifest_bytes = (project / SHARED_MANIFEST).read_bytes()
    manifest = json.loads(manifest_bytes)
    if not isinstance(manifest, dict) or manifest.get('format') != 'ipm_cr_shared_dependencies_v1':
        raise ValueError('Unsupported shared dependency manifest format')
    entries = manifest.get('files')
    if not isinstance(entries, list) or not entries:
        raise ValueError('Shared dependencies must be a nonempty list')
    files, seen = [], set()
    for entry in entries:
        if not isinstance(entry, str) or not entry:
            raise ValueError('Shared dependency paths must be nonempty strings')
        relative = PurePosixPath(entry)
        if relative.is_absolute() or '..' in relative.parts or relative.as_posix() != entry:
            raise ValueError('Use canonical project-relative dependency paths: ' + entry)
        source = (project / entry).resolve()
        if not source.is_relative_to(project) or not source.is_file():
            raise ValueError('Shared dependency is missing or outside the project: ' + entry)
        if source in seen:
            raise ValueError('Duplicate shared dependency: ' + entry)
        seen.add(source)
        files.append(Path(entry))
    return files, hashlib.sha256(manifest_bytes).hexdigest()


def quote(path):
    return str(path).replace("'", "''")


def available_memory():
    output = subprocess.run(['vm_stat'], check=True, capture_output=True, text=True).stdout
    page = int(re.search(r'page size of (\d+)', output).group(1))
    counts = dict((k, int(v)) for k, v in re.findall(r'^(Pages [^:]+):\s+(\d+)\.', output, re.M))
    return page * sum(counts.get(k, 0) for k in ('Pages free', 'Pages inactive', 'Pages speculative'))


def prepare(folder, args):
    project = HERE.parent
    shared, manifest_hash = shared_dependencies(project)
    shared_inputs = [SHARED_MANIFEST, *shared]
    shared_hashes = {str(path): sha(project / path) for path in shared_inputs}
    if shared_hashes[str(SHARED_MANIFEST)] != manifest_hash:
        raise RuntimeError('Shared dependency manifest changed while being read')
    folder.mkdir(parents=True, exist_ok=False)
    runtime = folder / 'runtime'
    (runtime / 'test/tests').mkdir(parents=True)
    for path in shared_inputs:
        target = runtime / path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(project / path, target)
    shutil.copytree(HERE / '+ipm_cr', runtime / 'test/+ipm_cr')
    for name in ('run_ipm_cr_explicit.m',):
        source = HERE / 'tests' / name
        if source.is_file():
            shutil.copy2(source, runtime / 'test/tests' / name)
    (runtime / AUTHORITY).parent.mkdir(parents=True)
    shutil.copy2(project / AUTHORITY, runtime / AUTHORITY)
    shutil.copy2(__file__, folder / 'supervisor.py')
    for path, expected in shared_hashes.items():
        if sha(project / path) != expected or sha(runtime / path) != expected:
            raise RuntimeError('Shared input changed while preparing the runtime: ' + path)
    inventory = [dict(path=str(p.relative_to(folder)), sha256=sha(p))
                 for p in sorted(runtime.rglob('*')) if p.is_file()]
    options = json.loads(args.options.read_text()) if args.options else {}
    atomic(folder / 'explicit_options.json', options)
    inventory.extend(dict(path=p.name, sha256=sha(p))
                     for p in (folder / 'supervisor.py', folder / 'explicit_options.json'))
    atomic(folder / 'source_inventory.json', inventory)
    deadline = dt.datetime.now(dt.timezone.utc) + dt.timedelta(hours=args.hours)
    atomic(folder / 'mission.json', dict(
        kind='ipm_cr_actual_uniform_core_annular_explicit_experiment_v1', status='prepared',
        targetTau=args.target, createdUtc=now(), deadlineUtc=deadline.isoformat(),
        segmentSeconds=args.segment_seconds, threads=1, rssLimitBytes=2 * 1024**3,
        minimumFreeDiskBytes=8 * 1024**3, sourceFrozen=True,
        sharedDependencyManifest=dict(path=str(SHARED_MANIFEST),
                                      sha256=shared_hashes[str(SHARED_MANIFEST)], fileCount=len(shared)),
        timeIntegrator=options.get('timeIntegrator', 'rk4'), timeControl='all_active_both_patches_all_stages_physical_CFL',
        implicitStages=0, polarOriginPresent=False, trajectoryAccuracyCertified=False,
        options=options, segments=[]))


def dispatch(folder):
    with (folder / 'supervisor.log').open('ab') as log:
        child = subprocess.Popen([sys.executable, str(folder / 'supervisor.py'), '_supervise',
                                  '--output', str(folder)], stdout=log, stderr=subprocess.STDOUT,
                                 start_new_session=True)
    print(json.dumps(dict(output=str(folder), supervisorPid=child.pid)), flush=True)


def stop_child(child, native):
    if child and child.poll() is None:
        (native / 'STOP').touch()
        try:
            child.wait(timeout=45)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()


def supervise(folder):
    def interrupted(_signum, _frame):
        raise InterruptedError('Supervisor interrupted')
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    with (folder / 'owner.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        mission = json.loads((folder / 'mission.json').read_text())
        for item in json.loads((folder / 'source_inventory.json').read_text()):
            assert sha(folder / item['path']) == item['sha256'], item['path']
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
                index = len(mission['segments']) + 1
                script = folder / 'segment.m'
                script.write_text(
                    "restoredefaultpath;maxNumCompThreads(1);\n"
                    "addpath('" + quote(folder / 'runtime/test') + "','" + quote(folder / 'runtime/test/tests') + "','" + quote(folder / 'runtime/time_evolution') + "');\n"
                    "opts=struct('outputDir','" + quote(native) + "','tauEnd'," + repr(mission['targetTau']) + ","
                    "'wallSeconds'," + repr(mission['segmentSeconds']) + ",'checkpointEveryTau',5,'checkpointWallSeconds',60,'logEveryTau',5);\n"
                    "opts.explicitOptions=jsondecode(fileread('" + quote(folder / 'explicit_options.json') + "'));\n"
                    "summary=run_ipm_cr_explicit(opts);disp(summary.status);\n")
                entry = dict(index=index, startedUtc=now(), maximumRssBytes=0)
                mission['segments'].append(entry)
                mission.update(status='running', supervisorPid=os.getpid(), updatedUtc=now())
                atomic(folder / 'mission.json', mission)
                segment_started = time.monotonic()
                reason = None
                with (folder / ('segment_%05d.log' % index)).open('xb') as log:
                    child = subprocess.Popen([MATLAB, '-singleCompThread', '-batch', "run('" + quote(script) + "')"],
                                             stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                                             env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))
                    mission['matlabPid'] = child.pid
                    while child.poll() is None:
                        sample = subprocess.run(['ps', '-p', str(child.pid), '-o', 'rss='], capture_output=True, text=True)
                        entry['maximumRssBytes'] = max(entry['maximumRssBytes'], int(sample.stdout.strip() or 0) * 1024)
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
            mission.update(updatedUtc=now(), supervisorWallSeconds=time.monotonic() - started)
            atomic(folder / 'mission.json', mission)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('start', 'resume', 'status', '_supervise'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--options', type=Path)
    parser.add_argument('--target', type=float, default=200)
    parser.add_argument('--hours', type=float, default=50/60)
    parser.add_argument('--segment-seconds', type=float, default=600)
    parser.add_argument('--prepare-only', action='store_true',
                        help='Prepare and verify a new frozen runtime without dispatching MATLAB')
    args = parser.parse_args()
    if args.prepare_only and args.action != 'start':
        parser.error('--prepare-only is only valid with start')
    folder = args.output.resolve()
    if args.action == 'start':
        assert args.target > 0 and args.hours > 0 and args.segment_seconds > 0
        prepare(folder, args)
        if args.prepare_only:
            print(json.dumps(dict(status='prepared', output=str(folder), dispatched=False)), flush=True)
        else:
            dispatch(folder)
    elif args.action == 'resume':
        assert not (folder / 'STOP').exists() and not (folder / 'native/STOP').exists()
        dispatch(folder)
    elif args.action == '_supervise':
        supervise(folder)
    else:
        mission = json.loads((folder / 'mission.json').read_text())
        print(json.dumps(mission, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
