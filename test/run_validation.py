#!/usr/bin/env python3
"""Run the bounded, single-thread ipm_cr elliptic qualification suite."""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import signal
import subprocess
import time


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def shared_dependencies(project):
    """Resolve the declared shared inputs without scanning legacy packages."""
    project = project.resolve()
    manifest_bytes = (project / 'test/SHARED_DEPENDENCIES.json').read_bytes()
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
        files.append(project / entry)
    return files, hashlib.sha256(manifest_bytes).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matlab', default='/Applications/MATLAB_R2026a.app/bin/matlab')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--timeout', type=int, default=180)
    parser.add_argument('--prepare-only', action='store_true',
                        help='Write the validation input receipt without launching MATLAB')
    args = parser.parse_args()
    folder = Path(__file__).resolve().parent
    manifest = folder / 'SHARED_DEPENDENCIES.json'
    originals, manifest_hash = shared_dependencies(folder.parent)
    stamp = dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    output = (args.output or folder / 'results' / stamp).resolve()
    output.mkdir(parents=True, exist_ok=False)
    before = {str(p): digest(p) for p in originals}
    inputs = list((folder / '+ipm_cr').rglob('*.m'))
    inputs += list((folder / 'tests').rglob('*.m')) + [folder / 'main_cr.m', Path(__file__).resolve(), manifest]
    receipt = dict(status='prepared', startedUtc=dt.datetime.now(dt.timezone.utc).isoformat(),
                   threads=1, hardWallSeconds=args.timeout,
                   inputHashes={str(p.relative_to(folder)): digest(p) for p in inputs},
                   productionSourceHashes=before, trajectoryAdvanced=False,
                   sharedDependencyManifest=dict(path='test/SHARED_DEPENDENCIES.json',
                                                 sha256=manifest_hash, fileCount=len(originals)))
    if receipt['inputHashes']['SHARED_DEPENDENCIES.json'] != manifest_hash:
        raise RuntimeError('Shared dependency manifest changed while being read')

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
    if args.prepare_only:
        print(json.dumps(dict(status='prepared', output=str(output), matlabLaunched=False)))
        return
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
