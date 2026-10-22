#!/usr/bin/env python3
"""Bounded, serial MATLAB validation for the independent polar package."""
from pathlib import Path
import argparse
import fcntl
import hashlib
import json
import os
import re
import signal
import shutil
import subprocess
import time

HERE = Path(__file__).resolve().parent
MATLAB = Path('/Applications/MATLAB_R2026a.app/bin/matlab')
RUNNER_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def atomic(path, value):
    temp = path.with_suffix('.writing')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')
    os.replace(temp, path)


def quote(path):
    return str(path).replace("'", "''")


def stop_owned(child):
    if child.poll() is None:
        os.killpg(child.pid, signal.SIGTERM)
        try:
            child.wait(timeout=15)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGKILL)
            child.wait()


def interrupted(signum, _frame):
    raise SystemExit(128 + signum)


def run_matlab(root, name, body, seconds=480, rss_limit=int(2.25 * 1024**3)):
    """Own exactly one process group; keep costs even when validation fails."""
    script = root / (name + '.m')
    script.write_text("addpath('" + quote(HERE) + "');\n" + body + '\n')
    log = root / (name + '.log')
    record = dict(name=name, status='starting', script=str(script), log=str(log),
                  scriptSha256=hashlib.sha256(script.read_bytes()).hexdigest(),
                  runnerSourceSha256=RUNNER_SHA256,
                  secondsLimit=seconds, sampledRssLimit=rss_limit, samples=[])
    receipt = root / (name + '_execution.json')
    inventory = subprocess.run(['ps', '-axo', 'pid=,rss=,command='], capture_output=True, text=True, check=True)
    other = []
    for line in inventory.stdout.splitlines():
        fields = line.split(None, 2)
        if len(fields) == 3 and int(fields[1]) > 128*1024 and re.search(
                r'(^|[/\s])MATLAB(?:_mac[a-z0-9]+)?(?:\s|$)', fields[2]):
            other.append(line.strip())
    memory = subprocess.run(['vm_stat'], capture_output=True, text=True, check=True).stdout
    page_size = int(re.search(r'page size of (\d+)', memory).group(1))
    pages = dict((key, int(value)) for key, value in re.findall(r'^(Pages [^:]+):\s+(\d+)\.', memory, re.M))
    available = page_size * sum(pages.get(k, 0) for k in ('Pages free', 'Pages inactive', 'Pages speculative'))
    record['preflight'] = dict(otherComputeMatlab=other, reclaimableEstimateBytes=available,
                               vmStat=memory, noResourceWaiting=True)
    large_case = '1281' in name or name.startswith('independent')
    minimum_available = int((2.5 if large_case else 1) * 1024**3)
    if other or available < minimum_available:
        record.update(status='resource_not_admitted', requiredEstimateBytes=minimum_available)
        atomic(receipt, record)
        raise RuntimeError(f'{name}: resource preflight failed; no process started and no resource waiter created')
    atomic(receipt, record)
    started = time.monotonic()
    with log.open('wb') as output:
        child = subprocess.Popen([str(MATLAB), '-singleCompThread', '-batch',
                                  "run('" + quote(script) + "')"],
                                 stdout=output, stderr=subprocess.STDOUT,
                                 start_new_session=True)
        reason = None
        try:
            record.update(status='running', pid=child.pid)
            atomic(receipt, record)
            while child.poll() is None:
                # Sample the already running owned job; never wait for resources.
                sample = subprocess.run(['ps', '-p', str(child.pid), '-o', 'rss=,comm='],
                                        text=True, capture_output=True)
                rss = None
                line = sample.stdout.strip()
                if line:
                    rss = int(line.split(maxsplit=1)[0]) * 1024
                record['samples'].append(dict(seconds=time.monotonic()-started, rssBytes=rss))
                atomic(receipt, record)
                if time.monotonic() - started > seconds:
                    reason = 'wall_timeout'
                elif rss is not None and rss > rss_limit:
                    reason = 'sampled_rss_limit'
                if reason:
                    stop_owned(child)
                    break
                time.sleep(1)
        except BaseException:
            stop_owned(child)
            record.update(status='interrupted', returncode=child.returncode,
                          elapsedSeconds=time.monotonic()-started)
            atomic(receipt, record)
            raise
        code = child.wait()
    record.update(status=reason or ('completed' if code == 0 else 'failed'),
                  returncode=code, elapsedSeconds=time.monotonic()-started,
                  maximumSampledRssBytes=max((s['rssBytes'] or 0 for s in record['samples']), default=0))
    atomic(receipt, record)
    print(json.dumps({k: record[k] for k in ('name', 'status', 'elapsedSeconds', 'maximumSampledRssBytes')}, ensure_ascii=False), flush=True)
    if code or reason:
        raise RuntimeError(f'{name}: {record["status"]}; inspect {log}')
    return record


def continue_window(root, target, seconds):
    """Long physical validation, still bounded by the first .005 spatial window."""
    assert 1e-4 < target <= .005 and 0 < seconds <= 24*3600
    mission = json.loads((root/'mission.json').read_text())
    assert mission['status'] == 'short_common_endpoint_passed' and mission['pairPassed'] is True
    selected = Path(mission.get('acceptedPairReport', root/'pair.json'))
    pair = json.loads(selected.read_text())
    assert pair['passed'] is True and pair['canonicalTime'] == 1e-4
    comparison_source = pair['comparisonSource']
    expected_comparator = HERE/'validation_ipm_c'/'validate_pair_factored.m'
    assert Path(comparison_source['path']).resolve() == expected_comparator.resolve()
    assert hashlib.sha256(expected_comparator.read_bytes()).hexdigest() == comparison_source['sha256']
    if 'acceptedPairReportSha256' in mission:
        assert hashlib.sha256(selected.read_bytes()).hexdigest() == mission['acceptedPairReportSha256']
    milestone = root/'milestones'/'tau_0p0001'
    milestone.mkdir(parents=True, exist_ok=False)
    archive = []
    for nr in (641, 1281):
        folder = root/f'native_{nr}'
        pointer = json.loads((folder/'latest.json').read_text())
        for source, label in [(Path(pointer['path']), 'checkpoint.mat'),
                              (folder/'latest.json', 'pointer.json'),
                              (folder/'monitor.csv', 'monitor.csv'),
                              (folder/'report.json', 'report.json'),
                              (folder/'manifest.json', 'manifest.json')]:
            dest = milestone/f'native_{nr}_{label}'
            shutil.copy2(source, dest)
            archive.append(dict(path=str(dest), sha256=hashlib.sha256(dest.read_bytes()).hexdigest(), source=str(source)))
    shutil.copy2(root/'pair.json', milestone/'pair.json')
    shutil.copy2(selected, milestone/'accepted_pair.json')
    shutil.copy2(expected_comparator, milestone/'validate_pair_factored.m')
    atomic(milestone/'inventory.json', archive)
    mission.update(status='continuing_native_spatial_window', currentCommonEndpoint=target,
                   supervisorPid=os.getpid(), wallSecondsBudget=seconds,
                   longTimeQualified=False, completedSegments=[])
    atomic(root/'mission.json', mission)
    started = time.monotonic()
    deadline = started + seconds
    try:
        for nr in (641, 1281):
            folder = root/f'native_{nr}'
            while time.monotonic()-started < seconds:
                if (root/'STOP').exists():
                    mission['status'] = 'stopped'
                    return
                current = json.loads((folder/'report.json').read_text())
                if current['tau'] >= target:
                    break
                if current['status'] not in ('target_reached', 'segment_budget', 'wall_budget'):
                    raise RuntimeError(f'{nr}: previous numerical failure or stop; no automatic restart')
                index = len(mission['completedSegments'])+1
                name = f'continue_{nr}_{index:04d}'
                scalar = root/(name+'_report.json')
                remaining = deadline-time.monotonic()
                if remaining < 30:
                    mission['status'] = 'long_validation_budget_exhausted'
                    return
                soft_seconds = min(300, max(1, remaining-15))
                body = ("p=ipm_c.output.loadCheckpoint('"+quote(folder)+"');c=p.cfg;"
                        "c.targetTau="+repr(target)+";c.wallSeconds="+repr(soft_seconds)+";"
                        "r=ipm_c.solve(c,'"+quote(folder)+"');"
                        "ipm_c.output.atomicJson('"+quote(scalar)+"',r);disp(r);")
                record = run_matlab(root, name, body, seconds=min(480, remaining))
                row = json.loads(scalar.read_text())
                if row['status'] not in ('target_reached', 'segment_budget', 'wall_budget'):
                    raise RuntimeError(f'{name}: {row["status"]}; no numerical retry')
                mission['completedSegments'].append(dict(name=name, tau=row['tau'],
                                                        seconds=record['elapsedSeconds'],
                                                        rhsCallsCompleted=row['rhsCallsCompleted']))
                atomic(root/'mission.json', mission)
            if json.loads((folder/'report.json').read_text())['tau'] < target:
                mission['status'] = 'long_validation_budget_exhausted'
                return
        pair_path = root/'pair_next_window.json'
        direct_path = root/'pair_next_window_direct.json'
        assert hashlib.sha256(expected_comparator.read_bytes()).hexdigest() == comparison_source['sha256']
        remaining = deadline-time.monotonic()
        if remaining < 30:
            mission['status'] = 'long_validation_budget_exhausted'
            return
        run_matlab(root, 'independent_next_window',
                   "addpath('"+quote(HERE/'validation_ipm_c')+"');"
                   "d=ipm_c.validatePair('"+quote(root/'native_641')+"','"+quote(root/'native_1281')+"','"+quote(direct_path)+"');"
                   "r=validate_pair_factored('"+quote(root/'native_641')+"','"+quote(root/'native_1281')+"','"+quote(pair_path)+"','"+quote(direct_path)+"');disp(r.passed);", seconds=min(480, remaining))
        pair = json.loads(pair_path.read_text())
        mission.update(status='next_common_endpoint_passed' if pair['passed'] else 'spatial_comparison_failed',
                       nextPairPassed=pair['passed'], validatedThroughTau=target if pair['passed'] else 1e-4)
    except BaseException as error:
        if time.monotonic() >= deadline:
            mission.update(status='long_validation_budget_exhausted', error=repr(error))
            return
        mission.update(status='continuation_failed', error=repr(error))
        raise
    finally:
        mission['elapsedSeconds'] = time.monotonic()-started
        atomic(root/'mission.json', mission)


def main():
    signal.signal(signal.SIGTERM, interrupted)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['components', 'qualify', 'continue', 'status'])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--target', type=float, default=.00490287607550439)
    parser.add_argument('--seconds', type=float, default=8*3600)
    args = parser.parse_args()
    root = args.output.resolve()
    if args.mode == 'status':
        for file in sorted(root.glob('*_execution.json')):
            row = json.loads(file.read_text())
            print(json.dumps({k: row.get(k) for k in ('name', 'status', 'pid', 'elapsedSeconds', 'maximumSampledRssBytes')}))
        for file in sorted(root.glob('*/report.json')):
            row = json.loads(file.read_text())
            print(json.dumps(dict(path=str(file), **{k: row.get(k) for k in ('status', 'tau', 'acceptedTotal', 'fullQ')})))
        return
    if args.mode == 'continue':
        with (HERE / '.ipm_c_runner.lock').open('a+') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            continue_window(root, args.target, args.seconds)
        return
    root.mkdir(parents=True, exist_ok=False)
    with (HERE / '.ipm_c_runner.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        mission = dict(kind='ipm_c_independent_native_qualification_v1', targetTau=200,
                       currentCommonEndpoint=1e-4, status='components',
                       longTimeQualified=False, notificationAutomationCreated=False)
        atomic(root / 'mission.json', mission)
        try:
            run_matlab(root, 'components',
                       "addpath('"+quote(HERE/'tests_ipm_c')+"');r=run_components(struct('throwOnFailure',false));"
                       "ipm_c.output.atomicJson('"+quote(root/'components.json')+"',r);assert(r.passed);")
            if args.mode == 'components':
                mission['status'] = 'components_passed'
            else:
                # A native fresh t0 chain at each resolution, with the same sinh map.
                for nr, nt, dr in [(641, 513, .004979728894707017),
                                   (1281, 1025, .0024893583183658034)]:
                    folder = root / f'native_{nr}'
                    for segment in range(1, 4):
                        report_path = folder / 'report.json'
                        if report_path.exists():
                            report = json.loads(report_path.read_text())
                            if report['status'] == 'target_reached':
                                break
                            if report['status'] not in ('segment_budget', 'wall_budget'):
                                raise RuntimeError(f'{nr}: numerical failure or stop; no automatic retry')
                        body = ("c=ipm_c.config();c.grid.nr="+str(nr)+";c.grid.ntheta="+str(nt)+";"
                                "c.grid.firstSpacing="+repr(dr)+";c.targetTau=1e-4;"
                                "r=ipm_c.solve(c,'"+quote(folder)+"');disp(r);")
                        run_matlab(root, f'native_{nr}_segment{segment}', body)
                    report = json.loads((folder/'report.json').read_text())
                    if report['status'] != 'target_reached':
                        raise RuntimeError(f'{nr}: common endpoint not reached within three bounded segments')
                    before = (folder/'latest.json').read_bytes()
                    run_matlab(root, f'restart_{nr}',
                               "p=ipm_c.output.loadCheckpoint('"+quote(folder)+"');"
                               "r=ipm_c.solve(p.cfg,'"+quote(folder)+"');assert(strcmp(r.status,'target_reached'));")
                    assert before == (folder/'latest.json').read_bytes(), 'Restart changed the committed cursor'
                run_matlab(root, 'independent_pair',
                           "r=ipm_c.validatePair('"+quote(root/'native_641')+"','"+quote(root/'native_1281')+"','"+quote(root/'pair.json')+"');disp(r.passed);")
                pair = json.loads((root/'pair.json').read_text())
                factored_path = root/'pair_factored.json'
                run_matlab(root, 'independent_pair_factored',
                           "addpath('"+quote(HERE/'validation_ipm_c')+"');"
                           "r=validate_pair_factored('"+quote(root/'native_641')+"','"+quote(root/'native_1281')+"','"+quote(factored_path)+"','"+quote(root/'pair.json')+"');disp(r.passed);")
                factored = json.loads(factored_path.read_text())
                mission['status'] = 'short_common_endpoint_passed' if factored['passed'] else 'spatial_comparison_failed'
                mission['directPairPassed'] = pair['passed']
                mission['pairPassed'] = factored['passed']
                mission['acceptedPairReport'] = str(factored_path)
                mission['acceptedPairReportSha256'] = hashlib.sha256(factored_path.read_bytes()).hexdigest()
                mission['nextCommonEndpoint'] = .00490287607550439
        except BaseException as error:
            mission.update(status='failed', error=repr(error))
            raise
        finally:
            atomic(root/'mission.json', mission)


if __name__ == '__main__':
    main()
