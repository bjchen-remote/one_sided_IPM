#!/usr/bin/env python3
"""Run named ipm_cr component suites; expensive order studies are explicit."""
from pathlib import Path
import argparse
import datetime as dt
import json
import os
import subprocess
import time
from common import (atomic, initial_input, matlab_executable, now, quote, relative_file,
                    sha, snapshot, source_files, stop_child, verify_snapshot)

ROOT = Path(__file__).resolve().parent.parent

# Component reports intentionally retain their original qualification scope.
# Completing an order measurement does not assert that fourth order was seen.
SUCCESS_RULES = {
    'test_ipm_cr': 'r.passed',
    'test_runtime': 'r.passed',
    'test_ipm_cr_connection_operators': "r.passedStructuralAndPolynomialChecks&&r.sourceHashesUnchanged&&strcmp(r.status,'completed_static_connection_measurements')",
    'test_ipm_cr_background': 'r.allRequestedChecksPassed&&r.sourceHashesUnchanged',
    'test_ipm_cr_superposition': 'r.passed',
    'test_ipm_cr_tail_background': 'r.allRequestedChecksPassed&&r.sourcesUnchanged',
    'test_ipm_cr_fourth_order': "r.sourcesUnchanged&&any(strcmp(r.status,{'completed_requested_fourth_order_screens','completed_measurements_qualification_incomplete'}))",
}


def suite_plan(suite, include_order):
    elliptic = [('elliptic', 'test_ipm_cr', "struct('includeIsolatedRefinement',true)")]
    runtime = [('runtime', 'test_runtime', '')]
    components = [
        ('connection', 'test_ipm_cr_connection_operators', ''),
        ('background', 'test_ipm_cr_background', ''),
        ('superposition', 'test_ipm_cr_superposition', ''),
        ('tail', 'test_ipm_cr_tail_background', ''),
    ]
    plan = {'elliptic': elliptic, 'runtime': runtime, 'full': elliptic + components + runtime}[suite]
    if include_order:
        plan = plan + [('order', 'test_ipm_cr_fourth_order', '')]
    return plan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite', choices=('elliptic', 'runtime', 'full'), default='elliptic')
    parser.add_argument('--include-order', action='store_true', help='Also run the expensive manufactured temporal/spatial study')
    parser.add_argument('--matlab', help='MATLAB executable; otherwise discover it on this host')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--timeout', type=float, default=180)
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    if not args.timeout > 0:
        parser.error('--timeout must be positive')
    plan = suite_plan(args.suite, args.include_order)
    executable = args.matlab if args.prepare_only else matlab_executable(args.matlab)
    originals = source_files(ROOT)
    authority = initial_input(ROOT)
    helpers = [Path('tools/common.py'), Path('tools/run_validation.py')]
    tests = [relative_file(ROOT, 'tests/' + function + '.m') for _, function, _ in plan]
    # Selected tests may call other maintained test helpers (for example the
    # native build receipt reader), so record all test sources, never old trees.
    tests += [p.relative_to(ROOT) for p in (ROOT / 'tests').rglob('*.m')]
    original_hashes = snapshot(ROOT, originals + [authority])
    input_hashes = snapshot(ROOT, tests + helpers)
    stamp = dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    output = (args.output or ROOT / 'results' / ('validation_' + stamp)).resolve()
    output.mkdir(parents=True, exist_ok=False)
    lines = ["restoredefaultpath;maxNumCompThreads(1);", "addpath('" + quote(ROOT / 'src') + "','" + quote(ROOT / 'tests') + "');", "components=struct();"]
    for label, function, options in plan:
        target = output / label
        target.mkdir()
        suffix = ',' + options if options else ''
        expression = SUCCESS_RULES[function]
        lines += ["r=" + function + "('" + quote(target) + "'" + suffix + ");",
                  "entry=struct('requiredChecksSatisfied',logical(" + expression + "));",
                  "if isfield(r,'status'),entry.reportStatus=r.status;else,entry.reportStatus='returned';end"]
        if label == 'order':
            lines += ["entry.acceptanceMeaning='measurements completed; accuracy qualification is separate';",
                      "entry.fourthOrderQualified=r.globalFourthOrderQualified;"]
        else:
            lines += ["entry.acceptanceMeaning='the requested component checks passed';"]
        lines += ["components." + label + "=entry;ipm_cr.io.json('" + quote(output / 'components.json') + "',components);",
                  "assert(entry.requiredChecksSatisfied,'ipm_cr:ValidationChecks','Required checks failed: " + label + "');"]
    script = output / 'validate.m'
    script.write_text('\n'.join(lines) + '\n')
    command = [executable or 'matlab', '-singleCompThread', '-batch', "run('" + quote(script) + "')"]
    receipt = dict(kind='ipm_cr_validation_v2', status='prepared', startedUtc=now(),
                   suite=args.suite, includedOrderStudy=args.include_order,
                   selectedComponents=[label for label, _, _ in plan],
                   componentAcceptanceRules={label: SUCCESS_RULES[function] for label, function, _ in plan},
                   threads=1, hardWallSeconds=args.timeout, inputHashes=input_hashes,
                   productionSourceHashes=original_hashes, command=command,
                   validationScriptSha256=sha(script), trajectoryAccuracyCertified=False,
                   scope='Selected component and runtime checks; no unrequested resolution study or tau200 qualification.')

    def persist():
        atomic(output / 'execution.json', receipt)

    persist()
    if args.prepare_only:
        print(json.dumps(dict(status='prepared', suite=args.suite, output=str(output), matlabLaunched=False)))
        return
    start = time.monotonic()
    child = None
    try:
        with (output / 'console.log').open('xb') as log:
            child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                                     env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))
            receipt.update(status='running', pid=child.pid)
            persist()
            try:
                code = child.wait(timeout=args.timeout)
            except subprocess.TimeoutExpired:
                stop_child(child)
                raise RuntimeError('MATLAB exceeded bounded validation wall time')
        receipt['exitCode'] = code
        verify_snapshot(ROOT, original_hashes)
        verify_snapshot(ROOT, input_hashes)
        if source_files(ROOT) != originals:
            raise RuntimeError('Source file set changed while validation ran')
        receipt.update(productionSourceUnchanged=True, testInputsUnchanged=True)
        if code:
            raise RuntimeError('MATLAB failed; inspect console.log')
        receipt['components'] = json.loads((output / 'components.json').read_text())
        receipt['status'] = 'completed'
    except BaseException as error:
        receipt.update(status='failed', error=f'{type(error).__name__}: {error}')
        raise
    finally:
        stop_child(child)
        receipt.update(seconds=time.monotonic() - start, endedUtc=now())
        persist()
    print(json.dumps(dict(status=receipt['status'], suite=args.suite, output=str(output), seconds=receipt['seconds'])))


if __name__ == '__main__':
    main()
