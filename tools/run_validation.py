#!/usr/bin/env python3
"""Run the single validation registry; research/history are explicit selections."""
from pathlib import Path
import argparse
import copy
import datetime as dt
import json
import math
import os
import re
import subprocess
import time
from common import (atomic, initial_input, matlab_executable, mex_extension, now,
                    quote, relative_file, sha, snapshot, source_files,
                    stop_child, verify_snapshot)

ROOT = Path(__file__).resolve().parent.parent
REGISTRY = ROOT / 'tools/validation.json'
CALLS = {'output', 'output_options', 'options_only', 'output_baseline', 'output_repetitions'}


def registry():
    value = json.loads(REGISTRY.read_text())
    if value.get('format') != 'ipm_cr_validation_registry_v1':
        raise ValueError('Unsupported validation registry')
    entries = {}
    for item in value['tests']:
        key = item['id']
        if not re.fullmatch(r'[a-z][a-z0-9_]*', key) or key in entries:
            raise ValueError('Invalid or duplicate validation id: ' + key)
        if (item['group'] not in value['groups'] or item['call'] not in CALLS
                or not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', item['function'])
                or not isinstance(item['options'], dict) or item['timeoutSeconds'] <= 0):
            raise ValueError('Invalid validation entry: ' + key)
        for field in ('success', 'qualification'):
            if field in item and (not isinstance(item[field], str) or '\n' in item[field]):
                raise ValueError('Invalid registry acceptance expression: ' + key)
        entries[key] = item
    for name in value['suites']:
        expand_suite(value, entries, name)
    return value, entries


def expand_suite(value, entries, name, seen=()):
    if name not in value['suites'] or name in seen:
        raise ValueError('Unknown or recursive validation suite: ' + name)
    result = []
    for key in value['suites'][name]:
        if key.startswith('@'):
            result += expand_suite(value, entries, key[1:], seen + (name,))
        elif key in entries:
            result.append(key)
        else:
            raise ValueError('Unknown validation id in suite: ' + key)
    return list(dict.fromkeys(result))


def substitute(value, inputs):
    if isinstance(value, dict):
        if set(value) == {'$input'}:
            return str(inputs[value['$input']])
        return {key: substitute(item, inputs) for key, item in value.items()}
    if isinstance(value, list):
        return [substitute(item, inputs) for item in value]
    return value


def uses_input(value):
    if isinstance(value, dict):
        return '$input' in value or any(uses_input(item) for item in value.values())
    return isinstance(value, list) and any(uses_input(item) for item in value)


def input_files(path, kind):
    """Missing optional research data is a skip; malformed existing data is an error."""
    if path is None:
        return None
    path = path.expanduser().resolve()
    if kind == 'file':
        return [path] if path.is_file() else None
    cursor = path / ('level_2/latest.json' if kind == 'study' else 'latest.json')
    if not cursor.is_file():
        return None
    pointer = json.loads(cursor.read_text())
    recorded = Path(pointer['path'])
    checkpoint = recorded if recorded.is_absolute() else cursor.parent / recorded
    if not checkpoint.is_file():
        return None
    if sha(checkpoint) != pointer['sha256']:
        raise ValueError('External checkpoint SHA mismatch: ' + str(checkpoint))
    return [cursor, checkpoint]


def build_receipt(binary, requested=None):
    if not binary.is_file():
        if requested:
            raise ValueError('Requested build evidence requires the current platform native binary')
        return None
    if requested:
        candidates = [requested.expanduser().resolve()]
    else:
        candidates = sorted(set((ROOT / 'results').glob('native_build*/native_build.json')) |
                            set((ROOT / 'results').glob('validation_*/native_build/native_build.json')), reverse=True)
    source_sha, binary_sha = sha(ROOT / 'src/+ipm_cr/+op/wenoMex.cpp'), sha(binary)
    for path in candidates:
        if path is None or not path.is_file():
            continue
        try:
            report = json.loads(path.read_text())
        except (OSError, ValueError):
            if requested:
                raise ValueError('Requested native build report is unreadable: ' + str(path))
            continue
        if not isinstance(report, dict):
            continue
        if (report.get('status') == 'completed' and report.get('sourceSha256') == source_sha
                and report.get('binarySha256') == binary_sha):
            return path.resolve()
    if requested:
        raise ValueError('Requested build report does not bind the current source and binary')
    return None


def main():
    value, entries = registry()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite', choices=tuple(value['suites']))
    parser.add_argument('--test', action='append', default=[], choices=tuple(entries), help='Select a registered test; repeat to combine')
    parser.add_argument('--list', action='store_true', help='Print the registry without MATLAB or source/native preflight')
    parser.add_argument('--include-order', action='store_true', help='Compatibility option: append the explicit order study')
    parser.add_argument('--matlab', help='MATLAB executable; otherwise discover it on this host')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--timeout', type=float, help='Hard wall time for the entire MATLAB process; default is the sum of selected entry budgets')
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--options', type=Path, help='JSON object keyed by validation id, containing per-test option overrides')
    parser.add_argument('--baseline', type=Path, help='Independent frozen baseline MAT for refactor')
    parser.add_argument('--checkpoint', type=Path, help='Accepted checkpoint folder for explicit endpoint/history checks')
    parser.add_argument('--baseline-study', type=Path, help='Old study containing level_2/latest.json')
    parser.add_argument('--build-report', type=Path, help='Existing source/binary SHA-matching native build receipt')
    parser.add_argument('--no-build-native', action='store_true', help='Require a matching existing native build instead of compiling it when needed')
    args = parser.parse_args()
    if args.list:
        print(json.dumps(value, ensure_ascii=False, indent=2))
        return
    if args.timeout is not None and (not math.isfinite(args.timeout) or args.timeout <= 0):
        parser.error('--timeout must be finite and positive')
    selected = expand_suite(value, entries, args.suite) if args.suite else []
    selected += args.test
    if not selected:
        selected = expand_suite(value, entries, 'quick')
    if args.include_order:
        selected.append('order')
    selected = list(dict.fromkeys(selected))
    overrides = json.loads(args.options.read_text()) if args.options else {}
    if not isinstance(overrides, dict) or set(overrides) - set(entries):
        parser.error('--options must be an object keyed by registered validation ids')
    if any(not isinstance(item, dict) for item in overrides.values()):
        parser.error('Each per-test options value must be an object')
    for key, options in overrides.items():
        if options and entries[key]['call'] in {'output', 'output_baseline'}:
            parser.error('This test has no options argument: ' + key)
        bound_options = {name for name, option in entries[key]['options'].items() if uses_input(option)}
        if set(options) & bound_options:
            parser.error('Use the named external input flag instead of overriding bound options for ' + key)
        if 'outputDir' in options or (key == 'native' and 'buildReport' in options):
            parser.error('Output and native build evidence are owned by the validation runner: ' + key)
    inputs = dict(baseline=args.baseline, checkpoint=args.checkpoint, baseline_study=args.baseline_study)
    external = {}
    plan, skipped = [], {}
    for key in selected:
        item = copy.deepcopy(entries[key])
        relative_file(ROOT, 'tests/' + item['function'] + '.m')
        missing = []
        for name, kind in item.get('inputs', {}).items():
            files = input_files(inputs[name], kind)
            if files is None:
                missing.append(name)
            else:
                for file in files:
                    external[str(file)] = sha(file)
                inputs[name] = inputs[name].expanduser().resolve()
        if missing:
            skipped[key] = dict(status='skipped', requiredChecksSatisfied=False,
                                qualificationSatisfied=None, missingInputs=missing,
                                reason='Required external data were not supplied or are unavailable',
                                acceptanceMeaning='Not executed; no passed or accuracy claim')
            continue
        item['options'] = substitute(item['options'], inputs)
        item['options'].update(overrides.get(key, {}))
        plan.append(item)
    needs_native = any(item['requiresNative'] for item in plan)
    binary = ROOT / ('src/+ipm_cr/+op/wenoMex' + mex_extension())
    existing_receipt = build_receipt(binary, args.build_report) if needs_native else None
    build_native = needs_native and existing_receipt is None
    if build_native and args.no_build_native:
        parser.error('A source/binary SHA-matching native build is required; build it first or omit --no-build-native')
    originals = source_files(ROOT)
    authority = initial_input(ROOT)
    original_hashes = snapshot(ROOT, originals + [authority])
    guard_hashes = dict(original_hashes)
    if build_native:
        guard_hashes.pop(binary.relative_to(ROOT).as_posix(), None)
    helpers = [Path('tools/common.py'), Path('tools/run_validation.py'), Path('tools/validation.json')]
    if build_native:
        helpers.append(Path('tools/build_ipm_cr_native.m'))
    tests = [p.relative_to(ROOT) for p in (ROOT / 'tests').rglob('*.m')]
    input_hashes = snapshot(ROOT, tests + helpers)
    if args.options:
        external[str(args.options.resolve())] = sha(args.options)
    if existing_receipt:
        external[str(existing_receipt)] = sha(existing_receipt)
    stamp = dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    output = (args.output or ROOT / 'results' / ('validation_' + stamp)).resolve()
    output.mkdir(parents=True, exist_ok=False)
    options_by_id = {}
    for item in plan:
        target = output / item['id']
        if item['call'] != 'output_repetitions':
            target.mkdir()
        options = item['options']
        if item['id'] == 'native':
            options['buildReport'] = str(output / 'native_build/native_build.json' if build_native else existing_receipt)
        if item['call'] == 'options_only':
            options['outputDir'] = str(target)
            # The endpoint inspector owns creation of its fresh directory.
            target.rmdir()
        options_by_id[item['id']] = options
    options_file = output / 'test_options.json'
    atomic(options_file, options_by_id)
    components_file = output / 'components.json'
    atomic(components_file, skipped)
    lines = ["restoredefaultpath;maxNumCompThreads(1);",
             "addpath('" + quote(ROOT / 'src') + "','" + quote(ROOT / 'tests') + "','" + quote(ROOT / 'tools') + "');",
             "components=jsondecode(fileread('" + quote(components_file) + "'));",
             "testOptions=jsondecode(fileread('" + quote(options_file) + "'));"]
    if build_native:
        lines += ["build_ipm_cr_native('" + quote(output / 'native_build') + "');"]
    for item in plan:
        key, function = item['id'], item['function']
        target = "'" + quote(output / key) + "'"
        options = 'testOptions.' + key
        calls = {'output': target, 'output_options': target + ',' + options,
                 'options_only': options,
                 'output_baseline': target + ",'" + quote(inputs.get('baseline') or '') + "'",
                 'output_repetitions': target + ',' + options + '.repetitions'}
        lines += ["entry=struct('status','running','requiredChecksSatisfied',false,'qualificationSatisfied',[]);",
                  'components.' + key + '=entry;ipm_cr.io.json(\'' + quote(components_file) + "',components);",
                  'componentTimer=tic;try',
                  'r=' + function + '(' + calls[item['call']] + ');',
                  'entry.requiredChecksSatisfied=logical(' + item['success'] + ');',
                  "if isfield(r,'status'),entry.reportStatus=r.status;else,entry.reportStatus='returned';end",
                  "entry.acceptanceMeaning='" + quote(item['meaning']) + "';"]
        if 'qualification' in item:
            lines += ['entry.qualificationSatisfied=logical(' + item['qualification'] + ');']
        lines += ["if ~entry.requiredChecksSatisfied,entry.status='failed_checks';",
                  "elseif ~isempty(entry.qualificationSatisfied)&&~entry.qualificationSatisfied,entry.status='completed_measurements_unqualified';",
                  "else,entry.status='passed_requested_checks';end",
                  'entry.seconds=toc(componentTimer);components.' + key + '=entry;',
                  "ipm_cr.io.json('" + quote(components_file) + "',components);",
                  "assert(entry.requiredChecksSatisfied,'ipm_cr:ValidationChecks','Required checks failed: " + key + "');",
                  'catch problem',
                  "entry.status='failed';entry.errorIdentifier=problem.identifier;entry.errorMessage=problem.message;entry.seconds=toc(componentTimer);",
                  'components.' + key + '=entry;',
                  "ipm_cr.io.json('" + quote(components_file) + "',components);rethrow(problem);end"]
    script = output / 'validate.m'
    script.write_text('\n'.join(lines) + '\n')
    executable = args.matlab if args.prepare_only or not plan else matlab_executable(args.matlab)
    command = [executable or 'matlab', '-singleCompThread', '-batch', "run('" + quote(script) + "')"] if plan else None
    timeout = args.timeout or sum(item['timeoutSeconds'] for item in plan) + (180 if build_native else 0)
    receipt = dict(kind='ipm_cr_validation_v3', status='prepared', startedUtc=now(),
                   suite=args.suite or ('selected' if args.test else 'quick'),
                   registrySha256=sha(REGISTRY), selectedComponents=selected,
                   componentPlan=plan, components=skipped, threads=1,
                   hardWallSeconds=timeout, timeoutPolicy='One process hard bound; default is sum of selected entry budgets plus native build allowance',
                   inputHashes=input_hashes, productionSourceHashesBefore=original_hashes,
                   guardedSourceHashes=guard_hashes, externalInputHashes=external,
                   nativeBuildPlanned=build_native, nativeBuildPerformed=False,
                   nativeBuildReceipt=str(existing_receipt) if existing_receipt else None,
                   command=command, validationScriptSha256=sha(script), testOptionsSha256=sha(options_file),
                   trajectoryAccuracyCertified=False,
                   scope='Requested assertions, measurement completion and scientific qualification remain separate; skipped inputs are not passed.')

    def persist():
        atomic(output / 'execution.json', receipt)

    persist()
    if args.prepare_only:
        print(json.dumps(dict(status='prepared', suite=receipt['suite'], output=str(output),
                              selectedComponents=selected, skippedComponents=list(skipped),
                              nativeBuildPlanned=build_native, matlabLaunched=False)))
        return
    if not plan:
        receipt.update(status='skipped', endedUtc=now(), seconds=0, matlabLaunched=False)
        persist()
        print(json.dumps(dict(status='skipped', output=str(output), skippedComponents=list(skipped), matlabLaunched=False)))
        return
    start = time.monotonic()
    child = None
    try:
        with (output / 'console.log').open('xb') as log:
            child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                                     env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'))
            receipt.update(status='running', pid=child.pid, matlabLaunched=True)
            persist()
            try:
                code = child.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                stop_child(child)
                raise RuntimeError('MATLAB exceeded bounded validation wall time')
        receipt['exitCode'] = code
        verify_snapshot(ROOT, guard_hashes)
        verify_snapshot(ROOT, input_hashes)
        if sha(script) != receipt['validationScriptSha256'] or sha(options_file) != receipt['testOptionsSha256']:
            raise RuntimeError('Generated validation script or test options changed during execution')
        for path, expected in external.items():
            if sha(Path(path)) != expected:
                raise RuntimeError('External test input changed: ' + path)
        expected_sources = set(originals) | ({binary.relative_to(ROOT)} if build_native else set())
        if set(source_files(ROOT)) != expected_sources:
            raise RuntimeError('Unexpected production source file-set change during validation')
        receipt.update(productionSourceUnchanged=True, testInputsUnchanged=True,
                       productionSourceHashesAfter=snapshot(ROOT, source_files(ROOT) + [authority]))
        if build_native:
            new_receipt = output / 'native_build/native_build.json'
            if build_receipt(binary, new_receipt) is None:
                raise RuntimeError('Native build did not produce a verified receipt')
            receipt.update(nativeBuildPerformed=True, nativeBuildReceipt=str(new_receipt),
                           nativeBuildReceiptSha256=sha(new_receipt))
        if code:
            raise RuntimeError('MATLAB failed; inspect console.log and per-component status')
        receipt['components'] = json.loads(components_file.read_text())
        for item in plan:
            entry = receipt['components'].get(item['id'], {})
            if (entry.get('requiredChecksSatisfied') is not True or entry.get('status') not in
                    {'passed_requested_checks', 'completed_measurements_unqualified'}):
                raise RuntimeError('MATLAB returned without completed required checks: ' + item['id'])
        receipt['unqualifiedComponents'] = [item['id'] for item in plan
                                           if receipt['components'][item['id']].get('qualificationSatisfied') is False]
        receipt['requestedChecksPassed'] = True
        receipt['status'] = 'completed_with_skips' if skipped else 'completed'
    except BaseException as error:
        receipt.update(status='failed', error=f'{type(error).__name__}: {error}')
        if components_file.is_file():
            receipt['components'] = json.loads(components_file.read_text())
        raise
    finally:
        stop_child(child)
        receipt.update(seconds=time.monotonic() - start, endedUtc=now())
        persist()
    print(json.dumps(dict(status=receipt['status'], suite=receipt['suite'], output=str(output), seconds=receipt['seconds'])))


if __name__ == '__main__':
    main()
