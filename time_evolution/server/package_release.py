#!/usr/bin/env python3
"""Create a complete portable release directory without running numerical tests.

ZIP creation is opt-in. Validation evidence is copied when explicitly supplied
or available beside the source, and is never inferred from packaging success.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import zipfile
from datetime import datetime, timezone

SOURCE_ROOT = Path(__file__).resolve().parent.parent
RELEASE_NAME = 'one_sided_ipm_wall_integral_20261007'
SERVER_ROOT = '/data/user/hd58131/ipm/' + RELEASE_NAME


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_directory', nargs='?', type=Path,
                        default=SOURCE_ROOT.parent / 'releases')
    parser.add_argument('--zip', action='store_true',
                        help='Also create a ZIP; the default is a complete directory only.')
    parser.add_argument('--validation-status', type=Path,
                        help='Copy this saved validation-status JSON; no tests are run.')
    args = parser.parse_args()
    parent = args.output_directory.expanduser().resolve()
    stage = parent / RELEASE_NAME
    archive = parent / (RELEASE_NAME + '.zip')
    checksum = parent / (archive.name + '.sha256')
    protected = (stage, archive, checksum) if args.zip else (stage,)
    if any(path.exists() for path in protected):
        raise SystemExit('Refusing to overwrite an existing release artifact.')
    entry = SOURCE_ROOT / 'server/main_server.m'
    expected = "releaseRoot='" + SERVER_ROOT + "';"
    if expected not in entry.read_text().splitlines():
        raise SystemExit('The fixed server path and release name disagree.')
    status_source = args.validation_status or SOURCE_ROOT / 'VALIDATION_STATUS.json'
    if args.validation_status and not status_source.is_file():
        raise SystemExit('The supplied validation-status file does not exist.')
    validation_status = None
    if status_source.is_file():
        validation_status = json.loads(status_source.read_text())
        if not isinstance(validation_status, dict):
            raise SystemExit('Validation status must be a JSON object.')
    parent.mkdir(parents=True, exist_ok=True)
    stage.mkdir()
    # Keep the portable MATLAB fallback; platform-specific MEX is optional.
    shutil.copytree(SOURCE_ROOT / '+ipm', stage / '+ipm',
                    ignore=shutil.ignore_patterns('*.mex*', '.DS_Store', '__pycache__'))
    shutil.copytree(SOURCE_ROOT / 'tests', stage / 'tests',
                    ignore=shutil.ignore_patterns('.DS_Store', '__pycache__'))
    shutil.copy2(entry, stage / 'main.m')
    shutil.copy2(SOURCE_ROOT / 'server/launch.sh', stage / 'launch.sh')
    shutil.copy2(SOURCE_ROOT / 'README_SERVER_ZH.md', stage / 'README_SERVER_ZH.md')
    shutil.copy2(SOURCE_ROOT / 'ARCHITECTURE.md', stage / 'ARCHITECTURE.md')
    readme = (SOURCE_ROOT / 'README.md').read_text().replace(
        '[server/main_server.m](server/main_server.m)', '[main.m](main.m)')
    (stage / 'README.md').write_text(readme)
    # Evidence and runtime provenance are package data, not runtime imports.
    for name in ('validation_evidence', 'validation_controls'):
        evidence = SOURCE_ROOT / name
        if evidence.is_dir():
            shutil.copytree(evidence, stage / name)
    for name in ('TEST_FIXTURE_PROVENANCE.json', 'CANDIDATE20_MANIFEST.json',
                 'PARENT19_MANIFEST.json'):
        if (SOURCE_ROOT / name).is_file():
            shutil.copy2(SOURCE_ROOT / name, stage / name)
    runtime_rows = [
        {'path': p.relative_to(stage).as_posix(), 'sha256': digest(p)}
        for p in sorted((stage / '+ipm').rglob('*')) if p.is_file()
    ]
    runtime_text = ''.join(r['sha256'] + '  ' + r['path'] + '\n' for r in runtime_rows)
    (stage / 'PACKAGE_RUNTIME_MANIFEST.json').write_text(json.dumps({
        'kind': 'actual_portable_packaged_runtime_identity',
        'platformMexExcluded': True,
        'runtimeMatlabFileCount': sum(r['path'].endswith('.m') for r in runtime_rows),
        'runtimeTreeSha256': hashlib.sha256(runtime_text.encode()).hexdigest(),
        'rows': runtime_rows,
    }, indent=2) + '\n')
    if validation_status is not None:
        shutil.copy2(status_source, stage / 'VALIDATION_STATUS.json')
        evidence_provenance = str(status_source.resolve())
    else:
        # Unknown is distinct from "not run": tests may exist outside this
        # source directory. A maintainer must supply/update real evidence.
        validation_status = {
            'schemaVersion': 1,
            'status': 'validation_evidence_not_supplied_to_packager',
            'createdUtc': datetime.now(timezone.utc).isoformat(),
            'numericalTests': {'executed': None, 'passed': None},
            'performanceMeasurements': {'executed': None, 'speedup': None},
            'note': 'Packaging does not execute or certify tests. Supply saved '
                    'evidence with --validation-status or update this file '
                    'from actual reports before claiming validation.',
        }
        (stage / 'VALIDATION_STATUS.json').write_text(
            json.dumps(validation_status, ensure_ascii=False, indent=2) + '\n')
        evidence_provenance = 'not_supplied_unknown_test_state'
    info = '\n'.join([
        'release=' + RELEASE_NAME,
        'server_release_root=' + SERVER_ROOT,
        'created_utc=' + datetime.now(timezone.utc).isoformat(),
        'matlab_server_target=R2022b',
        'matlab_local_validation_reference=R2026a_see_validation_status',
        'integral=int_{-L}^{L}R_X(X,0)dX',
        'numerical_validation=see_VALIDATION_STATUS.json',
        'validation_evidence_source=' + evidence_provenance,
        'packaging_performs_numerical_tests=false',
        'packaging_certifies_copied_evidence=false',
        'local_execution_authorization=not_controlled_by_packager',
        'platform_mex=excluded_portable_matlab_fallback',
        'prior_release=ipm_long_time_server_20260915_v4',
        '',
    ])
    (stage / 'PACKAGE_INFO.txt').write_text(info)
    files = sorted(p for p in stage.rglob('*') if p.is_file())
    (stage / 'SHA256SUMS').write_text(''.join(
        f'{digest(p)}  {p.relative_to(stage).as_posix()}\n' for p in files))
    if args.zip:
        with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as bundle:
            for path in sorted(p for p in stage.rglob('*') if p.is_file()):
                bundle.write(path, path.relative_to(parent).as_posix())
        checksum.write_text(f'{digest(archive)}  {archive.name}\n')
    print(stage)
    if args.zip:
        print(archive)


if __name__ == '__main__':
    main()
