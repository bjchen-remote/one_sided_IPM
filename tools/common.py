"""Shared immutable-input and process helpers for the MATLAB tools."""
from pathlib import Path, PurePosixPath
import datetime as dt
import hashlib
import json
import os
import platform
import re
import shutil
import signal
import subprocess


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def atomic(path, value):
    temporary = path.with_name(path.name + '.writing')
    with temporary.open('w') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def relative_file(root, name):
    """Require a canonical relative file whose symlink target stays in root."""
    if not isinstance(name, str) or not name:
        raise ValueError('Input file names must be nonempty strings')
    relative = PurePosixPath(name)
    if relative.is_absolute() or '\\' in name or '..' in relative.parts or relative.as_posix() != name:
        raise ValueError('Use canonical project-relative input names: ' + name)
    source = root / name
    if not source.resolve().is_relative_to(root.resolve()) or not source.is_file():
        raise ValueError('Input is missing or outside the project: ' + name)
    return Path(name)


def initial_input(root):
    descriptor = root / 'data/initial.json'
    meta = json.loads(descriptor.read_bytes())
    if not isinstance(meta, dict) or meta.get('format') != 'ipm_cr_initial_v1':
        raise ValueError('Unsupported initial-input manifest')
    source = relative_file(root, meta.get('file'))
    if sha(root / source) != meta.get('sha256'):
        raise ValueError('Initial input SHA differs from data/initial.json')
    return source


def mex_extension():
    if platform.system() == 'Darwin':
        return '.mexmaca64' if platform.machine() == 'arm64' else '.mexmaci64'
    if platform.system() == 'Linux':
        return '.mexa64'
    raise RuntimeError('These process tools support macOS and Linux')


def source_files(root):
    package = root / 'src/+ipm_cr'
    if not package.is_dir():
        raise ValueError('Missing src/+ipm_cr package')
    suffixes = {'.m', '.c', '.cpp', mex_extension()}
    files = [relative_file(root, p.relative_to(root).as_posix())
             for p in (root / 'src').rglob('*') if p.is_file() and p.suffix in suffixes]
    files.append(relative_file(root, 'data/initial.json'))
    return sorted(set(files))


def snapshot(root, files):
    return {p.as_posix(): sha(root / p) for p in sorted(set(files))}


def verify_snapshot(root, expected):
    for relative, expected_sha in expected.items():
        if sha(root / relative) != expected_sha:
            raise RuntimeError('Input changed: ' + relative)


def freeze(root, target, files):
    """Copy one checked input set and detect edits/additions during copying."""
    files = sorted(set(files))
    before_sources = source_files(root)
    if not set(before_sources).issubset(files):
        raise RuntimeError('Source file set changed while selecting frozen inputs')
    before = snapshot(root, files)
    target.mkdir(parents=True, exist_ok=False)
    for relative in files:
        destination = target / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / relative, destination)
    verify_snapshot(root, before)
    verify_snapshot(target, before)
    if source_files(root) != before_sources:
        raise RuntimeError('Source file set changed while preparing the runtime')
    return before


def matlab_executable(requested=None):
    candidates = [requested] if requested else [shutil.which('matlab')]
    if not requested and platform.system() == 'Darwin':
        candidates += [str(p) for p in sorted(Path('/Applications').glob('MATLAB_R*.app/bin/matlab'), reverse=True)]
    for candidate in candidates:
        if not candidate:
            continue
        found = shutil.which(candidate)
        path = Path(found or candidate).expanduser()
        if path.is_file() and os.access(path, os.X_OK):
            return str(path.resolve())
    raise RuntimeError('MATLAB executable unavailable; pass --matlab /absolute/path/to/matlab')


def require_kernel(root, options):
    reconstruction = options.get('transport', {}).get('reconstruction', 'linear5')
    if reconstruction != 'legacy_weno5':
        path = root / ('src/+ipm_cr/+op/wenoMex' + mex_extension())
        if not path.is_file():
            raise RuntimeError('Build the current platform wenoMex kernel before starting a trajectory')


def quote(path):
    return str(path).replace("'", "''")


def available_memory():
    if platform.system() == 'Linux':
        for line in Path('/proc/meminfo').read_text().splitlines():
            if line.startswith('MemAvailable:'):
                return int(line.split()[1]) * 1024
        raise RuntimeError('Linux MemAvailable is unavailable')
    if platform.system() == 'Darwin':
        output = subprocess.run(['vm_stat'], check=True, capture_output=True, text=True).stdout
        page = int(re.search(r'page size of (\d+)', output).group(1))
        counts = {key: int(value) for key, value in re.findall(r'^(Pages [^:]+):\s+(\d+)\.', output, re.M)}
        return page * sum(counts.get(key, 0) for key in ('Pages free', 'Pages inactive', 'Pages speculative'))
    raise RuntimeError('Resource monitoring supports macOS and Linux')


def process_rss(pid):
    sample = subprocess.run(['ps', '-p', str(pid), '-o', 'rss='], capture_output=True, text=True)
    return int(sample.stdout.strip() or 0) * 1024


def stop_child(child, native=None):
    if child is None or child.poll() is not None:
        return
    if native is not None:
        (native / 'STOP').touch()
    try:
        child.wait(timeout=45 if native is not None else 2)
    except subprocess.TimeoutExpired:
        os.killpg(child.pid, signal.SIGTERM)
        try:
            child.wait(timeout=10)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGKILL)
            child.wait(timeout=10)
