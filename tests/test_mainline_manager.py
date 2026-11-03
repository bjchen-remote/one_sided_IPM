"""Portable frozen-input/seed regressions; no MATLAB or PDE execution."""
from pathlib import Path
from types import SimpleNamespace
import copy
import json
import shutil
import signal
import sys
import tempfile
import unittest
from unittest.mock import patch

TOOLS = Path(__file__).resolve().parents[1] / 'tools'
sys.path.insert(0, str(TOOLS))
import common
import run_explicit_experiment as manager


class MainlineManagerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)
        self.root = self.base / 'project'
        for name in ('src/+ipm_cr', 'data', 'tools'):
            (self.root / name).mkdir(parents=True)
        (self.root / 'src/+ipm_cr/init.m').write_text('function init; end\n')
        self.kernel = self.root / ('src/+ipm_cr/+op/wenoMex' + common.mex_extension())
        self.kernel.parent.mkdir(parents=True)
        self.kernel.write_bytes(b'opaque platform-kernel fixture; never executed')
        (self.root / 'main.m').write_text('function main; end\n')
        initial = self.root / 'data/initial.mat'
        initial.write_bytes(b'opaque authority fixture; no MATLAB data are loaded')
        common.atomic(self.root / 'data/initial.json', dict(
            format='ipm_cr_initial_v1', file='data/initial.mat', sha256=common.sha(initial)))
        common.atomic(self.root / 'data/mesh.json', dict(format='geometry fixture'))
        for name in ('common.py', 'run_explicit_experiment.py'):
            shutil.copy2(TOOLS / name, self.root / 'tools' / name)
        self.root_patch = patch.object(manager, 'ROOT', self.root)
        self.root_patch.start()
        self.addCleanup(self.root_patch.stop)
        self.addCleanup(self.temporary.cleanup)

    def prepare(self, seed=None, name='run'):
        options = dict(format='ipm_cr_case_v1', run=dict(tauEnd=400))
        if seed is not None:
            options['run']['seed'] = seed
        case = self.base / (name + '.json')
        common.atomic(case, options)
        args = SimpleNamespace(case=case, options=None, prepare_only=False,
                               hours=1, target=400, segment_seconds=30, matlab=None)
        folder = self.base / name
        manager.prepare(folder, args)
        return folder

    def verify(self, folder, mission=None):
        if mission is None:
            mission = manager.read_options(folder / 'mission.json')
        frozen = manager.read_options(folder / 'source_inventory.json')
        expected = manager.read_options(folder / 'options_identity.json')['sha256']
        manager.verify_frozen(folder, mission, frozen, expected)
        return mission

    def original_seed(self):
        path = self.base / "original τ390 seed's data.mat"
        path.write_bytes(b'accepted checkpoint fixture at actual tau390')
        return path, dict(file=str(path), sha256=common.sha(path), tau=390)

    def test_linear_default_and_weno_modes_need_a_native_kernel(self):
        self.kernel.unlink()
        for mode in ('central4', 'central6', 'legacy_weno5'):
            common.require_kernel(self.root, dict(transport=dict(reconstruction=mode)))
        for mode in (None, 'linear5', 'smooth_weno5'):
            options = {} if mode is None else dict(transport=dict(reconstruction=mode))
            with self.assertRaises(RuntimeError):
                common.require_kernel(self.root, options)
        self.kernel.write_bytes(b'opaque platform-kernel fixture; never executed')
        common.require_kernel(self.root, {})

    def test_mesh_is_frozen_and_new_or_changed_sources_are_rejected(self):
        folder = self.prepare()
        self.verify(folder)
        frozen = manager.read_options(folder / 'source_inventory.json')
        self.assertIn('data/mesh.json', frozen)
        self.assertIn('main.m', frozen)
        extra = folder / 'runtime/src/+ipm_cr/unrecorded.m'
        extra.write_text('function unrecorded; end\n')
        with self.assertRaises(RuntimeError):
            self.verify(folder)
        extra.unlink()
        mesh = folder / 'runtime/data/mesh.json'
        mesh.write_text('{}')
        with self.assertRaises(RuntimeError):
            self.verify(folder)

    def test_seed_copy_survives_original_removal_and_directory_move(self):
        original, seed = self.original_seed()
        folder = self.prepare(seed)
        mission = self.verify(folder)
        frozen_seed = mission['options']['run']['seed']
        self.assertEqual(frozen_seed, dict(file='input_seed.mat', sha256=seed['sha256'], tau=390))
        self.assertEqual(mission['seedInput']['originalFile'], str(original.resolve()))
        self.assertEqual(common.sha(folder / 'input_seed.mat'), seed['sha256'])
        original.unlink()
        moved = self.base / 'moved run'
        shutil.move(folder, moved)
        self.verify(moved)

    def test_changed_seed_and_relabelled_tau_are_rejected(self):
        _, seed = self.original_seed()
        folder = self.prepare(seed)
        mission = self.verify(folder)
        changed = copy.deepcopy(mission)
        changed['options']['run']['seed']['tau'] = 400
        changed['seedInput']['tau'] = 400
        with self.assertRaises(RuntimeError):
            self.verify(folder, changed)
        (folder / 'input_seed.mat').write_bytes(b'changed checkpoint')
        with self.assertRaises(RuntimeError):
            self.verify(folder)

    def test_seed_copy_detects_source_race(self):
        original, seed = self.original_seed()
        copy2 = shutil.copy2

        def racing_copy(source, target, *args, **kwargs):
            result = copy2(source, target, *args, **kwargs)
            if Path(source).resolve() == original.resolve():
                original.write_bytes(b'source changed during copy')
            return result

        with patch.object(manager.shutil, 'copy2', racing_copy):
            with self.assertRaises(RuntimeError):
                self.prepare(seed)

    def test_seed_schema_and_sha_are_checked_before_use(self):
        _, seed = self.original_seed()
        for field, value in (('tau', -1), ('tau', True), ('tau', float('inf')),
                             ('sha256', 'f' * 63), ('file', '')):
            bad = dict(seed)
            bad[field] = value
            with self.assertRaises(ValueError):
                manager.seed_input(dict(run=dict(seed=bad)), True)
        bad = dict(seed)
        bad['sha256'] = '0' * 64
        with self.assertRaises(RuntimeError):
            self.prepare(bad)

    def test_stopped_resume_checks_frozen_inputs_without_launch(self):
        _, seed = self.original_seed()
        folder = self.prepare(seed)
        before = self.verify(folder)
        (folder / 'STOP').touch()
        previous = {sig: signal.getsignal(sig) for sig in (signal.SIGINT, signal.SIGTERM)}
        try:
            with patch.object(manager.sys, 'platform', 'test'), \
                    patch.object(manager.subprocess, 'Popen', side_effect=AssertionError('No process expected')):
                manager.supervise(folder, 'unused-matlab')
        finally:
            for sig, handler in previous.items():
                signal.signal(sig, handler)
        after = self.verify(folder)
        self.assertEqual(after['status'], 'stopped')
        self.assertEqual(after['deadlineUtc'], before['deadlineUtc'])
        self.assertEqual(after['segments'], [])
        (folder / 'input_seed.mat').write_bytes(b'corrupt after accepted stop')
        try:
            with self.assertRaises(RuntimeError):
                manager.supervise(folder, 'unused-matlab')
        finally:
            for sig, handler in previous.items():
                signal.signal(sig, handler)


if __name__ == '__main__':
    unittest.main()
