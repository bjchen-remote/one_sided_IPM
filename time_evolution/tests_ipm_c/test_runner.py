"""Actual disposable processes exercise watchdog cleanup, without MATLAB."""
from pathlib import Path
import importlib.util
import os
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

SOURCE = Path(__file__).resolve().parents[1] / 'run_ipm_c.py'
spec = importlib.util.spec_from_file_location('ipm_c_runner_under_test', SOURCE)
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class RunnerLifecycle(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.fake = self.root/'fake_matlab'
        self.fake.write_text(f'#!{sys.executable}\nimport os,time\n'
                             f'open({str(self.root/"child.pid")!r},"w").write(str(os.getpid()))\n'
                             'time.sleep(30)\n')
        self.fake.chmod(0o700)
        self.real_run = subprocess.run

    def tearDown(self):
        self.temp.cleanup()

    def preflight(self, args, **kwargs):
        if args[:2] == ['ps', '-axo']:
            return subprocess.CompletedProcess(args, 0, '', '')
        if args == ['vm_stat']:
            return subprocess.CompletedProcess(args, 0,
                'Mach Virtual Memory Statistics: (page size of 16384 bytes)\nPages free: 300000.\n', '')
        return self.real_run(args, **kwargs)

    def assert_gone(self, pid):
        with self.assertRaises(ProcessLookupError):
            os.kill(pid, 0)

    def test_wall_deadline_retires_real_owned_group(self):
        with patch.object(runner, 'MATLAB', self.fake), patch.object(runner.subprocess, 'run', self.preflight):
            with self.assertRaises(RuntimeError):
                runner.run_matlab(self.root, 'timeout', '', seconds=.05)
        pid = int((self.root/'child.pid').read_text())
        self.assert_gone(pid)

    def test_receipt_failure_after_spawn_retires_child(self):
        atomic = runner.atomic
        calls = 0
        pids = []
        popen = subprocess.Popen

        def capture(*args, **kwargs):
            child = popen(*args, **kwargs)
            pids.append(child.pid)
            return child

        def fail_once(path, value):
            nonlocal calls
            calls += 1
            if calls == 2:
                raise OSError('injected receipt failure after spawn')
            atomic(path, value)

        with patch.object(runner, 'MATLAB', self.fake), patch.object(runner.subprocess, 'run', self.preflight), \
             patch.object(runner.subprocess, 'Popen', capture), patch.object(runner, 'atomic', fail_once):
            with self.assertRaises(OSError):
                runner.run_matlab(self.root, 'receipt_fault', '')
        self.assertEqual(len(pids), 1)
        self.assert_gone(pids[0])

    def test_sigterm_of_supervisor_retires_separate_session(self):
        launcher = self.root/'launcher.py'
        launcher.write_text(
            'import importlib.util,signal,subprocess\nfrom pathlib import Path\n'
            f's=importlib.util.spec_from_file_location("r",{str(SOURCE)!r});r=importlib.util.module_from_spec(s);s.loader.exec_module(r)\n'
            f'r.MATLAB=Path({str(self.fake)!r})\nreal=subprocess.run\n'
            'def probe(a,**k):\n'
            ' if a[:2]==["ps","-axo"]:return subprocess.CompletedProcess(a,0,"","")\n'
            ' if a==["vm_stat"]:return subprocess.CompletedProcess(a,0,"Mach Virtual Memory Statistics: (page size of 16384 bytes)\\nPages free: 300000.\\n","")\n'
            ' return real(a,**k)\n'
            'r.subprocess.run=probe\nsignal.signal(signal.SIGTERM,r.interrupted)\n'
            f'r.run_matlab(Path({str(self.root)!r}),"sigterm","")\n')
        child = subprocess.Popen([sys.executable, str(launcher)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        try:
            deadline = time.monotonic()+10
            while not (self.root/'child.pid').exists() and time.monotonic()<deadline:
                if child.poll() is not None:
                    self.fail('Supervisor exited before disposable child was ready')
                time.sleep(.02)
            pid = int((self.root/'child.pid').read_text())
            child.terminate()
            child.wait(timeout=10)
            self.assert_gone(pid)
        finally:
            if child.poll() is None:
                child.kill()
                child.wait()


if __name__ == '__main__':
    unittest.main()
