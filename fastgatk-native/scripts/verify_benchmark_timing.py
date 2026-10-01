#!/usr/bin/env python3
"""Defend timeout cancellation of descendants, not just the timed parent."""
from __future__ import annotations

import os
from pathlib import Path
import signal
import sys
import tempfile
import time
import unittest

from benchmark_lib import ROOT, run_timed


class TimingCancellation(unittest.TestCase):
    def test_timeout_kills_sigterm_ignoring_grandchild(self):
        (ROOT / 'work').mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='timing-cancellation-', dir=ROOT / 'work') as temporary:
            work = Path(temporary)
            pid_file = work / 'grandchild.pid'
            child = ('import os,signal,sys,time; from pathlib import Path; '
                     'signal.signal(signal.SIGTERM,signal.SIG_IGN); '
                     'Path(sys.argv[1]).write_text(str(os.getpid())); time.sleep(60)')
            parent = ('import subprocess,sys,time; from pathlib import Path; '
                      'subprocess.Popen([sys.executable,"-c",sys.argv[1],sys.argv[2]]); '
                      '\nwhile not Path(sys.argv[2]).exists(): time.sleep(.005)\n'
                      'print("descendant-ready",flush=True); time.sleep(60)')
            pid = None
            try:
                result = run_timed([sys.executable, '-c', parent, child, str(pid_file)], work, 'timeout',
                                   check=False, timeout=2.0, terminate_grace_seconds=0.1)
                self.assertEqual(result['stdout'], 'descendant-ready\n')
                self.assertTrue(result['timeout'])
                pid = int(pid_file.read_text())
                stat = Path(f'/proc/{pid}/stat')
                deadline = time.monotonic() + 3
                while stat.exists():
                    try:
                        state = stat.read_text().rsplit(')', 1)[1].split()[0]
                    except FileNotFoundError:
                        break
                    if state == 'Z':
                        break
                    self.assertLess(time.monotonic(), deadline,
                                    'a live grandchild survived cancellation after its parent exited')
                    time.sleep(.01)
            finally:
                if pid is None and pid_file.exists():
                    pid = int(pid_file.read_text())
                if pid:
                    command = Path(f'/proc/{pid}/cmdline')
                    if command.exists() and str(pid_file).encode() in command.read_bytes():
                        try:
                            os.kill(pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass


if __name__ == '__main__':
    unittest.main()
