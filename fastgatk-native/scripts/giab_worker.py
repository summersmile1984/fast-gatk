#!/usr/bin/env python3
"""Run exactly one timed workload inside its fresh Docker resource envelope."""
from __future__ import annotations

import json
import os
from pathlib import Path
import re
import sys
import threading
import time

from benchmark_lib import run_timed
from giab_environment import CGROUP_FILES, save_json


def cgroup_snapshot() -> dict:
    root = Path('/sys/fs/cgroup')
    return {name: (root / name).read_text().strip() for name in CGROUP_FILES}


def counters(text: str) -> dict:
    return {fields[0]: int(fields[1]) for line in text.splitlines()
            if len(fields := line.split()) == 2}


def device_counters(text: str) -> dict:
    devices = {}
    for line in text.splitlines():
        fields = line.split()
        devices[fields[0]] = {key: int(value) for key, value in (item.split('=', 1) for item in fields[1:])}
    return devices


def directory_bytes(path: Path) -> int:
    size = 0
    for directory, _, files in os.walk(path):
        for name in files:
            try:
                entry = Path(directory, name)
                if not entry.is_symlink():
                    size += entry.stat().st_size
            except FileNotFoundError:
                pass
    return size


def output_bytes(command: list[str]) -> int:
    paths = set()
    for index, arg in enumerate(command[:-1]):
        if arg in ('-O', '-o', '--output', '--reports-prefix', '-M', '--metrics-file', '--telemetry', '--output-manifest'):
            output = Path(command[index + 1])
            if output.parent.is_dir():
                paths.update(p for p in output.parent.glob(output.name + '*') if p.is_file())
            if output.suffix == '.bam':
                paths.add(output.with_suffix('.bai'))
    return sum(path.stat().st_size for path in paths if path.is_file())


def run_request(request: dict) -> dict:
    stage = Path(request['stage_dir'])
    scratch = Path(request['scratch'])
    stage.mkdir(parents=True, exist_ok=True)
    scratch.mkdir(parents=True, exist_ok=True)
    environment = {key: value for key, value in os.environ.items()
                   if not key.startswith(('OMP_', 'SLURM_')) and
                   key not in ('JAVA_TOOL_OPTIONS', '_JAVA_OPTIONS', 'JDK_JAVA_OPTIONS', 'PYTHONHOME', 'PYTHONPATH')}
    environment.update(request['env'])
    before = cgroup_snapshot()
    cpus = request['cpus']
    effective = set()
    for entry in before['cpuset.cpus.effective'].split(','):
        bounds = list(map(int, entry.split('-')))
        effective.update(range(bounds[0], bounds[-1] + 1))
    quota, period = map(int, before['cpu.max'].split())
    if (sorted(os.sched_getaffinity(0)) != cpus or effective != set(cpus) or quota != len(cpus) * period
            or int(before['memory.max']) != request['memory_bytes'] or int(before['memory.swap.max']) != 0):
        raise RuntimeError('observed workload affinity/cgroup CPU/memory/swap differs from approved envelope')
    stop = threading.Event()
    scratch_peak = [directory_bytes(scratch)]
    sampling_errors = []

    def monitor():
        while not stop.wait(1.0):
            try:
                scratch_peak[0] = max(scratch_peak[0], directory_bytes(scratch))
            except OSError as error:
                sampling_errors.append(str(error))

    sampler = threading.Thread(target=monitor, daemon=True)
    window_start = time.monotonic()
    sampler.start()
    try:
        result = run_timed(request['command'], stage, 'workload', check=False, env=environment,
                           capture_output=False, timeout=request['timeout_seconds'], cwd=scratch,
                           time_command=request['time_command'])
    finally:
        stop.set()
        sampler.join()
    if not result['timing_measured']:
        for metric in ('max_rss_kb', 'filesystem_inputs', 'filesystem_outputs'):
            result[metric] = None
    scratch_peak[0] = max(scratch_peak[0], directory_bytes(scratch))
    after = cgroup_snapshot()
    window_seconds = time.monotonic() - window_start
    first_events, last_events = counters(before['memory.events']), counters(after['memory.events'])
    memory_events = {key: value - first_events.get(key, 0) for key, value in last_events.items()}
    first_cpu, last_cpu = counters(before['cpu.stat']), counters(after['cpu.stat'])
    cpu = {key: value - first_cpu.get(key, 0) for key, value in last_cpu.items()}
    before_io, after_io = device_counters(before['io.stat']), device_counters(after['io.stat'])
    io = {device: {key: value - before_io.get(device, {}).get(key, 0) for key, value in values.items()}
          for device, values in after_io.items()}
    status = 'completed' if result['returncode'] == 0 else 'failed'
    if result['timeout']:
        status = 'timeout'
    if memory_events.get('oom_kill', 0) or memory_events.get('oom_group_kill', 0):
        status = 'oom'
        result['signal'] = 9
    error_kinds = set()
    if status != 'completed':
        # Scan bounded lines, never load full workload logs into the supervisor.
        with Path(result['stderr_path']).open(errors='replace') as stream:
            while text := stream.read(65536):
                if 'RESOURCE_EXHAUSTED' in text:
                    error_kinds.add('RESOURCE_EXHAUSTED')
                if 'No space left on device' in text or 'ENOSPC' in text:
                    error_kinds.add('ENOSPC')
        if status == 'failed' and error_kinds:
            status = 'resource_exhausted'
    # GNU time itself may exit normally after its child died from a signal.
    if result['signal'] is None and Path(result['timing_path']).is_file():
        signal_match = re.search(r'Command terminated by signal (\d+)', Path(result['timing_path']).read_text())
        if signal_match:
            result['signal'] = int(signal_match[1])
    result.update(execution_status=status,
                  resource_status='limit_hit' if status in ('oom', 'resource_exhausted') or memory_events.get('max', 0) else 'enforced_within_limit',
                  argv=request['original_command'], env=environment,
                  cgroup_before=before, cgroup_after=after, memory_events_delta=memory_events,
                  cpu_stat_delta=cpu, io_stat_delta=io, memory_peak_bytes=int(after['memory.peak']),
                  cgroup_cpu_core_seconds=cpu.get('usage_usec', 0) / 1e6,
                  cgroup_measurement_window_seconds=window_seconds,
                  allocated_cpu_utilization=cpu.get('usage_usec', 0) / 1e6 / window_seconds / len(cpus),
                  scratch_peak_sampled_bytes=scratch_peak[0], scratch_sample_interval_seconds=1.0,
                  scratch_peak_semantics='sampled logical file bytes; transient between-sample peaks may be missed',
                  scratch_sampling_errors=sampling_errors,
                  output_bytes=output_bytes(request['original_command']), resource_error_kinds=sorted(error_kinds),
                  measurement_semantics={'max_rss_kb': 'GNU time workload/process-tree maximum RSS, Linux KiB',
                                         'memory_peak_bytes': 'fresh cgroup peak including descendants and charged page cache',
                                         'filesystem_inputs_outputs': 'GNU time block counters, NOT bytes',
                                         'io_stat_delta': 'per-device cgroup I/O bytes and operation counters'})
    save_json(stage / 'worker-result.json', result)
    return result


def main() -> int:
    request = json.loads(Path(sys.argv[1]).read_text())
    try:
        result = run_request(request)
    except Exception as error:
        save_json(Path(request['stage_dir']) / 'worker-error.json',
                  {'error': str(error), 'type': type(error).__name__})
        raise
    code = result['returncode']
    return code if code >= 0 else 128 - code


if __name__ == '__main__':
    raise SystemExit(main())
