#!/usr/bin/env python3
"""One fresh Docker cgroup per GIAB workload; collect evidence before removal."""
from __future__ import annotations

from datetime import datetime, timezone
import fcntl
import json
import os
from pathlib import Path
import subprocess
import time
import uuid

from giab_assets import load_lock
from giab_environment import ROOT, docker_limits, filesystem_type, logged, save_json, sha256


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _worker_runtime(data_root: Path, environment: dict) -> dict:
    """Export the pinned native image's Python 3; leave hap.py's Python 2 intact."""
    source_image = environment['envelope']['image']['id']
    directory = data_root / 'tools' / ('worker-python-' + source_image.split(':')[-1][:16])
    directory.parent.mkdir(parents=True, exist_ok=True)
    with directory.with_suffix('.lock').open('a') as lease:
        fcntl.flock(lease, fcntl.LOCK_EX)
        manifest_path = directory / 'runtime-manifest.json'
        if not directory.exists():
            identity = uuid.uuid4().hex
            staging = directory.with_name(directory.name + '.partial-' + identity)
            logs = directory.parent / 'worker-python-preparation' / identity
            (staging / 'usr/bin').mkdir(parents=True)
            (staging / 'usr/lib').mkdir()
            container = 'giab-python-export-' + identity
            try:
                logged(['docker', 'create', '--name', container, '--entrypoint', '/bin/true', source_image],
                       logs, 'create', timeout=60)
                probe = logged(['docker', 'run', '--rm', '--network', 'none', source_image, 'python3',
                                '-c', 'import json,sys; print(json.dumps({"version":sys.version,"minor":"%s.%s"%sys.version_info[:2],"prefix":sys.prefix}))'],
                               logs, 'version', timeout=60)
                version = json.loads(Path(probe['stdout']).read_text())
                if version['prefix'] != '/usr' or not version['minor'].startswith('3.'):
                    raise RuntimeError('Pinned worker image has an unsupported Python layout')
                executable = 'usr/bin/python' + version['minor']
                library = 'usr/lib/python' + version['minor']
                for name, relative in (('executable', executable), ('stdlib', library)):
                    logged(['docker', 'cp', container + ':/' + relative, str(staging / Path(relative).parent)],
                           logs, 'copy-' + name, timeout=120)
                files = {}
                for path in sorted((staging / 'usr').rglob('*')):
                    if path.is_symlink():
                        files[str(path.relative_to(staging))] = {'symlink': os.readlink(path)}
                    elif path.is_file():
                        files[str(path.relative_to(staging))] = {'sha256': sha256(path), 'bytes': path.stat().st_size}
                manifest = {'source_image': source_image, 'version': version,
                            'executable': str(directory / executable), 'python_home': str(directory / 'usr'),
                            'files': files, 'preparation_logs': str(logs),
                            'policy': 'Pinned Python 3 supervisor only; Python 2 comparator and its libraries are unchanged'}
                save_json(staging / 'runtime-manifest.json', manifest)
                staging.rename(directory)
            finally:
                logged(['docker', 'rm', '--force', container], logs, 'remove', timeout=60, check=False)
        if not manifest_path.is_file():
            raise RuntimeError(f'Unpublished worker runtime; preserving evidence: {directory}')
        manifest = json.loads(manifest_path.read_text())
        if manifest.get('source_image') != source_image or not manifest.get('files'):
            raise RuntimeError('Worker runtime provenance does not match the pinned native image')
        for relative, expected in manifest['files'].items():
            path = directory / relative
            if 'symlink' in expected:
                valid = path.is_symlink() and os.readlink(path) == expected['symlink']
            else:
                valid = path.is_file() and not path.is_symlink() and sha256(path) == expected['sha256']
            if not valid:
                raise RuntimeError(f'Worker runtime artifact changed: {path}')
        return {**manifest, 'manifest_path': str(manifest_path)}


def execute_stage(command: list[str], *, data_root: Path, stage_dir: Path,
                  scratch: Path, environment: dict, timeout_seconds: int, name: str,
                  mounts: list[Path] | None = None, image: str | None = None,
                  extra_env: dict[str, str] | None = None) -> dict:
    data_root, stage_dir, scratch = (Path(path).resolve() for path in (data_root, stage_dir, scratch))
    stage_dir.mkdir(parents=True, exist_ok=True)
    scratch.mkdir(parents=True, exist_ok=True)
    if any((stage_dir / filename).exists() for filename in ('stage.json', 'worker-result.json', 'request.json')):
        raise FileExistsError(f'stage directory was already attempted; preserving previous evidence: {stage_dir}')
    result = {'execution_status': 'failed', 'resource_status': 'not_measured', 'name': name,
              'argv': command, 'started_at': _now(), 'returncode': None, 'signal': None,
              'timeout': False, 'metadata_path': str(stage_dir / 'stage.json'),
              'output_bytes': None, 'scratch_peak_sampled_bytes': None, 'memory_peak_bytes': None,
              'max_rss_kb': None, 'cpu_seconds': None, 'user_seconds': None, 'system_seconds': None,
              'stdout_path': str(stage_dir / 'workload.stdout.log'),
              'stderr_path': str(stage_dir / 'workload.stderr.log')}
    save_json(stage_dir / 'stage.json', result)
    cpus, memory = environment.get('cpus', []), environment.get('memory_bytes')
    if environment.get('status') != 'eligible':
        result['error'] = 'execution environment preflight is not eligible'
        save_json(stage_dir / 'stage.json', result)
        return result
    if (len(cpus), memory) != (16, 64 * 1024**3) and environment.get('experiment_kind') != 'resource_smoke':
        result['error'] = 'production configuration requires exactly 16 logical CPUs and 64 GiB'
        save_json(stage_dir / 'stage.json', result)
        return result
    if not cpus or len(set(cpus)) != len(cpus) or not isinstance(memory, int) or memory <= 0 or timeout_seconds <= 0:
        raise ValueError('invalid resource/timeout envelope')
    if filesystem_type(scratch) == 'tmpfs' or filesystem_type(data_root) == 'tmpfs':
        result['error'] = 'GIAB data/scratch cannot use tmpfs'
        save_json(stage_dir / 'stage.json', result)
        return result
    runtime = environment['runtime']
    loader = [runtime['loader'], '--library-path', runtime['directory']]
    native_paths = {entry['path'] for entry in environment['build']['binaries'].values()}
    actual = [*loader, *command] if command[0] in native_paths | {environment['oracle']['java']} else command
    worker = ['python3']
    worker_options = []
    selected_image = environment['envelope']['image']['id']
    if image is not None:
        try:
            if image != load_lock()['images']['happy']:
                raise ValueError('Only the content-pinned hap.py evaluator image is supported')
            python = _worker_runtime(data_root, environment)
            selected_image = image
            worker = [*loader, python['executable'], '-B', '-S']
            worker_options = ['--env', 'PYTHONHOME=' + python['python_home'], '--env', 'PYTHONPATH=']
            result['worker_runtime'] = python['manifest_path']
        except Exception as error:
            result['error'] = str(error)
            save_json(stage_dir / 'stage.json', result)
            return result
    env = {'OMP_NUM_THREADS': str(len(cpus)), 'OMP_PROC_BIND': 'false', 'TMPDIR': str(scratch),
           'HOME': str(scratch), 'LC_ALL': 'C.UTF-8', 'LANG': 'C.UTF-8'}
    if extra_env:
        env.update(extra_env)
    request = {'command': actual, 'original_command': command, 'stage_dir': str(stage_dir),
               'scratch': str(scratch), 'timeout_seconds': timeout_seconds, 'cpus': cpus,
               'memory_bytes': memory, 'env': env, 'time_command': [*loader, str(Path(runtime['directory']) / 'time')]}
    save_json(stage_dir / 'request.json', request)
    container = 'giab-stage-' + uuid.uuid4().hex
    # Parent output can lie outside data-root; mount both without changing argv paths.
    mount_paths = list(dict.fromkeys([data_root, stage_dir.parent, scratch,
                                     *(Path(path).resolve() for path in mounts or [])]))
    docker = ['docker', 'run', '--name', container, *docker_limits(cpus, memory),
              '--user', f'{os.getuid()}:{os.getgid()}', '--mount', f'type=bind,src={ROOT},dst={ROOT},readonly',
              *(arg for path in mount_paths for arg in ('--mount', f'type=bind,src={path},dst={path}')),
              *worker_options, '--entrypoint', worker[0], selected_image, *worker[1:],
              str(Path(__file__).with_name('giab_worker.py')), str(stage_dir / 'request.json')]
    result.update(env=env, docker_argv=docker, container_name=container)
    supervisor_start = time.monotonic()
    interruption = None
    with (stage_dir / 'supervisor.stdout.log').open('wb') as stdout, (stage_dir / 'supervisor.stderr.log').open('wb') as stderr:
        process = subprocess.Popen(docker, stdout=stdout, stderr=stderr, start_new_session=True)
        try:
            client_code = process.wait(timeout=timeout_seconds + 60)
        except (subprocess.TimeoutExpired, KeyboardInterrupt) as error:
            interruption = error
            result['timeout'] = isinstance(error, subprocess.TimeoutExpired)
            # Docker clients are not workload process groups. Kill the owned cgroup.
            killed = subprocess.run(['docker', 'kill', container], capture_output=True, text=True)
            result['supervisor_kill'] = {'returncode': killed.returncode, 'stdout': killed.stdout, 'stderr': killed.stderr}
            try:
                client_code = process.wait(timeout=15)
            except subprocess.TimeoutExpired:
                process.kill()
                client_code = process.wait()
        finally:
            inspected = subprocess.run(['docker', 'inspect', container], capture_output=True, text=True)
            if inspected.returncode == 0:
                info = json.loads(inspected.stdout)[0]
                result['docker_state'] = info['State']
                save_json(stage_dir / 'container-inspect.json', info)
            else:
                result['docker_inspect_error'] = inspected.stderr
    result['docker_client_returncode'] = client_code
    worker = stage_dir / 'worker-result.json'
    if worker.is_file():
        result.update(json.loads(worker.read_text()))
    else:
        result.update(returncode=client_code, error='worker did not publish final evidence', resource_status='not_measured')
        worker_error = stage_dir / 'worker-error.json'
        if worker_error.is_file():
            result['worker_error'] = json.loads(worker_error.read_text())
    state = result.get('docker_state', {})
    if state.get('OOMKilled'):
        result.update(execution_status='oom', resource_status='limit_hit', signal=9)
    elif interruption is not None:
        result.update(execution_status='timeout' if isinstance(interruption, subprocess.TimeoutExpired) else 'failed',
                      timeout=isinstance(interruption, subprocess.TimeoutExpired))
    elif client_code != 0 and result.get('execution_status') == 'completed':
        result.update(execution_status='failed', error='Docker/worker exit disagrees with workload success')
    if not worker.is_file() and result.get('execution_status') not in ('oom', 'timeout'):
        result['execution_status'] = 'failed'
    result['supervisor_wall_seconds'] = time.monotonic() - supervisor_start
    result.setdefault('seconds', result['supervisor_wall_seconds'])
    result['ended_at'] = _now()
    result['partial'] = result['execution_status'] != 'completed'
    save_json(stage_dir / 'stage.json', result)
    removed = subprocess.run(['docker', 'rm', '--force', container], capture_output=True, text=True)
    result['container_removal'] = {'returncode': removed.returncode, 'stdout': removed.stdout, 'stderr': removed.stderr}
    save_json(stage_dir / 'stage.json', result)
    if isinstance(interruption, KeyboardInterrupt):
        raise interruption
    return result
