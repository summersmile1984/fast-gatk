#!/usr/bin/env python3
"""GIAB-specific tool provenance and enforced Docker execution prerequisites."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import tarfile
import time
import uuid
import urllib.request

from oracle_toolchain import require_toolchain, repo_root

ROOT = repo_root()
GIB = 1024 ** 3
TOOLS = ('mark-duplicates', 'bqsr', 'apply-bqsr', 'hc-call', 'genotype-gvcf')
CGROUP_FILES = ('cpu.max', 'cpuset.cpus.effective', 'memory.max', 'memory.swap.max',
                'memory.peak', 'cpu.stat', 'memory.events', 'io.stat')


def save_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.part')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temporary.replace(path)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def logged(argv: list[str], work: Path, name: str, *, stdout: Path | None = None,
           timeout: float = 600, check: bool = True) -> dict:
    """Preparation commands: disk logs; wall is supervisor wall, not workload RSS."""
    work.mkdir(parents=True, exist_ok=True)
    out = stdout or work / (name + '.stdout.log')
    err = work / (name + '.stderr.log')
    out.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    start = time.monotonic()
    timed_out = False
    with out.open('wb') as o, err.open('wb') as e:
        process = subprocess.Popen(argv, stdout=o, stderr=e, start_new_session=True)
        try:
            rc = process.wait(timeout=timeout)
        except (subprocess.TimeoutExpired, KeyboardInterrupt) as failure:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                pass
            # A surviving grandchild must not outlive the command group.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            rc = process.wait()
            timed_out = isinstance(failure, subprocess.TimeoutExpired)
            if not timed_out:
                raise
    result = {'argv': argv, 'started_unix': started, 'ended_unix': time.time(),
              'supervisor_wall_seconds': time.monotonic() - start,
              'returncode': rc, 'timeout': timed_out, 'stdout': str(out), 'stderr': str(err)}
    save_json(work / (name + '.command.json'), result)
    if check and (rc or timed_out):
        raise RuntimeError(f'{name}: exit={rc}, timeout={timed_out}; stderr={err}; stdout={out}')
    return result


def log_text(result: dict, key: str = 'stdout', limit: int = 128 * 1024) -> str:
    with Path(result[key]).open('r', errors='replace') as stream:
        return stream.read(limit)


def filesystem_type(path: Path) -> str:
    actual = str(path.resolve())
    matches = []
    for line in Path('/proc/self/mountinfo').read_text().splitlines():
        left, right = line.split(' - ', 1)
        mount = re.sub(r'\\([0-7]{3})', lambda m: chr(int(m[1], 8)), left.split()[4])
        if actual == mount or actual.startswith(mount.rstrip('/') + '/'):
            matches.append((len(mount), right.split()[0]))
    if not matches:
        raise RuntimeError(f'cannot determine filesystem for {path}')
    return max(matches)[1]


def docker_limits(cpus: list[int], memory_bytes: int) -> list[str]:
    return ['--cpuset-cpus', ','.join(map(str, cpus)), '--cpus', str(len(cpus)),
            '--memory', str(memory_bytes), '--memory-swap', str(memory_bytes),
            '--network', 'none', '--cap-drop', 'ALL', '--security-opt', 'no-new-privileges']


def inspect_image(image: str, work: Path, name: str) -> dict:
    result = logged(['docker', 'image', 'inspect', image], work, name)
    info = json.loads(log_text(result))[0]
    if info['Os'] != 'linux' or info['Architecture'] != 'amd64':
        raise RuntimeError(f'{image}: expected linux/amd64')
    if '@sha256:' in image and not any(
            entry.split('@')[-1] == image.split('@')[-1] for entry in info.get('RepoDigests', [])):
        raise RuntimeError(f'{image}: local image lacks the requested repository digest')
    return {'requested': image, 'id': info['Id'], 'repo_digests': info.get('RepoDigests', []),
            'created': info['Created'], 'os': info['Os'], 'architecture': info['Architecture'],
            'inspection': result}


class PinnedTools:
    """Only asset preparation tools, mounted at identical absolute paths."""
    def __init__(self, data_root: Path, images: dict, cpus: list[int], memory_bytes: int,
                 *, mounts: list[Path] | None = None):
        self.root = data_root.resolve()
        self.images = images
        self.cpus = cpus
        self.memory_bytes = memory_bytes
        self.logs = self.root / 'preparation-logs'
        self.scratch = self.root / 'scratch' / 'preparation'
        self.scratch.mkdir(parents=True, exist_ok=True)
        self.mounts = list(dict.fromkeys([self.root, *(Path(path).resolve() for path in mounts or [])]))

    def __call__(self, tool: str, args: list[str], name: str, *, stdout: Path | None = None) -> dict:
        if tool not in ('samtools', 'bcftools'):
            raise ValueError(f'not an asset preparation tool: {tool}')
        container = 'giab-prepare-' + uuid.uuid4().hex
        command = ['docker', 'run', '--name', container, *docker_limits(self.cpus, self.memory_bytes),
                   '--user', f'{os.getuid()}:{os.getgid()}',
                   *(arg for path in self.mounts for arg in ('--mount', f'type=bind,src={path},dst={path}')),
                   '--env', f'TMPDIR={self.scratch}', '--entrypoint', tool, self.images[tool], *map(str, args)]
        try:
            return logged(command, self.logs, name, stdout=stdout, timeout=24 * 3600)
        finally:
            # Killing a docker client does not kill the workload. Explicitly remove this owned container.
            state = subprocess.run(['docker', 'inspect', container], capture_output=True, text=True)
            if state.returncode == 0:
                save_json(self.logs / (name + '.container.json'), json.loads(state.stdout)[0])
                subprocess.run(['docker', 'rm', '--force', container], stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL, check=True)


def source_fingerprint() -> dict:
    files = []
    for component in ('fastgatk-native', 'fastgatk-core', 'fastgatk-kernels', 'fastgatk-runtime'):
        directory = ROOT / component
        candidates = [directory / 'CMakeLists.txt']
        for name in ('src', 'include'):
            candidates.extend(p for p in (directory / name).rglob('*') if p.is_file())
        for path in sorted(candidates):
            files.append({'path': str(path.relative_to(ROOT)), 'sha256': sha256(path)})
    return {'files': files, 'sha256': hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()}


def build_native(data_root: Path, work: Path) -> dict:
    build = ROOT / 'fastgatk-native/build-giab-openmp'
    cmake = ROOT / 'third_party/toolchains/cmake-4.3.4-linux-x86_64/bin/cmake'
    before = source_fingerprint()
    configure = [str(cmake), '-S', str(ROOT / 'fastgatk-native'), '-B', str(build),
                 '-DCMAKE_BUILD_TYPE=Release', '-DFASTGATK_KOKKOS_BACKEND=OPENMP',
                 '-DKokkos_ENABLE_TESTS=OFF',
                 '-DFAST_GATK_HTSLIB_ROOT=' + str(ROOT / 'third_party/htslib-build/htslib-src'),
                 '-DFAST_GATK_ZLIB_ROOT=' + str(ROOT / 'third_party/htslib-build/install')]
    configured = logged(configure, work, 'native-configure', timeout=300)
    built = logged([str(cmake), '--build', str(build), '--parallel', '16', '--target',
                    *('fastgatk-' + tool for tool in TOOLS)], work, 'native-build', timeout=3600)
    after = source_fingerprint()
    if before['sha256'] != after['sha256']:
        raise RuntimeError('native source changed during build; refusing mixed provenance')
    cache = (build / 'CMakeCache.txt').read_text()
    for required in ('CMAKE_BUILD_TYPE:STRING=Release', 'FASTGATK_KOKKOS_BACKEND:STRING=OPENMP',
                     'Kokkos_ENABLE_OPENMP:BOOL=ON'):
        if required not in cache:
            raise RuntimeError('native build cache missing ' + required)
    git = {}
    for key, args in [('head', ['rev-parse', 'HEAD']), ('status', ['status', '--porcelain=v1'])]:
        result = logged(['git', '-C', str(ROOT), *args], work, 'source-' + key, check=False)
        git[key] = {'returncode': result['returncode'], 'value': log_text(result)}
    result = {'build_directory': str(build), 'source': after, 'git': git,
              'cmake_sha256': sha256(cmake), 'cache_sha256': sha256(build / 'CMakeCache.txt'),
              'cache': cache, 'configure': configured, 'build': built,
              'binaries': {tool: {'path': str(build / ('fastgatk-' + tool)),
                                  'sha256': sha256(build / ('fastgatk-' + tool))} for tool in TOOLS}}
    save_json(data_root / 'build-manifest.json', result)
    return result


def runtime_bundle(data_root: Path) -> dict:
    """Freeze native loader/libs; do not change the production image or install a JVM."""
    runtime = data_root / 'tools' / 'host-runtime'
    runtime.mkdir(parents=True, exist_ok=True)
    names = ('ld-linux-x86-64.so.2', 'libc.so.6', 'libm.so.6', 'libgcc_s.so.1',
             'libstdc++.so.6', 'libgomp.so.1')
    artifacts = {}
    for name in names:
        source = Path('/usr/lib/x86_64-linux-gnu') / name
        target = runtime / name
        digest = sha256(source)
        if not target.is_file() or sha256(target) != digest:
            shutil.copyfile(source, target)
            target.chmod(0o755)
        artifacts[name] = {'path': str(target), 'sha256': digest, 'source': str(source.resolve())}
    target = runtime / 'time'
    shutil.copyfile('/usr/bin/time', target)
    target.chmod(0o755)
    artifacts['time'] = {'path': str(target), 'sha256': sha256(target), 'source': '/usr/bin/time'}
    return {'directory': str(runtime), 'loader': str(runtime / names[0]), 'artifacts': artifacts}


def preflight_environment(data_root: Path, images: dict, *, cpus: int = 16,
                          memory_gib: int = 64, native_image: str = 'fastgatk-native:ci') -> dict:
    root = data_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    work = root / 'preflight' / ('environment-' + uuid.uuid4().hex)
    work.mkdir(parents=True)
    result = {'schema_version': 1, 'status': 'blocked', 'errors': [], 'data_root': str(root),
              'input_status': 'pending', 'accuracy_status': 'not_evaluated',
              'resource_status': 'not_measured', 'evidence': str(work)}
    def attempt(name, function):
        try:
            result[name] = function()
        except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
            result['errors'].append({'gate': name, 'error': str(error)})
        save_json(root / 'environment-manifest.json', result)
    if (cpus, memory_gib) != (16, 64):
        raise ValueError('approved baseline requires exactly --cpus 16 --memory-gib 64')
    allowed = sorted(os.sched_getaffinity(0))
    available_kb = int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines()
                            if line.startswith('MemAvailable:')))
    result['host'] = {'allowed_cpus': allowed, 'memory_available_bytes': available_kb * 1024,
                      'disk_free_bytes': shutil.disk_usage(root).free, 'filesystem': filesystem_type(root)}
    if len(allowed) < cpus or available_kb * 1024 < memory_gib * GIB:
        result['errors'].append({'gate': 'host_capacity', 'error': 'less than approved CPU/memory budget'})
    if result['host']['filesystem'] == 'tmpfs' or result['host']['disk_free_bytes'] < 600 * GIB:
        result['errors'].append({'gate': 'disk_reservation', 'error': 'requires 600 GiB free on non-tmpfs'})
    result['cpus'] = allowed[:cpus]
    result['memory_bytes'] = memory_gib * GIB
    def attest_happy_source():
        # This pinned image lost git-describe metadata during Docker COPY.
        # Recover its source layer rather than inventing a runtime version.
        url = 'https://codeload.github.com/Illumina/hap.py/tar.gz/refs/tags/v0.3.12'
        expected = '1ba8d0245b8c7d6a68ca5b78ab8d0d87104e9a3b36c86169319c5add2395b6a5'
        archive = root / 'preflight/happy-v0.3.12-source.tar.gz'
        if not archive.is_file():
            part = archive.with_suffix(archive.suffix + '.part')
            with urllib.request.urlopen(url, timeout=120) as response, part.open('wb') as output:
                shutil.copyfileobj(response, output, 1024 * 1024)
            if sha256(part) != expected:
                part.unlink()
                raise RuntimeError('official hap.py v0.3.12 source archive SHA-256 mismatch')
            part.replace(archive)
        if sha256(archive) != expected:
            raise RuntimeError('cached hap.py v0.3.12 source archive SHA-256 mismatch')
        image_tar = work / 'happy-image.tar'
        logged(['docker', 'image', 'save', images['happy'], '-o', str(image_tar)],
               work, 'happy-source-image', timeout=300)
        def hashes(tar, prefix):
            values = {}
            for member in tar:
                if member.isfile() and member.name.startswith(prefix):
                    with tar.extractfile(member) as stream:
                        hasher = hashlib.sha256()
                        for block in iter(lambda: stream.read(1024 * 1024), b''):
                            hasher.update(block)
                        values[member.name[len(prefix):]] = hasher.hexdigest()
            return values
        with tarfile.open(archive, mode='r|gz') as release:
            official = hashes(release, 'hap.py-0.3.12/')
        with tarfile.open(image_tar) as image:
            manifest = json.load(image.extractfile('manifest.json'))[0]
            config = json.load(image.extractfile(manifest['Config']))
            layers = iter(manifest['Layers'])
            source = None
            for entry in config['history']:
                layer = None if entry.get('empty_layer') else next(layers)
                command = entry.get('created_by', '')
                if 'COPY dir:' in command and 'in /opt/hap.py-source/' in command:
                    source = layer
            if source is None:
                raise RuntimeError('pinned hap.py image has no recoverable source COPY layer')
            with tarfile.open(fileobj=image.extractfile(source), mode='r|*') as copied:
                shipped = hashes(copied, 'opt/hap.py-source/')
        missing = sorted(set(official) - set(shipped))
        different = [p for p in official if p in shipped and official[p] != shipped[p]]
        extra = [p for p in set(shipped) - set(official) if not Path(p).name.startswith('.wh.')]
        proof = {'version': 'v0.3.12', 'method': 'full Docker source layer versus official release',
                 'source_url': url, 'source_sha256': expected, 'source_layer': source,
                 'files': official, 'missing': missing, 'different': different, 'additional': extra,
                 'image_history': config['history'], 'image': images['happy']}
        save_json(work / 'happy-source-attestation.json', proof)
        if not official or missing or different or extra:
            raise RuntimeError('pinned hap.py source differs from official v0.3.12 release')
        image_tar.unlink()
        return proof
    def image_gates():
        found = {}
        for name, image in images.items():
            if '@sha256:' not in image:
                raise RuntimeError(f'{name}: image must be digest pinned')
            present = subprocess.run(['docker', 'image', 'inspect', image],
                                     stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
            if present.returncode:
                logged(['docker', 'pull', image], work, name + '-pull', timeout=1200)
            found[name] = inspect_image(image, work, name + '-inspect')
        return found
    attempt('images', image_gates)
    def versions():
        commands = {'samtools': ('samtools', ['--version']), 'bcftools': ('bcftools', ['--version']),
                    'happy': ('/opt/hap.py/bin/hap.py', ['--version']),
                    'rtg': ('/opt/hap.py/libexec/rtg-tools-install/rtg', ['version'])}
        found = {}
        for name, (binary, args) in commands.items():
            image = images['happy' if name == 'rtg' else name]
            run = logged(['docker', 'run', '--rm', *docker_limits(result['cpus'], result['memory_bytes']),
                          '--entrypoint', binary, image, *args], work, name + '-version', timeout=120)
            text = log_text(run) + log_text(run, 'stderr')
            found[name] = {'output': text, 'command': run}
        save_json(work / 'actual-versions.json', found)
        for name in ('samtools', 'bcftools'):
            if not found[name]['output'].startswith(name + ' 1.21\n'):
                raise RuntimeError(f'{name}: expected actual version 1.21; see {work}/actual-versions.json')
        found['happy']['source_attestation'] = attest_happy_source()
        found['happy']['runtime_version_identifies_release'] = bool(
            re.search(r'\b0\.3\.12\b', found['happy']['output']))
        save_json(work / 'actual-versions.json', found)
        return found
    if 'images' in result:
        attempt('versions', versions)
    def oracle():
        java, gatk = require_toolchain(ROOT)
        expected_java = ROOT / 'third_party/jdk17/bin/java'
        expected_gatk = ROOT / 'third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar'
        if Path(java).resolve() != expected_java.resolve() or gatk.resolve() != expected_gatk.resolve():
            raise RuntimeError('GIAB disallows JAVA/GATK_JAR overrides outside the repository pinned toolchain')
        jrun = logged([java, '-version'], work, 'jdk-version', timeout=30)
        grun = logged([java, '-jar', str(gatk), '--version'], work, 'gatk-version', timeout=60)
        if not re.search(r'version "17\.', log_text(jrun, 'stderr')) or 'v4.6.2.0' not in log_text(grun):
            raise RuntimeError('pinned oracle actual versions do not match JDK17/GATK 4.6.2.0')
        return {'java': java, 'java_sha256': sha256(Path(java)), 'gatk': str(gatk),
                'gatk_sha256': sha256(gatk), 'jdk_release': expected_java.parents[1].joinpath('release').read_text(),
                'jdk_version': jrun, 'gatk_version': grun}
    attempt('oracle', oracle)
    def envelope():
        image = inspect_image(native_image, work, 'envelope-image')
        probe = "import os,json,pathlib; p=pathlib.Path('/sys/fs/cgroup'); print(json.dumps({'affinity':sorted(os.sched_getaffinity(0)),'cgroup':{n:(p/n).read_text() for n in " + repr(CGROUP_FILES) + "}}))"
        run = logged(['docker', 'run', '--rm', *docker_limits(result['cpus'], result['memory_bytes']),
                      '--entrypoint', 'python3', image['id'], '-c', probe], work, 'cgroup-probe')
        observed = json.loads(log_text(run))
        cg = observed['cgroup']
        if observed['affinity'] != result['cpus'] or int(cg['memory.max']) != result['memory_bytes'] or int(cg['memory.swap.max']) != 0:
            raise RuntimeError('observed cgroup memory/swap/affinity differs from requested limits')
        quota, period = map(int, cg['cpu.max'].split())
        effective = set()
        for part in cg['cpuset.cpus.effective'].strip().split(','):
            bounds = list(map(int, part.split('-')))
            effective.update(range(bounds[0], bounds[-1] + 1))
        if quota != cpus * period or effective != set(result['cpus']):
            raise RuntimeError('CPU quota/cpuset was not enforced')
        result['resource_status'] = 'enforced_within_limit'
        return {'image': image, 'observed': observed, 'command': run,
                'note': 'limit configuration observed; destructive OOM/timeout smoke is a separate gate'}
    attempt('envelope', envelope)
    attempt('build', lambda: build_native(root, work))
    attempt('runtime', lambda: runtime_bundle(root))
    if all(key in result for key in ('build', 'runtime', 'envelope')):
        def native_startup():
            rt = result['runtime']
            commands = {}
            for tool, binary in result['build']['binaries'].items():
                command = ['docker', 'run', '--rm', *docker_limits(result['cpus'], result['memory_bytes']),
                           '--mount', f'type=bind,src={ROOT},dst={ROOT},readonly',
                           '--mount', f'type=bind,src={root},dst={root}',
                           '--env', 'OMP_NUM_THREADS=16', '--env', 'OMP_PROC_BIND=false',
                           '--entrypoint', rt['loader'], result['envelope']['image']['id'],
                           '--library-path', rt['directory'], binary['path'], '--help']
                commands[tool] = logged(command, work, tool + '-startup', timeout=60)
            return commands
        attempt('native_startup', native_startup)
    result['status'] = 'eligible' if not result['errors'] else 'blocked'
    save_json(root / 'environment-manifest.json', result)
    return result
