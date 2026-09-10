#!/usr/bin/env python3
"""Exercise the scheduler-neutral scatter/gather failure contract.

The production SLURM executor and Nextflow are intentionally not required for
this test: a pair of tiny executables emulate a native shard and GatherVcfs.
Each fails once after writing a partial artifact, then succeeds.  The runner
must retry in a private attempt directory, publish the manifest last, and
resume without invoking either executable when the complete bundle is present.
The test also proves that retry exhaustion leaves no visible partial bundle.
"""

from __future__ import annotations

import json
import os
import stat
import subprocess
import tempfile
from pathlib import Path


def executable(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def run(
    command: list[str],
    environment: dict[str, str],
    *,
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command, env=environment, text=True, capture_output=True, check=check
    )


def manifest(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    runner = root / "fastgatk-native/workflow/scatter_gather_smoke.sh"
    slurm_wrapper = root / "fastgatk-native/workflow/slurm_smoke.sh"
    if not runner.is_file() or not slurm_wrapper.is_file():
        raise SystemExit("scatter-gather or SLURM wrapper is missing")

    with tempfile.TemporaryDirectory(prefix="fastgatk-scheduler-retry-") as directory:
        work = Path(directory)
        source = work / "reads.fixture"
        source.write_text("synthetic input\n", encoding="utf-8")
        shard_state = work / "shard-state"
        gather_state = work / "gather-state"
        invocation_log = work / "invocations.jsonl"
        sbatch_log = work / "sbatch.jsonl"

        shard = executable(
            work / "fake-shard",
            r'''#!/usr/bin/env python3
import json
import fcntl
import os
import pathlib
import sys

def value(name):
    for index, argument in enumerate(sys.argv[1:]):
        if argument == name:
            return sys.argv[index + 2]
        if argument.startswith(name + "="):
            return argument.split("=", 1)[1]
    raise SystemExit("missing " + name)

out = pathlib.Path(value("-O"))
manifest = pathlib.Path(value("--output-manifest"))
state = pathlib.Path(os.environ["FAKE_FASTGATK_STATE"])
# The production runner intentionally launches shards concurrently.  Serialize
# only this test double's shared fail-once counter so the oracle remains
# deterministic under real process scheduling (the workflow itself stays
# parallel).
state.parent.mkdir(parents=True, exist_ok=True)
with state.open("a+", encoding="utf-8") as stream:
    fcntl.flock(stream.fileno(), fcntl.LOCK_EX)
    stream.seek(0)
    count = int(stream.read() or "0")
    stream.seek(0)
    stream.truncate()
    stream.write(str(count + 1))
    stream.flush()
    os.fsync(stream.fileno())
    fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
with pathlib.Path(os.environ["FAKE_FASTGATK_LOG"]).open("a", encoding="utf-8") as log:
    log.write(json.dumps({"kind": "shard", "out": str(out), "attempt": count + 1}) + "\n")
if os.environ.get("FAKE_FASTGATK_FAIL_IF_CALLED") == "1":
    raise SystemExit("resume invoked shard")
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("partial shard\n", encoding="utf-8")
if count == 0:
    manifest.write_text("{\"schema_version\":1,\"outputs\":[{\"complete\":false}]}\n", encoding="utf-8")
    raise SystemExit(17)
out.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", encoding="utf-8")
manifest.write_text(json.dumps({
    "schema_version": 1, "tool": "fake-shard", "status": "smoke",
    "primary_output": str(out), "primary_output_kind": "vcf",
    "input": value("-I"), "interval": value("-L"),
    "outputs": [{"path": str(out), "kind": "vcf", "complete": True}],
}, sort_keys=True) + "\n", encoding="utf-8")
''',
        )
        gather = executable(
            work / "fake-gather",
            r'''#!/usr/bin/env python3
import json
import os
import pathlib
import sys

def value(name):
    for index, argument in enumerate(sys.argv[1:]):
        if argument == name:
            return sys.argv[index + 2]
        if argument.startswith(name + "="):
            return argument.split("=", 1)[1]
    raise SystemExit("missing " + name)

out = pathlib.Path(value("-O"))
manifest = pathlib.Path(value("--output-manifest"))
state = pathlib.Path(os.environ["FAKE_FASTGATK_GATHER_STATE"])
count = int(state.read_text() or "0") if state.exists() else 0
state.write_text(str(count + 1), encoding="utf-8")
with pathlib.Path(os.environ["FAKE_FASTGATK_LOG"]).open("a", encoding="utf-8") as log:
    log.write(json.dumps({"kind": "gather", "out": str(out), "attempt": count + 1}) + "\n")
if os.environ.get("FAKE_FASTGATK_FAIL_IF_CALLED") == "1":
    raise SystemExit("resume invoked gather")
out.parent.mkdir(parents=True, exist_ok=True)
out.write_bytes(b"compressed-ish gathered output\n")
out.with_name(out.name + ".tbi").write_bytes(b"tabix index\n")
if count == 0:
    manifest.write_text("{\"schema_version\":1,\"outputs\":[{\"complete\":false}]}\n", encoding="utf-8")
    raise SystemExit(19)
manifest.write_text(json.dumps({
    "schema_version": 1, "tool": "fake-gather", "status": "prototype",
    "primary_output": str(out), "primary_output_kind": "vcf",
    "outputs": [
        {"path": str(out), "kind": "vcf", "complete": True},
        {"path": str(out) + ".tbi", "kind": "vcf-index", "complete": True},
    ],
}, sort_keys=True) + "\n", encoding="utf-8")
''',
        )
        fake_bin = work / "bin"
        fake_bin.mkdir()
        sbatch = executable(
            fake_bin / "sbatch",
            r'''#!/usr/bin/env python3
import json
import os
import pathlib
import subprocess
import sys

args = sys.argv[1:]
pathlib.Path(os.environ["FAKE_SBATCH_LOG"]).open("a", encoding="utf-8").write(
    json.dumps(args) + "\n")
wrap = None
for index, argument in enumerate(args):
    if argument == "--wrap" and index + 1 < len(args):
        wrap = args[index + 1]
    elif argument.startswith("--wrap="):
        wrap = argument.split("=", 1)[1]
if wrap is None:
    raise SystemExit("fake sbatch did not receive --wrap")
environment = os.environ.copy()
environment["SLURM_JOB_ID"] = "fake-job"
environment["SLURM_STEP_ID"] = "0"
raise SystemExit(subprocess.run(["bash", "-c", wrap], env=environment).returncode)
''',
        )

        output = work / "results"
        env = os.environ.copy()
        env.update(
            {
                "FAKE_FASTGATK_STATE": str(shard_state),
                "FAKE_FASTGATK_GATHER_STATE": str(gather_state),
                "FAKE_FASTGATK_LOG": str(invocation_log),
                "FAKE_SBATCH_LOG": str(sbatch_log),
                "FASTGATK_SCATTER_RETRIES": "1",
                "FASTGATK_USE_SBATCH": "1",
                "PATH": str(fake_bin) + os.pathsep + os.environ.get("PATH", ""),
            }
        )
        command = [
            str(runner),
            "-I",
            str(source),
            "-L",
            "chr1:1-2,chr1:3-4",
            "--binary",
            str(shard),
            "--gather-binary",
            str(gather),
            "--outdir",
            str(output),
            "--threads",
            "1",
            "--retries",
            "1",
            "--slurm-wrapper",
            str(slurm_wrapper),
        ]

        first = run(command, env)
        summary = json.loads(first.stdout.strip().splitlines()[-1])
        assert summary["status"] == "pass", first.stdout
        assert summary["resume"] == {"resumed_shards": 0, "computed_shards": 2}
        assert summary["retry"]["retried_shards"] == 1, summary
        assert summary["retry"]["gather_retries"] == 1
        assert summary["retry"]["gather_state"] == "computed"
        assert summary["compatibility"]["atomic_attempt_publish"] is True
        assert summary["compatibility"]["failure_retry"] is True
        sbatch_invocations = [
            json.loads(line) for line in sbatch_log.read_text(encoding="utf-8").splitlines()
        ]
        assert len(sbatch_invocations) == 5, sbatch_invocations
        assert all("--wait" in args and any(arg.startswith("--cpus-per-task=") for arg in args)
                   for args in sbatch_invocations)
        assert (output / "gathered.vcf.gz").is_file()
        assert (output / "gathered.vcf.gz.tbi").is_file()
        assert (output / "gathered.vcf.gz.manifest.json").is_file()
        assert manifest(output / "gathered.vcf.gz.manifest.json")["compatibility"][
            "atomic_attempt_publish"
        ] is True
        assert not list(output.glob(".fastgatk-*.attempt.*"))
        assert not list(output.glob("*.tmp"))

        before = invocation_log.read_text(encoding="utf-8")
        env["FAKE_FASTGATK_FAIL_IF_CALLED"] = "1"
        resumed = run(command, env)
        resumed_summary = json.loads(resumed.stdout.strip().splitlines()[-1])
        assert resumed_summary["resume"] == {
            "resumed_shards": 2,
            "computed_shards": 0,
        }
        assert resumed_summary["retry"]["gather_state"] == "resumed"
        assert invocation_log.read_text(encoding="utf-8") == before

        # Retry exhaustion is fail-closed.  A fresh output directory and a
        # permanently failing shard must not expose calls-*.vcf, a manifest,
        # or attempt-directory residue to a downstream gather.
        failure_output = work / "failure-results"
        failure_env = dict(env)
        failure_env.pop("FAKE_FASTGATK_FAIL_IF_CALLED", None)
        failure_env["FAKE_FASTGATK_STATE"] = str(work / "always-fail-shard-state")
        always_fail = executable(
            work / "always-fail-shard",
            r'''#!/usr/bin/env python3
import pathlib
import sys
out = pathlib.Path(sys.argv[sys.argv.index("-O") + 1])
manifest = pathlib.Path(sys.argv[sys.argv.index("--output-manifest") + 1])
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("partial\n", encoding="utf-8")
manifest.write_text("{\"schema_version\":1,\"outputs\":[{\"complete\":false}]}\n", encoding="utf-8")
raise SystemExit(23)
''',
        )
        failure_command = command.copy()
        failure_command[failure_command.index("--binary") + 1] = str(always_fail)
        failure_command[failure_command.index("--outdir") + 1] = str(failure_output)
        failure_command[failure_command.index("--retries") + 1] = "1"
        failed = run(failure_command, failure_env, check=False)
        assert failed.returncode != 0
        if failure_output.exists():
            assert not list(failure_output.glob("calls-*.vcf"))
            assert not list(failure_output.glob(".fastgatk-*.attempt.*"))
        assert not (failure_output / "scatter-gather.manifest.json").exists()

    print(
        json.dumps(
            {
                "status": "pass",
                "retry_after_partial_failure": True,
                "gather_retry_after_partial_failure": True,
                "resume_without_reexecution": True,
                "retry_exhaustion_fail_closed": True,
                "real_slurm_multi_node": False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
