#!/usr/bin/env python3
"""Verify SLURM resource translation and the no-nested-job boundary.

The wrapper is intentionally tested with a fake ``sbatch`` and a tiny
environment-capturing executable.  This keeps the contract runnable on a
developer workstation while checking the exact argv/resource boundary that a
real Nextflow SLURM process relies on.
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


def clean_slurm_environment() -> dict[str, str]:
    environment = os.environ.copy()
    for key in (
        "SLURM_JOB_ID",
        "SLURM_STEP_ID",
        "SLURM_CPUS_PER_TASK",
        "SLURM_MEM_PER_NODE",
        "SLURM_MEM_PER_CPU",
        "SLURM_TMPDIR",
        "SLURM_JOB_GRES",
        "CUDA_VISIBLE_DEVICES",
    ):
        environment.pop(key, None)
    return environment


def run_capture(wrapper: Path, tool: Path, capture: Path, environment: dict[str, str]) -> dict:
    environment = dict(environment)
    environment["FASTGATK_CAPTURE"] = str(capture)
    result = subprocess.run(
        [str(wrapper), str(tool), "--sentinel", "value with spaces"],
        check=True,
        text=True,
        capture_output=True,
        env=environment,
    )
    payload = json.loads(capture.read_text(encoding="utf-8"))
    assert payload["argv"] == ["--sentinel", "value with spaces"], result.stdout
    return payload


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    wrapper = root / "fastgatk-native/workflow/slurm_smoke.sh"
    if not wrapper.is_file():
        raise SystemExit(f"SLURM wrapper not found: {wrapper}")

    with tempfile.TemporaryDirectory(prefix="fastgatk-slurm-wrapper-") as directory:
        work = Path(directory)
        tool = executable(
            work / "capture-tool",
            """#!/usr/bin/env python3
import json
import os
import pathlib
import sys
payload = {
    "argv": sys.argv[1:],
    "env": {key: os.environ.get(key) for key in (
        "SLURM_JOB_ID", "SLURM_STEP_ID", "SLURM_CPUS_PER_TASK",
        "SLURM_MEM_PER_NODE", "SLURM_MEM_PER_CPU", "SLURM_TMPDIR",
        "SLURM_JOB_GRES", "CUDA_VISIBLE_DEVICES")},
}
pathlib.Path(os.environ["FASTGATK_CAPTURE"]).write_text(
    json.dumps(payload, sort_keys=True), encoding="utf-8")
""",
        )
        fake_bin = work / "bin"
        fake_bin.mkdir()
        sbatch_capture = work / "sbatch.json"
        executable(
            fake_bin / "sbatch",
            """#!/usr/bin/env python3
import json
import os
import pathlib
import subprocess
import sys
args = sys.argv[1:]
pathlib.Path(os.environ["FASTGATK_SBATCH_CAPTURE"]).write_text(
    json.dumps(args), encoding="utf-8")
wrap = None
for index, arg in enumerate(args):
    if arg == "--wrap" and index + 1 < len(args):
        wrap = args[index + 1]
    elif arg.startswith("--wrap="):
        wrap = arg.split("=", 1)[1]
if wrap is None:
    raise SystemExit("fake sbatch did not receive --wrap")
completed = subprocess.run(["bash", "-c", wrap], env=os.environ.copy())
raise SystemExit(completed.returncode)
""",
        )

        submit_env = clean_slurm_environment()
        submit_env.update(
            {
                "PATH": f"{fake_bin}{os.pathsep}{submit_env.get('PATH', '')}",
                "FASTGATK_USE_SBATCH": "1",
                "FASTGATK_SBATCH_CAPTURE": str(sbatch_capture),
                "FASTGATK_SLURM_CPUS_PER_TASK": "11",
                "FASTGATK_SLURM_MEMORY": "12G",
                "FASTGATK_SLURM_GRES": "gpu:1",
                "FASTGATK_SLURM_TMPDIR": str(work / "scratch"),
                "FASTGATK_SLURM_JOB_NAME": "unit-wrapper",
                "FASTGATK_SLURM_PARTITION": "accelerated",
                "FASTGATK_SLURM_ACCOUNT": "biology",
                "FASTGATK_SLURM_TIME": "00:10:00",
            }
        )
        submit_payload = run_capture(wrapper, tool, work / "submit-tool.json", submit_env)
        submit_args = json.loads(sbatch_capture.read_text(encoding="utf-8"))
        for expected in (
            "--wait",
            "--job-name=unit-wrapper",
            "--cpus-per-task=11",
            "--mem=12G",
            "--gres=gpu:1",
            "--partition=accelerated",
            "--account=biology",
            "--time=00:10:00",
        ):
            assert expected in submit_args, (expected, submit_args)
        assert submit_payload["env"]["SLURM_CPUS_PER_TASK"] == "11"
        assert submit_payload["env"]["SLURM_MEM_PER_NODE"] == "12G"
        assert submit_payload["env"]["SLURM_MEM_PER_CPU"] is None
        assert submit_payload["env"]["SLURM_TMPDIR"] == str(work / "scratch")
        assert submit_payload["env"]["SLURM_JOB_GRES"] == "gpu:1"

        per_cpu_env = clean_slurm_environment()
        per_cpu_capture = work / "per-cpu-tool.json"
        per_cpu_env.update(
            {
                "PATH": f"{fake_bin}{os.pathsep}{per_cpu_env.get('PATH', '')}",
                "FASTGATK_USE_SBATCH": "1",
                "FASTGATK_SBATCH_CAPTURE": str(work / "per-cpu-sbatch.json"),
                "SLURM_MEM_PER_CPU": "3G",
            }
        )
        per_cpu_payload = run_capture(wrapper, tool, per_cpu_capture, per_cpu_env)
        per_cpu_args = json.loads(
            (work / "per-cpu-sbatch.json").read_text(encoding="utf-8")
        )
        assert "--mem-per-cpu=3G" in per_cpu_args
        assert per_cpu_payload["env"]["SLURM_MEM_PER_CPU"] == "3G"

        # The GPU-count shorthand is translated only when an explicit GRES
        # string is absent.  An explicit CUDA mapping is passed through for
        # site wrappers/fake sbatch, while real SLURM may populate it itself.
        gpu_env = clean_slurm_environment()
        gpu_capture = work / "gpu-tool.json"
        gpu_sbatch_capture = work / "gpu-sbatch.json"
        gpu_env.update(
            {
                "PATH": f"{fake_bin}{os.pathsep}{gpu_env.get('PATH', '')}",
                "FASTGATK_USE_SBATCH": "1",
                "FASTGATK_SBATCH_CAPTURE": str(gpu_sbatch_capture),
                "FASTGATK_GPU_COUNT": "2",
                "FASTGATK_GPU_TYPE": "a100",
                "FASTGATK_CUDA_VISIBLE_DEVICES": "3,5",
            }
        )
        gpu_payload = run_capture(wrapper, tool, gpu_capture, gpu_env)
        gpu_args = json.loads(gpu_sbatch_capture.read_text(encoding="utf-8"))
        assert "--gpus=a100:2" in gpu_args
        assert gpu_payload["env"]["CUDA_VISIBLE_DEVICES"] == "3,5"

        invalid_gpu_env = dict(gpu_env)
        invalid_gpu_env["FASTGATK_GPU_COUNT"] = "zero"
        invalid = subprocess.run(
            [str(wrapper), str(tool), "--sentinel", "value with spaces"],
            check=False, text=True, capture_output=True, env=invalid_gpu_env,
        )
        assert invalid.returncode == 2
        assert "FASTGATK_SLURM_GPU_COUNT" in invalid.stderr

        nested_env = dict(submit_env)
        nested_env.update(
            {
                "SLURM_JOB_ID": "12345",
                "SLURM_STEP_ID": "0",
                "SLURM_CPUS_PER_TASK": "7",
                "SLURM_MEM_PER_NODE": "8192",
                "SLURM_TMPDIR": str(work / "allocated-tmp"),
                "SLURM_JOB_GRES": "gpu:2",
            }
        )
        nested_capture = work / "nested-tool.json"
        nested_payload = run_capture(wrapper, tool, nested_capture, nested_env)
        assert not (work / "nested-sbatch.json").exists()
        assert nested_payload["env"] == {
            "SLURM_JOB_ID": "12345",
            "SLURM_STEP_ID": "0",
            "SLURM_CPUS_PER_TASK": "7",
            "SLURM_MEM_PER_NODE": "8192",
            "SLURM_MEM_PER_CPU": None,
            "SLURM_TMPDIR": str(work / "allocated-tmp"),
            "SLURM_JOB_GRES": "gpu:2",
            "CUDA_VISIBLE_DEVICES": None,
        }

        local_env = clean_slurm_environment()
        local_env.update(
            {
                "FASTGATK_USE_SBATCH": "0",
                "FASTGATK_SLURM_CPUS_PER_TASK": "5",
                "FASTGATK_SLURM_MEMORY": "6G",
                "FASTGATK_SLURM_GRES": "gpu:local",
                "FASTGATK_SLURM_TMPDIR": str(work / "local-tmp"),
                "PATH": f"{fake_bin}{os.pathsep}{local_env.get('PATH', '')}",
                "FASTGATK_SBATCH_CAPTURE": str(work / "must-not-submit.json"),
            }
        )
        local_payload = run_capture(wrapper, tool, work / "local-tool.json", local_env)
        assert local_payload["env"] == {
            "SLURM_JOB_ID": "local-fastgatk-smoke",
            "SLURM_STEP_ID": None,
            "SLURM_CPUS_PER_TASK": "5",
            "SLURM_MEM_PER_NODE": "6G",
            "SLURM_MEM_PER_CPU": None,
            "SLURM_TMPDIR": str(work / "local-tmp"),
            "SLURM_JOB_GRES": "gpu:local",
            "CUDA_VISIBLE_DEVICES": None,
        }

    print(
        json.dumps(
            {
                "status": "pass",
                "submit_resource_mapping": True,
                "existing_allocation_no_nested_sbatch": True,
                "local_resource_environment": True,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
