#!/usr/bin/env python3
"""Exercise scheduler-aware batch selection and fail-closed staging."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    # CTest injects the backend-specific target path.  Keep the OpenMP path as
    # the standalone default so this remains a useful local smoke command.
    binary = Path(os.environ.get(
        "FASTGATK_RESOURCE_LIMITS_BINARY",
        str(root / "fastgatk-native/build/fastgatk-hc-call"),
    ))
    if not binary.is_file():
        raise SystemExit(f"resource-limit contract binary not found: {binary}")
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    with tempfile.TemporaryDirectory(prefix="fastgatk-resource-") as directory:
        work = Path(directory)
        env = os.environ.copy()
        env.update({"SLURM_MEM_PER_NODE": "1M", "SLURM_CPUS_PER_TASK": "2",
                    "FASTGATK_DEVICE_MEMORY_BYTES": "4096",
                    "FASTGATK_DEVICE_FREE_BYTES": "2048",
                    "FASTGATK_LOCAL_SSD": "1",
                    "FASTGATK_REMOTE_INPUT": "1"})
        summary_path = work / "summary.json"
        result = subprocess.run([
            str(binary), "-I", str(bam), "-L", "17:69000-69100",
            "-O", str(work / "calls.json"), "--batch-records", "4096",
            "--threads", "8",
        ], check=True, text=True, capture_output=True, env=env)
        summary_path.write_text(result.stdout, encoding="utf-8")
        summary = json.loads(result.stdout.strip().splitlines()[-1])
        assert summary["requested_batch_records"] == 4096
        assert 0 < summary["effective_batch_records"] < 4096
        assert summary["effective_threads"] == 2
        resources = summary["resources"]
        assert resources["host_hard_bytes"] == 1024 * 1024
        # Keep launcher values auditable, but a host-only Kokkos binary must
        # not use a CUDA allocation as a device staging budget.
        assert resources["reported_device_memory_bytes"] == 4096
        assert resources["reported_device_free_bytes"] == 2048
        if resources["device_backend_compiled"]:
            assert resources["device_assignment_available"] is True
            assert resources["device_telemetry_rejected"] is False
            assert resources["device_memory_bytes"] == 4096
            assert resources["device_free_bytes"] == 2048
        else:
            assert resources["device_assignment_available"] is False
            assert resources["device_telemetry_rejected"] is True
            assert resources["device_memory_bytes"] == 0
            assert resources["device_free_bytes"] == 0
        assert resources["local_ssd"] is True
        assert resources["remote_input"] is True

        tiny_env = dict(env)
        tiny_env["SLURM_MEM_PER_NODE"] = "1K"
        rejected = subprocess.run([
            str(binary), "-I", str(bam), "-L", "17:69000-69100",
            "-O", str(work / "rejected.json"), "--batch-records", "4096",
        ], text=True, capture_output=True, env=tiny_env)
        assert rejected.returncode != 0
        assert "RESOURCE_EXHAUSTED" in rejected.stderr
    print(json.dumps({"status": "pass", "adaptive_batch": summary["effective_batch_records"],
                      "resource_error": "RESOURCE_EXHAUSTED"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
