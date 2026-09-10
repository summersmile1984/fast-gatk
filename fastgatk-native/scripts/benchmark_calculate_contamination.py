#!/usr/bin/env python3
"""Small end-to-end benchmark for the Kokkos ContaminationModel path."""

from __future__ import annotations

import json
import os
import statistics
import subprocess
import tempfile
import time
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_CONTAMINATION_BINARY",
        str(root / "fastgatk-native/build/fastgatk-calculate-contamination")))
    fixture = root / (
        "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/"
        "calculatecontamination/NA12891_0.08_NA12892_0.92.table"
    )
    if not binary.is_file() or not fixture.is_file():
        raise SystemExit("native binary and bundled fixture are required")
    samples: list[float] = []
    contamination = None
    with tempfile.TemporaryDirectory(prefix="fastgatk-contamination-benchmark-") as directory:
        work = Path(directory)
        for iteration in range(7):
            output = work / f"contamination-{iteration}.table"
            manifest = work / f"contamination-{iteration}.manifest.json"
            started = time.perf_counter()
            result = subprocess.run(
                [str(binary), "-I", str(fixture), "-O", str(output),
                 "--output-manifest", str(manifest), "--threads", "1"],
                text=True, capture_output=True, check=False,
                env={**os.environ, "OMP_PROC_BIND": "false"},
            )
            elapsed = time.perf_counter() - started
            if result.returncode != 0:
                raise SystemExit(result.stderr)
            if iteration >= 2:
                samples.append(elapsed)
            fields = output.read_text(encoding="utf-8").splitlines()[-1].split("\t")
            contamination = float(fields[1])
        telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
        assert telemetry["segmenter_kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect")
        assert telemetry["segmenter_kernel_execution_policy"] == "MDRangePolicy"
    samples.sort()
    p95_index = min(len(samples) - 1, int(0.95 * len(samples)))
    print(json.dumps({
        "schema_version": 1,
        "tool": "CalculateContamination",
        "backend": "Kokkos",
        "status": "pass",
        "fixture": fixture.name,
        "repetitions": len(samples),
        "p50_seconds": statistics.median(samples),
        "p95_seconds": samples[p95_index],
        "contamination": contamination,
        "plan_reuse": True,
        "kernel": {
            "execution_space": telemetry["kernel_execution_space"],
            "policy": telemetry["kernel_execution_policy"],
            "batches": telemetry["kernel_batches"],
            "observations": telemetry["kernel_observations"],
            "prepare_seconds": telemetry["kernel_prepare_seconds"],
            "execute_seconds": telemetry["kernel_execute_seconds"],
            "plan_allocations": telemetry.get("likelihood_plan_allocations", 0),
            "plan_reuses": telemetry.get("likelihood_plan_reuses", 0),
        },
        "segmenter_kernel": {
            "execution_space": telemetry.get("segmenter_kernel_execution_space", ""),
            "policy": telemetry.get("segmenter_kernel_execution_policy", "MDRangePolicy"),
            "batches": telemetry.get("segmenter_kernel_batches", 0),
            "observations": telemetry.get("segmenter_kernel_observations", 0),
            "prepare_seconds": telemetry.get("segmenter_kernel_prepare_seconds", 0.0),
            "execute_seconds": telemetry.get("segmenter_kernel_execute_seconds", 0.0),
        },
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
