#!/usr/bin/env python3
"""Contract checks for the GATK CallCopyRatioSegments arithmetic path."""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = pathlib.Path(os.environ.get(
    "FASTGATK_CALL_COPY_RATIO_BINARY",
    str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build"))
        / "fastgatk-call-copy-ratio-segments"),
))


def write_segments(path: pathlib.Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:GATKCopyNumber\tSM:SAMPLE\n")
        handle.write("CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO\n")
        values = [0.0, -0.15200309344504997, 0.13750352374993502, 0.5, 1.5, -1.5]
        for index, value in enumerate(values):
            start = index * 100 + 1
            handle.write(f"chr1\t{start}\t{start + 99}\t10\t{value}\n")


def output_calls(path: pathlib.Path) -> list[str]:
    for line in path.read_text(encoding="utf-8").splitlines():
        if line and not line.startswith("@") and not line.startswith("CONTIG"):
            return [field for field in line.split("\t")]
    raise AssertionError("missing called segment rows")


def main() -> None:
    if not BINARY.exists():
        raise SystemExit(f"missing native binary: {BINARY}")
    with tempfile.TemporaryDirectory(prefix="fastgatk-call-copy-ratio-") as temporary:
        work = pathlib.Path(temporary)
        input_path = work / "segments.tsv"
        output_path = work / "called.tsv"
        legacy_path = work / "called.seg"
        manifest_path = work / "called.json"
        write_segments(input_path)
        result = subprocess.run(
            [str(BINARY), "-I", str(input_path), "-O", str(output_path),
             "--legacy-output", str(legacy_path), "--output-manifest", str(manifest_path)],
            text=True, capture_output=True, check=True,
        )
        rows = [line.split("\t") for line in output_path.read_text(encoding="utf-8").splitlines()
                if line and not line.startswith("@") and not line.startswith("CONTIG")]
        assert [row[5] for row in rows] == ["0", "0", "0", "+", "+", "-"], rows
        assert legacy_path.read_text(encoding="utf-8").splitlines()[0].startswith("Sample\tChromosome")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        assert manifest["tool"] == "CallCopyRatioSegments"
        assert manifest["determinism"] == "strict"
        assert manifest["execution_space"] in {"OpenMP", "Serial"}
        assert manifest["telemetry"]["copy_ratio_kernel_execution_space"] in {"OpenMP", "Serial"}
        assert manifest["telemetry"]["copy_ratio_kernel_simd_width"] >= 1
        assert manifest["telemetry"]["copy_ratio_kernel_simd_groups"] >= 1
        assert manifest["telemetry"]["copy_ratio_kernel_execute_seconds"] >= 0.0
        assert manifest["telemetry"]["kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect"
        )
        assert manifest["telemetry"]["kernel_execution_policy"] == "RangePolicy"
        assert manifest["telemetry"]["kernel_batches"] == 1
        assert manifest["telemetry"]["kernel_observations"] == 6
        assert manifest["telemetry"]["kernel_prepare_seconds"] >= 0.0
        assert manifest["segments"] == 6
        assert manifest["telemetry"]["output_bytes"] == output_path.stat().st_size
        assert manifest["telemetry"]["legacy_output_bytes"] == legacy_path.stat().st_size
        assert manifest["telemetry"]["wall_seconds"] >= 0.0
        assert all(item["complete"] for item in manifest["outputs"])
        assert json.loads(result.stdout)["segments"] == 6

        bad = subprocess.run(
            [str(BINARY), "-I", str(input_path), "-O", str(work / "bad.tsv"),
             "--calling-z-score", "0"], text=True, capture_output=True,
        )
        assert bad.returncode != 0
        assert "Calling z-score threshold must be positive" in bad.stderr
        print(json.dumps({
            "status": "pass",
            "segments": len(rows),
            "calls": [row[5] for row in rows],
            "legacy_output": True,
            "invalid_threshold_rejected": True,
        }, sort_keys=True))


if __name__ == "__main__":
    main()
