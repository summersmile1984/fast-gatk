#!/usr/bin/env python3
"""File-boundary benchmark for the native orientation-prior hand-off."""

from __future__ import annotations

import json
import io
import os
import subprocess
import tarfile
import tempfile
import time
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_LEARN_ORIENTATION_BINARY",
        str(root / "fastgatk-native/build/fastgatk-learn-read-orientation-model")))
    records = int(os.environ.get("FASTGATK_ORIENTATION_BENCH_RECORDS", "100000"))
    if records < 1:
        raise SystemExit("FASTGATK_ORIENTATION_BENCH_RECORDS must be positive")
    with tempfile.TemporaryDirectory(prefix="fastgatk-orientation-benchmark-") as directory:
        work = Path(directory)
        input_path = work / "input.f1r2.tsv"
        output_path = work / "artifact-prior.tar.gz"
        manifest_path = work / "artifact-prior.manifest.json"
        with input_path.open("w", encoding="utf-8") as stream:
            stream.write("# FASTGATK-MUTECT2-F1R2 v2\n# contig\tposition\talt\tf1r2\tr1f2\ttotal\n")
            for index in range(records):
                alt = "ACGT"[index % 4]
                f1r2 = (index % 17) + 1
                r1f2 = (index % 11) + 1
                stream.write(f"17\t{index + 1}\t{alt}\t{f1r2}\t{r1f2}\t{f1r2 + r1f2}\n")
        env = os.environ.copy()
        env.setdefault("OMP_PROC_BIND", "true")
        env.setdefault("OMP_PLACES", "threads")
        start = time.perf_counter()
        completed = subprocess.run([
            str(native), "-I", str(input_path), "-O", str(output_path),
            "--sample", "BENCH", "--threads", "4", "--output-manifest", str(manifest_path)
        ], check=True, text=True, capture_output=True, env=env)
        wall_seconds = time.perf_counter() - start
        summary = json.loads(completed.stdout.splitlines()[-1])
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        # Also exercise the real three-member CollectF1R2Counts contract.  The
        # tiny synthetic archive keeps the benchmark self-contained while
        # forcing the standard parser and Java-compatible Beta-Binomial EM.
        standard_archive = work / "standard-collect-f1r2.tar.gz"
        ref_metrics = (
            "## htsjdk.samtools.metrics.StringHeader\n# BENCH\n\n"
            "## HISTOGRAM\tjava.lang.Integer\n"
            "depth\tACT\tAGT\n"
            "1\t100\t100\n2\t80\t80\n3\t40\t40\n"
        ).encode("utf-8")
        alt_metrics = (
            "## htsjdk.samtools.metrics.StringHeader\n# BENCH\n\n"
            "## HISTOGRAM\tjava.lang.Integer\n"
            "depth\tACT_C_F1R2\tACT_C_F2R1\tAGT_G_F2R1\tAGT_G_F1R2\n"
            "1\t20\t5\t5\t20\n2\t0\t0\t0\t0\n3\t0\t0\t0\t0\n"
        ).encode("utf-8")
        alt_table = (
            "#<METADATA>SAMPLE=BENCH\n"
            "context\tref_count\talt_count\tref_f1r2\talt_f1r2\tdepth\talt\n"
            "ACT\t30\t20\t15\t12\t50\tC\n"
        ).encode("utf-8")
        with tarfile.open(standard_archive, "w:gz") as archive:
            for name, payload in (("BENCH.ref_histogram", ref_metrics),
                                  ("BENCH.alt_histogram", alt_metrics),
                                  ("BENCH.alt_table", alt_table)):
                info = tarfile.TarInfo(name)
                info.size = len(payload)
                archive.addfile(info, fileobj=io.BytesIO(payload))
        standard_output = work / "standard-prior.tar.gz"
        standard_manifest_path = work / "standard-prior.manifest.json"
        standard_start = time.perf_counter()
        standard_completed = subprocess.run([
            str(native), "-I", str(standard_archive), "-O", str(standard_output),
            "--sample", "BENCH", "--threads", "4",
            "--output-manifest", str(standard_manifest_path),
        ], check=True, text=True, capture_output=True, env=env)
        standard_wall_seconds = time.perf_counter() - standard_start
        standard_summary = json.loads(standard_completed.stdout.splitlines()[-1])
        standard_manifest = json.loads(standard_manifest_path.read_text(encoding="utf-8"))
        assert standard_manifest["compatibility"]["full_gatk_em"] is True
        assert standard_manifest["compatibility"]["em_kokkos"] is True
        standard_telemetry = standard_manifest["telemetry"]
        assert standard_telemetry["em_kokkos_batches"] > 0
        assert standard_telemetry["em_kokkos_observations"] > 0
        assert standard_telemetry["kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect"
        )
        assert standard_telemetry["kernel_execution_policy"] == "RangePolicy"
        assert standard_telemetry["kernel_batches"] > 0
        assert standard_telemetry["kernel_observations"] > 0
        print(json.dumps({
            "status": "pass",
            "tool": "LearnReadOrientationModel",
            "execution_space": manifest["execution_space"],
            "records": records,
            "wall_seconds": wall_seconds,
            "records_per_second": records / wall_seconds,
            "output_bytes": output_path.stat().st_size,
            "f1r2_alt_support": summary["f1r2_alt_support"],
            "r1f2_alt_support": summary["r1f2_alt_support"],
            "kernel": {
                "lifecycle": manifest["telemetry"]["kernel_lifecycle"],
                "execution_policy": manifest["telemetry"]["kernel_execution_policy"],
                "batches": manifest["telemetry"]["kernel_batches"],
                "observations": manifest["telemetry"]["kernel_observations"],
                "prepare_seconds": manifest["telemetry"]["kernel_prepare_seconds"],
                "execute_seconds": manifest["telemetry"]["kernel_execute_seconds"],
            },
            "standard_em_wall_seconds": standard_wall_seconds,
            "standard_em_contexts_with_data": standard_summary["contexts_with_data"],
            "standard_em_records": standard_summary["records"],
            "standard_em_execution_space": standard_telemetry["em_execution_space"],
            "standard_em_kokkos_batches": standard_telemetry["em_kokkos_batches"],
            "standard_em_kokkos_observations": standard_telemetry["em_kokkos_observations"],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
