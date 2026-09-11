#!/usr/bin/env python3
"""Contract/oracle test for CountBasesInReference's Kokkos count path."""

from __future__ import annotations

import json
import gzip
import os
import pathlib
import subprocess
import tempfile
import oracle_guard


def write_reference(path: pathlib.Path) -> None:
    records = [("chr1", "AACCGGTTNN"), ("chr2", "TTAA")]
    offset = 0
    fai: list[tuple[str, int, int, int, int]] = []
    with path.open("wb") as stream:
        for name, sequence in records:
            header = f">{name}\n".encode("ascii")
            stream.write(header)
            offset += len(header)
            body = (sequence + "\n").encode("ascii")
            stream.write(body)
            fai.append((name, len(sequence), offset, len(sequence), len(sequence) + 1))
            offset += len(body)
    path.with_suffix(path.suffix + ".fai").write_text(
        "".join(f"{name}\t{length}\t{seq_offset}\t{line_bases}\t{line_width}\n"
                for name, length, seq_offset, line_bases, line_width in fai),
        encoding="ascii",
    )
    path.with_suffix(".dict").write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n@SQ\tSN:chr2\tLN:4\n",
        encoding="ascii",
    )


def main() -> int:
    root = pathlib.Path(__file__).resolve().parents[2]
    native = pathlib.Path(os.environ.get(
        "FASTGATK_COUNT_BASES_BINARY",
        str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
            / "fastgatk-count-bases-in-reference"),
    ))
    if not native.exists():
        raise SystemExit(f"missing native binary: {native}")
    with tempfile.TemporaryDirectory(prefix="fastgatk-count-bases-") as directory:
        work = pathlib.Path(directory)
        reference = work / "reference.fasta"
        write_reference(reference)
        output = work / "counts.txt"
        manifest = work / "counts.manifest.json"
        result = subprocess.run([
            str(native), "-R", str(reference), "-L", "chr1:2-7", "-L", "chr2",
            "-O", str(output), "--output-manifest", str(manifest), "--threads", "2",
        ], check=True, text=True, capture_output=True)
        expected = "A : 3\nC : 2\nG : 2\nT : 3\n"
        if output.read_text(encoding="ascii") != expected:
            raise AssertionError({"expected": expected, "actual": output.read_text(encoding="ascii")})
        if result.stdout != expected:
            raise AssertionError({"stdout": result.stdout, "expected": expected})
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        telemetry = metadata["telemetry"]
        assert metadata["regions"] == 2 and metadata["bases"] == 10
        assert telemetry["kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect"
        )
        assert telemetry["kernel_execution_space"] in {"OpenMP", "Serial"}
        assert telemetry["kernel_execution_policy"] == "RangePolicy"
        assert telemetry["kernel_records"] == 10
        assert telemetry["kernel_prepare_seconds"] >= 0.0
        assert telemetry["kernel_execute_seconds"] >= 0.0

        interval_file = work / "targets.interval_list"
        interval_file.write_text("@HD\tVN:1.6\nchr1\t1\t2\nchr1\t2\t4\n",
                                 encoding="ascii")
        interval_output = work / "interval-counts.txt"
        interval_result = subprocess.run([
            str(native), "--reference", str(reference), "--intervals", str(interval_file),
            "--output", str(interval_output),
        ], text=True, capture_output=True, check=True)
        # Adjacent/overlapping interval-list rows are merged before fetching;
        # chr1:1-4 therefore contributes AACC exactly once.
        assert interval_output.read_text(encoding="ascii") == "A : 2\nC : 2\n"
        assert interval_result.stdout == "A : 2\nC : 2\n"

        bed_file = work / "targets.bed"
        bed_file.write_text("track name=targets\nchr2\t1\t4\n", encoding="ascii")
        bed_result = subprocess.run([
            str(native), "-R", str(reference), "-L", str(bed_file),
        ], text=True, capture_output=True, check=True)
        assert bed_result.stdout == "A : 2\nT : 1\n"

        compressed_interval_file = work / "targets.interval_list.gz"
        with gzip.open(compressed_interval_file, "wt", encoding="ascii") as stream:
            stream.write("@HD\tVN:1.6\nchr1\t5\t7\n")
        compressed_interval_result = subprocess.run([
            str(native), "--reference", str(reference), "--intervals", str(compressed_interval_file),
        ], text=True, capture_output=True, check=True)
        assert compressed_interval_result.stdout == "G : 2\nT : 1\n"

        compressed_bed_file = work / "targets.bed.gz"
        with gzip.open(compressed_bed_file, "wt", encoding="ascii") as stream:
            stream.write("track name=targets\nchr2\t0\t2\n")
        compressed_bed_result = subprocess.run([
            str(native), "-R", str(reference), "-L", str(compressed_bed_file),
        ], text=True, capture_output=True, check=True)
        assert compressed_bed_result.stdout == "T : 2\n"

        malformed_interval_file = work / "malformed.interval_list"
        malformed_interval_file.write_text("chr1\t1\n", encoding="ascii")
        malformed_result = subprocess.run([
            str(native), "-R", str(reference), "-L", str(malformed_interval_file),
        ], text=True, capture_output=True, check=False)
        assert malformed_result.returncode != 0
        assert "malformed interval" in malformed_result.stderr

        java = root / "third_party/jdk17/bin/java"
        jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        java_oracle = False
        if oracle_guard.oracle_ready('verify_count_bases_in_reference.py', java, jar):
            java_result = subprocess.run([
                str(java), "-jar", str(jar), "CountBasesInReference", "-R", str(reference),
                "-L", "chr1:2-7", "-L", "chr2",
            ], text=True, capture_output=True)
            if java_result.returncode == 0:
                java_lines = [line for line in java_result.stdout.splitlines() if " : " in line]
                if "\n".join(java_lines) + "\n" != expected:
                    raise AssertionError({"java": java_lines, "expected": expected})
                java_oracle = True
        print(json.dumps({"status": "pass", "regions": 2, "bases": 10,
                          "java_oracle": java_oracle}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
