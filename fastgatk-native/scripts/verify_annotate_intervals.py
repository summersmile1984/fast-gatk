#!/usr/bin/env python3
"""Contract/oracle test for native AnnotateIntervals GC annotations."""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile


def write_reference(path: pathlib.Path) -> None:
    records = [("chr1", "ACGTACGT" "GGCC" "NNNN" "ATAT"), ("chr2", "NNNN")]
    offset = 0
    fai = []
    with path.open("wb") as output:
        for name, sequence in records:
            header = f">{name}\n".encode()
            output.write(header)
            offset += len(header)
            body = (sequence + "\n").encode()
            output.write(body)
            fai.append((name, len(sequence), offset, len(sequence), len(sequence) + 1))
            offset += len(body)
    path.with_suffix(path.suffix + ".fai").write_text(
        "".join(f"{name}\t{length}\t{seq_offset}\t{line_bases}\t{line_width}\n"
                for name, length, seq_offset, line_bases, line_width in fai),
        encoding="utf-8",
    )
    path.with_suffix(".dict").write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:20\n@SQ\tSN:chr2\tLN:4\n",
        encoding="utf-8",
    )


def rows(path: pathlib.Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("@") and not line.startswith("CONTIG")]


def main() -> int:
    root = pathlib.Path(__file__).resolve().parents[2]
    native = pathlib.Path(os.environ.get(
        "FASTGATK_ANNOTATE_INTERVALS_BINARY",
        str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
            / "fastgatk-annotate-intervals"),
    ))
    if not native.exists():
        raise SystemExit(f"missing native binary: {native}")
    with tempfile.TemporaryDirectory(prefix="fastgatk-annotate-intervals-") as directory:
        work = pathlib.Path(directory)
        reference = work / "reference.fasta"
        write_reference(reference)
        intervals = work / "targets.interval_list"
        intervals.write_text("@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:20\n@SQ\tSN:chr2\tLN:4\n"
                            "chr1\t1\t8\t+\tbin1\nchr1\t9\t16\t+\tbin2\nchr2\t1\t4\t+\tbin3\n", encoding="utf-8")
        native_output = work / "native.tsv"
        manifest = work / "native.manifest.json"
        native_result = subprocess.run([
            str(native), "-R", str(reference), "-L", str(intervals), "-O", str(native_output),
            "--interval-merging-rule", "OVERLAPPING_ONLY", "--feature-query-lookahead", "0",
            "--output-manifest", str(manifest),
        ], check=True, text=True, capture_output=True)
        summary = json.loads(native_result.stdout.splitlines()[-1])
        expected = ["chr1\t1\t8\t0.500000", "chr1\t9\t16\t1.000000", "chr2\t1\t4\tNaN"]
        if rows(native_output) != expected:
            raise AssertionError({"expected": expected, "actual": rows(native_output)})
        if summary["intervals"] != 3 or not manifest.is_file():
            raise AssertionError(summary)
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        telemetry = metadata["telemetry"]
        assert telemetry["kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect"
        )
        assert telemetry["kernel_execution_space"] in {"OpenMP", "Serial"}
        assert telemetry["kernel_execution_policy"] == "RangePolicy"
        assert telemetry["track_index_strategy"] == "per-contig-binary-search"
        assert telemetry["kernel_batches"] == 1
        assert telemetry["kernel_records"] == 3
        assert metadata["feature_query_lookahead"] == 0
        assert telemetry["feature_query_lookahead"] == 0
        assert telemetry["kernel_prepare_seconds"] >= 0.0
        assert telemetry["kernel_execute_seconds"] >= 0.0

        java = root / "third_party/jdk17/bin/java"
        jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        java_checked = False
        if java.exists() and jar.exists():
            java_output = work / "java.tsv"
            java_result = subprocess.run([
                str(java), "-jar", str(jar), "AnnotateIntervals", "-R", str(reference),
                "-L", str(intervals), "--interval-merging-rule", "OVERLAPPING_ONLY",
                "--feature-query-lookahead", "0", "-O", str(java_output),
            ], text=True, capture_output=True)
            if java_result.returncode != 0:
                raise AssertionError("Java GATK AnnotateIntervals failed:\n" + java_result.stderr[-4000:])
            java_rows = rows(java_output)
            if java_rows != expected:
                raise AssertionError({"expected": expected, "java": java_rows})
            java_checked = True

        # Optional BED tracks exercise the same htsjdk BEDCodec coordinate
        # convention (start+1 / end-1), score-column handling, NaN→unity,
        # and length-weighted averages used by AnnotateIntervals.
        mappability = work / "mappability.bed"
        mappability.write_text(
            "chr1\t0\t5\tmap-unscored\n"
            "chr1\t10\t20\tmap-scored\t0.5\n"
            "chr2\t0\t4\tmap-nan\tNaN\n", encoding="utf-8",
        )
        segmental = work / "segmental.bed"
        segmental.write_text("chr1\t5\t15\tseg\t0.25\n", encoding="utf-8")
        track_native = work / "tracks-native.tsv"
        track_manifest = work / "tracks-native.manifest.json"
        track_result = subprocess.run([
            str(native), "-R", str(reference), "-L", str(intervals),
            "--mappability-track", str(mappability),
            "--segmental-duplication-track", str(segmental),
            "--interval-merging-rule", "OVERLAPPING_ONLY", "-O", str(track_native),
            "--output-manifest", str(track_manifest),
        ], check=True, text=True, capture_output=True)
        track_expected = [
            "chr1\t1\t8\t0.500000\t0.500000\t0.093750",
            "chr1\t9\t16\t1.000000\t0.375000\t0.187500",
            "chr2\t1\t4\tNaN\t0.750000\t0.000000",
        ]
        if rows(track_native) != track_expected:
            raise AssertionError({"expected": track_expected, "actual": rows(track_native)})
        track_metadata = json.loads(track_manifest.read_text(encoding="utf-8"))
        assert track_metadata["mappability_segments"] == 3
        assert track_metadata["segmental_duplication_segments"] == 1
        assert track_metadata["telemetry"]["track_index_strategy"] == "per-contig-binary-search"
        if java.exists() and jar.exists():
            for track in (mappability, segmental):
                index_result = subprocess.run([
                    str(java), "-jar", str(jar), "IndexFeatureFile", "-I", str(track),
                ], text=True, capture_output=True)
                if index_result.returncode != 0:
                    raise AssertionError("Java IndexFeatureFile failed:\n" + index_result.stderr[-4000:])
            track_java = work / "tracks-java.tsv"
            java_track_result = subprocess.run([
                str(java), "-jar", str(jar), "AnnotateIntervals", "-R", str(reference),
                "-L", str(intervals), "--mappability-track", str(mappability),
                "--segmental-duplication-track", str(segmental),
                "--interval-merging-rule", "OVERLAPPING_ONLY", "-O", str(track_java),
            ], text=True, capture_output=True)
            if java_track_result.returncode != 0:
                raise AssertionError("Java GATK AnnotateIntervals tracks failed:\n" +
                                     java_track_result.stderr[-4000:])
            if track_native.read_text(encoding="utf-8").splitlines() != \
                    track_java.read_text(encoding="utf-8").splitlines():
                raise AssertionError("native/Java AnnotateIntervals track output differs")

        # Copy-number tools require OVERLAPPING_ONLY but still apply -XL
        # subtraction.  Padding is deliberately fail-closed: GATK's
        # CopyNumberArgumentValidationUtils rejects every non-zero value.
        boundary_cases = [
            ("overlap_only", ["-L", "chr1:1-4", "-L", "chr1:5-8",
                              "--interval-merging-rule", "OVERLAPPING_ONLY"],
             ["chr1\t1\t4\t0.500000", "chr1\t5\t8\t0.500000"]),
            ("exclude", ["-L", "chr1:1-12", "-XL", "chr1:5-8",
                         "--interval-merging-rule", "OVERLAPPING_ONLY"],
             ["chr1\t1\t4\t0.500000", "chr1\t9\t12\t1.000000"]),
        ]
        boundary_java_checked = False
        for label, arguments, expected_rows in boundary_cases:
            native_case = work / f"native-{label}.tsv"
            native_manifest = work / f"native-{label}.manifest.json"
            native_run = subprocess.run([
                str(native), "-R", str(reference), *arguments, "-O", str(native_case),
                "--output-manifest", str(native_manifest),
            ], text=True, capture_output=True)
            if native_run.returncode != 0:
                raise AssertionError(f"native AnnotateIntervals {label} failed:\n" + native_run.stderr[-4000:])
            native_rows = rows(native_case)
            if native_rows != expected_rows:
                raise AssertionError({"case": label, "expected": expected_rows, "native": native_rows})
            native_meta = json.loads(native_manifest.read_text(encoding="utf-8"))
            if label == "exclude":
                assert native_meta["excluded_intervals"] == 1
            if java.exists() and jar.exists():
                java_case = work / f"java-{label}.tsv"
                java_run = subprocess.run([
                    str(java), "-Xmx1g", "-jar", str(jar), "AnnotateIntervals", "-R", str(reference),
                    *arguments, "-O", str(java_case),
                ], text=True, capture_output=True)
                if java_run.returncode != 0:
                    raise AssertionError(f"Java GATK AnnotateIntervals {label} failed:\n" + java_run.stderr[-4000:])
                java_rows = rows(java_case)
                if java_rows != expected_rows or native_rows != java_rows:
                    raise AssertionError({"case": label, "expected": expected_rows,
                                          "java": java_rows, "native": native_rows})
                boundary_java_checked = True
        for option in ("--interval-padding", "--interval-exclusion-padding"):
            native_reject = subprocess.run([
                str(native), "-R", str(reference), "-L", "chr1:1-8",
                "--interval-merging-rule", "OVERLAPPING_ONLY", option, "1",
                "-O", str(work / f"native-reject-{option[11:]}.tsv"),
            ], text=True, capture_output=True)
            expected_error = (
                "Interval exclusion padding must be set to 0."
                if option == "--interval-exclusion-padding"
                else "Interval padding must be set to 0."
            )
            if native_reject.returncode == 0 or expected_error not in native_reject.stderr:
                raise AssertionError({"option": option, "native": native_reject.stderr})
            if java.exists() and jar.exists():
                java_reject = subprocess.run([
                    str(java), "-Xmx1g", "-jar", str(jar), "AnnotateIntervals", "-R", str(reference),
                    "-L", "chr1:1-8", "--interval-merging-rule", "OVERLAPPING_ONLY", option, "1",
                    "-O", str(work / f"java-reject-{option[11:]}.tsv"),
                ], text=True, capture_output=True)
                if java_reject.returncode == 0 or expected_error not in java_reject.stderr:
                    raise AssertionError({"option": option, "java": java_reject.stderr[-4000:]})
        print(json.dumps({"status": "pass", "intervals": 3, "java_oracle": java_checked,
                          "boundary_java_oracle": boundary_java_checked,
                          "native_rows": len(rows(native_output))}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
