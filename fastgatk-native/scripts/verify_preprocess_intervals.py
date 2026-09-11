#!/usr/bin/env python3
"""Contract/oracle test for native PreprocessIntervals."""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile
import oracle_guard


def write_reference(path: pathlib.Path) -> None:
    records = [("chr1", "ACGTACGTGGCCNNNNATAT"), ("chr2", "NNNN")]
    offset = 0
    fai: list[tuple[str, int, int, int, int]] = []
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
            if line and not line.startswith("@")]


def main() -> int:
    root = pathlib.Path(__file__).resolve().parents[2]
    native = pathlib.Path(os.environ.get(
        "FASTGATK_PREPROCESS_INTERVALS_BINARY",
        str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
            / "fastgatk-preprocess-intervals"),
    ))
    if not native.exists():
        raise SystemExit(f"missing native binary: {native}")
    with tempfile.TemporaryDirectory(prefix="fastgatk-preprocess-intervals-") as directory:
        work = pathlib.Path(directory)
        reference = work / "reference.fasta"
        write_reference(reference)
        intervals = work / "targets.interval_list"
        intervals.write_text(
            "@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:20\n@SQ\tSN:chr2\tLN:4\n"
            "chr1\t1\t20\t+\twhole\nchr2\t1\t4\t+\tn_only\n",
            encoding="utf-8",
        )
        native_output = work / "native.interval_list"
        manifest = work / "native.manifest.json"
        native_result = subprocess.run([
            str(native), "-R", str(reference), "-L", str(intervals),
            "--bin-length", "5", "--padding", "0",
            "--interval-merging-rule", "OVERLAPPING_ONLY", "-O", str(native_output),
            "--output-manifest", str(manifest),
        ], check=True, text=True, capture_output=True)
        summary = json.loads(native_result.stdout.splitlines()[-1])
        native_rows = rows(native_output)
        expected = ["chr1\t1\t5\t+\t.", "chr1\t6\t10\t+\t.", "chr1\t11\t15\t+\t.", "chr1\t16\t20\t+\t."]
        if native_rows != expected:
            raise AssertionError({"expected": expected, "actual": native_rows})
        if summary["output_intervals"] != 4 or not manifest.is_file():
            raise AssertionError(summary)
        telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
        assert telemetry["kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect"
        )
        assert telemetry["kernel_execution_space"] in {"OpenMP", "Serial"}
        assert telemetry["kernel_execution_policy"] == "RangePolicy"
        assert telemetry["kernel_batches"] >= 1
        assert telemetry["kernel_records"] == summary["unfiltered_bins"]
        assert telemetry["kernel_prepare_seconds"] >= 0.0
        assert telemetry["kernel_execute_seconds"] >= 0.0

        # GATK's interval argument collection subtracts -XL regions before
        # padding/binning.  PreprocessIntervals intentionally requires
        # exclusion padding to remain zero; exercise the split-at-boundary
        # behavior and the manifest accounting so a no-op parser cannot pass.
        exclusions = work / "excluded.interval_list"
        exclusions.write_text(
            "@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:20\n@SQ\tSN:chr2\tLN:4\n"
            "chr1\t6\t10\t+\tmasked\n",
            encoding="utf-8",
        )
        excluded_output = work / "native-excluded.interval_list"
        excluded_manifest = work / "native-excluded.manifest.json"
        excluded_result = subprocess.run([
            str(native), "-R", str(reference), "-L", str(intervals),
            "-XL", str(exclusions), "--bin-length", "0", "--padding", "0",
            "--interval-exclusion-padding", "0",
            "--interval-merging-rule", "OVERLAPPING_ONLY", "-O", str(excluded_output),
            "--output-manifest", str(excluded_manifest),
        ], check=True, text=True, capture_output=True)
        excluded_summary = json.loads(excluded_result.stdout.splitlines()[-1])
        excluded_rows = rows(excluded_output)
        expected_excluded = ["chr1\t1\t5\t+\t.", "chr1\t11\t20\t+\t."]
        if excluded_rows != expected_excluded:
            raise AssertionError({"expected": expected_excluded, "actual": excluded_rows})
        excluded_metadata = json.loads(excluded_manifest.read_text(encoding="utf-8"))
        assert excluded_summary["excluded_intervals"] == 1
        # chr2 remains as an all-N interval and is removed later by the
        # reference-bin filter, so the post-subtraction interval count is 3
        # while only two intervals are emitted.
        assert excluded_summary["intervals_after_exclusion"] == 3
        assert excluded_metadata["exclusion_padding"] == 0

        # The generic interval argument default is ALL, but GATK's
        # PreprocessIntervals validation rejects both omitted/default ALL and
        # explicit ALL. Keep native replacement fail-closed at this boundary.
        merge_rejection = True
        for merge_args in ([], ["--interval-merging-rule", "ALL"]):
            rejected = subprocess.run([
                str(native), "-R", str(reference), "-L", str(intervals),
                "--bin-length", "5", "--padding", "0", *merge_args,
                "-O", str(work / ("rejected-default.interval_list" if not merge_args
                                  else "rejected-all.interval_list")),
            ], text=True, capture_output=True)
            if rejected.returncode == 0 or "Interval merging rule must be set to OVERLAPPING_ONLY." not in rejected.stderr:
                raise AssertionError({"merge_args": merge_args, "returncode": rejected.returncode,
                                      "stderr": rejected.stderr})

        java = root / "third_party/jdk17/bin/java"
        jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        java_checked = False
        if oracle_guard.oracle_ready('verify_preprocess_intervals.py', java, jar):
            java_output = work / "java.interval_list"
            java_result = subprocess.run([
                str(java), "-jar", str(jar), "PreprocessIntervals", "-R", str(reference),
                "-L", str(intervals), "--bin-length", "5", "--padding", "0",
                "--interval-merging-rule", "OVERLAPPING_ONLY", "-O", str(java_output),
            ], text=True, capture_output=True)
            if java_result.returncode != 0:
                raise AssertionError("Java GATK PreprocessIntervals failed:\n" + java_result.stderr[-4000:])
            java_rows = rows(java_output)
            if java_rows != expected:
                raise AssertionError({"expected": expected, "java": java_rows})
            java_excluded_output = work / "java-excluded.interval_list"
            java_excluded_result = subprocess.run([
                str(java), "-jar", str(jar), "PreprocessIntervals", "-R", str(reference),
                "-L", str(intervals), "-XL", str(exclusions), "--bin-length", "0",
                "--padding", "0", "--interval-exclusion-padding", "0",
                "--interval-merging-rule", "OVERLAPPING_ONLY", "-O", str(java_excluded_output),
            ], text=True, capture_output=True)
            if java_excluded_result.returncode != 0:
                raise AssertionError("Java GATK PreprocessIntervals -XL failed:\n" + java_excluded_result.stderr[-4000:])
            java_excluded_rows = rows(java_excluded_output)
            if java_excluded_rows != expected_excluded:
                raise AssertionError({"expected": expected_excluded, "java": java_excluded_rows})
            for merge_args in ([], ["--interval-merging-rule", "ALL"]):
                java_rejected = subprocess.run([
                    str(java), "-jar", str(jar), "PreprocessIntervals", "-R", str(reference),
                    "-L", str(intervals), "--bin-length", "5", "--padding", "0", *merge_args,
                    "-O", str(work / ("java-rejected-default.interval_list" if not merge_args
                                      else "java-rejected-all.interval_list")),
                ], text=True, capture_output=True)
                if java_rejected.returncode == 0 or "Interval merging rule must be set to OVERLAPPING_ONLY." not in (
                        java_rejected.stdout + java_rejected.stderr):
                    raise AssertionError({"java_merge_args": merge_args, "returncode": java_rejected.returncode,
                                          "stdout": java_rejected.stdout[-1000:], "stderr": java_rejected.stderr[-1000:]})
            java_checked = True
        print(json.dumps({"status": "pass", "input_intervals": 2,
                          "output_intervals": len(native_rows), "java_oracle": java_checked,
                          "interval_merging_rule_rejection": merge_rejection}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
