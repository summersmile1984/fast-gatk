#!/usr/bin/env python3
"""Contract/oracle test for annotation-based native FilterIntervals."""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile
import oracle_guard


def data_rows(path: pathlib.Path) -> list[tuple[str, int, int]]:
    rows: list[tuple[str, int, int]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        rows.append((fields[0], int(fields[1]), int(fields[2])))
    return rows


def main() -> int:
    root = pathlib.Path(__file__).resolve().parents[2]
    native = pathlib.Path(os.environ.get(
        "FASTGATK_FILTER_INTERVALS_BINARY",
        str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
            / "fastgatk-filter-intervals"),
    ))
    collect = pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build")) / "fastgatk-collect-read-counts"
    if not native.exists():
        raise SystemExit(f"missing native binary: {native}")
    with tempfile.TemporaryDirectory(prefix="fastgatk-filter-intervals-") as directory:
        work = pathlib.Path(directory)
        intervals = work / "targets.interval_list"
        intervals.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"
            "chr1\t1\t100\t+\ti1\nchr1\t201\t300\t+\ti2\n"
            "chr1\t401\t500\t+\ti3\nchr1\t601\t700\t+\ti4\n"
            "chr1\t801\t900\t+\ti5\n",
            encoding="utf-8",
        )
        annotations = work / "annotated.tsv"
        annotations.write_text(
            "@HD\tVN:1.5\n@SQ\tSN:chr1\tLN:1000\n"
            "CONTIG\tSTART\tEND\tGC_CONTENT\n"
            "chr1\t1\t100\t0.500000\n"
            "chr1\t201\t300\t0.050000\n"
            "chr1\t401\t500\t0.800000\n"
            "chr1\t601\t700\t0.950000\n"
            "chr1\t801\t900\t0.700000\n",
            encoding="utf-8",
        )
        output = work / "filtered.interval_list"
        manifest = work / "filtered.manifest.json"
        completed = subprocess.run([
            str(native), "-L", str(intervals), "--annotated-intervals", str(annotations),
            "-XL", "chr1:801-900", "--minimum-gc-content", "0.1",
            "--maximum-gc-content", "0.9", "--interval-merging-rule", "OVERLAPPING_ONLY", "-O", str(output),
            "--output-manifest", str(manifest),
        ], check=True, text=True, capture_output=True)
        summary = json.loads(completed.stdout.splitlines()[-1])
        expected = [("chr1", 1, 100), ("chr1", 401, 500)]
        native_rows = data_rows(output)
        if native_rows != expected:
            raise AssertionError({"expected": expected, "native": native_rows})
        if summary["output_intervals"] != len(expected) or not manifest.is_file():
            raise AssertionError(summary)
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        telemetry = metadata["telemetry"]
        assert telemetry["kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect"
        )
        assert telemetry["kernel_execution_space"] in {"OpenMP", "Serial"}
        assert telemetry["kernel_execution_policy"] == "RangePolicy"
        assert telemetry["kernel_batches"] == 1
        assert telemetry["kernel_records"] == 5
        assert telemetry["kernel_prepare_seconds"] >= 0.0
        assert telemetry["kernel_execute_seconds"] >= 0.0

        # GATK's short aliases are part of the command contract.  FilterIntervals
        # rejects non-zero include/exclude padding during argument validation;
        # accepting it would silently change annotation-row matching.
        alias_output = work / "alias.interval_list"
        alias_manifest = work / "alias.manifest.json"
        subprocess.run([
            str(native), "-L", str(intervals), "--annotated-intervals", str(annotations),
            "-XL", "chr1:801-900", "--minimum-gc-content", "0.1",
            "--maximum-gc-content", "0.9", "-imr", "OVERLAPPING_ONLY",
            "-isr", "UNION", "-O", str(alias_output), "--output-manifest", str(alias_manifest),
        ], check=True, text=True, capture_output=True)
        if data_rows(alias_output) != expected:
            raise AssertionError({"expected_alias": expected, "native_alias": data_rows(alias_output)})
        alias_metadata = json.loads(alias_manifest.read_text(encoding="utf-8"))
        if alias_metadata["interval_merging_rule"] != "OVERLAPPING_ONLY" or \
                alias_metadata["interval_set_rule"] != "UNION":
            raise AssertionError({"alias_manifest": alias_metadata})
        for padding_flag, message in (("-ip", "Interval padding must be set to 0."),
                                      ("-ixp", "Interval exclusion padding must be set to 0.")):
            rejected = subprocess.run([
                str(native), "-L", str(intervals), "--annotated-intervals", str(annotations),
                padding_flag, "1", "-O", str(work / f"rejected-{padding_flag[1:]}.list"),
            ], text=True, capture_output=True)
            if rejected.returncode == 0 or message not in rejected.stderr:
                raise AssertionError({"padding_flag": padding_flag, "returncode": rejected.returncode,
                                      "stderr": rejected.stderr})
        for merge_args in ([], ["-imr", "ALL"]):
            rejected_merge = subprocess.run([
                str(native), "-L", str(intervals), "--annotated-intervals", str(annotations),
                *merge_args, "-O", str(work / ("rejected-default.list" if not merge_args else "rejected-all.list")),
            ], text=True, capture_output=True)
            if rejected_merge.returncode == 0 or \
                    "Interval merging rule must be set to OVERLAPPING_ONLY." not in rejected_merge.stderr:
                raise AssertionError({"merge_args": merge_args, "returncode": rejected_merge.returncode,
                                      "stderr": rejected_merge.stderr})

        # Although the generic GATK argument help lists INTERSECTION, the
        # FilterIntervals tool validates the shared interval collection and
        # rejects it.  Keep that fail-closed behavior in the native path.
        left_selector = work / "left.interval_list"
        left_selector.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"
            "chr1\t1\t100\t+\tleft1\n"
            "chr1\t401\t500\t+\tleft2\n"
            "chr1\t801\t900\t+\tleft3\n", encoding="utf-8")
        right_selector = work / "right.interval_list"
        right_selector.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"
            "chr1\t401\t500\t+\tright1\n"
            "chr1\t801\t900\t+\tright2\n", encoding="utf-8")
        intersection_result = subprocess.run([
            str(native), "-L", str(left_selector), "-L", str(right_selector),
            "--interval-set-rule", "INTERSECTION", "--annotated-intervals", str(annotations),
            "--minimum-gc-content", "0.1", "--maximum-gc-content", "0.9",
            "--interval-merging-rule", "OVERLAPPING_ONLY", "-O", str(work / "intersection.list"),
        ], text=True, capture_output=True)
        if intersection_result.returncode == 0 or \
                "Interval set rule must be set to UNION." not in intersection_result.stderr:
            raise AssertionError({"intersection_returncode": intersection_result.returncode,
                                  "stderr": intersection_result.stderr})

        # Sequence dictionaries are not necessarily lexically sorted (chr2
        # precedes chr10 here).  UNION must preserve dictionary order, matching
        # Picard/GATK output rather than Python/C++ lexical order.
        multi_intervals = work / "multi.interval_list"
        multi_intervals.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr2\tLN:1000\n@SQ\tSN:chr10\tLN:1000\n"
            "chr10\t1\t100\t+\ta10\nchr10\t201\t300\t+\tb10\n"
            "chr2\t1\t100\t+\ta2\nchr2\t201\t300\t+\tb2\n", encoding="utf-8")
        multi_right = work / "multi-right.interval_list"
        multi_right.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr2\tLN:1000\n@SQ\tSN:chr10\tLN:1000\n"
            "chr2\t1\t100\t+\tr2a\nchr2\t201\t300\t+\tr2b\n"
            "chr10\t1\t100\t+\tr10a\nchr10\t201\t300\t+\tr10b\n", encoding="utf-8")
        multi_annotations = work / "multi-annotated.tsv"
        multi_annotations.write_text(
            "@HD\tVN:1.5\n@SQ\tSN:chr2\tLN:1000\n@SQ\tSN:chr10\tLN:1000\n"
            "CONTIG\tSTART\tEND\tGC_CONTENT\n"
            "chr2\t1\t100\t0.5\nchr2\t201\t300\t0.5\n"
            "chr10\t1\t100\t0.5\nchr10\t201\t300\t0.5\n", encoding="utf-8")
        multi_output = work / "multi-output.interval_list"
        subprocess.run([
            str(native), "-L", str(multi_intervals), "-L", str(multi_right),
            "-imr", "OVERLAPPING_ONLY", "-isr", "UNION",
            "--annotated-intervals", str(multi_annotations),
            "-O", str(multi_output),
        ], check=True, text=True, capture_output=True)
        expected_multi = [("chr2", 1, 100), ("chr2", 201, 300),
                          ("chr10", 1, 100), ("chr10", 201, 300)]
        if data_rows(multi_output) != expected_multi:
            raise AssertionError({"expected_multi": expected_multi,
                                  "native_multi": data_rows(multi_output)})

        java = root / "third_party/jdk17/bin/java"
        jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        java_checked = False
        if oracle_guard.oracle_ready('verify_filter_intervals.py', java, jar):
            java_output = work / "java.interval_list"
            result = subprocess.run([
                str(java), "-jar", str(jar), "FilterIntervals", "-L", str(intervals),
                "--annotated-intervals", str(annotations), "-XL", "chr1:801-900",
                "--minimum-gc-content", "0.1", "--maximum-gc-content", "0.9",
                "--interval-merging-rule", "OVERLAPPING_ONLY",
                "-O", str(java_output),
            ], text=True, capture_output=True)
            if result.returncode != 0:
                raise AssertionError("Java GATK FilterIntervals failed:\n" + result.stderr[-4000:])
            java_rows = data_rows(java_output)
            if java_rows != expected:
                raise AssertionError({"expected": expected, "java": java_rows})
            multi_java_output = work / "multi-java.interval_list"
            multi_java = subprocess.run([
                str(java), "-jar", str(jar), "FilterIntervals", "-L", str(multi_intervals),
                "-L", str(multi_right), "-imr", "OVERLAPPING_ONLY", "-isr", "UNION",
                "--annotated-intervals", str(multi_annotations), "-O", str(multi_java_output),
            ], text=True, capture_output=True)
            if multi_java.returncode != 0:
                raise AssertionError("Java GATK multi-contig FilterIntervals failed:\n" + multi_java.stderr[-4000:])
            if data_rows(multi_java_output) != expected_multi:
                raise AssertionError({"expected_multi": expected_multi,
                                      "java_multi": data_rows(multi_java_output)})
            java_checked = True

        count_one = work / "counts-1.tsv"
        count_two = work / "counts-2.tsv"
        count_header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:S\tSM:S\nCONTIG\tSTART\tEND\tCOUNT\n"
        count_rows_one = "".join(
            f"chr1\t{index * 100 + 1}\t{index * 100 + 100}\t{value}\n"
            for index, value in enumerate([1, 10, 10, 10, 10, 100])
        )
        count_rows_two = "".join(
            f"chr1\t{index * 100 + 1}\t{index * 100 + 100}\t{value}\n"
            for index, value in enumerate([1, 10, 10, 10, 10, 100])
        )
        count_one.write_text(count_header + count_rows_one, encoding="utf-8")
        count_two.write_text(count_header + count_rows_two, encoding="utf-8")
        count_intervals = work / "count-targets.interval_list"
        count_intervals.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n" + "".join(
                f"chr1\t{index * 100 + 1}\t{index * 100 + 100}\t+\ti{index}\n"
                for index in range(6)
            ), encoding="utf-8")
        count_output = work / "count-filtered.interval_list"
        count_result = subprocess.run([
            str(native), "-L", str(count_intervals), "-I", str(count_one), "-I", str(count_two),
            "--interval-merging-rule", "OVERLAPPING_ONLY",
            "--low-count-filter-count-threshold", "5",
            "--low-count-filter-percentage-of-samples", "50",
            "--extreme-count-filter-minimum-percentile", "1",
            "--extreme-count-filter-maximum-percentile", "99",
            "--extreme-count-filter-percentage-of-samples", "90",
            "-O", str(count_output),
        ], check=True, text=True, capture_output=True)
        count_summary = json.loads(count_result.stdout.splitlines()[-1])
        expected_count = [("chr1", 101, 200), ("chr1", 201, 300),
                          ("chr1", 301, 400), ("chr1", 401, 500),
                          ("chr1", 501, 600)]
        if data_rows(count_output) != expected_count:
            raise AssertionError({"expected_count": expected_count, "native_count": data_rows(count_output)})

        count_java_checked = False
        if oracle_guard.oracle_ready('verify_filter_intervals.py', java, jar):
            count_java_output = work / "java-count-filtered.interval_list"
            count_java_result = subprocess.run([
                str(java), "-jar", str(jar), "FilterIntervals", "-L", str(count_intervals),
                "-I", str(count_one), "-I", str(count_two),
                "--interval-merging-rule", "OVERLAPPING_ONLY",
                "--low-count-filter-count-threshold", "5",
                "--low-count-filter-percentage-of-samples", "50",
                "--extreme-count-filter-minimum-percentile", "1",
                "--extreme-count-filter-maximum-percentile", "99",
                "--extreme-count-filter-percentage-of-samples", "90",
                "-O", str(count_java_output),
            ], text=True, capture_output=True)
            if count_java_result.returncode != 0:
                raise AssertionError("Java GATK count FilterIntervals failed:\n" + count_java_result.stderr[-4000:])
            if data_rows(count_java_output) != expected_count:
                raise AssertionError({"expected_count": expected_count, "java_count": data_rows(count_java_output)})
            count_java_checked = True
        # Exercise the native HDF5 SimpleCountCollection reader through the
        # same count-filter path.  The values are intentionally unconstrained
        # here; the TSV case above remains the numerical GATK percentile oracle.
        hdf5_intervals = work / "hdf5-targets.interval_list"
        hdf5_intervals.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"
            "chr1\t1\t100\t+\th1\nchr1\t201\t300\t+\th2\n",
            encoding="utf-8")
        sam = work / "hdf5-reads.sam"
        sam.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n" +
            "r0\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI\n"
            "r1\t0\tchr1\t201\t60\t1M\t*\t0\t0\tA\tI\n", encoding="utf-8")
        hdf5_count = work / "counts.h5"
        subprocess.run([
            str(collect), "-I", str(sam), "-L", str(hdf5_intervals), "-O", str(hdf5_count),
            "--format", "HDF5", "--sample", "HDF5_SAMPLE", "--minimum-mapping-quality", "0",
        ], check=True, text=True, capture_output=True)
        hdf5_output = work / "hdf5-filtered.interval_list"
        hdf5_result = subprocess.run([
            str(native), "-L", str(hdf5_intervals), "-I", str(hdf5_count),
            "--interval-merging-rule", "OVERLAPPING_ONLY",
            "--low-count-filter-count-threshold", "0",
            "--low-count-filter-percentage-of-samples", "100",
            "--extreme-count-filter-minimum-percentile", "0",
            "--extreme-count-filter-maximum-percentile", "100",
            "--extreme-count-filter-percentage-of-samples", "100",
            "-O", str(hdf5_output),
        ], text=True, capture_output=True)
        if hdf5_result.returncode != 0:
            raise AssertionError("native HDF5 FilterIntervals failed:\n" + hdf5_result.stderr[-4000:])
        hdf5_summary = json.loads(hdf5_result.stdout.splitlines()[-1])
        if data_rows(hdf5_output) != [("chr1", 1, 100), ("chr1", 201, 300)]:
            raise AssertionError({"hdf5_rows": data_rows(hdf5_output)})
        print(json.dumps({"status": "pass", "input_intervals": 5, "count_input_intervals": 6,
                          "output_intervals": len(native_rows), "java_oracle": java_checked,
                          "short_aliases": True, "padding_rejection": True,
                          "merge_validation": True, "set_rule_validation": True,
                          "dictionary_order_union": True,
                          "count_output_intervals": count_summary["output_intervals"],
                          "count_java_oracle": count_java_checked,
                          "hdf5_output_intervals": hdf5_summary["output_intervals"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
