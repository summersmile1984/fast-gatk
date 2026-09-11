#!/usr/bin/env python3
"""Contract and Java-reader check for the native GATK HDF5 PoN writer."""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile
import oracle_guard


def write_counts(path: pathlib.Path, sample: str, values: list[int]) -> None:
    lines = ["@HD\tVN:1.6", "@SQ\tSN:chr1\tLN:100000", f"@RG\tID:{sample}\tSM:{sample}",
             "CONTIG\tSTART\tEND\tCOUNT"]
    for index, value in enumerate(values):
        start = index * 100 + 1
        lines.append(f"chr1\t{start}\t{start + 99}\t{value}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    root = pathlib.Path(__file__).resolve().parents[2]
    binary = pathlib.Path(os.environ.get(
        "FASTGATK_CREATE_PON_BINARY",
        str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
            / "fastgatk-create-read-count-panel-of-normals")))
    if not binary.exists():
        raise SystemExit(f"missing native binary: {binary}")
    with tempfile.TemporaryDirectory(prefix="fastgatk-pon-contract-") as directory:
        work = pathlib.Path(directory)
        values = [
            [100 + (i % 7) for i in range(12)],
            [96 + ((i * 3) % 11) for i in range(12)],
            [104 + ((i * 5) % 13) for i in range(12)],
        ]
        inputs = []
        for index, row in enumerate(values):
            path = work / f"normal-{index + 1}.tsv"
            write_counts(path, f"NORMAL_{index + 1}", row)
            inputs.append(path)
        annotated = work / "annotated-intervals.tsv"
        annotated_lines = ["@HD\tVN:1.5", "@SQ\tSN:chr1\tLN:100000", "CONTIG\tSTART\tEND\tGC_CONTENT"]
        for index in range(12):
            start = index * 100 + 1
            annotated_lines.append(f"chr1\t{start}\t{start + 99}\t{(index % 5) / 4.0:.6f}")
        annotated.write_text("\n".join(annotated_lines) + "\n", encoding="utf-8")
        output = work / "pon.hdf5"
        manifest = work / "pon.manifest.json"
        command = [str(binary)]
        for path in inputs:
            command += ["-I", str(path)]
        command += ["-O", str(output), "--output-manifest", str(manifest),
                    "--minimum-interval-median-percentile", "0",
                    "--maximum-zeros-in-sample-percentage", "100",
                    "--maximum-zeros-in-interval-percentage", "100",
                    "--extreme-sample-median-percentile", "0",
                    "--number-of-eigensamples", "2", "--maximum-chunk-size", "20",
                    "--annotated-intervals", str(annotated),
                    "--threads", "2"]
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        summary = json.loads(result.stdout.splitlines()[-1])
        if summary["panel_samples"] != 3 or summary["panel_intervals"] != 12 or summary["eigensamples"] != 2 or not summary["gc_content"]:
            raise AssertionError(summary)
        if not output.is_file() or output.stat().st_size == 0:
            raise AssertionError("native PoN output is empty")
        manifest_data = json.loads(manifest.read_text(encoding="utf-8"))
        if manifest_data["format"] != "HDF5-SVD-ReadCountPanelOfNormals-v7" or not manifest_data["gc_content"]:
            raise AssertionError(manifest_data)
        telemetry = manifest_data["telemetry"]
        if telemetry["kernel_lifecycle"] != (
                "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect"):
            raise AssertionError(telemetry)
        if telemetry["kernel_execution_space"] not in {"OpenMP", "Serial"}:
            raise AssertionError(telemetry)
        if telemetry["kernel_execution_policy"] != "MDRangePolicy+RangePolicy":
            raise AssertionError(telemetry)
        if telemetry["kernel_batches"] < 1 or telemetry["kernel_observations"] < 1:
            raise AssertionError(telemetry)
        if telemetry["kernel_prepare_seconds"] < 0.0 or telemetry["kernel_execute_seconds"] < 0.0:
            raise AssertionError(telemetry)

        native_standardized = work / "native-standardized.tsv"
        native_denoised = work / "native-denoised.tsv"
        native_manifest = work / "native-denoise.manifest.json"
        native_result = subprocess.run([
            str(root / "fastgatk-native/build/fastgatk-denoise-read-counts"),
            "-I", str(inputs[0]), "--count-panel-of-normals", str(output),
            "--standardized-copy-ratios", str(native_standardized),
            "--denoised-copy-ratios", str(native_denoised),
            "--number-of-eigensamples", "2", "--output-manifest", str(native_manifest),
        ], check=True, text=True, capture_output=True)
        native_summary = json.loads(native_result.stdout.splitlines()[-1])
        if native_summary["panel_samples"] != 3 or native_summary["intervals"] != 12 or not native_summary["svd"]:
            raise AssertionError(native_summary)
        if not native_standardized.is_file() or not native_denoised.is_file():
            raise AssertionError("native DenoiseReadCounts did not produce both outputs")

        java = root / "third_party/jdk17/bin/java"
        jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        java_checked = False
        if oracle_guard.oracle_ready('verify_create_read_count_panel_of_normals.py', java, jar):
            standardized = work / "java-standardized.tsv"
            denoised = work / "java-denoised.tsv"
            java_command = [str(java), "-jar", str(jar), "DenoiseReadCounts",
                            "-I", str(inputs[0]), "--count-panel-of-normals", str(output),
                            "--standardized-copy-ratios", str(standardized),
                            "--denoised-copy-ratios", str(denoised),
                            "--number-of-eigensamples", "2"]
            java_result = subprocess.run(java_command, text=True, capture_output=True)
            if java_result.returncode != 0:
                raise AssertionError("Java GATK could not read native HDF5 PoN:\n" + java_result.stderr[-4000:])
            java_checked = standardized.is_file() and denoised.is_file() and standardized.stat().st_size > 0 and denoised.stat().st_size > 0
            if not java_checked:
                raise AssertionError("Java DenoiseReadCounts did not produce both outputs")
            java_pon = work / "java-pon.hdf5"
            java_create = subprocess.run([
                str(java), "-jar", str(jar), "CreateReadCountPanelOfNormals",
                *sum((["-I", str(path)] for path in inputs), []),
                "-O", str(java_pon), "--annotated-intervals", str(annotated),
                "--minimum-interval-median-percentile", "0",
                "--maximum-zeros-in-sample-percentage", "100",
                "--maximum-zeros-in-interval-percentage", "100",
                "--extreme-sample-median-percentile", "0",
                "--number-of-eigensamples", "2", "--maximum-chunk-size", "20",
            ], text=True, capture_output=True)
            if java_create.returncode != 0:
                raise AssertionError("Java GATK CreateReadCountPanelOfNormals failed:\n" + java_create.stderr[-4000:])
            java_standardized = work / "java-standardized-from-java-pon.tsv"
            java_denoised = work / "java-denoised-from-java-pon.tsv"
            java_denoise = subprocess.run([
                str(java), "-jar", str(jar), "DenoiseReadCounts", "-I", str(inputs[0]),
                "--count-panel-of-normals", str(java_pon),
                "--standardized-copy-ratios", str(java_standardized),
                "--denoised-copy-ratios", str(java_denoised), "--number-of-eigensamples", "2",
            ], text=True, capture_output=True)
            if java_denoise.returncode != 0:
                raise AssertionError("Java GATK DenoiseReadCounts with Java PoN failed:\n" + java_denoise.stderr[-4000:])
            native_from_java = work / "native-standardized-from-java-pon.tsv"
            native_from_java_manifest = work / "native-from-java.manifest.json"
            native_from_java_result = subprocess.run([
                str(root / "fastgatk-native/build/fastgatk-denoise-read-counts"),
                "-I", str(inputs[0]), "--count-panel-of-normals", str(java_pon),
                "--standardized-copy-ratios", str(native_from_java), "--denoised-copy-ratios",
                str(work / "native-denoised-from-java-pon.tsv"), "--number-of-eigensamples", "2",
                "--output-manifest", str(native_from_java_manifest),
            ], text=True, capture_output=True)
            if native_from_java_result.returncode != 0:
                raise AssertionError("native DenoiseReadCounts could not read Java PoN:\n" + native_from_java_result.stderr[-4000:])
            def data_values(path: pathlib.Path) -> list[float]:
                values = []
                for line in path.read_text(encoding="utf-8").splitlines():
                    if not line or line.startswith("@") or line.startswith("CONTIG"):
                        continue
                    values.append(float(line.split("\t")[-1]))
                return values
            native_values = data_values(native_standardized)
            java_values = data_values(java_standardized)
            native_from_java_values = data_values(native_from_java)
            if len(native_from_java_values) != len(java_values) or any(abs(a - b) > 2.0e-5 for a, b in zip(native_from_java_values, java_values)):
                raise AssertionError({"native_from_java_pon": native_from_java_values, "java_standardized": java_values,
                                      "native_manifest": json.loads(native_from_java_manifest.read_text(encoding="utf-8"))})
        print(json.dumps({"status": "pass", "input_samples": 3, "intervals": 12,
                          "eigensamples": 2, "java_reader": java_checked,
                          "native_reader": True}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
