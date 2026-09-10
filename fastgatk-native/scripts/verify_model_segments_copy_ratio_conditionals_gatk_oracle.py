#!/usr/bin/env python3
"""Pinned GATK oracle for the CopyRatioModeller conditional-model slice.

The native path intentionally does not claim Java's release-specific MCMC
draws are bit-identical.  This gate therefore compares the stable modeling
partition/schema with GATK 4.6.2.0 and checks the native conditional state on a
fixture containing one clear copy-ratio outlier.  It catches regressions where
the Gibbs conditional prepass is skipped, reports non-finite global
parameters, or lets the outlier drag the segment posterior mean to the wrong
mode.
"""

from __future__ import annotations

import json
import math
import os
import pathlib
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
BUILD = pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native/build"))
NATIVE = pathlib.Path(os.environ.get(
    "FASTGATK_MODEL_SEGMENTS_BINARY", str(BUILD / "fastgatk-model-segments")))
JAVA = ROOT / "third_party/jdk17/bin/java"
JAR = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=True)


def write_copy_ratios(path: pathlib.Path) -> None:
    lines = [
        "@HD\tVN:1.6",
        "@SQ\tSN:chr1\tLN:1000",
        "@RG\tID:GATKCopyNumber\tSM:SAMPLE",
        "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO",
    ]
    # The fourth value in the first supplied segment is an intentional
    # outlier.  GATK's model assigns by midpoint and estimates a global
    # variance/outlier probability before sampling segment means.
    for index, value in enumerate((0.00, 0.01, 0.02, 4.00, 0.80, 0.82, 0.81, 0.79)):
        lines.append(f"chr1\t{index * 10 + 1}\t{index * 10 + 10}\t{value}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_segments(path: pathlib.Path) -> None:
    path.write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"
        "chr1\t1\t40\t+\tseg1\nchr1\t41\t80\t+\tseg2\n",
        encoding="utf-8",
    )


def write_allelic_counts(path: pathlib.Path) -> None:
    lines = [
        "@HD\tVN:1.6",
        "@SQ\tSN:chr1\tLN:1000",
        "@RG\tID:GATKAllelicCounts\tSM:SAMPLE",
        "CONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE",
    ]
    counts = ((50, 50), (45, 55), (48, 52), (80, 20),
              (50, 50), (51, 49), (49, 51), (50, 50))
    for index, (ref, alt) in enumerate(counts):
        lines.append(f"chr1\t{index * 10 + 5}\t{ref}\t{alt}\tA\tC")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def data_rows(path: pathlib.Path) -> list[list[str]]:
    return [
        line.split("\t")
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("@") and not line.startswith("CONTIG")
    ]


def parameter_deciles(path: pathlib.Path, parameter_name: str) -> list[float]:
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if fields and fields[0] == parameter_name:
            return [float(value) for value in fields[1:]]
    raise AssertionError(f"missing parameter {parameter_name} in {path}")


def main() -> int:
    for path in (NATIVE, JAVA, JAR):
        if not path.exists():
            raise SystemExit(f"missing ModelSegments conditional oracle input: {path}")

    with tempfile.TemporaryDirectory(prefix="fastgatk-model-segments-conditionals-") as directory:
        work = pathlib.Path(directory)
        copy_ratios = work / "copy-ratios.tsv"
        segments = work / "segments.interval_list"
        allelic_counts = work / "allelic-counts.tsv"
        write_copy_ratios(copy_ratios)
        write_segments(segments)
        write_allelic_counts(allelic_counts)

        native_dir = work / "native"
        java_dir = work / "java"
        native_dir.mkdir()
        java_dir.mkdir()
        native_manifest = work / "native.json"
        native_result = run([
            str(NATIVE), "-I", str(copy_ratios), "--segments", str(segments),
            "--output-prefix", "sample", "-O", str(native_dir),
            "--output-manifest", str(native_manifest), "--allelic-counts", str(allelic_counts),
            "--mode", "BOTH",
            "--number-of-samples-copy-ratio", "12",
            "--number-of-burn-in-samples-copy-ratio", "6",
        ])
        java_result = run([
            str(JAVA), "-jar", str(JAR), "ModelSegments",
            "--denoised-copy-ratios", str(copy_ratios), "--segments", str(segments),
            "--allelic-counts", str(allelic_counts),
            "--output-prefix", "sample", "-O", str(java_dir),
            "--number-of-samples-copy-ratio", "12",
            "--number-of-burn-in-samples-copy-ratio", "6",
            "--number-of-samples-allele-fraction", "2",
            "--number-of-burn-in-samples-allele-fraction", "1",
            "--maximum-number-of-smoothing-iterations", "0",
            "--number-of-smoothing-iterations-per-fit", "0",
        ])

        native_prefix = native_dir / "sample"
        java_prefix = java_dir / "sample"
        native_final = data_rows(native_prefix.with_suffix(".modelFinal.seg"))
        java_final = data_rows(java_prefix.with_suffix(".modelFinal.seg"))
        expected_intervals = [["chr1", "1", "40"], ["chr1", "41", "80"]]
        assert [row[:3] for row in native_final] == expected_intervals, native_final
        assert [row[:3] for row in java_final] == expected_intervals, java_final
        assert [int(row[3]) for row in native_final] == [4, 4], native_final
        assert [int(row[3]) for row in java_final] == [4, 4], java_final
        assert [int(row[4]) for row in native_final] == [4, 4], native_final
        assert [int(row[4]) for row in java_final] == [4, 4], java_final

        manifest = json.loads(native_manifest.read_text(encoding="utf-8"))
        assert manifest["posterior_sampler"] == "deterministic-gibbs-conditional-random-walk-v2", manifest
        # GATK's family-specific number-of-samples includes burn-in.  Native's
        # compact shared-chain control retains num_samples *after* burn-in.
        # Lock both interpretations in provenance rather than implying the
        # projected chain is the Java Gibbs/minibatch-slice chain.
        assert manifest["mcmc_controls_provenance"] == "gatk-family-max-projection", manifest
        assert manifest["copy_ratio_requested_total_samples"] == 12, manifest
        assert manifest["copy_ratio_requested_burn_in_samples"] == 6, manifest
        assert manifest["copy_ratio_requested_retained_samples"] == 6, manifest
        assert manifest["allele_fraction_requested_total_samples"] == 100, manifest
        assert manifest["allele_fraction_requested_burn_in_samples"] == 50, manifest
        assert manifest["allele_fraction_requested_retained_samples"] == 50, manifest
        assert manifest["native_chain_burn_in_iterations"] == 50, manifest
        assert manifest["native_chain_retained_samples"] == 100, manifest
        assert manifest["native_chain_total_iterations"] == 150, manifest
        assert manifest["posterior_rng"] == "splitmix64-box-muller-fixed-seed-1216", manifest
        assert manifest["posterior_report_semantics"] == "bounded-native-random-walk-quantiles", manifest
        assert manifest["parameter_report_semantics"] == (
            "conditional-point-estimate-repeated-deciles"), manifest
        assert manifest["parameter_report_distinct_draws"] == 1, manifest
        assert manifest["model_begin_final_segment_reports_distinct"] is False, manifest
        assert manifest["model_begin_final_parameter_reports_distinct"] is False, manifest
        assert manifest["gatk_full_gibbs_slice_equivalent"] is False, manifest
        assert manifest["java_rng_bit_identity"] is False, manifest
        assert manifest["copy_ratio_conditional_model"] == "deterministic-gibbs-responsibility-v1", manifest
        assert manifest["copy_ratio_conditional_iterations"] >= 4, manifest
        assert manifest["allele_fraction_conditional_model"] == (
            "deterministic-binomial-responsibility-v1"), manifest
        assert manifest["allele_fraction_conditional_iterations"] >= 4, manifest
        assert math.isfinite(manifest["allele_fraction_mean_bias_conditional"])
        assert 0.0 <= manifest["allele_fraction_bias_variance_conditional"] <= 0.5
        assert 0.0 <= manifest["allele_fraction_outlier_probability_conditional"] <= 0.15
        assert math.isfinite(manifest["copy_ratio_variance_conditional"])
        assert manifest["copy_ratio_variance_conditional"] >= 1.0e-6
        assert 0.001 <= manifest["copy_ratio_outlier_probability_conditional"] <= 0.999
        # The outlier is assigned to an explicit uniform component; the
        # segment posterior remains near the inlier mode rather than 4.0.
        assert float(native_final[0][7]) < 1.0, native_final
        assert all(math.isfinite(float(row[index])) for row in native_final for index in (8, 9, 10)), native_final

        native_parameter = native_prefix.with_suffix(".modelFinal.cr.param").read_text(encoding="utf-8")
        native_allele_parameter = native_prefix.with_suffix(".modelFinal.af.param").read_text(encoding="utf-8")
        java_parameter = java_prefix.with_suffix(".modelFinal.cr.param").read_text(encoding="utf-8")
        java_allele_parameter = java_prefix.with_suffix(".modelFinal.af.param").read_text(encoding="utf-8")
        assert "VARIANCE\t" in native_parameter and "OUTLIER_PROBABILITY\t" in native_parameter
        assert "VARIANCE\t" in java_parameter and "OUTLIER_PROBABILITY\t" in java_parameter
        assert "MEAN_BIAS\t" in native_allele_parameter and "BIAS_VARIANCE\t" in native_allele_parameter
        assert "OUTLIER_PROBABILITY\t" in native_allele_parameter
        assert "MEAN_BIAS\t" in java_allele_parameter and "BIAS_VARIANCE\t" in java_allele_parameter
        assert "OUTLIER_PROBABILITY\t" in java_allele_parameter
        native_variance_deciles = parameter_deciles(
            native_prefix.with_suffix(".modelFinal.cr.param"), "VARIANCE")
        java_variance_deciles = parameter_deciles(
            java_prefix.with_suffix(".modelFinal.cr.param"), "VARIANCE")
        assert len(native_variance_deciles) == len(java_variance_deciles) == 9
        assert len(set(native_variance_deciles)) == 1, native_variance_deciles
        # Six post-burn-in Java draws are enough for this noisy/outlier fixture
        # to make the real posterior report non-degenerate.  This strict pinned
        # contrast prevents native point summaries from being mislabeled as
        # Java posterior deciles.
        assert len(set(java_variance_deciles)) > 1, java_variance_deciles
        assert java_result.returncode == 0 and native_result.returncode == 0

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "conditional_model": manifest["copy_ratio_conditional_model"],
            "conditional_iterations": manifest["copy_ratio_conditional_iterations"],
            "variance_finite": True,
            "outlier_probability_finite": True,
            "outlier_responsibility_active": True,
            "allele_fraction_conditional": True,
            "partition_and_schema_exact": True,
            "posterior_report_provenance_explicit": True,
            "java_variance_report_non_degenerate": True,
            "native_parameter_report_point_estimate": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
