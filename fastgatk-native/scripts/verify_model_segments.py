#!/usr/bin/env python3
"""Contract checks for the TSV ModelSegments segmentation boundary."""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = pathlib.Path(os.environ.get(
    "FASTGATK_MODEL_SEGMENTS_BINARY",
    str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native" / "build"))
        / "fastgatk-model-segments")))
JAVA = ROOT / "third_party" / "jdk17" / "bin" / "java"
GATK_JAR = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"


def write_copy_ratios(path: pathlib.Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:GATKCopyNumber\tSM:SAMPLE\n")
        handle.write("CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO\n")
        values = [0.00, 0.02, 0.03, 0.80, 0.82, 0.81, -0.50]
        for index, value in enumerate(values):
            start = index * 10 + 1
            handle.write(f"chr1\t{start}\t{start + 9}\t1\t{value}\n")


def write_allelic(path: pathlib.Path) -> None:
    path.write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:GATKCopyNumber\tSM:SAMPLE\n"
        "CONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE\n"
        + "\n".join(f"chr1\t{index * 10 + 1}\t10\t10\tA\tC" for index in range(7))
        + "\n",
        encoding="utf-8",
    )


def rows(path: pathlib.Path) -> list[list[str]]:
    return [
        line.split("\t")
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("@") and not line.startswith("CONTIG")
    ]


def main() -> None:
    if not BINARY.exists():
        raise SystemExit(f"missing native binary: {BINARY}")
    with tempfile.TemporaryDirectory(prefix="fastgatk-model-segments-") as temporary:
        work = pathlib.Path(temporary)
        input_path = work / "denoised.tsv"
        allelic_path = work / "allelic.tsv"
        output_path = work / "segments.tsv"
        manifest_path = work / "segments.json"
        write_copy_ratios(input_path)
        write_allelic(allelic_path)
        result = subprocess.run(
            [str(BINARY), "--denoised-copy-ratios", str(input_path),
             "--allelic-counts", str(allelic_path), "--output-prefix", str(work / "sample"),
             "--output-manifest", str(manifest_path), "--change-point-threshold", "0.15"],
            text=True, capture_output=True, check=True,
        )
        assert output_path.exists() is False  # output-prefix derives sample.modelFinal.segments.tsv
        derived = work / "sample.modelFinal.segments.tsv"
        segment_rows = rows(derived)
        assert [int(row[3]) for row in segment_rows] == [3, 3, 1], segment_rows
        assert [row[0:3] for row in segment_rows] == [
            ["chr1", "1", "30"], ["chr1", "31", "60"], ["chr1", "61", "70"]
        ], segment_rows
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        assert manifest["tool"] == "ModelSegments"
        assert manifest["status"] == "prototype"
        assert manifest["input_points"] == 7
        assert manifest["segments"] == 3
        assert manifest["allelic_loci"] == 7
        assert manifest["allelic_signal_used"] is True
        assert manifest["execution_policy"] == "RangePolicy"
        assert manifest["determinism"] == "strict"
        telemetry = manifest["telemetry"]
        assert telemetry["output_bytes"] > 0
        assert telemetry["sidecar_bytes"] > 0
        assert telemetry["wall_seconds"] > 0
        assert telemetry["kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect")
        assert telemetry["kernel_execution_policy"] == "RangePolicy"
        assert telemetry["execution_space"] in ("OpenMP", "Serial")
        assert telemetry["kernel_batches"] >= 1
        assert telemetry["kernel_observations"] >= 1
        assert telemetry["kernel_prepare_seconds"] >= 0.0
        assert telemetry["kernel_execute_seconds"] >= 0.0
        header = next(line for line in derived.read_text(encoding="utf-8").splitlines()
                      if line.startswith("CONTIG\t"))
        assert "NUM_POINTS_ALLELE_FRACTION" in header
        assert len(segment_rows[0]) == 12, segment_rows
        for suffix in (".modelFinal.seg", ".modelBegin.seg", ".cr.seg", ".cr.igv.seg", ".af.igv.seg",
                       ".modelBegin.cr.param", ".modelBegin.af.param",
                       ".modelFinal.cr.param", ".modelFinal.af.param"):
            assert (work / f"sample{suffix}").is_file(), suffix
        assert len(manifest["sidecars"]) == 9
        for suffix in (".modelBegin.cr.param", ".modelBegin.af.param",
                       ".modelFinal.cr.param", ".modelFinal.af.param"):
            parameter_lines = (work / f"sample{suffix}").read_text(encoding="utf-8").splitlines()
            assert any(line.startswith("PARAMETER_NAME\t") for line in parameter_lines)
        assert json.loads(result.stdout)["segments"] == 3

        # Allelic imbalance alone must be able to split a copy-ratio-flat
        # interval.  This exercises the optional Kokkos allele-fraction signal
        # rather than merely validating that loci exist.
        flat_input = work / "flat.denoised.tsv"
        flat_input.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"
            "CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO\n"
            + "\n".join(f"chr1\t{index * 10 + 1}\t{index * 10 + 10}\t1\t0.0" for index in range(6))
            + "\n", encoding="utf-8")
        flat_allelic = work / "flat.allelic.tsv"
        flat_allelic.write_text(
            "@HD\tVN:1.6\nCONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE\n"
            + "\n".join(
                f"chr1\t{index * 10 + 1}\t10\t{10 if index < 3 else 8}\tA\tC"
                for index in range(6)) + "\n", encoding="utf-8")
        flat_output = work / "flat.modelFinal.segments.tsv"
        flat_manifest = work / "flat.json"
        subprocess.run(
            [str(BINARY), "-I", str(flat_input), "--allelic-counts", str(flat_allelic),
             "-O", str(flat_output), "--output-manifest", str(flat_manifest),
             "--change-point-threshold", "0.15", "--allelic-change-point-threshold", "0.05"],
            text=True, capture_output=True, check=True,
        )
        flat_rows = rows(flat_output)
        assert [int(row[3]) for row in flat_rows] == [3, 3], flat_rows
        assert json.loads(flat_manifest.read_text(encoding="utf-8"))["allelic_signal_used"] is True

        # Explicit kernel controls select the deterministic Kokkos
        # KernelSegmenter path.  The two-level synthetic profile should yield
        # at most the requested three segments and must be repeatable.
        kernel_output = work / "kernel.modelFinal.segments.tsv"
        kernel_manifest = work / "kernel.json"
        kernel_args = [str(BINARY), "-I", str(input_path), "--allelic-counts", str(allelic_path),
                       "-O", str(kernel_output), "--output-manifest", str(kernel_manifest),
                       "--maximum-number-of-segments-per-chromosome", "3", "--window-size", "2",
                       "--kernel-approximation-dimension", "7", "--number-of-changepoints-penalty-factor", "0"]
        first_kernel = subprocess.run(kernel_args, text=True, capture_output=True, check=True)
        first_rows = rows(kernel_output)
        assert 1 <= len(first_rows) <= 3, first_rows
        kernel_manifest_data = json.loads(kernel_manifest.read_text(encoding="utf-8"))
        assert kernel_manifest_data["segmentation_method"] == "KernelSegmenter"
        assert kernel_manifest_data["execution_policy"] == "RangePolicy"
        first_bytes = kernel_output.read_bytes()
        subprocess.run(kernel_args, text=True, capture_output=True, check=True)
        assert kernel_output.read_bytes() == first_bytes
        assert json.loads(first_kernel.stdout)["segmentation_method"] == "KernelSegmenter"

        output_dir = work / "gatk-output"
        directory_result = subprocess.run(
            [str(BINARY), "--denoised-copy-ratios", str(input_path),
             "--output-prefix", "sample", "-O", str(output_dir),
             # The native threshold path is an explicit compatibility
             # extension; GATK's default without --segments is KernelSegmenter.
             "--change-point-threshold", "0.15"],
            text=True, capture_output=True, check=True,
        )
        assert (output_dir / "sample.modelFinal.segments.tsv").is_file()
        assert (output_dir / "sample.modelFinal.seg").is_file()
        assert (output_dir / "sample.cr.seg").is_file()
        assert json.loads(directory_result.stdout)["segments"] == 3

        # Matched-normal genotyping follows GATK's beta-integrated
        # homozygous-vs-heterozygous log-ratio rule.  Homozygous normal sites
        # are removed from both the normal and case het sidecars, while
        # interval-contained heterozygous sites are retained.
        normal_allelic = work / "normal.allelic.tsv"
        normal_allelic.write_text(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:GATKCopyNumber\tSM:NORMAL\n"
            "CONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE\n"
            + "\n".join(
                f"chr1\t{index * 10 + 1}\t10\t{0 if index in (2, 5) else (8 if index == 3 else 10)}\tA\tC"
                for index in range(7)) + "\n", encoding="utf-8")
        normal_manifest = work / "normal.json"
        normal_args = [str(BINARY), "-I", str(input_path), "--allelic-counts", str(allelic_path),
                       "--normal-allelic-counts", str(normal_allelic), "--output-prefix", "normal",
                       "-O", str(work), "--output-manifest", str(normal_manifest),
                       "--minimum-total-allele-count-normal", "10"]
        subprocess.run(normal_args, text=True, capture_output=True, check=True)
        normal_data = json.loads(normal_manifest.read_text(encoding="utf-8"))
        assert normal_data["matched_normal_mode"] is True
        assert normal_data["heterozygous_allelic_loci"] == 5, normal_data
        assert normal_data["genotyping_execution_policy"] == "RangePolicy"
        assert normal_data["genotyping_execution_space"] in ("OpenMP", "Serial")
        normal_hets = rows(work / "normal.hets.tsv")
        normal_hets_sidecar = rows(work / "normal.hets.normal.tsv")
        assert len(normal_hets) == 5 and len(normal_hets_sidecar) == 5
        assert all(int(row[3]) > 0 for row in normal_hets_sidecar)

        # Compare the heterozygous-site hand-off itself with the pinned GATK
        # 4.6.2.0 implementation.  The full Java MCMC model remains a
        # separate fallback, but this oracle proves the input filtering and
        # sidecar contract rather than only comparing a synthetic count.
        assert JAVA.exists() and GATK_JAR.exists()
        java_copy = work / "java.denoised.tsv"
        java_copy_lines = []
        for line in input_path.read_text(encoding="utf-8").splitlines():
            if line.startswith("CONTIG\t"):
                java_copy_lines.append("CONTIG\tSTART\tEND\tLOG2_COPY_RATIO")
            elif line and not line.startswith("@"):
                fields = line.split("\t")
                java_copy_lines.append("\t".join(fields[:3] + [fields[-1]]))
            else:
                java_copy_lines.append(line)
        java_copy.write_text("\n".join(java_copy_lines) + "\n", encoding="utf-8")
        java_output = work / "java-output"
        java_result = subprocess.run(
            [str(JAVA), "-jar", str(GATK_JAR), "ModelSegments",
             "--denoised-copy-ratios", str(java_copy), "--allelic-counts", str(allelic_path),
             "--normal-allelic-counts", str(normal_allelic),
             "--minimum-total-allele-count-normal", "10",
             "--number-of-samples-copy-ratio", "2",
             "--number-of-burn-in-samples-copy-ratio", "1",
             "--number-of-samples-allele-fraction", "2",
             "--number-of-burn-in-samples-allele-fraction", "1",
             "--number-of-smoothing-iterations-per-fit", "0",
             "--output-prefix", "normal", "-O", str(java_output)],
            text=True, capture_output=True, check=False,
        )
        assert java_result.returncode == 0, java_result.stderr
        assert (work / "normal.hets.tsv").read_text(encoding="utf-8") == \
            (java_output / "normal.hets.tsv").read_text(encoding="utf-8")
        assert (work / "normal.hets.normal.tsv").read_text(encoding="utf-8") == \
            (java_output / "normal.hets.normal.tsv").read_text(encoding="utf-8")

        # Repeatable denoised/allelic inputs select the bounded native
        # multisample segmentation path and emit a Picard interval list.
        multi_manifest = work / "multi.json"
        multi_args = [str(BINARY), "--denoised-copy-ratios", str(input_path),
                      "--denoised-copy-ratios", str(input_path),
                      "--allelic-counts", str(allelic_path),
                      "--allelic-counts", str(allelic_path),
                      "--output-prefix", "multi", "-O", str(work),
                      "--change-point-threshold", "0.15",
                      "--output-manifest", str(multi_manifest)]
        multi_result = subprocess.run(multi_args, text=True, capture_output=True, check=True)
        multi_output = work / "multi.interval_list"
        assert multi_output.is_file()
        multi_intervals = [line for line in multi_output.read_text(encoding="utf-8").splitlines()
                           if line and not line.startswith("@")]
        assert len(multi_intervals) == 3, multi_intervals
        multi_manifest_data = json.loads(multi_manifest.read_text(encoding="utf-8"))
        assert multi_manifest_data["multisample"] is True
        assert multi_manifest_data["sample_count"] == 2
        assert multi_manifest_data["common_points"] == 7
        assert len(multi_manifest_data["inputs"]) == 2
        assert multi_manifest_data["determinism"] == "strict"
        multi_telemetry = multi_manifest_data["telemetry"]
        assert multi_telemetry["output_bytes"] > 0
        assert multi_telemetry["wall_seconds"] > 0
        assert multi_telemetry["kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect")
        assert multi_telemetry["kernel_batches"] >= 1
        assert multi_telemetry["kernel_observations"] >= 1
        assert json.loads(multi_result.stdout)["multisample"] is True

        # The bounded native posterior sampler is deterministic across runs
        # and consumes the GATK-compatible probability controls.  It is not a
        # claim of Java GATK's release-specific posterior bit identity; the
        # manifest makes the native sampler explicit for workflow auditing.
        posterior_output = work / "posterior.modelFinal.segments.tsv"
        posterior_manifest = work / "posterior.json"
        posterior_args = [str(BINARY), "-I", str(input_path),
                          "--allelic-counts", str(allelic_path), "-O", str(posterior_output),
                          "--output-manifest", str(posterior_manifest), "--mode", "BOTH",
                          "--change-point-threshold", "0.15",
                          "--num-burn-in-iterations", "16", "--num-samples", "64"]
        posterior_result = subprocess.run(posterior_args, text=True, capture_output=True, check=True)
        posterior_bytes = posterior_output.read_bytes()
        posterior_rows = rows(posterior_output)
        assert len(posterior_rows) == 3 and all(len(row) == 12 for row in posterior_rows)
        for row in posterior_rows:
            assert all(float(row[index]) == float(row[index]) for index in range(5, 12))
        posterior_manifest_data = json.loads(posterior_manifest.read_text(encoding="utf-8"))
        assert posterior_manifest_data["mcmc"] is True
        assert posterior_manifest_data["posterior_sampler"] == "deterministic-gibbs-conditional-random-walk-v2"
        assert posterior_manifest_data["copy_ratio_conditional_model"] == "deterministic-gibbs-responsibility-v1"
        assert posterior_manifest_data["copy_ratio_conditional_iterations"] >= 4
        assert posterior_manifest_data["copy_ratio_variance_conditional"] >= 1.0e-6
        assert 0.0 < posterior_manifest_data["copy_ratio_outlier_probability_conditional"] < 1.0
        assert posterior_manifest_data["allele_fraction_conditional_model"] == (
            "deterministic-binomial-responsibility-v1")
        assert posterior_manifest_data["allele_fraction_conditional_iterations"] >= 4
        assert posterior_manifest_data["allele_fraction_mean_bias_conditional"] >= 0.0
        assert 0.0 <= posterior_manifest_data["allele_fraction_bias_variance_conditional"] <= 0.5
        assert 0.0 <= posterior_manifest_data["allele_fraction_outlier_probability_conditional"] <= 0.15
        assert posterior_manifest_data["posterior_sample_block_size"] == 4096
        assert posterior_manifest_data["posterior_reservoir_capacity"] == 2048
        assert posterior_manifest_data["num_burn_in_iterations"] == 16
        assert posterior_manifest_data["num_samples"] == 64
        assert posterior_manifest_data["mcmc_parameter_surface"] == "gatk-compatible-bounded-shared-chain"
        assert posterior_manifest_data["mcmc_controls_provenance"] == "native-generic"
        assert posterior_manifest_data["native_chain_burn_in_iterations"] == 16
        assert posterior_manifest_data["native_chain_retained_samples"] == 64
        assert posterior_manifest_data["native_chain_total_iterations"] == 80
        assert posterior_manifest_data["posterior_rng"] == "splitmix64-box-muller-fixed-seed-1216"
        assert posterior_manifest_data["posterior_report_semantics"] == (
            "bounded-native-random-walk-quantiles")
        assert posterior_manifest_data["parameter_report_semantics"] == (
            "conditional-point-estimate-repeated-deciles")
        assert posterior_manifest_data["parameter_report_distinct_draws"] == 1
        assert posterior_manifest_data["model_begin_final_segment_reports_distinct"] is False
        assert posterior_manifest_data["model_begin_final_parameter_reports_distinct"] is False
        assert posterior_manifest_data["gatk_full_gibbs_slice_equivalent"] is False
        assert posterior_manifest_data["java_rng_bit_identity"] is False
        assert posterior_manifest_data["minor_allele_fraction_prior_alpha"] == 25.0
        assert posterior_manifest_data["smoothing_applied"] is False
        posterior_telemetry = posterior_manifest_data["telemetry"]
        assert posterior_telemetry["kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect")
        assert posterior_telemetry["kernel_execution_policy"] == "RangePolicy"
        assert posterior_telemetry["kernel_batches"] >= 1
        assert posterior_telemetry["kernel_observations"] >= 1
        subprocess.run(posterior_args, text=True, capture_output=True, check=True)
        assert posterior_output.read_bytes() == posterior_bytes
        assert json.loads(posterior_result.stdout)["mcmc"] is True

        # Standard GATK ModelSegments modeling names are accepted alongside
        # the native compact controls and are recorded without silently
        # dropping per-model values.
        gatk_surface_output = work / "gatk-surface.modelFinal.segments.tsv"
        gatk_surface_manifest = work / "gatk-surface.json"
        subprocess.run(
            [str(BINARY), "-I", str(input_path), "--allelic-counts", str(allelic_path),
             "-O", str(gatk_surface_output), "--output-manifest", str(gatk_surface_manifest),
             "--change-point-threshold", "0.15",
             "--number-of-samples-copy-ratio", "8",
             "--number-of-burn-in-samples-copy-ratio", "3",
             "--number-of-samples-allele-fraction", "6",
             "--number-of-burn-in-samples-allele-fraction", "2",
             "--minor-allele-fraction-prior-alpha", "30",
             "--maximum-number-of-smoothing-iterations", "0",
             "--number-of-smoothing-iterations-per-fit", "0",
             "--smoothing-credible-interval-threshold-copy-ratio", "1.5",
             "--smoothing-credible-interval-threshold-allele-fraction", "1.25"],
            text=True, capture_output=True, check=True,
        )
        surface_manifest = json.loads(gatk_surface_manifest.read_text(encoding="utf-8"))
        assert surface_manifest["num_samples"] == 8
        assert surface_manifest["num_burn_in_iterations"] == 3
        assert surface_manifest["number_of_samples_copy_ratio"] == 8
        assert surface_manifest["number_of_samples_allele_fraction"] == 6
        assert surface_manifest["minor_allele_fraction_prior_alpha"] == 30.0
        assert surface_manifest["smoothing_applied"] is False
        print(json.dumps({
            "status": "pass",
            "segments": 3,
            "points": 7,
            "allelic_validation": True,
            "kernel_segmentation": True,
            "multisample_segmentation": True,
            "probabilistic_controls": True,
            "posterior_repeatable": True,
            "gatk_modeling_surface": True,
        }, sort_keys=True))


if __name__ == "__main__":
    main()
