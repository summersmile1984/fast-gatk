#!/usr/bin/env python3
"""Benchmark the native FilterMutectCalls VCF/orientation-prior boundary."""

from __future__ import annotations

import gzip
import io
import json
import os
import subprocess
import tarfile
import tempfile
import time
from pathlib import Path


def make_prior(path: Path, sample: str) -> None:
    columns = [
        "context", "rev_comp", "f1r2_a", "f1r2_c", "f1r2_g", "f1r2_t",
        "f2r1_a", "f2r1_c", "f2r1_g", "f2r1_t", "hom_ref", "germline_het",
        "somatic_het", "hom_var", "num_examples", "num_alt_examples",
    ]
    complements = str.maketrans("ACGT", "TGCA")
    rows = []
    for left in "ACGT":
        for middle in "ACGT":
            for right in "ACGT":
                context = left + middle + right
                masses = [0.0] * 12
                masses[1] = 0.90
                masses[5] = 0.01
                masses[8:] = [0.01, 0.01, 0.03, 0.01]
                normalizer = sum(masses)
                rows.append("\t".join([
                    context, context.translate(complements)[::-1],
                    *[f"{value / normalizer:.12g}" for value in masses], "10", "10",
                ]))
    text = f"#<METADATA>SAMPLE={sample}\n" + "\t".join(columns) + "\n" + "\n".join(rows) + "\n"
    with tarfile.open(path, "w:gz") as archive:
        member = tarfile.TarInfo(f"{sample}.orientation_priors")
        member.size = len(text.encode("utf-8"))
        archive.addfile(member, io.BytesIO(text.encode("utf-8")))


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_FILTER_MUTECT_BINARY",
        str(root / "fastgatk-native/build/fastgatk-filter-mutect-calls"),
    ))
    with tempfile.TemporaryDirectory(prefix="fastgatk-filter-mutect-benchmark-") as directory:
        work = Path(directory)
        reference = work / "reference.fa"
        reference.write_text(">chr1\nAACGTT\n", encoding="utf-8")
        (work / "reference.fa.fai").write_text("chr1\t6\t6\t6\t7\n", encoding="utf-8")
        prior_tumor1 = work / "TUMOR1.orientation_priors.tar.gz"
        prior_tumor2 = work / "TUMOR2.orientation_priors.tar.gz"
        make_prior(prior_tumor1, "TUMOR1")
        make_prior(prior_tumor2, "TUMOR2")
        header = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=6>
##normal_sample=NORMAL
##INFO=<ID=TLOD,Number=A,Type=Float,Description=Tumor log odds>
##INFO=<ID=POPAF,Number=A,Type=Float,Description=Population negative log10 allele frequency>
##INFO=<ID=AF,Number=A,Type=Float,Description=Allele fraction>
##INFO=<ID=F1R2,Number=A,Type=Integer,Description=Forward support>
##INFO=<ID=R1F2,Number=A,Type=Integer,Description=Reverse support>
##INFO=<ID=NALOD,Number=A,Type=Float,Description=Negative log10 normal artifact odds>
##INFO=<ID=RPA,Number=R,Type=Integer,Description=Repeats per allele>
##INFO=<ID=RU,Number=1,Type=String,Description=Repeat unit>
##INFO=<ID=PARTIFACT,Number=1,Type=Float,Description=Artifact posterior>
##INFO=<ID=AS_SB_TABLE,Number=A,Type=String,Description=Allele-specific strand table>
##INFO=<ID=AS_UNIQ_ALT_READ_COUNT,Number=A,Type=String,Description=Unique alternate reads>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=F1R2,Number=A,Type=Integer,Description=Forward orientation alternate support>
##FORMAT=<ID=F2R1,Number=A,Type=Integer,Description=Reverse orientation alternate support>
##FORMAT=<ID=PGT,Number=1,Type=String,Description=Phased genotype>
##FORMAT=<ID=PID,Number=1,Type=String,Description=Phasing identifier>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR1\tTUMOR2\tNORMAL
"""
        body = "".join(
            f"chr1\t{i + 2}\t.\tA\tC\t.\tPASS\tTLOD=20;POPAF=2;NALOD={2 if i % 2 == 0 else -3};RPA=10,9;RU=A;PARTIFACT={0.8 if i % 2 == 0 else 0.2};AS_SB_TABLE=20,20|{100 if i % 2 == 0 else 1},{1 if i % 2 == 0 else 100};AS_UNIQ_ALT_READ_COUNT={100 if i % 2 == 0 else 1}\tGT:AD:F1R2:F2R1:PGT:PID"
            f"\t0/1:20,{100 if i % 2 == 0 else 1}:{100 if i % 2 == 0 else 1}:{0 if i % 2 == 0 else 0}:0|1:100_A_C"
            f"\t0/1:20,{1 if i % 2 == 0 else 100}:{1 if i % 2 == 0 else 50}:{1 if i % 2 == 0 else 50}:0|1:100_A_C"
            "\t0/0:100,0:0:0:.:.\n"
            for i in range(128)
        )
        input_path = work / "input.vcf.gz"
        with gzip.open(input_path, "wt", encoding="utf-8") as stream:
            stream.write(header + body)
        env = os.environ.copy()
        env.setdefault("OMP_PROC_BIND", "true")
        env.setdefault("OMP_PLACES", "threads")
        durations = []
        summary = None
        telemetry = None
        repeats = 3
        for repeat in range(repeats):
            output = work / f"filtered-{repeat}.vcf.gz"
            manifest = work / f"filtered-{repeat}.manifest.json"
            started = time.perf_counter()
            summary = json.loads(subprocess.check_output([
                str(binary), "-R", str(reference), "-V", str(input_path), "-O", str(output),
                "--orientation-bias-artifact-priors", str(prior_tumor1),
                "--orientation-bias-artifact-priors", str(prior_tumor2),
                "--max-orientation-artifact-probability", "0.5",
                "--contamination-estimate", "0.1",
                "--output-manifest", str(manifest),
            ], text=True, env=env).splitlines()[-1])
            telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
            durations.append(time.perf_counter() - started)
        assert summary is not None and summary["input_records"] == 128
        assert telemetry is not None and telemetry["haplotype_records"] > 0

        # Exercise the multi-ALT batch path separately.  The ordinary workload
        # is intentionally single-ALT and therefore uses the low-overhead
        # scalar fast path; this fixture proves that two concrete ALTs share
        # one tumor-evidence launch per record.
        multi_header = header.replace("length=6", "length=200")
        multi_body = "".join(
            f"chr1\t{index + 2}\t.\tA\tC,G\t.\tPASS\t"
            "TLOD=20,18;POPAF=2,2;AF=0.5,0.3;F1R2=20,20;R1F2=20,20;"
            "AS_SB_TABLE=20,20|20,20|20,20;AS_UNIQ_ALT_READ_COUNT=40,30\t"
            "GT:AD:F1R2:F2R1:PGT:PID"
            "\t0/1:20,40,20:20,20:20,20:0|1:200_A_C"
            "\t0/0:100,0,0:0,0:0,0:.:."
            "\t0/0:100,0,0:0,0:0,0:.:.\n"
            for index in range(64)
        )
        multi_input = work / "multi-alt-input.vcf.gz"
        with gzip.open(multi_input, "wt", encoding="utf-8") as stream:
            stream.write(multi_header + multi_body)
        multi_output = work / "multi-alt-output.vcf.gz"
        multi_manifest = work / "multi-alt-output.manifest.json"
        multi_summary = json.loads(subprocess.check_output([
            str(binary), "-V", str(multi_input), "-O", str(multi_output),
            "--output-manifest", str(multi_manifest),
        ], text=True, env=env).splitlines()[-1])
        multi_alt_telemetry = json.loads(multi_manifest.read_text(encoding="utf-8"))["telemetry"]
        assert multi_summary["input_records"] == 64
        assert multi_alt_telemetry["kernel_execution_policy"] == "RangePolicy"
        assert multi_alt_telemetry["kernel_batches"] > 0
        assert multi_alt_telemetry["kernel_observations"] >= 2 * multi_summary["input_records"]

        # A separate empirical-model workload measures the streaming GATK
        # SomaticClusteringModel path (quantile peak split + BIC + EM), rather
        # than only the ordinary file-boundary filters above.
        model_values = [5, 10, 15, 20, 25] + [60] * 90
        model_header = header.replace("length=6", "length=200")
        model_body = "".join(
            f"chr1\t{index}\t.\tA\tG\t.\tPASS\tTLOD=8;AF={af / 100:.2f};F1R2=20;R1F2=20\tGT:AD"
            f"\t0/1:{100 - af},{af}\t0/0:100,0\t0/0:100,0\n"
            for index, af in enumerate(model_values, start=1)
        )
        model_input = work / "empirical-model-input.vcf.gz"
        model_stats_input = work / "empirical-model-input.stats"
        with gzip.open(model_input, "wt", encoding="utf-8") as stream:
            stream.write(model_header + model_body)
        model_stats_input.write_text("statistic\tvalue\ncallable\t1000000\n", encoding="utf-8")
        model_durations = []
        model_stats = None
        model_repeats = 3
        for repeat in range(model_repeats):
            model_output = work / f"empirical-model-{repeat}.vcf.gz"
            model_stats_path = work / f"empirical-model-{repeat}.stats.json"
            started = time.perf_counter()
            subprocess.check_output([
                str(binary), "-V", str(model_input), "-O", str(model_output),
                "--stats", str(model_stats_input), "--filtering-stats", str(model_stats_path),
                "--threshold-strategy", "CONSTANT", "--initial-threshold", "0.1",
            ], text=True, env=env)
            model_stats = json.loads(model_stats_path.read_text(encoding="utf-8"))
            model_durations.append(time.perf_counter() - started)
        assert model_stats is not None and model_stats["empirical_somatic_model_records"] == 95
        assert model_stats["empirical_cluster_count"] >= 3
        model_durations.sort()
        model_p50 = model_durations[len(model_durations) // 2]
        durations.sort()
        p50 = durations[len(durations) // 2]
        print(json.dumps({
            "status": "pass", "tool": "FilterMutectCalls", "backend": "Kokkos",
            "repeats": repeats, "records": summary["input_records"],
            "pass_records": summary["pass_records"], "wall_seconds": sum(durations) / repeats,
            "p50_wall_seconds": p50, "records_per_second": summary["input_records"] / p50,
            "orientation_depth_disagreements": telemetry["orientation_depth_disagreements"] if telemetry else 0,
            "kernel_execution_space": telemetry.get("kernel_execution_space", "") if telemetry else "",
            "kernel_execution_policy": telemetry.get("kernel_execution_policy", "") if telemetry else "",
            "kernel_batches": telemetry.get("kernel_batches", 0) if telemetry else 0,
            "kernel_observations": telemetry.get("kernel_observations", 0) if telemetry else 0,
            "kernel_batch_calls": telemetry.get("kernel_batch_calls", 0) if telemetry else 0,
            "kernel_batch_observations": telemetry.get("kernel_batch_observations", 0) if telemetry else 0,
            "kernel_prepare_seconds": telemetry.get("kernel_prepare_seconds", 0.0) if telemetry else 0.0,
            "kernel_execute_seconds": telemetry.get("kernel_execute_seconds", 0.0) if telemetry else 0.0,
            "multi_alt_kernel_batches": multi_alt_telemetry.get("kernel_batches", 0),
            "multi_alt_kernel_observations": multi_alt_telemetry.get("kernel_observations", 0),
            "multi_alt_batch_calls": multi_alt_telemetry.get("kernel_batch_calls", 0),
            "multi_alt_batch_observations": multi_alt_telemetry.get("kernel_batch_observations", 0),
            "allele_specific_strand_records": telemetry["allele_specific_strand_records"] if telemetry else 0,
            "allele_specific_unique_records": telemetry["allele_specific_unique_records"] if telemetry else 0,
            "normal_artifact_records": telemetry["normal_artifact_records"] if telemetry else 0,
            "contamination_posterior_records": telemetry.get("contamination_posterior_records", 0) if telemetry else 0,
            "contamination_posterior_alleles": telemetry.get("contamination_posterior_alleles", 0) if telemetry else 0,
            "weak_evidence_records": telemetry.get("weak_evidence_records", 0) if telemetry else 0,
            "slippage_records": telemetry.get("slippage_records", 0) if telemetry else 0,
            "haplotype_records": telemetry["haplotype_records"] if telemetry else 0,
            "empirical_model_records": model_stats["empirical_somatic_model_records"],
            "empirical_cluster_count": model_stats["empirical_cluster_count"],
            "empirical_model_p50_wall_seconds": model_p50,
            "empirical_model_records_per_second": model_stats["empirical_somatic_model_records"] / model_p50,
        }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
