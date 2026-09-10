#!/usr/bin/env python3
"""Exercise the native BaseRecalibrator -> ApplyBQSR compatibility loop."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from collections import Counter
from pathlib import Path


def run_json(command: list[str], extra_env: dict[str, str] | None = None) -> dict:
    env = os.environ.copy()
    env.setdefault("OMP_PROC_BIND", "true")
    env.setdefault("OMP_PLACES", "threads")
    if extra_env:
        env.update(extra_env)
    result = subprocess.run(command, check=True, text=True, capture_output=True, env=env)
    return json.loads(result.stdout.strip().splitlines()[-1])


def table_rows(text: str, table: str) -> list[list[str]]:
    marker = f"#:GATKTable:{table}:"
    lines = text.splitlines()
    try:
        start = next(i for i, line in enumerate(lines) if line.startswith(marker)) + 2
    except StopIteration:
        return []
    rows: list[list[str]] = []
    for line in lines[start:]:
        if line.startswith("#:GATKTable:"):
            break
        if line.strip():
            rows.append(line.split())
    return rows


def run_command(command: list[str]) -> None:
    env = os.environ.copy()
    env.setdefault("OMP_PROC_BIND", "true")
    env.setdefault("OMP_PLACES", "threads")
    subprocess.run(command, check=True, text=True, capture_output=True, env=env)


def sam_record_map(path: Path) -> dict[str, tuple[str, str, str, str]]:
    """Read stable SAM record fields after Java GATK PrintReads conversion."""
    result: dict[str, tuple[str, str, str, str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            raise AssertionError(f"malformed SAM record in {path}: {line[:120]}")
        # ApplyBQSR must preserve alignment and sequence fields; only QUAL is
        # recalibrated.  Compare the entire stable record core rather than a
        # quality sum so a compensating per-base error cannot pass silently.
        result[fields[0]] = (fields[1], fields[5], fields[9], fields[10])
    return result


def sam_record_lines(path: Path) -> list[str]:
    """Return records in writer order, including optional tags such as OQ."""
    return [line for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("@")] 


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    # Use the materialized chr17 fixture rather than the repository's LFS
    # pointer for human_g1k_v37.20.21.fasta; BQSR's reference-base and
    # low-quality skip semantics must be exercised, not silently treated as
    # an unavailable reference.
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    smoke = native / "fastgatk-hc-smoke"
    with tempfile.TemporaryDirectory(prefix="fastgatk-bqsr-") as directory:
        work = Path(directory)
        report = work / "recal.tsv"
        report_manifest = work / "recal.tsv.manifest.json"
        covariate_table = work / "recal.tsv.covariates.tsv"
        checkpoint = work / "recal.checkpoint"
        checkpoint_report = work / "checkpoint-recal.tsv"
        checkpoint_manifest = work / "checkpoint-recal.tsv.manifest.json"
        resumed_report = work / "resumed-recal.tsv"
        resumed_manifest = work / "resumed-recal.tsv.manifest.json"
        known_sites = work / "known-sites.vcf"
        known_sites_second = work / "known-sites-second.vcf"
        known_report = work / "known-recal.tsv"
        multi_known_report = work / "multi-known-recal.tsv"
        multi_known_manifest = work / "multi-known-recal.tsv.manifest.json"
        full_report = work / "full-recal.tsv"
        indel_report = work / "indel-recal.tsv"
        known_manifest = work / "known-recal.tsv.manifest.json"
        merged_report = work / "merged-recal.tsv"
        merged_manifest = work / "merged-recal.tsv.manifest.json"
        recal_bam = work / "recal.bam"
        apply_manifest = work / "recal.bam.manifest.json"
        gatk_no_index_bam = work / "gatk-no-index.bam"
        native_no_index_bam = work / "native-no-index.bam"
        gatk_no_index_sam = work / "gatk-no-index.sam"
        native_no_index_sam = work / "native-no-index.sam"
        apply_checkpoint = work / "apply.checkpoint"
        apply_checkpoint_bam = work / "recal-checkpoint.bam"
        apply_checkpoint_manifest = work / "recal-checkpoint.bam.manifest.json"
        apply_resumed_bam = work / "recal-resumed.bam"
        apply_resumed_manifest = work / "recal-resumed.bam.manifest.json"
        report_without_sidecar = work / "gatk-like-recal.tsv"
        recal_bam_without_sidecar = work / "recal-no-sidecar.bam"

        base = run_json([
            str(native / "fastgatk-bqsr"), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-69100", "--batch-records", "3", "-O", str(report),
            "--output-manifest", str(report_manifest),
        ])
        assert base["tool"] == "BaseRecalibrator" and base["records"] == 7
        assert base["covariate_observations"] > 0
        report_text = report.read_text(encoding="utf-8")
        assert report_text.startswith("#:GATKReport.v1.1:5\n")
        assert "#:GATKTable:RecalTable0:" in report_text
        quantized = table_rows(report_text, "Quantized")
        rows = table_rows(report_text, "RecalTable1")
        assert len(quantized) == 94 and any(int(row[1]) > 0 for row in quantized)
        assert rows and any(row[2] == "M" and int(row[4]) > 0 for row in rows)
        baseline_observations = sum(int(row[4]) for row in rows if row[2] == "M")
        covariate_rows = [line for line in covariate_table.read_text(encoding="utf-8").splitlines()
                          if line and not line.startswith("#")]
        assert covariate_rows and len(covariate_rows) == base["covariate_rows"]
        metadata = json.loads(report_manifest.read_text(encoding="utf-8"))
        assert metadata["execution_space"] in {"OpenMP", "Serial", "Host"}
        assert metadata["determinism"] == "strict"
        assert metadata["primary_output_kind"] == "recalibration-report"
        assert metadata["compatibility"]["integer_table"] is True
        assert metadata["compatibility"]["covariate_table"] is True
        assert metadata["compatibility"]["hierarchical_covariate_model"] is True
        assert metadata["compatibility"]["covariate_integer_reduction"] is True
        assert metadata["telemetry"]["kernel_execution_space"]
        assert metadata["telemetry"]["kernel_execution_policy"] == "TeamPolicy"
        assert metadata["telemetry"]["kernel_execute_seconds"] >= 0.0
        assert metadata["telemetry"]["covariate_kernel_execution_space"] in {
            "OpenMP", "Serial", "Host"
        }
        assert metadata["telemetry"]["covariate_kernel_execution_policy"] == "TeamPolicy"
        assert metadata["telemetry"]["covariate_kernel_observations"] > 0
        assert metadata["telemetry"]["covariate_kernel_observations"] == base[
            "covariate_observations"
        ]
        assert metadata["telemetry"]["covariate_kernel_prepare_seconds"] >= 0.0
        assert metadata["telemetry"]["covariate_kernel_execute_seconds"] >= 0.0
        assert metadata["telemetry"]["covariate_kernel_workspace_bytes"] > 0
        assert metadata["telemetry"]["covariate_kernel_team_local_histogram"] is True
        assert metadata["telemetry"]["pipeline_lifecycle"] == (
            "Host decode->bounded queue->Kokkos compute->encode->sink"
        )
        assert metadata["telemetry"]["pipeline_decoded_items"] == metadata["telemetry"]["pipeline_computed_items"]
        assert metadata["telemetry"]["pipeline_computed_items"] == metadata["telemetry"]["pipeline_encoded_items"]
        assert metadata["telemetry"]["pipeline_decoded_items"] >= 1
        assert metadata["telemetry"]["pipeline_peak_decoded_bytes"] > 0
        assert metadata["telemetry"]["pipeline_peak_computed_bytes"] > 0
        assert metadata["telemetry"]["pipeline_peak_encoded_bytes"] > 0
        assert metadata["telemetry"]["output_bytes"] == report.stat().st_size
        assert metadata["telemetry"]["covariate_bytes"] == covariate_table.stat().st_size
        assert metadata["telemetry"]["wall_seconds"] >= 0.0
        assert all(item["complete"] for item in metadata["outputs"])
        assert any(item["kind"] == "bqsr-covariate-table" for item in metadata["outputs"])
        assert metadata["compatibility"]["resource_adaptive_batch"] is True
        assert metadata["telemetry"]["requested_batch_records"] == 3
        assert metadata["telemetry"]["initial_batch_records"] == 3
        assert metadata["telemetry"]["effective_batch_records"] == 3

        # The first decode must honor a small scheduler allocation before the
        # Kokkos count kernel sees any records.  This is intentionally a
        # reference-backed interval with a tiny number of reads so the oracle
        # checks scheduling telemetry without changing the report semantics.
        bounded_report = work / "bounded-recal.tsv"
        bounded_manifest = work / "bounded-recal.tsv.manifest.json"
        bounded = run_json([
            str(native / "fastgatk-bqsr"), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-69100", "--batch-records", "4096",
            "-O", str(bounded_report), "--output-manifest", str(bounded_manifest),
        ], {"SLURM_MEM_PER_NODE": "1M", "SLURM_CPUS_PER_TASK": "1"})
        bounded_metadata = json.loads(bounded_manifest.read_text(encoding="utf-8"))
        assert bounded["records"] == base["records"]
        assert bounded_metadata["telemetry"]["requested_batch_records"] == 4096
        assert 0 < bounded_metadata["telemetry"]["initial_batch_records"] < 4096
        assert 0 < bounded_metadata["telemetry"]["effective_batch_records"] < 4096
        assert bounded_metadata["telemetry"]["adaptive_batch_reductions"] >= 0

        # BaseRecalibrator checkpoints are published atomically at batch
        # boundaries. A second invocation validates the input/options
        # signature, skips the completed batches and must reproduce both the
        # GATKReport and covariate sidecar byte-for-byte.
        checkpoint_run = run_json([
            str(native / "fastgatk-bqsr"), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-69100", "--batch-records", "3", "--checkpoint", str(checkpoint),
            "-O", str(checkpoint_report), "--output-manifest", str(checkpoint_manifest),
        ])
        resumed_run = run_json([
            str(native / "fastgatk-bqsr"), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-69100", "--batch-records", "3", "--resume-checkpoint", str(checkpoint),
            "-O", str(resumed_report), "--output-manifest", str(resumed_manifest),
        ])
        assert checkpoint_run["records"] == resumed_run["records"] == 7
        assert checkpoint_report.read_bytes() == resumed_report.read_bytes()
        assert Path(f"{checkpoint_report}.covariates.tsv").read_bytes() == Path(
            f"{resumed_report}.covariates.tsv").read_bytes()
        checkpoint_metadata = json.loads(checkpoint_manifest.read_text(encoding="utf-8"))
        resumed_metadata = json.loads(resumed_manifest.read_text(encoding="utf-8"))
        assert checkpoint_metadata["telemetry"]["checkpoint_enabled"] is True
        assert checkpoint_metadata["telemetry"]["checkpoint_resumed"] is False
        assert resumed_metadata["telemetry"]["checkpoint_resumed"] is True
        assert resumed_metadata["telemetry"]["checkpoint_records"] == 7
        assert resumed_metadata["compatibility"]["checkpoint_resume"] is True
        mismatched = subprocess.run([
            str(native / "fastgatk-bqsr"), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-69100", "--batch-records", "4", "--resume-checkpoint", str(checkpoint),
            "-O", str(work / "mismatched.tsv"),
        ], text=True, capture_output=True, check=False)
        assert mismatched.returncode != 0
        assert "batch size mismatch" in mismatched.stderr

        known_sites.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=17,length=1000000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "17\t69067\t.\tT\tG\t.\tPASS\t.\n", encoding="utf-8")
        known_sites_second.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=17,length=1000000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            # A two-base REF span exercises VariantContext/Locatable masking
            # beyond the VCF record's first coordinate.
            "17\t69068\t.\tCA\tC\t.\tPASS\t.\n", encoding="utf-8")
        known = run_json([
            str(native / "fastgatk-bqsr"), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-69100", "--known-sites", str(known_sites),
            "-O", str(known_report), "--output-manifest", str(known_manifest),
        ])
        known_rows = table_rows(known_report.read_text(encoding="utf-8"), "RecalTable1")
        known_observations = sum(int(row[4]) for row in known_rows if row[2] == "M")
        assert known["tool"] == "BaseRecalibrator"
        assert 0 < known_observations < baseline_observations
        known_metadata = json.loads(known_manifest.read_text(encoding="utf-8"))
        assert known_metadata["compatibility"]["known_sites_filter"] is True
        multi_known = run_json([
            str(native / "fastgatk-bqsr"), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-69100", "--known-sites", str(known_sites),
            "--known-sites", str(known_sites_second), "-O", str(multi_known_report),
            "--output-manifest", str(multi_known_manifest),
        ])
        multi_known_metadata = json.loads(multi_known_manifest.read_text(encoding="utf-8"))
        assert multi_known["tool"] == "BaseRecalibrator"
        assert multi_known_metadata["compatibility"]["known_sites_filter"] is True
        assert multi_known_metadata["compatibility"]["known_sites_count"] == 2
        multi_known_rows = table_rows(multi_known_report.read_text(encoding="utf-8"), "RecalTable1")
        assert sum(int(row[4]) for row in multi_known_rows if row[2] == "M") < known_observations

        # Keep a full-file BaseRecalibrator run in the compatibility matrix;
        # the interval smoke above catches coordinate selection while this
        # path exercises EOF batching and report aggregation across all reads.
        full = run_json([
            str(native / "fastgatk-bqsr"), "-I", str(bam), "-R", str(reference),
            "--known-sites", str(known_sites), "-O", str(full_report),
        ])
        # BaseRecalibrator's standard filters remove 120 records from this
        # fixture (mapping-quality, duplicate/QC and malformed-read cases).
        assert full["records"] == 373 and full["bases"] == 28348

        # GATK's optional indel model emits I/D rows at the default indel
        # quality for every eligible read offset.  CIGAR events contribute the
        # error count; the context/cycle covariates use the independent indel
        # context size and the four-base cycle cushion.
        indel = run_json([
            str(native / "fastgatk-bqsr"), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-69100", "--known-sites", str(known_sites),
            "--compute-indel-bqsr-tables", "--indels-context-size", "4",
            "-O", str(indel_report),
        ])
        indel_rows = table_rows(indel_report.read_text(encoding="utf-8"), "RecalTable1")
        event_counts = {event: sum(int(row[4]) for row in indel_rows if row[2] == event)
                        for event in ("M", "I", "D")}
        assert event_counts["M"] > 0
        assert event_counts["I"] == event_counts["D"] == event_counts["M"]
        indel_metadata = json.loads(Path(f"{indel_report}.manifest.json").read_text(encoding="utf-8"))
        assert indel_metadata["telemetry"]["compute_indel_bqsr_tables"] is True
        assert indel_metadata["telemetry"]["indels_context_size"] == 4

        gathered = run_json([
            str(native / "fastgatk-gather-bqsr-reports"), "-I", str(report),
            "--input", str(report), "-O", str(merged_report),
            "--output-manifest", str(merged_manifest),
        ])
        assert gathered["tool"] == "GatherBQSRReports" and gathered["input_reports"] == 2
        merged_rows = table_rows(merged_report.read_text(encoding="utf-8"), "RecalTable1")
        original_counts = [int(row[4]) for row in rows if row[2] == "M"]
        merged_counts = [int(row[4]) for row in merged_rows if row[2] == "M"]
        assert merged_counts == [count * 2 for count in original_counts]
        merged_metadata = json.loads(merged_manifest.read_text(encoding="utf-8"))
        assert merged_metadata["execution_space"] == "Host"
        assert merged_metadata["determinism"] == "strict"
        assert merged_metadata["primary_output"] == str(merged_report)
        assert merged_metadata["primary_output_kind"] == "recalibration-report"
        assert merged_metadata["compatibility"]["deterministic_merge"] is True
        assert merged_metadata["compatibility"]["hierarchical_covariate_model"] is True
        assert merged_metadata["telemetry"]["input_reports"] == 2
        assert merged_metadata["telemetry"]["output_bytes"] == merged_report.stat().st_size
        assert merged_metadata["telemetry"]["covariate_bytes"] == Path(
            f"{merged_report}.covariates.tsv").stat().st_size
        assert merged_metadata["telemetry"]["wall_seconds"] >= 0.0
        assert all(item["complete"] for item in merged_metadata["outputs"])
        merged_covariates = [line for line in (merged_report.with_name(merged_report.name + ".covariates.tsv"))
                             .read_text(encoding="utf-8").splitlines()
                             if line and not line.startswith("#")]
        assert len(merged_covariates) == len(covariate_rows)

        apply = run_json([
            str(native / "fastgatk-apply-bqsr"), "-I", str(bam), "-R", str(reference),
            "--bqsr-recal-file", str(report), "-O", str(recal_bam),
            "--output-manifest", str(apply_manifest),
        ])
        assert apply["tool"] == "ApplyBQSR" and apply["records"] == 493
        assert recal_bam.exists() and recal_bam.stat().st_size > 0
        assert Path(f"{recal_bam}.bai").exists()
        applied_metadata = json.loads(apply_manifest.read_text(encoding="utf-8"))
        assert applied_metadata["execution_space"] in {"OpenMP", "Serial", "Host"}
        assert applied_metadata["determinism"] == "strict"
        assert applied_metadata["primary_output_kind"] == "bam-or-cram"
        assert applied_metadata["compatibility"]["report_recal_table_fallback"] is True
        assert applied_metadata["compatibility"]["hts_index"] is True
        assert applied_metadata["telemetry"]["kernel_execution_space"]
        assert applied_metadata["telemetry"]["kernel_execution_policy"] == "RangePolicy"
        assert applied_metadata["telemetry"]["kernel_execute_seconds"] >= 0.0
        assert applied_metadata["telemetry"]["pipeline_lifecycle"] == (
            "Host decode->bounded queue->Kokkos compute->encode->sink"
        )
        pipeline = applied_metadata["telemetry"]
        assert pipeline["pipeline_decoded_items"] == pipeline["pipeline_computed_items"]
        assert pipeline["pipeline_computed_items"] == pipeline["pipeline_encoded_items"]
        assert pipeline["pipeline_decoded_items"] >= 1
        assert pipeline["pipeline_peak_decoded_bytes"] > 0
        assert pipeline["pipeline_peak_computed_bytes"] > 0
        assert pipeline["pipeline_peak_encoded_bytes"] > 0
        assert applied_metadata["telemetry"]["output_bytes"] == recal_bam.stat().st_size
        assert applied_metadata["telemetry"]["index_bytes"] == Path(f"{recal_bam}.bai").stat().st_size
        assert applied_metadata["telemetry"]["wall_seconds"] >= 0.0
        assert all(item["complete"] for item in applied_metadata["outputs"])
        assert any(item["kind"] == "bam-or-cram-index" for item in applied_metadata["outputs"])
        assert applied_metadata["compatibility"]["resource_adaptive_batch"] is True
        assert applied_metadata["telemetry"]["requested_batch_records"] == 4096
        assert applied_metadata["telemetry"]["initial_batch_records"] == 4096

        bounded_apply_bam = work / "bounded-apply.bam"
        bounded_apply_manifest = work / "bounded-apply.bam.manifest.json"
        bounded_apply = run_json([
            str(native / "fastgatk-apply-bqsr"), "-I", str(bam), "-R", str(reference),
            "--bqsr-recal-file", str(report), "--batch-records", "4096",
            "-O", str(bounded_apply_bam), "--output-manifest", str(bounded_apply_manifest),
        ], {"SLURM_MEM_PER_NODE": "1M", "SLURM_CPUS_PER_TASK": "1"})
        bounded_apply_metadata = json.loads(bounded_apply_manifest.read_text(encoding="utf-8"))
        assert bounded_apply["records"] == apply["records"]
        assert bounded_apply_metadata["telemetry"]["requested_batch_records"] == 4096
        assert 0 < bounded_apply_metadata["telemetry"]["initial_batch_records"] < 4096
        assert 0 < bounded_apply_metadata["telemetry"]["effective_batch_records"] < 4096

        # ApplyBQSR resumes by validating and reusing the durable output prefix
        # recorded by the checkpoint, then atomically replacing that path with
        # the completed stream and rebuilding its index.
        checkpoint_apply = run_json([
            str(native / "fastgatk-apply-bqsr"), "-I", str(bam), "-R", str(reference),
            "--bqsr-recal-file", str(report), "--batch-records", "17",
            "--checkpoint", str(apply_checkpoint), "-O", str(apply_checkpoint_bam),
            "--output-manifest", str(apply_checkpoint_manifest),
        ])
        checkpoint_bytes = apply_checkpoint_bam.read_bytes()
        checkpoint_index_bytes = Path(f"{apply_checkpoint_bam}.bai").read_bytes()
        resumed_apply = run_json([
            str(native / "fastgatk-apply-bqsr"), "-I", str(bam), "-R", str(reference),
            "--bqsr-recal-file", str(report), "--batch-records", "17",
            "--resume-checkpoint", str(apply_checkpoint), "-O", str(apply_checkpoint_bam),
            "--output-manifest", str(apply_resumed_manifest),
        ])
        assert checkpoint_apply["records"] == resumed_apply["records"] == 493
        assert apply_checkpoint_bam.read_bytes() == checkpoint_bytes
        assert Path(f"{apply_checkpoint_bam}.bai").read_bytes() == checkpoint_index_bytes
        checkpoint_apply_metadata = json.loads(apply_checkpoint_manifest.read_text(encoding="utf-8"))
        resumed_apply_metadata = json.loads(apply_resumed_manifest.read_text(encoding="utf-8"))
        assert checkpoint_apply_metadata["telemetry"]["checkpoint_enabled"] is True
        assert checkpoint_apply_metadata["telemetry"]["checkpoint_resumed"] is False
        assert resumed_apply_metadata["telemetry"]["checkpoint_resumed"] is True
        assert resumed_apply_metadata["telemetry"]["checkpoint_records"] == 493
        assert resumed_apply_metadata["telemetry"]["checkpoint_prefix_reused"] is True
        before = run_json([str(smoke), "-I", str(bam)])
        after = run_json([str(smoke), "-I", str(recal_bam)])
        assert before["reads"] == after["reads"] == 493
        assert before["bases"] == after["bases"] == 37427
        assert before["mapq_sum"] == after["mapq_sum"]
        assert before["quality_sum"] != after["quality_sum"]

        # A report produced by the Java GATK does not carry our optional
        # sidecar.  RecalTable1/2 must still provide deterministic fallback
        # adjustments instead of silently reverting to an all-zero model.
        shutil.copyfile(report, report_without_sidecar)
        apply_without_sidecar = run_json([
            str(native / "fastgatk-apply-bqsr"), "-I", str(bam), "-R", str(reference),
            "--bqsr-recal-file", str(report_without_sidecar),
            "-O", str(recal_bam_without_sidecar),
        ])
        assert apply_without_sidecar["tool"] == "ApplyBQSR"
        fallback = run_json([str(smoke), "-I", str(recal_bam_without_sidecar)])
        assert fallback["reads"] == 493 and fallback["bases"] == 37427
        assert fallback["quality_sum"] != before["quality_sum"]

        prior_bam = work / "recal-prior30.bam"
        prior_manifest = work / "recal-prior30.bam.manifest.json"
        prior_apply = run_json([
            str(native / "fastgatk-apply-bqsr"), "-I", str(bam), "-R", str(reference),
            "--bqsr-recal-file", str(report), "--global-qscore-prior", "30",
            "-O", str(prior_bam), "--output-manifest", str(prior_manifest),
        ])
        prior_after = run_json([str(smoke), "-I", str(prior_bam)])
        assert prior_apply["tool"] == "ApplyBQSR"
        assert json.loads(prior_manifest.read_text(encoding="utf-8"))["telemetry"]["global_qscore_prior"] == 30.0
        assert prior_after["quality_sum"] != after["quality_sum"]

        # When the pinned Java GATK fixture is available, make the report and
        # ApplyBQSR compatibility claim executable rather than relying only on
        # the native round-trip.  The native report is compared table-by-table
        # and the Java report is then consumed directly by the native reader.
        java = root / "third_party/jdk17/bin/java"
        gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        oracle = {"status": "skip", "reason": "bundled GATK/JDK oracle not present"}
        if java.exists() and gatk_jar.exists():
            gatk_report = work / "oracle-gatk-recal.tsv"
            gatk_full_report = work / "oracle-gatk-full-recal.tsv"
            gatk_indel_report = work / "oracle-gatk-indel-recal.tsv"
            gatk_gather_report = work / "oracle-gatk-gather-recal.tsv"
            native_gather_report = work / "oracle-native-gather-recal.tsv"
            gatk_bam = work / "oracle-gatk.bam"
            native_oracle_bam = work / "oracle-native.bam"
            gatk_sam = work / "oracle-gatk.sam"
            native_sam = work / "oracle-native.sam"
            gatk_oq_bam = work / "oracle-gatk-oq.bam"
            native_oq_bam = work / "oracle-native-oq.bam"
            gatk_oq_sam = work / "oracle-gatk-oq.sam"
            native_oq_sam = work / "oracle-native-oq.sam"
            source_sam = work / "oracle-source.sam"
            no_oq_sam = work / "oracle-source-no-oq.sam"
            gatk_emit_bam = work / "oracle-gatk-emit-oq.bam"
            native_emit_bam = work / "oracle-native-emit-oq.bam"
            gatk_emit_sam = work / "oracle-gatk-emit-oq.sam"
            native_emit_sam = work / "oracle-native-emit-oq.sam"
            gatk_static_bam = work / "oracle-gatk-static.bam"
            native_static_bam = work / "oracle-native-static.bam"
            gatk_static_sam = work / "oracle-gatk-static.sam"
            native_static_sam = work / "oracle-native-static.sam"
            gatk_static_down_bam = work / "oracle-gatk-static-down.bam"
            native_static_down_bam = work / "oracle-native-static-down.bam"
            gatk_static_down_sam = work / "oracle-gatk-static-down.sam"
            native_static_down_sam = work / "oracle-native-static-down.sam"
            gatk_dynamic_bam = work / "oracle-gatk-dynamic.bam"
            native_dynamic_bam = work / "oracle-native-dynamic.bam"
            gatk_dynamic_sam = work / "oracle-gatk-dynamic.sam"
            native_dynamic_sam = work / "oracle-native-dynamic.sam"
            missing_rg_report = work / "oracle-missing-rg-recal.tsv"
            gatk_missing_rg_bam = work / "oracle-gatk-missing-rg.bam"
            native_missing_rg_bam = work / "oracle-native-missing-rg.bam"
            gatk_missing_rg_sam = work / "oracle-gatk-missing-rg.sam"
            native_missing_rg_sam = work / "oracle-native-missing-rg.sam"
            run_command([
                str(java), "-jar", str(gatk_jar), "IndexFeatureFile", "-I", str(known_sites),
            ])
            run_command([
                str(java), "-jar", str(gatk_jar), "IndexFeatureFile", "-I", str(known_sites_second),
            ])
            run_command([
                str(java), "-jar", str(gatk_jar), "BaseRecalibrator", "-R", str(reference),
                "-I", str(bam), "-L", "17:69000-69100", "--known-sites", str(known_sites),
                "-O", str(gatk_report),
            ])
            run_command([
                str(java), "-jar", str(gatk_jar), "BaseRecalibrator", "-R", str(reference),
                "-I", str(bam), "--known-sites", str(known_sites),
                "-O", str(gatk_full_report),
            ])
            gatk_multi_report = work / "oracle-gatk-multi-known-recal.tsv"
            run_command([
                str(java), "-jar", str(gatk_jar), "BaseRecalibrator", "-R", str(reference),
                "-I", str(bam), "-L", "17:69000-69100", "--known-sites", str(known_sites),
                "--known-sites", str(known_sites_second), "-O", str(gatk_multi_report),
            ])
            run_command([
                str(java), "-jar", str(gatk_jar), "BaseRecalibrator", "-R", str(reference),
                "-I", str(bam), "-L", "17:69000-69100", "--known-sites", str(known_sites),
                "--compute-indel-bqsr-tables", "--indels-context-size", "4",
                "-O", str(gatk_indel_report),
            ])
            # GatherBQSRReports must consume a report emitted by Java GATK
            # directly; no FASTGATK covariate sidecar is present.  Compare the
            # merged hierarchical tables against the pinned Java gatherer.
            run_command([
                str(java), "-jar", str(gatk_jar), "GatherBQSRReports",
                "-I", str(gatk_report), "-I", str(gatk_report),
                "-O", str(gatk_gather_report),
            ])
            native_gather = run_json([
                str(native / "fastgatk-gather-bqsr-reports"),
                "-I", str(gatk_report), "--input", str(gatk_report),
                "-O", str(native_gather_report),
            ])
            native_gather_rows = {
                table: Counter(tuple(row) for row in table_rows(native_gather_report.read_text(encoding="utf-8"), table))
                for table in ("Arguments", "Quantized", "RecalTable0", "RecalTable1", "RecalTable2")
            }
            gatk_gather_rows = {
                table: Counter(tuple(row) for row in table_rows(gatk_gather_report.read_text(encoding="utf-8"), table))
                for table in ("Arguments", "Quantized", "RecalTable0", "RecalTable1", "RecalTable2")
            }
            assert native_gather_rows == gatk_gather_rows
            native_rows = {
                table: Counter(tuple(row) for row in table_rows(known_report.read_text(encoding="utf-8"), table))
                for table in ("Quantized", "RecalTable0", "RecalTable1", "RecalTable2")
            }
            gatk_rows = {
                table: Counter(tuple(row) for row in table_rows(gatk_report.read_text(encoding="utf-8"), table))
                for table in ("Quantized", "RecalTable0", "RecalTable1", "RecalTable2")
            }
            assert native_rows == gatk_rows
            native_full_rows = {
                table: Counter(tuple(row) for row in table_rows(full_report.read_text(encoding="utf-8"), table))
                for table in ("Quantized", "RecalTable0", "RecalTable1", "RecalTable2")
            }
            gatk_full_rows = {
                table: Counter(tuple(row) for row in table_rows(gatk_full_report.read_text(encoding="utf-8"), table))
                for table in ("Quantized", "RecalTable0", "RecalTable1", "RecalTable2")
            }
            assert native_full_rows == gatk_full_rows
            native_indel_rows = {
                table: Counter(tuple(row) for row in table_rows(indel_report.read_text(encoding="utf-8"), table))
                for table in ("Quantized", "RecalTable0", "RecalTable1", "RecalTable2")
            }
            gatk_indel_rows = {
                table: Counter(tuple(row) for row in table_rows(gatk_indel_report.read_text(encoding="utf-8"), table))
                for table in ("Quantized", "RecalTable0", "RecalTable1", "RecalTable2")
            }
            assert native_indel_rows == gatk_indel_rows
            assert {
                table: Counter(tuple(row) for row in table_rows(multi_known_report.read_text(encoding="utf-8"), table))
                for table in ("Quantized", "RecalTable0", "RecalTable1", "RecalTable2")
            } == {
                table: Counter(tuple(row) for row in table_rows(gatk_multi_report.read_text(encoding="utf-8"), table))
                for table in ("Quantized", "RecalTable0", "RecalTable1", "RecalTable2")
            }
            run_command([
                str(java), "-jar", str(gatk_jar), "ApplyBQSR", "-R", str(reference),
                "-I", str(bam), "--bqsr-recal-file", str(gatk_report), "-O", str(gatk_bam),
            ])
            native_java_apply = run_json([
                str(native / "fastgatk-apply-bqsr"), "-I", str(bam), "-R", str(reference),
                "--bqsr-recal-file", str(gatk_report), "-O", str(native_oracle_bam),
            ])
            # Convert both BAMs through the same pinned Java reader so the
            # oracle compares alignment/sequence/QUAL fields record-by-record,
            # not merely the aggregate quality sum.
            run_command([
                str(java), "-jar", str(gatk_jar), "PrintReads", "-I", str(gatk_bam),
                "-O", str(gatk_sam), "--create-output-bam-index", "false",
            ])
            run_command([
                str(java), "-jar", str(gatk_jar), "PrintReads", "-I", str(native_oracle_bam),
                "-O", str(native_sam), "--create-output-bam-index", "false",
            ])
            gatk_records = sam_record_map(gatk_sam)
            native_records = sam_record_map(native_sam)
            assert gatk_records == native_records
            gatk_after = run_json([str(smoke), "-I", str(gatk_bam)])
            native_after = run_json([str(smoke), "-I", str(native_oracle_bam)])
            assert native_after["reads"] == gatk_after["reads"] == 493
            assert native_after["bases"] == gatk_after["bases"] == 37427
            assert native_after["mapq_sum"] == gatk_after["mapq_sum"]
            assert native_after["quality_sum"] == gatk_after["quality_sum"]

            # --create-output-bam-index is an optional boolean.  A false
            # value must suppress the sidecar in both Java and native paths;
            # silently creating an index changes Nextflow/Picard artifact
            # contracts even when the BAM records themselves are identical.
            run_command([
                str(java), "-jar", str(gatk_jar), "ApplyBQSR", "-R", str(reference),
                "-I", str(bam), "--bqsr-recal-file", str(gatk_report),
                "--create-output-bam-index", "false", "-O", str(gatk_no_index_bam),
            ])
            native_no_index_apply = run_json([
                str(native / "fastgatk-apply-bqsr"), "-I", str(bam), "-R", str(reference),
                "--bqsr-recal-file", str(gatk_report), "--create-output-bam-index=false",
                "-O", str(native_no_index_bam),
            ])
            assert not Path(f"{gatk_no_index_bam}.bai").exists()
            assert not Path(f"{native_no_index_bam}.bai").exists()
            no_index_metadata = json.loads(
                Path(f"{native_no_index_bam}.manifest.json").read_text(encoding="utf-8"))
            assert no_index_metadata["compatibility"]["create_output_bam_index"] is False
            assert no_index_metadata["telemetry"]["create_output_bam_index"] is False
            assert no_index_metadata["compatibility"]["hts_index"] is False
            assert all(item["kind"] != "bam-or-cram-index" for item in no_index_metadata["outputs"])
            run_command([
                str(java), "-jar", str(gatk_jar), "PrintReads", "-I", str(gatk_no_index_bam),
                "-O", str(gatk_no_index_sam), "--create-output-bam-index", "false",
            ])
            run_command([
                str(java), "-jar", str(gatk_jar), "PrintReads", "-I", str(native_no_index_bam),
                "-O", str(native_no_index_sam), "--create-output-bam-index", "false",
            ])
            assert sam_record_lines(gatk_no_index_sam) == sam_record_lines(native_no_index_sam)
            # --use-original-qualities must reset the model lookup to the
            # FASTQ-encoded OQ tag while preserving that tag in the output.
            # Compare complete record lines (not a QNAME map, since paired
            # reads legitimately share QNAMEs) after the same Java reader.
            run_command([
                str(java), "-jar", str(gatk_jar), "ApplyBQSR", "-R", str(reference),
                "-I", str(bam), "--bqsr-recal-file", str(gatk_report),
                "--use-original-qualities", "-O", str(gatk_oq_bam),
            ])
            native_oq_apply = run_json([
                str(native / "fastgatk-apply-bqsr"), "-I", str(bam), "-R", str(reference),
                "--bqsr-recal-file", str(gatk_report), "--use-original-qualities",
                "-O", str(native_oq_bam),
            ])
            run_command([
                str(java), "-jar", str(gatk_jar), "PrintReads", "-I", str(gatk_oq_bam),
                "-O", str(gatk_oq_sam), "--create-output-bam-index", "false",
            ])
            run_command([
                str(java), "-jar", str(gatk_jar), "PrintReads", "-I", str(native_oq_bam),
                "-O", str(native_oq_sam), "--create-output-bam-index", "false",
            ])
            gatk_oq_records = sam_record_lines(gatk_oq_sam)
            native_oq_records = sam_record_lines(native_oq_sam)
            assert len(gatk_oq_records) == len(native_oq_records) == 493
            assert gatk_oq_records == native_oq_records
            assert any("\tOQ:Z:" in record for record in native_oq_records)
            # Exercise the complementary --emit-original-quals path with a
            # real SAM input in which every OQ tag has been removed.  This
            # keeps the fixture self-contained (no pysam/samtools dependency)
            # while comparing HTSJDK and HTSlib output through PrintReads.
            run_command([
                str(java), "-jar", str(gatk_jar), "PrintReads", "-I", str(bam),
                "-O", str(source_sam), "--create-output-bam-index", "false",
            ])
            source_without_oq = []
            for line in source_sam.read_text(encoding="utf-8").splitlines():
                if line.startswith("@"):
                    source_without_oq.append(line)
                    continue
                fields = [field for field in line.split("\t") if not field.startswith("OQ:")]
                source_without_oq.append("\t".join(fields))
            no_oq_sam.write_text("\n".join(source_without_oq) + "\n", encoding="utf-8")
            assert all("\tOQ:" not in record for record in sam_record_lines(no_oq_sam))
            run_command([
                str(java), "-jar", str(gatk_jar), "ApplyBQSR", "-R", str(reference),
                "-I", str(no_oq_sam), "--bqsr-recal-file", str(gatk_report),
                "--emit-original-quals", "true", "-O", str(gatk_emit_bam),
            ])
            native_emit_apply = run_json([
                str(native / "fastgatk-apply-bqsr"), "-I", str(no_oq_sam), "-R", str(reference),
                "--bqsr-recal-file", str(gatk_report), "--emit-original-quals=true",
                "-O", str(native_emit_bam),
            ])
            run_command([
                str(java), "-jar", str(gatk_jar), "PrintReads", "-I", str(gatk_emit_bam),
                "-O", str(gatk_emit_sam), "--create-output-bam-index", "false",
            ])
            run_command([
                str(java), "-jar", str(gatk_jar), "PrintReads", "-I", str(native_emit_bam),
                "-O", str(native_emit_sam), "--create-output-bam-index", "false",
            ])
            gatk_emit_records = sam_record_lines(gatk_emit_sam)
            native_emit_records = sam_record_lines(native_emit_sam)
            assert len(gatk_emit_records) == len(native_emit_records) == 493
            assert gatk_emit_records == native_emit_records
            assert all("\tOQ:Z:" in record for record in native_emit_records)

            # Static quantization is separate from dynamic --quantize-quals.
            # GATK maps through probability-space nearest bins (or the lower
            # bin when --round-down-quantized is set); compare every record
            # against the pinned Java transformer through a common reader.
            run_command([
                str(java), "-jar", str(gatk_jar), "ApplyBQSR", "-R", str(reference),
                "-I", str(bam), "--bqsr-recal-file", str(gatk_report),
                "--static-quantized-quals", "10", "--static-quantized-quals", "20",
                "--static-quantized-quals", "30", "-O", str(gatk_static_bam),
            ])
            native_static_apply = run_json([
                str(native / "fastgatk-apply-bqsr"), "-I", str(bam), "-R", str(reference),
                "--bqsr-recal-file", str(gatk_report),
                "--static-quantized-quals", "10", "--static-quantized-quals", "20",
                "--static-quantized-quals", "30", "-O", str(native_static_bam),
            ])
            run_command([
                str(java), "-jar", str(gatk_jar), "PrintReads", "-I", str(gatk_static_bam),
                "-O", str(gatk_static_sam), "--create-output-bam-index", "false",
            ])
            run_command([
                str(java), "-jar", str(gatk_jar), "PrintReads", "-I", str(native_static_bam),
                "-O", str(native_static_sam), "--create-output-bam-index", "false",
            ])
            assert sam_record_lines(gatk_static_sam) == sam_record_lines(native_static_sam)

            run_command([
                str(java), "-jar", str(gatk_jar), "ApplyBQSR", "-R", str(reference),
                "-I", str(bam), "--bqsr-recal-file", str(gatk_report),
                "--static-quantized-quals", "10", "--static-quantized-quals", "20",
                "--static-quantized-quals", "30", "--round-down-quantized",
                "-O", str(gatk_static_down_bam),
            ])
            native_static_down_apply = run_json([
                str(native / "fastgatk-apply-bqsr"), "-I", str(bam), "-R", str(reference),
                "--bqsr-recal-file", str(gatk_report),
                "--static-quantized-quals", "10", "--static-quantized-quals", "20",
                "--static-quantized-quals", "30", "--round-down-quantized",
                "-O", str(native_static_down_bam),
            ])
            run_command([
                str(java), "-jar", str(gatk_jar), "PrintReads", "-I", str(gatk_static_down_bam),
                "-O", str(gatk_static_down_sam), "--create-output-bam-index", "false",
            ])
            run_command([
                str(java), "-jar", str(gatk_jar), "PrintReads", "-I", str(native_static_down_bam),
                "-O", str(native_static_down_sam), "--create-output-bam-index", "false",
            ])
            assert sam_record_lines(gatk_static_down_sam) == sam_record_lines(native_static_down_sam)

            # A positive --quantize-quals requests a fresh QualQuantizer run
            # over the report's empirical histogram (rather than blindly
            # reusing the report's precomputed level count).
            run_command([
                str(java), "-jar", str(gatk_jar), "ApplyBQSR", "-R", str(reference),
                "-I", str(bam), "--bqsr-recal-file", str(gatk_report),
                "--quantize-quals", "4", "-O", str(gatk_dynamic_bam),
            ])
            native_dynamic_apply = run_json([
                str(native / "fastgatk-apply-bqsr"), "-I", str(bam), "-R", str(reference),
                "--bqsr-recal-file", str(gatk_report), "--quantize-quals", "4",
                "-O", str(native_dynamic_bam),
            ])
            run_command([
                str(java), "-jar", str(gatk_jar), "PrintReads", "-I", str(gatk_dynamic_bam),
                "-O", str(gatk_dynamic_sam), "--create-output-bam-index", "false",
            ])
            run_command([
                str(java), "-jar", str(gatk_jar), "PrintReads", "-I", str(native_dynamic_bam),
                "-O", str(native_dynamic_sam), "--create-output-bam-index", "false",
            ])
            assert sam_record_lines(gatk_dynamic_sam) == sam_record_lines(native_dynamic_sam)

            # Dynamic and static quantization are mutually exclusive in GATK;
            # reject the combination rather than choosing one silently.
            static_mutex = subprocess.run([
                str(native / "fastgatk-apply-bqsr"), "-I", str(bam), "-R", str(reference),
                "--bqsr-recal-file", str(gatk_report), "--quantize-quals", "4",
                "--static-quantized-quals", "10", "-O", str(work / "mutex.bam"),
            ], text=True, capture_output=True, check=False)
            assert static_mutex.returncode != 0
            assert "cannot be combined" in static_mutex.stderr

            # A read group absent from the recalibration table is a hard error
            # by default.  With --allow-missing-read-group, GATK bypasses the
            # covariate model and maps the original qualities through the
            # selected quantizer; exercise both branches and compare output
            # records byte-for-byte after Java decoding.
            # Keep the read/header RG valid and remove that group from the
            # recalibration model instead (the real GATK integration test
            # filters a group out during BaseRecalibrator).
            missing_rg_report.write_text(
                gatk_report.read_text(encoding="utf-8").replace(
                    "809R9ABXX101220.5", "MISSING_RECAL_GROUP"),
                encoding="utf-8")
            missing_rg_rejected = subprocess.run([
                str(native / "fastgatk-apply-bqsr"), "-I", str(bam),
                "-R", str(reference), "--bqsr-recal-file", str(missing_rg_report),
                "-O", str(work / "missing-rg-rejected.bam"),
            ], text=True, capture_output=True, check=False)
            assert missing_rg_rejected.returncode != 0
            assert "read group not found" in missing_rg_rejected.stderr
            run_command([
                str(java), "-jar", str(gatk_jar), "ApplyBQSR", "-R", str(reference),
                "-I", str(bam), "--bqsr-recal-file", str(missing_rg_report),
                "--allow-missing-read-group", "-O", str(gatk_missing_rg_bam),
            ])
            native_missing_rg_apply = run_json([
                str(native / "fastgatk-apply-bqsr"), "-I", str(bam),
                "-R", str(reference), "--bqsr-recal-file", str(missing_rg_report),
                "--allow-missing-read-group", "-O", str(native_missing_rg_bam),
            ])
            run_command([
                str(java), "-jar", str(gatk_jar), "PrintReads", "-I", str(gatk_missing_rg_bam),
                "-O", str(gatk_missing_rg_sam), "--create-output-bam-index", "false",
            ])
            run_command([
                str(java), "-jar", str(gatk_jar), "PrintReads", "-I", str(native_missing_rg_bam),
                "-O", str(native_missing_rg_sam), "--create-output-bam-index", "false",
            ])
            assert sam_record_lines(gatk_missing_rg_sam) == sam_record_lines(native_missing_rg_sam)
            oracle = {
                "status": "pass",
                "report_tables_bit_identical": True,
                "full_report_tables_bit_identical": True,
                "indel_report_tables_bit_identical": True,
                "gather_report_tables_bit_identical": True,
                "gather_report_all_tables_bit_identical": True,
                "apply_record_core_bit_identical": True,
                "apply_records_compared": len(gatk_records),
                "apply_quality_sum": native_after["quality_sum"],
                "apply_quality_sum_gatk": gatk_after["quality_sum"],
                "native_apply_records": native_java_apply["records"],
                "apply_without_output_index_record_lines_bit_identical": True,
                "apply_without_output_index_records_compared": len(sam_record_lines(native_no_index_sam)),
                "native_apply_without_output_index_records": native_no_index_apply["records"],
                "apply_original_qualities_record_lines_bit_identical": True,
                "apply_original_qualities_records_compared": len(native_oq_records),
                "native_apply_original_qualities_records": native_oq_apply["records"],
                "apply_emit_original_qualities_record_lines_bit_identical": True,
                "apply_emit_original_qualities_records_compared": len(native_emit_records),
                "native_emit_original_qualities_records": native_emit_apply["records"],
                "apply_static_quantized_record_lines_bit_identical": True,
                "apply_static_quantized_records_compared": len(sam_record_lines(native_static_sam)),
                "apply_static_round_down_record_lines_bit_identical": True,
                "apply_static_round_down_records_compared": len(sam_record_lines(native_static_down_sam)),
                "native_apply_static_records": native_static_apply["records"],
                "native_apply_static_round_down_records": native_static_down_apply["records"],
                "apply_dynamic_quantized_record_lines_bit_identical": True,
                "apply_dynamic_quantized_records_compared": len(sam_record_lines(native_dynamic_sam)),
                "native_apply_dynamic_quantized_records": native_dynamic_apply["records"],
                "apply_missing_read_group_record_lines_bit_identical": True,
                "apply_missing_read_group_records_compared": len(sam_record_lines(native_missing_rg_sam)),
                "native_apply_missing_read_group_records": native_missing_rg_apply["records"],
            }
        elif os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise AssertionError("FASTGATK_REQUIRE_GATK_ORACLE=1 but bundled GATK/JDK oracle is missing")
        print(json.dumps({"status": "pass", "base": base, "apply": apply,
                          "quality_sum_before": before["quality_sum"],
                          "quality_sum_after": after["quality_sum"],
                          "apply_without_sidecar": apply_without_sidecar,
                          "quality_sum_without_sidecar": fallback["quality_sum"],
                          "gatk_oracle": oracle}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
