#!/usr/bin/env python3
"""Verify unified core + real HTSlib BAM/CRAM HaplotypeCaller smoke adapter."""
from __future__ import annotations

import json
import gzip
import os
import subprocess
import tempfile
from pathlib import Path


def run(binary: Path, input_path: Path, reference: Path | None = None,
        region: str | None = None) -> dict:
    command = ["taskset", "-c", "0-3", str(binary), f"--input={input_path}",
               "--batch-records=32", "--threads=4"]
    if reference:
        command.append(f"--reference={reference}")
    if region:
        command.append(f"--region={region}")
    env = os.environ.copy()
    env["OMP_PROC_BIND"] = "true"
    env["OMP_PLACES"] = "threads"
    output = subprocess.check_output(command, text=True, env=env)
    result = json.loads(output.strip().splitlines()[-1])
    assert result["status"] in {"smoke", "prototype"}
    assert result["htslib_backend"].startswith("HTSlib")
    assert result["execution_space"] in {"OpenMP", "Threads", "Serial"}
    assert result["execute_seconds"] > 0.0
    return result


def run_call(binary: Path, input_path: Path, reference: Path | None = None,
             region: str | None = None, batch_records: int = 32,
             extra: list[str] | None = None) -> dict:
    command = ["taskset", "-c", "0-3", str(binary), f"--input={input_path}",
               f"--batch-records={batch_records}", "--threads=4", "--max-candidates=16"]
    if reference:
        command.append(f"--reference={reference}")
    if region:
        command.append(f"--region={region}")
    if extra:
        command.extend(extra)
    env = os.environ.copy()
    env["OMP_PROC_BIND"] = "true"
    env["OMP_PLACES"] = "threads"
    output = subprocess.check_output(command, text=True, env=env)
    result = json.loads(output.strip().splitlines()[-1])
    assert result["status"] in {"smoke", "prototype"}
    assert result["pipeline"] == "assembly-likelihood-genotyping"
    assert result["htslib_backend"].startswith("HTSlib")
    assert result["candidate_sites"] <= 16
    assert result["variant_calls"] <= result["candidate_sites"]
    assert result["assembly_seconds"] > 0.0
    return result


def run_call_artifacts(binary: Path, input_path: Path) -> dict:
    """Exercise the GATK aliases and compressed VCF/manifest file contract."""
    workspace = Path(__file__).resolve().parents[2]
    with tempfile.TemporaryDirectory(prefix="fastgatk-call-artifacts-") as temp:
        root = Path(temp)
        vcf = root / "calls.vcf.gz"
        manifest = root / "calls.vcf.gz.manifest.json"
        command = ["taskset", "-c", "0-3", str(binary), "-I", str(input_path),
                   "-L", "17:69000-69100", "-O", str(vcf),
                   "--output-manifest", str(manifest), "--native-pair-hmm-threads", "4"]
        env = os.environ.copy()
        env["OMP_PROC_BIND"] = "true"
        env["OMP_PLACES"] = "threads"
        summary = json.loads(subprocess.check_output(command, text=True, env=env).splitlines()[-1])
        assert summary["status"] == "prototype"
        assert vcf.exists() and vcf.stat().st_size > 0
        assert Path(f"{vcf}.tbi").exists()
        assert manifest.exists()
        vcf_text = gzip.open(vcf, "rt", encoding="utf-8").read()
        assert "##fileformat=VCFv4.2" in vcf_text
        # hc-call uses a sample column in its prototype VCF.
        assert "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" not in vcf_text.splitlines(True)
        assert "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" in vcf_text
        # The native writer must preserve the sample identity from the BAM
        # @RG SM tag; falling back to a synthetic name would break direct
        # replacement in cohort/Nextflow channels.
        assert vcf_text.split("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", 1)[1].splitlines()[0] == "NA12878"
        metadata = json.loads(manifest.read_text())
        assert metadata["primary_output_kind"] == "vcf"
        assert all(item["complete"] for item in metadata["outputs"])
        assert metadata["compatibility"]["gatk_parameter_aliases"] is True
        assert metadata["compatibility"]["cigar_aware_projection"] is True
        assert metadata["compatibility"]["gvcf"] is False
        assert metadata["compatibility"]["vcf_index"] is True
        # GATK's explicit -ERC NONE must override a .g.vcf-looking output
        # name and remain an ordinary VCF.  This is distinct from the native
        # convenience inference used when no ERC mode was supplied.
        none_vcf = root / "calls.none.g.vcf.gz"
        none_manifest = root / "calls.none.g.vcf.gz.manifest.json"
        none_command = ["taskset", "-c", "0-3", str(binary), "-I", str(input_path),
                        "-R", str(workspace / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"),
                        "-L", "17:69000-69100", "-O", str(none_vcf), "-ERC", "NONE",
                        "--output-manifest", str(none_manifest), "--min-depth", "1",
                        "--min-alt-support", "1"]
        none_summary = json.loads(subprocess.check_output(none_command, text=True, env=env).splitlines()[-1])
        none_text = gzip.open(none_vcf, "rt", encoding="utf-8").read()
        assert "<NON_REF>" not in none_text
        assert Path(f"{none_vcf}.tbi").exists()
        none_metadata = json.loads(none_manifest.read_text())
        assert none_summary["gvcf_requested"] is False
        assert none_summary["gvcf_semantics"] == "not-requested"
        assert none_metadata["primary_output_kind"] == "vcf"
        assert none_metadata["compatibility"]["gvcf"] is False
        gvcf = root / "calls.g.vcf.gz"
        gvcf_manifest = root / "calls.g.vcf.gz.manifest.json"
        gvcf_command = ["taskset", "-c", "0-3", str(binary), "-I", str(input_path),
                        "-R", str(workspace / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"),
                        "-L", "17:69000-69100", "-O", str(gvcf), "-ERC", "GVCF",
                        "--output-manifest", str(gvcf_manifest), "--min-depth", "1",
                        "--min-alt-support", "1", "--native-pair-hmm-threads", "2"]
        gvcf_summary = json.loads(subprocess.check_output(gvcf_command, text=True, env=env).splitlines()[-1])
        gvcf_text = gzip.open(gvcf, "rt", encoding="utf-8").read()
        assert "##ALT=<ID=NON_REF" in gvcf_text
        assert "<NON_REF>" in gvcf_text
        assert Path(f"{gvcf}.tbi").exists()
        gvcf_metadata = json.loads(gvcf_manifest.read_text())
        assert gvcf_metadata["primary_output_kind"] == "gvcf"
        assert all(item["complete"] for item in gvcf_metadata["outputs"])
        assert gvcf_metadata["compatibility"]["gvcf"] is True
        assert gvcf_metadata["compatibility"]["vcf_index"] is True
        assert gvcf_metadata["compatibility"]["gvcf_semantics"] == "reference-blocks+candidate-sites"
        assert gvcf_metadata["telemetry"]["gvcf_reference_blocks"] > 0
        assert gvcf_summary["graph_reference_nodes"] > 0
        assert gvcf_summary["graph_reference_connected_nodes"] > 0
        assert gvcf_summary["graph_haplotype_sequence_count"] > 0
        assert gvcf_summary["read_filter_applied"] is True
        assert gvcf_summary["read_filter_min_mapping_quality"] == 20
        assert gvcf_summary["read_filter_exclude_mapping_quality_unavailable"] is True
        assert gvcf_summary["read_filter_exclude_mapping_quality_zero"] is False
        assert gvcf_summary["read_filter_exclude_duplicates"] is True
        assert gvcf_summary["read_filter_exclude_supplementary"] is False
        assert gvcf_summary["gvcf_interval_reference_blocks"] is True
        assert "END=" in gvcf_text
        # BP_RESOLUTION must not silently alias block mode: every simple
        # reference-confidence locus is a separate record with GATK's
        # GT:AD:DP:GQ:PL field order and no END/INFO block annotations.
        bp_gvcf = root / "calls.bp.g.vcf.gz"
        bp_manifest = root / "calls.bp.g.vcf.gz.manifest.json"
        bp_command = ["taskset", "-c", "0-3", str(binary), "-I", str(input_path),
                      "-R", str(workspace / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"),
                      "-L", "17:69000-69100", "-O", str(bp_gvcf), "-ERC", "BP_RESOLUTION",
                      "--output-manifest", str(bp_manifest), "--min-depth", "1",
                      "--min-alt-support", "1"]
        bp_summary = json.loads(subprocess.check_output(bp_command, text=True, env=env).splitlines()[-1])
        bp_text = gzip.open(bp_gvcf, "rt", encoding="utf-8").read()
        bp_records = [line.split("\t") for line in bp_text.splitlines()
                      if line and not line.startswith("#")]
        assert bp_records
        assert "##GVCFBlock" not in bp_text
        bp_reference_records = [row for row in bp_records if row[4] == "<NON_REF>"]
        assert bp_reference_records
        assert all(row[7] == "." and row[8] == "GT:AD:DP:GQ:PL"
                   for row in bp_reference_records)
        assert all("END=" not in row[7] for row in bp_reference_records)
        bp_candidates = [row for row in bp_records if row[4] != "<NON_REF>"]
        assert all("END=" not in row[7] for row in bp_candidates)
        bp_metadata = json.loads(bp_manifest.read_text())
        assert bp_metadata["compatibility"]["gvcf_semantics"] == \
            "base-pair-resolution+candidate-sites"
        assert bp_summary["gvcf_semantics"] == "base-pair-resolution+candidate-sites"
        return {"summary": summary, "vcf_bytes": vcf.stat().st_size,
                "gvcf_bytes": gvcf.stat().st_size,
                "bp_gvcf_bytes": bp_gvcf.stat().st_size,
                "bp_gvcf_records": len(bp_records)}


def run_projection_test(binary: Path, input_path: Path, region: str | None = None) -> dict:
    command = [str(binary), str(input_path)]
    if region:
        command.append(region)
    result = json.loads(subprocess.check_output(command, text=True).strip().splitlines()[-1])
    assert result["status"] == "pass"
    assert result["records"] > 0
    assert result["projected_bases"] <= result["bases"]
    assert result["cigar_ops"] > 0
    return result


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_HC_SMOKE_BINARY", root / "fastgatk-native/build/fastgatk-hc-smoke"))
    call_binary = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    projection_binary = Path(os.environ.get(
        "FASTGATK_HTS_READER_BINARY", root / "fastgatk-native/build/fastgatk-hts-reader-test"))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    cram = root / "gatk-source/src/test/resources/org/broadinstitute/hellbender/metrics/analysis/CollectQualityYieldMetrics/valid.cram"
    cram_ref = root / "gatk-source/src/test/resources/org/broadinstitute/hellbender/metrics/analysis/CollectQualityYieldMetrics/valid.fasta"
    full = run(binary, bam)
    region = run(binary, bam, region="17:69000-69100")
    cram_result = run(binary, cram, reference=cram_ref)
    projection_full = run_projection_test(projection_binary, bam)
    projection_region = run_projection_test(projection_binary, bam, "17:69000-69100")
    assert (full["reads"], full["bases"]) == (493, 37427)
    assert (region["reads"], region["bases"]) == (7, 532)
    assert (cram_result["reads"], cram_result["bases"]) == (8, 704)
    call_bam = run_call(call_binary, bam, batch_records=32)
    call_bam_rebatched = run_call(call_binary, bam, batch_records=493)
    call_cram = run_call(call_binary, cram, reference=cram_ref, batch_records=3)
    tuned_graph = run_call(
        call_binary, bam, batch_records=64,
        extra=["--kmer-size", "7", "--min-kmer-count", "1",
               "--max-num-haplotypes-in-population", "8", "--max-haplotype-depth", "64",
               "--dont-increase-kmer-sizes-for-cycles"])
    corrected_assembly = run_call(
        call_binary, bam, batch_records=64,
        extra=["--error-correct-reads",
               "--kmer-length-for-read-error-correction", "5",
               "--min-observations-for-kmer-to-be-solid", "2"])
    pileup_assembly = run_call(
        call_binary, bam, region="17:69000-69100", batch_records=32,
        extra=["--error-correction-log-odds", "3.5"])
    call_artifacts = run_call_artifacts(call_binary, bam)
    assert call_bam["reads"] == 493
    assert call_bam["observations"] == 37427
    assert call_bam["projected_observations"] <= projection_full["projected_bases"]
    assert call_bam["read_filter_applied"] is True
    # CIGAR-aware projection intentionally removes insertions/soft clips from
    # SNP pileup, and the deterministic read-filter mask removes flagged
    # records before every downstream kernel.
    # The default mask now matches GATK's MAPQ>=20 and duplicate-read policy;
    # the pre-assembly overlapping-mate quality correction is also enabled, so
    # keep the fixture count explicit so a filter/correction-order change is
    # visible.
    assert call_bam["candidate_sites"] == 2
    assert call_bam["overlapping_quality_correction_used"] is True
    assert call_bam["overlapping_quality_correction_metadata_available"] is True
    assert call_bam["overlapping_pairs"] == 37
    assert call_bam["overlapping_bases"] == 1002
    assert call_bam["overlapping_conflicting_bases"] == 13
    assert call_bam["graph_used"] is True
    # The no-reference smoke call intentionally exercises the explicit
    # diversity fallback.  The production graph path is reference-backed and
    # is checked below through the GVCF artifact call; without a FASTA there
    # is no coordinate scaffold, so a zero-node graph is an expected,
    # fail-closed telemetry state rather than a fabricated graph.
    assert call_bam["graph_reference_nodes"] == 0
    assert call_bam["graph_nodes"] >= 0 and call_bam["graph_edges"] >= 0
    assert call_bam["read_filter_used"] is True
    assert call_bam["pairhmm_used"] is False
    assert call_bam["pairhmm_skip_reason"] in {
        "no-reference", "no-candidates", "no-read-haplotype-overlap",
        "no-reference-backed-candidate"
    }
    # GATK's --disable-read-filter takes a class name; native applies the
    # supported class selectively instead of disabling the whole mask.
    selective = subprocess.check_output([
        str(call_binary), "-I", str(bam), "-O", "-",
        "--disable-read-filter", "NotDuplicateReadFilter"],
        text=True, env={**os.environ, "OMP_PROC_BIND": "true", "OMP_PLACES": "threads"}).splitlines()[-1]
    selective_summary = json.loads(selective)
    assert selective_summary["read_filter_exclude_duplicates"] is False
    assert selective_summary["read_filter_min_mapping_quality"] == 20

    # Explicitly adding a supported native filter must be reflected in the
    # effective mask; unknown filter classes fail closed instead of being
    # silently recorded and ignored.
    explicit = subprocess.check_output([
        str(call_binary), "-I", str(bam), "-O", "-", "--read-filter", "NotDuplicateReadFilter"
    ], text=True, env={**os.environ, "OMP_NUM_THREADS": "1"}).splitlines()[-1]
    explicit_summary = json.loads(explicit)
    assert explicit_summary["read_filter_exclude_duplicates"] is True
    good_cigar = subprocess.check_output([
        str(call_binary), "-I", str(bam), "-O", "-", "--read-filter", "GoodCigarReadFilter"
    ], text=True, env={**os.environ, "OMP_NUM_THREADS": "1"}).splitlines()[-1]
    good_cigar_summary = json.loads(good_cigar)
    assert good_cigar_summary["read_filter_require_good_cigar"] is True
    assert good_cigar_summary["read_filter_used"] is True
    length_filter = subprocess.run([
        str(call_binary), "-I", str(bam), "-O", "-", "--read-filter", "ReadLengthReadFilter",
        "--min-read-length", "80", "--max-read-length", "100"
    ], text=True, capture_output=True, env={**os.environ, "OMP_NUM_THREADS": "1"})
    assert length_filter.returncode == 0, (length_filter.stdout, length_filter.stderr)
    length_summary = json.loads(length_filter.stdout.strip().splitlines()[-1])
    assert length_summary["read_filter_require_read_length"] is True
    assert length_summary["read_filter_min_read_length"] == 80
    assert length_summary["read_filter_max_read_length"] == 100
    assert length_summary["filtered_reads"] == 493
    # --disable-tool-default-read-filters removes only the implicit GATK
    # policy.  The shared Kokkos mask remains active and explicit predicates
    # must survive regardless of argument ordering.
    no_defaults = subprocess.check_output([
        str(call_binary), "-I", str(bam), "-O", "-",
        "--disable-tool-default-read-filters"], text=True,
        env={**os.environ, "OMP_NUM_THREADS": "1"}).splitlines()[-1]
    no_defaults_summary = json.loads(no_defaults)
    for field in ("read_filter_min_mapping_quality", "read_filter_exclude_duplicates",
                  "read_filter_exclude_unmapped", "read_filter_exclude_secondary",
                  "read_filter_exclude_supplementary", "read_filter_exclude_qcfail"):
        assert no_defaults_summary[field] in (0, False), (field, no_defaults_summary[field])
    assert no_defaults_summary["read_filter_applied"] is True
    explicit_no_defaults = subprocess.check_output([
        str(call_binary), "-I", str(bam), "-O", "-",
        "--disable-tool-default-read-filters", "--read-filter", "NotDuplicateReadFilter",
        "--read-filter", "MappedReadFilter"], text=True,
        env={**os.environ, "OMP_NUM_THREADS": "1"}).splitlines()[-1]
    explicit_no_defaults_summary = json.loads(explicit_no_defaults)
    assert explicit_no_defaults_summary["read_filter_min_mapping_quality"] == 0
    assert explicit_no_defaults_summary["read_filter_exclude_duplicates"] is True
    assert explicit_no_defaults_summary["read_filter_exclude_unmapped"] is True
    reverse_no_defaults = subprocess.check_output([
        str(call_binary), "-I", str(bam), "-O", "-",
        "--read-filter", "NotDuplicateReadFilter", "--disable-tool-default-read-filters"],
        text=True, env={**os.environ, "OMP_NUM_THREADS": "1"}).splitlines()[-1]
    reverse_no_defaults_summary = json.loads(reverse_no_defaults)
    assert reverse_no_defaults_summary["read_filter_min_mapping_quality"] == 0
    assert reverse_no_defaults_summary["read_filter_exclude_duplicates"] is True
    unsupported = subprocess.run([
        str(call_binary), "-I", str(bam), "-O", "-", "--read-filter", "UnsupportedReadFilter"
    ], text=True, capture_output=True, env={**os.environ, "OMP_NUM_THREADS": "1"})
    assert unsupported.returncode != 0
    assert "UNSUPPORTED_PARAMETER" in unsupported.stderr
    assert call_bam["signature"] == call_bam_rebatched["signature"]
    # The GATK-compatible default read-filter mask removes one duplicate from
    # this CRAM fixture before HC batching (8 physical records -> 7 reads).
    assert call_cram["reads"] == 7
    # The permissive HC CLI defaults may expose a reference-backed candidate in
    # this CRAM fixture; in that case the production PairHMM path is expected.
    assert call_cram["pairhmm_skip_reason"] in {
        "no-reference", "no-candidates", "no-read-haplotype-overlap",
        "no-reference-backed-candidate", "executed"
    }
    assert call_cram["pairhmm_used"] == (call_cram["pairhmm_skip_reason"] == "executed")
    assert tuned_graph["graph_kmer_size"] == 7
    assert tuned_graph["graph_kmer_size_selected"] == 7
    assert tuned_graph["graph_kmer_iterations"] == 1
    assert tuned_graph["graph_dont_increase_kmer_sizes_for_cycles"] is True
    assert tuned_graph["graph_min_kmer_count"] == 1
    assert tuned_graph["graph_max_paths"] == 8
    assert tuned_graph["graph_max_depth"] == 64
    assert tuned_graph["graph_dangling_recovered_paths"] >= 0
    assert tuned_graph["graph_dangling_recovered_bases"] >= 0
    assert tuned_graph["graph_seqgraph_nodes"] >= 0
    assert tuned_graph["graph_seqgraph_edges"] >= 0
    assert tuned_graph["graph_non_unique_kmers"] >= 0
    assert corrected_assembly["error_correction_used"] is True
    assert corrected_assembly["error_correction_kmer_length"] == 5
    assert corrected_assembly["error_correction_min_solid_observations"] == 2
    assert corrected_assembly["error_correction_uncorrectable_kmers"] >= 0
    assert corrected_assembly["error_correction_execution_space"] in {"OpenMP", "Serial", "Threads"}
    assert pileup_assembly["error_correction_used"] is True
    assert pileup_assembly["error_correction_mode"] == "pileup"
    assert pileup_assembly["error_correction_pileup_loci"] >= 0
    assert pileup_assembly["error_correction_pileup_skipped_indel_adjacent_bases"] >= 0
    print(json.dumps({"status": "pass", "bam": full, "bam_region": region,
                      "cram": cram_result, "calling_bam": call_bam,
                      "calling_cram": call_cram, "calling_artifacts": call_artifacts,
                      "calling_pileup": pileup_assembly,
                      "projection": {"full": projection_full, "region": projection_region}}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
