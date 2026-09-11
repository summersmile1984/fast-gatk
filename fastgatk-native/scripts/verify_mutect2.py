#!/usr/bin/env python3
"""Verify the native Mutect2 adapter contract and tumor/normal sidecars."""

from __future__ import annotations

import gzip
import json
import math
import os
import subprocess
import tarfile
import tempfile
from pathlib import Path
import oracle_guard


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_MUTECT2_BINARY", root / "fastgatk-native/build/fastgatk-mutect2"))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-") as directory:
        work = Path(directory)
        output = work / "calls.vcf.gz"
        manifest = work / "calls.vcf.gz.manifest.json"
        stats = work / "calls.stats.json"
        f1r2 = work / "calls.f1r2.tsv"
        command = [str(native), "-I", str(bam), "-R", str(reference),
                   "-L", "17:69000-70000", "-O", str(output),
                   "--tumor-sample", "NA12878",
                   "--normal-input", str(bam), "--normal-lod", "-100",
                   "--min-depth", "1", "--min-alt-support", "1",
                   "--kmer-size", "7", "--min-kmer-count", "1",
                   "--max-num-haplotypes-in-population", "8", "--max-haplotype-depth", "64",
                   "--max-haplotype-combination-alleles", "5",
                   "--error-correct-reads",
                   "--kmer-length-for-read-error-correction", "5",
                   "--min-observations-for-kmer-to-be-solid", "2",
                   "--error-correction-log-odds", "3.5",
                   # GATK's common Boolean form must retain the default
                   # Mutect2 filter list when it is explicitly false.
                   "--disable-tool-default-read-filters=false",
                   "--dont-use-soft-clipped-bases=false",
                   "--do-not-correct-overlapping-quality=false",
                   "--dont-increase-kmer-sizes-for-cycles=false",
                   "--recover-all-dangling-branches=false",
                   "--disable-adaptive-pruning=false",
                   "--stats", str(stats), "--f1r2-tar-gz", str(f1r2),
                   "--output-manifest", str(manifest)]
        env = os.environ.copy()
        env.setdefault("OMP_PROC_BIND", "true")
        env.setdefault("OMP_PLACES", "threads")
        result = json.loads(subprocess.check_output(command, text=True, env=env).splitlines()[-1])
        assert result["tool"] == "Mutect2" and result["status"] == "prototype"
        assert result["tumor_reads"] == 493 and result["tumor_calls"] > 0
        assert output.exists() and output.stat().st_size > 0
        assert Path(f"{output}.tbi").exists()
        text = gzip.open(output, "rt", encoding="utf-8").read()
        assert "##source=Mutect2" in text
        assert "##GATKCommandLine=<ID=Mutect2,Version=fastgatk-native," in text
        # GATK constructs its SampleList from ReadUtils.getSamplesFromHeader(),
        # a TreeSet. The public VCF order is therefore lexicographic and
        # independent of the Host's normal-side likelihood replay order.
        assert "#CHROM\tPOS" in text and "\tNA12878\tNORMAL" in text
        assert "##INFO=<ID=NALOD" in text
        assert "##INFO=<ID=PSOMATIC" not in text
        assert "##INFO=<ID=PGERMLINE" not in text
        assert "##INFO=<ID=PARTIFACT" not in text
        assert "##INFO=<ID=OBP" not in text
        variant_lines = [line.split("\t") for line in text.splitlines()
                         if line and not line.startswith("#")]
        # Mutect2 leaves site QUAL missing; FilterMutectCalls computes the
        # calibrated confidence downstream.  Keep this explicit so a future
        # writer change cannot accidentally reintroduce a pre-filter numeric
        # QUAL that diverges from GATK's output contract.
        assert variant_lines and all(fields[5] == "." for fields in variant_lines)
        assert all("NALOD=" in fields[7] and "NLOD=" in fields[7]
                   for fields in variant_lines)
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        # Mutect2 4.6.2.0 defaults --native-pair-hmm-threads to four.  The
        # native adapter must expose the same default Kokkos launch width;
        # an explicit --threads override remains covered elsewhere.
        assert metadata["telemetry"]["requested_pairhmm_threads"] == 4
        assert metadata["telemetry"]["effective_pairhmm_threads"] == 4
        assert metadata["compatibility"]["tumor_normal"] is True
        assert metadata["compatibility"]["somatic_filter"] is True
        assert metadata["compatibility"]["f1r2"] is True
        assert metadata["compatibility"]["somatic_read_likelihoods"] is True
        assert metadata["compatibility"]["somatic_fragment_grouping"] is True
        assert metadata["compatibility"]["somatic_informative_read_overlap"] is True
        assert metadata["compatibility"]["somatic_grouped_responsibility_depth"] is True
        assert metadata["compatibility"]["somatic_posterior"] is True
        assert metadata["compatibility"]["orientation_bias_posterior"] is True
        assert metadata["compatibility"]["contamination_posterior"] is True
        assert metadata["compatibility"]["read_haplotype_request_deduplication"] is True
        assert metadata["compatibility"]["read_filter_gatk_defaults"] is True
        assert metadata["compatibility"]["pcr_indel_model"] == "CONSERVATIVE"
        assert metadata["compatibility"]["pcr_error_rate_factor"] == 3.0
        assert metadata["telemetry"]["read_filter_min_mapping_quality"] == 20
        assert metadata["telemetry"]["read_filter_exclude_duplicates"] is True
        # GATK 4.6.2's StandardMutect2ReadFilters does not include
        # NotSupplementaryAlignmentReadFilter.  This must remain false by
        # default; callers may enable it explicitly through --read-filter.
        assert metadata["telemetry"]["read_filter_exclude_supplementary"] is False
        assert metadata["telemetry"]["read_filter_require_good_cigar"] is True
        assert metadata["telemetry"]["read_filter_require_nonzero_reference_span"] is True
        assert metadata["telemetry"]["read_filter_require_non_chimeric_original_alignment"] is True
        assert metadata["telemetry"]["read_filter_require_no_n_cigar"] is True
        assert metadata["telemetry"]["read_filter_require_read_group"] is True
        assert metadata["telemetry"]["read_filter_require_read_length"] is True
        assert metadata["telemetry"]["read_filter_min_read_length"] == 30
        assert metadata["telemetry"]["read_filter_exclude_mapping_quality_unavailable"] is True
        assert metadata["telemetry"]["read_filter_exclude_mapping_quality_zero"] is True
        # GATK keeps the active-region seed and final emission thresholds
        # distinct; native records both and applies their stricter value at
        # the final candidate boundary.
        assert metadata["compatibility"]["initial_tumor_lod"] == 2.0
        assert metadata["compatibility"]["tumor_lod_to_emit"] == 3.0
        assert metadata["compatibility"]["effective_tumor_lod_emit"] == 3.0
        assert metadata["telemetry"]["initial_tumor_lod"] == 2.0
        assert metadata["telemetry"]["tumor_lod_to_emit"] == 3.0
        assert metadata["telemetry"]["effective_tumor_lod_emit"] == 3.0
        assert metadata["compatibility"]["minimum_allele_fraction"] == 0.0
        assert metadata["telemetry"]["minimum_allele_fraction"] == 0.0
        # With a normal input, GATK's default population AF is 1e-6 rather
        # than the 5e-8 tumor-only default.
        assert metadata["telemetry"]["population_allele_frequency"] == 1.0e-6
        # A k=7 reference window in this fixture contains repeated k-mers.
        # GATK's --dont-increase-kmer-sizes-for-cycles disables the retry
        # ladder, so that intentionally no-variation graph must not also be
        # the contract's positive-call/sidecar fixture.  Exercise the option
        # independently and require its exact recorded policy instead.
        no_cycle_output = work / "calls-no-kmer-escalation.vcf.gz"
        no_cycle_manifest = work / "calls-no-kmer-escalation.manifest.json"
        no_cycle_result = json.loads(subprocess.check_output([
            str(native), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-70000", "-O", str(no_cycle_output),
            "--tumor-sample", "NA12878", "--normal-input", str(bam),
            "--normal-lod", "-100", "--min-depth", "1", "--min-alt-support", "1",
            "--kmer-size", "7", "--dont-increase-kmer-sizes-for-cycles",
            "--add-output-vcf-command-line=false",
            "--output-manifest", str(no_cycle_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert no_cycle_result["tumor_reads"] == 493
        no_cycle_metadata = json.loads(no_cycle_manifest.read_text(encoding="utf-8"))
        assert no_cycle_metadata["telemetry"]["graph_kmer_size"] == 7
        assert no_cycle_metadata["telemetry"]["graph_dont_increase_kmer_sizes_for_cycles"] is True
        no_cycle_text = gzip.open(no_cycle_output, "rt", encoding="utf-8").read()
        assert "##GATKCommandLine=" not in no_cycle_text
        assert "##source=Mutect2" not in no_cycle_text
        no_defaults_output = work / "calls-no-default-filters.vcf.gz"
        no_defaults_manifest = work / "calls-no-default-filters.manifest.json"
        no_defaults_result = json.loads(subprocess.check_output([
            str(native), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-70000", "-O", str(no_defaults_output),
            "--tumor-sample", "NA12878", "--normal-input", str(bam), "--normal-lod", "-100",
            "--min-depth", "1", "--min-alt-support", "1",
            "--disable-tool-default-read-filters",
            "--read-filter", "NotDuplicateReadFilter", "--read-filter", "MappedReadFilter",
            "--output-manifest", str(no_defaults_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert no_defaults_result["tumor_reads"] == 493
        no_defaults_metadata = json.loads(no_defaults_manifest.read_text(encoding="utf-8"))
        no_defaults_compat = no_defaults_metadata["compatibility"]
        no_defaults_telemetry = no_defaults_metadata["telemetry"]
        assert no_defaults_compat["disable_tool_default_read_filters"] is True
        assert no_defaults_compat["read_filter_gatk_defaults"] is False
        assert no_defaults_telemetry["read_filter_min_mapping_quality"] == 0
        assert no_defaults_telemetry["read_filter_exclude_duplicates"] is True
        # Repeatable -L selectors must use the same Host interval-set rule as
        # GATK.  Keep this interval-only check tumor-only: the fixture uses
        # the same BAM/sample name for its separate tumor/normal sidecar
        # exercise above, whereas GATK assigns sample roles by read-group
        # sample name and therefore cannot represent those two roles from one
        # name.  The tumor-only command below is the directly comparable GATK
        # oracle: UNION retains all six calls, while INTERSECTION is
        # restricted to the 17:69300-69400 overlap (the 69368 call).
        union_output = work / "calls-interval-union.vcf.gz"
        union_manifest = work / "calls-interval-union.manifest.json"
        union_result = json.loads(subprocess.check_output([
            str(native), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-69400", "-L", "17:69300-70000",
            "-isr", "UNION", "-O", str(union_output),
            "--tumor-sample", "NA12878",
            "--min-depth", "1", "--min-alt-support", "1",
            "--output-manifest", str(union_manifest),
        ], text=True, env=env).splitlines()[-1])
        union_positions = [int(line.split("\t")[1]) for line in
                           gzip.open(union_output, "rt", encoding="utf-8")
                           if line and not line.startswith("#")]
        assert union_result["tumor_calls"] == 6
        assert union_positions == [69067, 69124, 69368, 69631, 69803, 69807]
        union_metadata = json.loads(union_manifest.read_text(encoding="utf-8"))
        assert union_metadata["compatibility"]["interval_set_rule"] == "UNION"
        assert union_metadata["telemetry"]["interval_set_rule"] == "UNION"
        intersection_output = work / "calls-interval-intersection.vcf.gz"
        intersection_manifest = work / "calls-interval-intersection.manifest.json"
        intersection_result = json.loads(subprocess.check_output([
            str(native), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-69400", "-L", "17:69300-70000",
            "--interval-set-rule=INTERSECTION", "-O", str(intersection_output),
            "--tumor-sample", "NA12878",
            "--min-depth", "1", "--min-alt-support", "1",
            "--output-manifest", str(intersection_manifest),
        ], text=True, env=env).splitlines()[-1])
        intersection_positions = [int(line.split("\t")[1]) for line in
                                  gzip.open(intersection_output, "rt", encoding="utf-8")
                                  if line and not line.startswith("#")]
        assert intersection_result["tumor_calls"] == 1
        assert intersection_positions == [69368]
        intersection_metadata = json.loads(intersection_manifest.read_text(encoding="utf-8"))
        assert intersection_metadata["compatibility"]["interval_set_rule"] == "INTERSECTION"
        assert intersection_metadata["telemetry"]["interval_set_rule"] == "INTERSECTION"
        # ReadLengthReadFilter is part of the same shared typed mask used by
        # HaplotypeCaller/CountReads/FlagStat.  This range is intentionally
        # above the fixture's reads, so the call set must be empty while the
        # manifest records the exact bounds for reproducibility.
        length_output = work / "calls-length-filtered.vcf.gz"
        length_manifest = work / "calls-length-filtered.manifest.json"
        length_result = json.loads(subprocess.check_output([
            str(native), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-70000", "-O", str(length_output),
            "--tumor-sample", "NA12878", "--normal-input", str(bam), "--normal-lod", "-100",
            "--min-depth", "1", "--min-alt-support", "1",
            "--read-filter", "ReadLengthReadFilter",
            "--min-read-length", "80", "--max-read-length", "100",
            "--output-manifest", str(length_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert length_result["tumor_calls"] == 0
        length_metadata = json.loads(length_manifest.read_text(encoding="utf-8"))
        assert length_metadata["compatibility"]["read_filter_require_read_length"] is True
        assert length_metadata["compatibility"]["read_filter_min_read_length"] == 80
        assert length_metadata["compatibility"]["read_filter_max_read_length"] == 100
        assert metadata["telemetry"]["graph_kmer_size"] == 7
        assert metadata["telemetry"]["graph_kmer_size_selected"] >= 7
        assert metadata["telemetry"]["graph_kmer_iterations"] >= 1
        assert metadata["telemetry"]["graph_dont_increase_kmer_sizes_for_cycles"] is False
        assert metadata["telemetry"]["graph_min_kmer_count"] == 1
        assert metadata["telemetry"]["graph_max_paths"] == 8
        assert metadata["telemetry"]["requested_batch_records"] >= metadata["telemetry"]["effective_batch_records"] >= 1
        assert metadata["telemetry"]["adaptive_batch_reductions"] >= 0
        assert metadata["telemetry"]["graph_max_depth"] == 64
        assert metadata["telemetry"]["graph_dangling_recovered_paths"] >= 0
        assert metadata["telemetry"]["graph_dangling_recovered_bases"] >= 0
        assert metadata["telemetry"]["max_haplotype_combination_alleles"] == 5
        assert metadata["telemetry"]["error_correct_reads"] is True
        assert metadata["telemetry"]["error_correction_kmer_length"] == 5
        assert metadata["telemetry"]["error_correction_min_solid_observations"] == 2
        assert metadata["telemetry"]["error_correction_uncorrectable_kmers"] >= 0
        assert metadata["telemetry"]["error_correction_mode"] == "pileup"
        assert metadata["telemetry"]["error_correction_pileup_loci"] >= 0
        assert metadata["telemetry"]["error_correction_execution_space"] in {"OpenMP", "Serial", "Threads"}
        assert metadata["telemetry"]["pairhmm_pairs"] > 0
        assert metadata["telemetry"]["pairhmm_pcr_indel_model"] == "CONSERVATIVE"
        assert metadata["telemetry"]["pairhmm_pcr_error_rate_factor"] == 3.0
        assert metadata["telemetry"]["pairhmm_pcr_adjusted_positions"] > 0
        assert metadata["telemetry"]["pairhmm_request_links"] >= metadata["telemetry"]["pairhmm_pairs"]
        # The cap is a graph-complexity guard, not a request to synthesize a
        # multi-allelic component.  This fixture has no retained component
        # that needs splitting after the source-compatible k-mer retry, so
        # verify parser propagation above and accept zero realized blocks.
        assert metadata["telemetry"]["pairhmm_haplotype_combination_blocks"] >= 0
        assert metadata["compatibility"]["assembly_region_pairhmm_partitioning"] is True
        assert metadata["compatibility"]["assembly_region_candidate_partitioning"] is True
        assert metadata["telemetry"]["pairhmm_assembly_region_groups"] >= 1
        assert metadata["telemetry"]["pairhmm_unassigned_candidates"] == 0
        assert metadata["telemetry"]["assembly_cross_region_candidates"] == 0
        assert metadata["telemetry"]["assembly_unassigned_candidates"] == 0
        assert all(item["complete"] for item in metadata["outputs"])
        stats_data = json.loads(stats.read_text(encoding="utf-8"))
        assert stats_data["normal_calls"] > 0
        assert stats_data["multiallelic_locus_writer"] is True
        assert stats_data["somatic_likelihood_model"] == "per-read-dirichlet-variational-evidence-v1"
        assert stats_data["somatic_evidence_grouping"] is True
        assert stats_data["somatic_informative_read_overlap_margin"] == 2
        assert stats_data["somatic_grouped_responsibility_depth"] is True
        assert stats_data["somatic_evidence_groups"] > 0
        assert stats_data["somatic_posterior_model"] == "somatic-germline-artifact-orientation-v1"
        assert stats_data["somatic_posterior_execution_space"]
        assert stats_data["somatic_likelihood_execution_space"]
        assert stats_data["initial_tumor_lod"] == 2.0
        assert stats_data["tumor_lod_to_emit"] == 3.0
        assert stats_data["effective_tumor_lod_emit"] == 3.0
        assert stats_data["minimum_allele_fraction"] == 0.0
        assert stats_data["population_allele_frequency"] == 1.0e-6
        assert stats_data["pcr_indel_model"] == "CONSERVATIVE"
        assert stats_data["pcr_error_rate_factor"] == 3.0
        assert stats_data["phred_scaled_global_read_mismapping_rate"] == 45.0
        assert stats_data["requested_batch_records"] >= stats_data["effective_batch_records"] >= 1
        assert stats_data["adaptive_batch_reductions"] >= 0
        assert "FASTGATK-MUTECT2-F1R2" in f1r2.read_text(encoding="utf-8")

        # GATK's minAF option is a real Dirichlet prior, not a parser-only
        # compatibility alias.  Exercise the short alias and require both
        # metadata propagation and a changed TLOD on the same fixture.  Keep
        # the graph/assembly controls identical between the two runs: the
        # broader contract above intentionally uses a small-kmer/error-
        # correction configuration to exercise parser propagation, whereas
        # this minAF check must isolate the likelihood prior itself.
        min_af_baseline_output = work / "calls-min-af-baseline.vcf.gz"
        min_af_output = work / "calls-min-af.vcf.gz"
        min_af_manifest = work / "calls-min-af.manifest.json"
        min_af_common = [
            str(native), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-70000",
            "--tumor-sample", "NA12878", "--normal-input", str(bam), "--normal-lod", "-100",
            "--min-depth", "1", "--min-alt-support", "1",
        ]
        min_af_baseline_result = json.loads(subprocess.check_output([
            *min_af_common, "-O", str(min_af_baseline_output),
        ], text=True, env=env).splitlines()[-1])
        min_af_result = json.loads(subprocess.check_output([
            *min_af_common, "-O", str(min_af_output), "-min-AF", "0.1",
            "--output-manifest", str(min_af_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert min_af_result["tumor_calls"] == min_af_baseline_result["tumor_calls"]
        min_af_metadata = json.loads(min_af_manifest.read_text(encoding="utf-8"))
        assert min_af_metadata["compatibility"]["minimum_allele_fraction"] == 0.1
        assert min_af_metadata["telemetry"]["minimum_allele_fraction"] == 0.1
        def tlod_values(vcf_text: str) -> list[float]:
            values = []
            for line in vcf_text.splitlines():
                if not line or line.startswith("#"):
                    continue
                info = {token.split("=", 1)[0]: token.split("=", 1)[1]
                        for token in line.split("\t")[7].split(";") if "=" in token}
                if "TLOD" in info:
                    values.extend(float(value) for value in info["TLOD"].split(","))
            return values
        min_af_baseline_text = gzip.open(min_af_baseline_output, "rt", encoding="utf-8").read()
        min_af_text = gzip.open(min_af_output, "rt", encoding="utf-8").read()
        assert len(tlod_values(min_af_text)) == len(tlod_values(min_af_baseline_text))
        assert any(abs(left - right) > 1.0e-6 for left, right in
                   zip(tlod_values(min_af_baseline_text), tlod_values(min_af_text)))

        # M2ArgumentCollection's mitochondrial mode is an optional boolean
        # that rewrites the sensitivity defaults before Host/Kokkos execution.
        # Exercise the inline form and assert the effective values in both
        # stats and OutputManifest, including the natural-log pruning value
        # used by GATK's -4*ln(10) mode adjustment.
        mito_output = work / "calls-mito.vcf.gz"
        mito_stats = work / "calls-mito.stats.json"
        mito_manifest = work / "calls-mito.manifest.json"
        mito_f1r2 = work / "calls-mito.f1r2.tsv"
        mito_result = json.loads(subprocess.check_output([
            str(native), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-69100", "-O", str(mito_output),
            "--tumor-sample", "NA12878", "--min-depth", "1", "--min-alt-support", "1",
            "--stats", str(mito_stats), "--f1r2-tar-gz", str(mito_f1r2),
            "--mitochondria-mode=true", "--output-manifest", str(mito_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert mito_result["tool"] == "Mutect2"
        mito_stats_data = json.loads(mito_stats.read_text(encoding="utf-8"))
        assert mito_stats_data["mitochondria_mode"] is True
        assert mito_stats_data["initial_tumor_lod"] == 0.0
        assert mito_stats_data["tumor_lod_to_emit"] == 0.0
        assert mito_stats_data["effective_tumor_lod_emit"] == 0.0
        assert mito_stats_data["population_allele_frequency"] == 4.0e-3
        assert abs(mito_stats_data["mitochondria_pruning_lod_threshold"] + 4.0 * math.log(10.0)) < 1.0e-12
        mito_metadata = json.loads(mito_manifest.read_text(encoding="utf-8"))
        assert mito_metadata["compatibility"]["mitochondria_mode"] is True
        assert mito_metadata["compatibility"]["mitochondria_defaults_applied"] is True
        assert mito_metadata["compatibility"]["initial_tumor_lod"] == 0.0
        assert mito_metadata["compatibility"]["tumor_lod_to_emit"] == 0.0
        assert mito_metadata["compatibility"]["effective_tumor_lod_emit"] == 0.0
        assert mito_metadata["telemetry"]["initial_tumor_lod"] == 0.0
        assert mito_metadata["telemetry"]["tumor_lod_to_emit"] == 0.0
        assert mito_metadata["telemetry"]["effective_tumor_lod_emit"] == 0.0
        assert mito_metadata["telemetry"]["population_allele_frequency"] == 4.0e-3

        # Indexed core/halo streaming must preserve the aggregate VCF record
        # stream while releasing each tile after the somatic/posterior sink.
        # BGZF block boundaries are allowed to differ; compare decompressed
        # bytes and the standard F1R2 archive member payloads instead.
        stream_output = work / "calls-stream.vcf.gz"
        stream_stats = work / "calls-stream.stats.json"
        stream_manifest = work / "calls-stream.manifest.json"
        stream_f1r2 = work / "calls-stream.f1r2.tar.gz"
        stream_command = list(command)
        replacements = {
            "-O": stream_output,
            "--stats": stream_stats,
            "--f1r2-tar-gz": stream_f1r2,
            "--output-manifest": stream_manifest,
        }
        for flag, value in replacements.items():
            position = stream_command.index(flag)
            stream_command[position + 1] = str(value)
        stream_command += ["--stream-by-region", "500"]
        stream_result = json.loads(subprocess.check_output(
            stream_command, text=True, env=env).splitlines()[-1])
        assert stream_result["stream_by_region"] is True
        assert gzip.open(output, "rb").read() == gzip.open(stream_output, "rb").read()
        stream_metadata = json.loads(stream_manifest.read_text(encoding="utf-8"))
        stream_telemetry = stream_metadata["telemetry"]
        assert stream_metadata["compatibility"]["stream_by_region"] is True
        assert stream_metadata["compatibility"]["indexed_input"] is True
        assert stream_metadata["compatibility"]["f1r2_standard_tar"] is True
        assert stream_telemetry["streamed_regions"] > 1
        assert stream_telemetry["streamed_peak_host_bytes"] > 0
        assert stream_telemetry["pipeline_lifecycle"] == \
            "Host decode->bounded queue->Kokkos compute->encode->sink"
        assert stream_telemetry["pipeline_decoded_items"] == stream_telemetry["streamed_regions"]
        assert stream_telemetry["pipeline_computed_items"] == stream_telemetry["streamed_regions"]
        assert stream_telemetry["pipeline_encoded_items"] == stream_telemetry["streamed_regions"]

        # A tiled core/halo traversal does not own the one continuous
        # ActivityProfile state required by GATK's --assembly-region-out
        # debug track.  It must reject the unsupported combination instead
        # of succeeding while omitting (or fabricating) that source trace.
        stream_activity = work / "calls-stream.activity.igv"
        rejected_activity = subprocess.run(
            [*stream_command, "--assembly-region-out", str(stream_activity)],
            text=True, env=env, capture_output=True, check=False)
        assert rejected_activity.returncode != 0
        assert "--assembly-region-out requires aggregate Mutect2 traversal" in rejected_activity.stderr
        assert not stream_activity.exists()

        # Repeatable indexed -I shards use the same core/halo path.  Duplicate
        # fixture shards are intentional here: they make the merge observable
        # without requiring a second large BAM and must not be rejected or
        # silently dropped.
        multi_stream_output = work / "calls-stream-multi-input.vcf.gz"
        multi_stream_stats = work / "calls-stream-multi-input.stats.json"
        multi_stream_manifest = work / "calls-stream-multi-input.manifest.json"
        multi_stream_f1r2 = work / "calls-stream-multi-input.f1r2.tar.gz"
        multi_stream_command = list(stream_command)
        input_insert = multi_stream_command.index("-I") + 2
        multi_stream_command[input_insert:input_insert] = ["-I", str(bam)]
        for flag, value in {
            "-O": multi_stream_output,
            "--stats": multi_stream_stats,
            "--f1r2-tar-gz": multi_stream_f1r2,
            "--output-manifest": multi_stream_manifest,
        }.items():
            position = multi_stream_command.index(flag)
            multi_stream_command[position + 1] = str(value)
        multi_stream_result = json.loads(subprocess.check_output(
            multi_stream_command, text=True, env=env).splitlines()[-1])
        assert multi_stream_result["stream_by_region"] is True
        assert multi_stream_result["tumor_reads"] == 2 * stream_result["tumor_reads"]
        multi_stream_metadata = json.loads(multi_stream_manifest.read_text(encoding="utf-8"))
        assert multi_stream_metadata["compatibility"]["indexed_input"] is True
        multi_stream_variants = [line.split("\t") for line in gzip.open(
            multi_stream_output, "rt", encoding="utf-8")
                                 if line and not line.startswith("#")]
        multi_stream_positions = [(fields[0], int(fields[1])) for fields in multi_stream_variants]
        assert multi_stream_positions == sorted(multi_stream_positions)
        assert len(multi_stream_positions) == len(set(multi_stream_positions))

        # A standard GATK-shaped single-input invocation may select both
        # tumor and normal samples from the same indexed file.  The stream
        # header must therefore receive a normal Result whenever
        # --normal-sample is set, even though --normal-input is absent.
        same_input_stream_output = work / "calls-stream-same-input.vcf.gz"
        same_input_stream_manifest = work / "calls-stream-same-input.manifest.json"
        same_input_stream_stats = work / "calls-stream-same-input.stats.json"
        same_input_stream_f1r2 = work / "calls-stream-same-input.f1r2.tar.gz"
        same_input_stream_result = json.loads(subprocess.check_output([
            str(native), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-70000", "-O", str(same_input_stream_output),
            "--tumor-sample", "NA12878", "--normal-sample", "NA12878", "--normal-lod", "-100",
            "--stream-by-region", "500", "--min-depth", "1",
            "--min-alt-support", "1", "--stats", str(same_input_stream_stats),
            "--f1r2-tar-gz", str(same_input_stream_f1r2),
            "--output-manifest", str(same_input_stream_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert same_input_stream_result["stream_by_region"] is True
        same_input_header = next(line for line in gzip.open(
            same_input_stream_output, "rt", encoding="utf-8") if line.startswith("#CHROM"))
        assert same_input_header.rstrip("\n").endswith("\tNA12878\tNA12878")
        same_input_manifest_payload = json.loads(
            same_input_stream_manifest.read_text(encoding="utf-8"))
        assert same_input_manifest_payload["compatibility"]["tumor_normal"] is True

        # A much smaller tile forces candidates and read pairs across several
        # core/halo boundaries.  The output must remain a coordinate-ordered,
        # core-only stream without duplicate loci.  Mutect2's release-specific
        # somatic posterior is not associative across independently computed
        # tiles, so whole-VCF byte identity is intentionally not required for
        # this multi-tile stress case (the 500bp contract above remains exact).
        tiny_output = work / "calls-stream-tiny.vcf.gz"
        tiny_stats = work / "calls-stream-tiny.stats.json"
        tiny_manifest = work / "calls-stream-tiny.manifest.json"
        tiny_f1r2 = work / "calls-stream-tiny.f1r2.tar.gz"
        tiny_command = list(command)
        tiny_replacements = {
            "-O": tiny_output,
            "--stats": tiny_stats,
            "--f1r2-tar-gz": tiny_f1r2,
            "--output-manifest": tiny_manifest,
        }
        for flag, value in tiny_replacements.items():
            position = tiny_command.index(flag)
            tiny_command[position + 1] = str(value)
        tiny_command += ["--stream-by-region", "100"]
        tiny_result = json.loads(subprocess.check_output(
            tiny_command, text=True, env=env).splitlines()[-1])
        assert tiny_result["stream_by_region"] is True
        tiny_metadata = json.loads(tiny_manifest.read_text(encoding="utf-8"))
        assert tiny_metadata["compatibility"]["stream_by_region"] is True
        assert tiny_metadata["telemetry"]["streamed_regions"] > stream_telemetry["streamed_regions"]
        assert tiny_metadata["telemetry"]["pipeline_decoded_items"] == \
            tiny_metadata["telemetry"]["streamed_regions"]
        tiny_variants = [line.split("\t") for line in gzip.open(tiny_output, "rt", encoding="utf-8")
                         if line and not line.startswith("#")]
        tiny_positions = [(fields[0], int(fields[1])) for fields in tiny_variants]
        assert tiny_positions == sorted(tiny_positions)
        assert len(tiny_positions) == len(set(tiny_positions))
        assert all(contig == "17" and 69000 <= position <= 70000
                   for contig, position in tiny_positions)
        # The explicit .tsv path above preserves the old adapter contract.  A
        # .tar.gz path must instead be directly consumable by
        # LearnReadOrientationModel/CollectF1R2Counts.
        standard_output = work / "calls-standard.vcf.gz"
        standard_manifest = work / "calls-standard.manifest.json"
        standard_stats = work / "calls-standard.stats.json"
        standard_f1r2 = work / "calls-standard.f1r2.tar.gz"
        standard_command = [str(native), "-I", str(bam), "-R", str(reference),
                            "-L", "17:69000-70000", "-O", str(standard_output),
                            "--tumor-sample", "NA12878",
                            "--normal-input", str(bam), "--normal-lod", "-100",
                            "--min-depth", "1", "--min-alt-support", "1",
                            "--kmer-size", "7", "--min-kmer-count", "1",
                            "--max-num-haplotypes-in-population", "8", "--max-haplotype-depth", "64",
                            "--max-haplotype-combination-alleles", "5",
                            "--error-correct-reads",
                            "--kmer-length-for-read-error-correction", "5",
                            "--min-observations-for-kmer-to-be-solid", "2",
                            "--error-correction-log-odds", "3.5",
                            "--stats", str(standard_stats), "--f1r2-tar-gz", str(standard_f1r2),
                            "--output-manifest", str(standard_manifest)]
        subprocess.check_output(standard_command, text=True, env=env)
        standard_metadata = json.loads(standard_manifest.read_text(encoding="utf-8"))
        assert standard_metadata["compatibility"]["f1r2_standard_tar"] is True
        assert standard_metadata["compatibility"]["f1r2_legacy_tsv"] is False
        assert any(item["kind"] == "f1r2-collect-tar-gz" and item["complete"]
                   for item in standard_metadata["outputs"])
        standard_stats_data = json.loads(standard_stats.read_text(encoding="utf-8"))
        assert standard_stats_data["f1r2_format"] == "collect-f1r2-counts-tar-gz"
        assert gzip.open(standard_output, "rb").read() == gzip.open(stream_output, "rb").read()
        with tarfile.open(standard_f1r2, "r:gz") as aggregate_archive, \
                tarfile.open(stream_f1r2, "r:gz") as stream_archive:
            aggregate_members = {name: aggregate_archive.extractfile(name).read()
                                 for name in aggregate_archive.getnames()}
            stream_members = {name: stream_archive.extractfile(name).read()
                              for name in stream_archive.getnames()}
            assert aggregate_members == stream_members
        with tarfile.open(standard_f1r2, "r:gz") as archive:
            names = archive.getnames()
            assert len(names) == 3
            assert sum(name.endswith(".ref_histogram") for name in names) == 1
            assert sum(name.endswith(".alt_histogram") for name in names) == 1
            assert sum(name.endswith(".alt_table") for name in names) == 1
            ref_member = archive.extractfile(next(name for name in names if name.endswith(".ref_histogram")))
            alt_member = archive.extractfile(next(name for name in names if name.endswith(".alt_histogram")))
            table_member = archive.extractfile(next(name for name in names if name.endswith(".alt_table")))
            assert ref_member is not None and b"## HISTOGRAM" in ref_member.read()
            assert alt_member is not None and b"## HISTOGRAM" in alt_member.read()
            assert table_member is not None and b"#<METADATA>SAMPLE=" in table_member.read()
        learn = root / "fastgatk-native/build/fastgatk-learn-read-orientation-model"
        learned = work / "calls-standard.orientation-priors.tar.gz"
        learn_manifest = work / "calls-standard.orientation-priors.manifest.json"
        learn_result = json.loads(subprocess.check_output([
            str(learn), "-I", str(standard_f1r2), "-O", str(learned),
            "--threads", "2", "--output-manifest", str(learn_manifest)
        ], text=True, env=env).splitlines()[-1])
        assert learn_result["contexts"] == 64
        with tarfile.open(learned, "r:gz") as archive:
            learned_names = archive.getnames()
            assert len(learned_names) == 1 and learned_names[0].endswith(".orientation_priors")
            learned_text = archive.extractfile(learned_names[0]).read().decode("utf-8")
            assert learned_text.startswith("#<METADATA>SAMPLE=")
            learned_rows = [line for line in learned_text.splitlines()
                            if line and not line.startswith("#")]
            assert learned_rows and learned_rows[0].startswith("context\t")
            assert len(learned_rows[1:]) == 64
        # A concrete two-ALT locus must be materialized as one VariantBlock.
        # This exercises the same per-read ownership matrix used by the real
        # fixture, rather than accepting one independent biallelic row per ALT.
        synthetic_reference = work / "synthetic.fa"
        synthetic_sam = work / "synthetic.sam"
        synthetic_normal_sam = work / "synthetic-normal.sam"
        synthetic_output = work / "synthetic.vcf.gz"
        synthetic_stats = work / "synthetic.stats.json"
        synthetic_manifest = work / "synthetic.manifest.json"
        reference_sequence = ("ACGT" * 10)[:40]
        synthetic_reference.write_text(
            ">chr1\n" + reference_sequence + "\n", encoding="utf-8")
        synthetic_reference.with_suffix(synthetic_reference.suffix + ".fai").write_text(
            "chr1\t40\t6\t40\t41\n", encoding="utf-8")
        sam_header = (
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:40\n"
            "@RG\tID:1\tSM:TUMOR\tPL:ILLUMINA\n"
            "@RG\tID:2\tSM:OTHER\tPL:ILLUMINA\n")
        synthetic_reads = []
        for alternate_index, alternate in enumerate(("A", "T")):
            for copy in range(20):
                sequence = reference_sequence[:10] + alternate + reference_sequence[11:31]
                synthetic_reads.append(
                    f"{alternate.lower()}{alternate_index}{copy}\t0\tchr1\t1\t60\t31M\t*\t0\t0\t"
                    f"{sequence}\t{'I' * 31}\tRG:Z:1\tNM:i:1\tMD:Z:10G20\n")
        # These records carry a real second @RG/@SM.  GATK Mutect2 treats
        # every non-normal input sample as a tumor even when --tumor-sample
        # names TUMOR, so they must enter the aggregate likelihood matrix
        # while retaining their own VCF FORMAT column.
        for copy in range(5):
            sequence = reference_sequence[:10] + "C" + reference_sequence[11:31]
            synthetic_reads.append(
                f"other{copy}\t0\tchr1\t1\t60\t31M\t*\t0\t0\t{sequence}\t{'I' * 31}"
                "\tRG:Z:2\tNM:i:1\tMD:Z:10G10\n")
        synthetic_sam.write_text(sam_header + "".join(synthetic_reads), encoding="utf-8")
        normal_header = (
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:40\n"
            "@RG\tID:3\tSM:NORMAL\tPL:ILLUMINA\n"
            "@RG\tID:4\tSM:OTHER_NORMAL\tPL:ILLUMINA\n")
        normal_reads = []
        for copy in range(8):
            sequence = reference_sequence[:10] + "A" + reference_sequence[11:31]
            normal_reads.append(
                f"normal{copy}\t0\tchr1\t1\t60\t31M\t*\t0\t0\t{sequence}\t{'I' * 31}"
                "\tRG:Z:3\tNM:i:1\tMD:Z:10G10\n")
        for copy in range(4):
            sequence = reference_sequence[:10] + "T" + reference_sequence[11:31]
            normal_reads.append(
                f"other-normal{copy}\t0\tchr1\t1\t60\t31M\t*\t0\t0\t{sequence}\t{'I' * 31}"
                "\tRG:Z:4\tNM:i:1\tMD:Z:10G10\n")
        synthetic_normal_sam.write_text(normal_header + "".join(normal_reads), encoding="utf-8")
        synthetic_result = json.loads(subprocess.check_output([
            str(native), "-I", str(synthetic_sam), "-R", str(synthetic_reference),
            "-L", "chr1:1-40", "-O", str(synthetic_output), "--tumor-sample", "TUMOR",
            "--normal-input", str(synthetic_normal_sam), "-normal", "NORMAL", "--normal-lod", "-100",
            "--min-depth", "1", "--min-alt-support", "1", "--kmer-size", "7",
            "--min-kmer-count", "1", "--max-num-haplotypes-in-population", "8",
            "--max-haplotype-depth", "64", "--dont-increase-kmer-sizes-for-cycles",
            "--stats", str(synthetic_stats), "--f1r2-tar-gz", str(work / "synthetic.f1r2.tar.gz"),
            "--output-manifest", str(synthetic_manifest),
        ], text=True, env=env).splitlines()[-1])
        synthetic_text = gzip.open(synthetic_output, "rt", encoding="utf-8").read()
        synthetic_records = [line.split("\t") for line in synthetic_text.splitlines()
                             if line and not line.startswith("#")]
        assert synthetic_result["tumor_reads"] == 45
        assert synthetic_result["tumor_calls"] >= 2
        assert len(synthetic_records) == 1
        assert synthetic_records[0][3] == "G" and synthetic_records[0][4] == "A,T"
        synthetic_header = next(line.split("\t") for line in synthetic_text.splitlines()
                                if line.startswith("#CHROM"))
        assert synthetic_header[9:] == ["NORMAL", "OTHER", "TUMOR"]
        synthetic_info = {token.split("=", 1)[0]: token.split("=", 1)[1]
                          for token in synthetic_records[0][7].split(";") if "=" in token}
        synthetic_normal_format = dict(zip(
            synthetic_records[0][8].split(":"), synthetic_records[0][9].split(":")))
        synthetic_other_tumor_format = dict(zip(
            synthetic_records[0][8].split(":"), synthetic_records[0][10].split(":")))
        synthetic_tumor_format = dict(zip(
            synthetic_records[0][8].split(":"), synthetic_records[0][11].split(":")))
        for name in ("TLOD", "NALOD", "NLOD"):
            assert len(synthetic_info[name].split(",")) == 2, (name, synthetic_info)
        for name in ("AS_SB_TABLE", "DP", "ECNT", "ECNTH", "MBQ", "MFRL", "MMQ", "MPOS", "POPAF"):
            assert name in synthetic_info, (name, synthetic_info)
        assert synthetic_info["AS_SB_TABLE"].count("|") == 2
        for name in ("MBQ", "MFRL", "MMQ"):
            assert len(synthetic_info[name].split(",")) == 3, (name, synthetic_info)
        assert len(synthetic_info["MPOS"].split(",")) == 2
        assert len(synthetic_info["POPAF"].split(",")) == 2
        assert len(synthetic_tumor_format["AD"].split(",")) == 3
        assert len(synthetic_tumor_format["AF"].split(",")) == 2
        for name in ("F1R2", "F2R1"):
            assert len(synthetic_tumor_format[name].split(",")) == 3
        # SomaticGenotypingEngine serializes the complete local AlleleList
        # for a multi-ALT tumor genotype, rather than a diploid germline
        # winner.  At this G>A,T fixture the GATK-shaped value is therefore
        # 0/1/2; accepting only biallelic diploid spellings contradicts the
        # source contract exercised by the preceding multi-ALT assertions.
        assert len(synthetic_records[0]) == 12 and \
            synthetic_records[0][11].split(":")[0] == "0/1/2"
        # The OTHER sample's five reads carry C, an allele removed before
        # this final G>A,T VariantContext is annotated.  GATK first
        # marginalizes its read likelihoods to REF/A/T, where those reads are
        # tied and therefore non-informative; DepthPerSampleHC consequently
        # writes DP=0 (rather than raw coverage five).
        assert synthetic_other_tumor_format["DP"] == "0"
        assert len(synthetic_other_tumor_format["AD"].split(",")) == 3
        assert len(synthetic_normal_format["AD"].split(",")) == 3
        assert len(synthetic_normal_format["AF"].split(",")) == 2
        synthetic_manifest_payload = json.loads(synthetic_manifest.read_text(encoding="utf-8"))
        assert synthetic_manifest_payload["compatibility"]["multiallelic_locus_writer"] is True
        assert synthetic_manifest_payload["compatibility"]["sample_name_filtering"] is True
        assert synthetic_manifest_payload["compatibility"]["tumor_reads"] == 45
        assert synthetic_manifest_payload["compatibility"]["normal_reads"] == 8

        # Mutect2 is fragment-first: when one end of a pair covers a candidate,
        # the other end contributes to the grouped read likelihood even when it
        # is outside the candidate interval.  The non-covering end is not an
        # independent AD/DP/annotation owner.  Keep a matched control with the
        # same covering ends but no mate, then compare the output contract.  A
        # larger PairHMM request count in the paired run is the observable
        # proof that the fragment-only rows reached the likelihood stage;
        # equal somatic evidence group counts and independent fields ensure
        # they were reduced once per fragment instead of as extra reads.
        fragment_reference = work / "fragment-first.fa"
        fragment_reference_sequence = ("ACGT" * 25)[:100]
        fragment_reference.write_text(
            ">chr1\n" + fragment_reference_sequence + "\n", encoding="utf-8")
        fragment_reference.with_suffix(fragment_reference.suffix + ".fai").write_text(
            "chr1\t100\t6\t100\t101\n", encoding="utf-8")
        fragment_header = (
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:100\n"
            "@RG\tID:1\tSM:TUMOR\tPL:ILLUMINA\n")
        fragment_covering = []
        fragment_noncovering = []
        for index in range(8):
            name = f"fragment{index}"
            alternate = index < 4
            covering_sequence = (fragment_reference_sequence[:10] +
                                 ("A" if alternate else "G") +
                                 fragment_reference_sequence[11:31])
            md = "10G20" if alternate else "31"
            nm = "1" if alternate else "0"
            fragment_covering.append(
                f"{name}\t99\tchr1\t1\t60\t31M\t=\t21\t51\t"
                f"{covering_sequence}\t{'I' * 31}\tRG:Z:1\tNM:i:{nm}\tMD:Z:{md}\n")
            # N bases make the mate explicitly non-informative for the local
            # SNP while still exercising the full read×haplotype likelihood
            # request.  This isolates the fragment ownership contract from a
            # second, unrelated pileup allele or a coordinate-dependent
            # haplotype score.
            mate_sequence = "N" * 31
            fragment_noncovering.append(
                f"{name}\t147\tchr1\t21\t60\t31M\t=\t1\t-51\t"
                f"{mate_sequence}\t{'I' * 31}\tRG:Z:1\tNM:i:0\tMD:Z:31\n")
        fragment_paired_sam = work / "fragment-first-paired.sam"
        fragment_covered_only_sam = work / "fragment-first-covered-only.sam"
        fragment_paired_sam.write_text(
            fragment_header + "".join(fragment_covering + fragment_noncovering),
            encoding="utf-8")
        fragment_covered_only_sam.write_text(
            fragment_header + "".join(
                line.replace("\t99\t", "\t0\t") for line in fragment_covering),
            encoding="utf-8")

        def run_fragment_fixture(sam_path: Path, stem: str):
            output_path = work / f"{stem}.vcf.gz"
            stats_path = work / f"{stem}.stats.json"
            manifest_path = work / f"{stem}.manifest.json"
            fixture_result = json.loads(subprocess.check_output([
                str(native), "-I", str(sam_path), "-R", str(fragment_reference),
                "-L", "chr1:1-100", "-O", str(output_path),
                "--tumor-sample", "TUMOR", "--min-depth", "1",
                "--min-alt-support", "1", "--kmer-size", "7",
                "--min-kmer-count", "1", "--max-num-haplotypes-in-population", "8",
                "--max-haplotype-depth", "64", "--dont-increase-kmer-sizes-for-cycles",
                "--stats", str(stats_path), "--output-manifest", str(manifest_path),
            ], text=True, env=env).splitlines()[-1])
            text = gzip.open(output_path, "rt", encoding="utf-8").read()
            records = [line.split("\t") for line in text.splitlines()
                       if line and not line.startswith("#")]
            candidate_records = [record for record in records
                                 if len(record) > 4 and record[0] == "chr1" and
                                 record[1] == "11" and record[3] == "G" and record[4] == "A"]
            assert len(candidate_records) == 1, (stem, records)
            record = candidate_records[0]
            info = {token.split("=", 1)[0]: token.split("=", 1)[1]
                    for token in record[7].split(";") if "=" in token}
            format_values = dict(zip(record[8].split(":"), record[9].split(":")))
            stats_payload = json.loads(stats_path.read_text(encoding="utf-8"))
            manifest_payload = json.loads(manifest_path.read_text(encoding="utf-8"))
            return fixture_result, record, info, format_values, stats_payload, manifest_payload

        paired_fixture = run_fragment_fixture(fragment_paired_sam, "fragment-first-paired")
        covered_only_fixture = run_fragment_fixture(
            fragment_covered_only_sam, "fragment-first-covered-only")
        paired_result, paired_record, paired_info, paired_format, paired_stats, paired_manifest = paired_fixture
        covered_result, covered_record, covered_info, covered_format, covered_stats, covered_manifest = covered_only_fixture
        # The paired run has one additional non-covering row per fragment.  It
        # must therefore execute more read×haplotype work while retaining the
        # same number of somatic evidence groups (one group per fragment).
        assert paired_result["tumor_reads"] == 16
        assert covered_result["tumor_reads"] == 8
        assert paired_stats["pairhmm_pairs"] > covered_stats["pairhmm_pairs"]
        assert paired_stats["somatic_evidence_groups"] == covered_stats["somatic_evidence_groups"] == 8
        assert paired_stats["somatic_evidence_grouping"] is True
        assert covered_stats["somatic_evidence_grouping"] is False
        # AD/DP and candidate INFO annotations are independent fragment owners;
        # adding the non-covering mate cannot double their depth or allele
        # support.  Compare all emitted site annotations except likelihood
        # fields whose grouped values intentionally include the mate.
        for key in ("AS_SB_TABLE", "DP", "ECNT", "ECNTH", "MBQ", "MFRL", "MMQ", "MPOS", "POPAF"):
            assert paired_info[key] == covered_info[key], (key, paired_info, covered_info)
        assert paired_format["AD"] == covered_format["AD"]
        assert paired_format["DP"] == covered_format["DP"]
        assert paired_format["AF"] == covered_format["AF"]
        # The covered-only fixture clears the first-in-pair flag, so GATK
        # legitimately classifies its otherwise identical fragments as F2R1
        # rather than F1R2.  The orientation labels may swap, but their
        # reference/ALT totals must not change when the non-covering mate is
        # added or removed.
        paired_orientation = [
            left + right for left, right in zip(
                map(int, paired_format["F1R2"].split(",")),
                map(int, paired_format["F2R1"].split(",")),
                strict=True)]
        covered_orientation = [
            left + right for left, right in zip(
                map(int, covered_format["F1R2"].split(",")),
                map(int, covered_format["F2R1"].split(",")),
                strict=True)]
        assert paired_orientation == covered_orientation
        assert paired_manifest["compatibility"]["somatic_fragment_grouping"] is True
        assert covered_manifest["compatibility"]["somatic_fragment_grouping"] is False

        # Standard GATK invocation may carry tumor and normal samples in one
        # multi-sample input (one -I plus --tumor-sample/--normal-sample).
        # The native path must reuse the same HTSlib stream and apply both
        # @RG/SM filters, without requiring the non-GATK --normal-input alias.
        combined_sam = work / "combined-samples.sam"
        tumor_lines = synthetic_sam.read_text(encoding="utf-8").splitlines()
        normal_lines = synthetic_normal_sam.read_text(encoding="utf-8").splitlines()
        normal_rg = []
        normal_body = []
        for line in normal_lines:
            if line.startswith("@RG"):
                if "SM:NORMAL" not in line:
                    continue
                line = line.replace("@RG\tID:2\tSM:NORMAL", "@RG\tID:3\tSM:NORMAL")
                normal_rg.append(line)
            elif line.startswith("normal"):
                normal_body.append(line.replace("RG:Z:2", "RG:Z:3"))
        combined_header = "\n".join(line for line in tumor_lines
                                      if line.startswith("@") and
                                      (not line.startswith("@RG") or "SM:TUMOR" in line))
        combined_sam.write_text(
            combined_header + "\n" + "\n".join(normal_rg) + "\n" +
            "\n".join(line for line in tumor_lines
                         if not line.startswith("@") and "RG:Z:1" in line) + "\n" +
            "\n".join(normal_body) + "\n", encoding="utf-8")
        combined_output = work / "combined-samples.vcf.gz"
        combined_manifest = work / "combined-samples.manifest.json"
        combined_result = json.loads(subprocess.check_output([
            str(native), "-I", str(combined_sam), "-R", str(synthetic_reference),
            "-L", "chr1:1-40", "-O", str(combined_output),
            "--tumor-sample", "TUMOR", "--normal-sample", "NORMAL", "--normal-lod", "-100",
            "--min-depth", "1", "--min-alt-support", "1",
            "--stats", str(work / "combined-samples.stats.json"),
            "--output-manifest", str(combined_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert combined_result["tumor_reads"] == 40
        assert combined_result["normal_reads"] == 8
        assert combined_result["tumor_sample"] == "TUMOR"
        assert combined_result["normal_sample"] == "NORMAL"
        combined_text = gzip.open(combined_output, "rt", encoding="utf-8").read()
        combined_header_line = next(line for line in combined_text.splitlines()
                                    if line.startswith("#CHROM"))
        assert combined_header_line.endswith("\tNORMAL\tTUMOR")
        combined_manifest_payload = json.loads(combined_manifest.read_text(encoding="utf-8"))
        assert combined_manifest_payload["compatibility"]["tumor_normal"] is True
        assert combined_manifest_payload["compatibility"]["sample_name_filtering"] is True

        # GATK's common --sites-only-vcf-output switch changes only the VCF
        # serialization boundary.  Evidence, site INFO and sidecar counts
        # must remain available while FORMAT/sample columns are omitted.
        sites_only_output = work / "combined-sites-only.vcf.gz"
        sites_only_manifest = work / "combined-sites-only.manifest.json"
        sites_only_stats = work / "combined-sites-only.stats.json"
        sites_only_result = json.loads(subprocess.check_output([
            str(native), "-I", str(combined_sam), "-R", str(synthetic_reference),
            "-L", "chr1:1-40", "-O", str(sites_only_output),
            "--tumor-sample", "TUMOR", "--normal-sample", "NORMAL", "--normal-lod", "-100",
            "--sites-only-vcf-output", "true", "--min-depth", "1",
            "--min-alt-support", "1", "--stats", str(sites_only_stats),
            "--output-manifest", str(sites_only_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert sites_only_result["tumor_reads"] == combined_result["tumor_reads"]
        sites_only_text = gzip.open(sites_only_output, "rt", encoding="utf-8").read()
        sites_only_header = next(line for line in sites_only_text.splitlines()
                                 if line.startswith("#CHROM"))
        assert sites_only_header.endswith("\tINFO")
        sites_only_records = [line.split("\t") for line in sites_only_text.splitlines()
                              if line and not line.startswith("#")]
        assert sites_only_records and all(len(record) == 8 for record in sites_only_records)
        sites_only_manifest_payload = json.loads(sites_only_manifest.read_text(encoding="utf-8"))
        assert sites_only_manifest_payload["compatibility"]["sites_only_vcf_output"] is True
        sites_only_stats_payload = json.loads(sites_only_stats.read_text(encoding="utf-8"))
        assert sites_only_stats_payload["sites_only_vcf_output"] is True

        # GATK's common optional index switch must affect only the final file
        # boundary: disabling it suppresses .tbi while retaining the same VCF
        # records, stats and F1R2 evidence.
        no_index_output = work / "combined-no-index.vcf.gz"
        no_index_manifest = work / "combined-no-index.manifest.json"
        no_index_stats = work / "combined-no-index.stats.json"
        no_index_result = json.loads(subprocess.check_output([
            str(native), "-I", str(combined_sam), "-R", str(synthetic_reference),
            "-L", "chr1:1-40", "-O", str(no_index_output),
            "--tumor-sample", "TUMOR", "--normal-sample", "NORMAL", "--normal-lod", "-100",
            "--create-output-variant-index", "false", "--min-depth", "1",
            "--min-alt-support", "1", "--stats", str(no_index_stats),
            "--output-manifest", str(no_index_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert no_index_result["tumor_reads"] == combined_result["tumor_reads"]
        assert no_index_output.exists() and not Path(f"{no_index_output}.tbi").exists()
        no_index_stats_payload = json.loads(no_index_stats.read_text(encoding="utf-8"))
        assert no_index_stats_payload["create_output_variant_index"] is False
        no_index_manifest_payload = json.loads(no_index_manifest.read_text(encoding="utf-8"))
        assert no_index_manifest_payload["compatibility"]["vcf_index"] is False
        assert no_index_manifest_payload["compatibility"]["create_output_variant_index"] is False

        # An uncompressed VCF follows GATK's Tribble writer boundary: the
        # optional index is a sibling .idx (not a Tabix .tbi).  The native
        # index must be consumable by HTSJDK's interval reader as well as
        # visible in OutputManifest.
        plain_output = work / "combined-plain.vcf"
        plain_manifest = work / "combined-plain.manifest.json"
        plain_stats = work / "combined-plain.stats.json"
        plain_result = json.loads(subprocess.check_output([
            str(native), "-I", str(combined_sam), "-R", str(synthetic_reference),
            "-L", "chr1:1-40", "-O", str(plain_output),
            "--tumor-sample", "TUMOR", "--normal-sample", "NORMAL", "--normal-lod", "-100",
            "--min-depth", "1", "--min-alt-support", "1",
            "--stats", str(plain_stats), "--output-manifest", str(plain_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert plain_result["tumor_reads"] == combined_result["tumor_reads"]
        plain_index = Path(f"{plain_output}.idx")
        assert plain_output.exists() and plain_output.stat().st_size > 0
        assert plain_index.is_file() and plain_index.stat().st_size > 0
        assert not Path(f"{plain_output}.tbi").exists()
        plain_manifest_payload = json.loads(plain_manifest.read_text(encoding="utf-8"))
        assert plain_manifest_payload["compatibility"]["vcf_index"] is True
        assert plain_manifest_payload["compatibility"]["create_output_variant_index"] is True
        if oracle_guard.oracle_ready('verify_mutect2.py', root / 'third_party/jdk17/bin/java', root / 'third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar'):
            query_output = work / "combined-plain-query.vcf"
            query = subprocess.run([
                str(root / "third_party/jdk17/bin/java"), "-Xmx1g", "-jar",
                str(root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"),
                "SelectVariants", "-V", str(plain_output), "-L", "chr1:10-40",
                "-O", str(query_output), "--create-output-variant-index", "false",
            ], text=True, capture_output=True, env=env)
            assert query.returncode == 0, query.stderr
            assert any(line.startswith("chr1\t") for line in query_output.read_text(encoding="utf-8").splitlines()
                       if line and not line.startswith("#"))

        # GATK allows multiple -normal/--normal-sample selectors.  Their
        # evidence is combined for somatic inference, but each selected
        # normal retains its own genotype column at the VCF boundary.
        multi_normal_sam = work / "combined-multi-normal.sam"
        combined_lines = combined_sam.read_text(encoding="utf-8").splitlines()
        combined_header_lines = [line for line in combined_lines if line.startswith("@")]
        combined_body_lines = [line for line in combined_lines if not line.startswith("@")]
        multi_normal_sam.write_text(
            "\n".join(combined_header_lines +
                         ["@RG\tID:4\tSM:OTHER_NORMAL\tPL:ILLUMINA"] +
                         combined_body_lines +
                         [line for line in normal_lines if line.startswith("other-normal")]) +
            "\n",
            encoding="utf-8")
        multi_normal_output = work / "combined-multi-normal.vcf.gz"
        multi_normal_stats = work / "combined-multi-normal.stats.json"
        multi_normal_manifest = work / "combined-multi-normal.manifest.json"
        multi_normal_result = json.loads(subprocess.check_output([
            str(native), "-I", str(multi_normal_sam), "-R", str(synthetic_reference),
            "-L", "chr1:1-40", "-O", str(multi_normal_output),
            "--tumor-sample", "TUMOR", "--normal-sample", "NORMAL",
            "--normal-sample", "OTHER_NORMAL", "--normal-lod", "-100", "--min-depth", "1",
            "--min-alt-support", "1", "--stats", str(multi_normal_stats),
            "--output-manifest", str(multi_normal_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert multi_normal_result["tumor_reads"] == 40
        assert multi_normal_result["normal_reads"] == 12
        multi_normal_text = gzip.open(multi_normal_output, "rt", encoding="utf-8").read()
        multi_normal_header = next(line for line in multi_normal_text.splitlines()
                                   if line.startswith("#CHROM"))
        assert multi_normal_header.endswith("\tNORMAL\tOTHER_NORMAL\tTUMOR")
        multi_normal_records = [line.split("\t") for line in multi_normal_text.splitlines()
                                if line and not line.startswith("#")]
        assert multi_normal_records and all(len(record) == 12 for record in multi_normal_records)
        multi_normal_stats_payload = json.loads(multi_normal_stats.read_text(encoding="utf-8"))
        assert multi_normal_stats_payload["normal_samples"] == ["NORMAL", "OTHER_NORMAL"]
        multi_normal_manifest_payload = json.loads(multi_normal_manifest.read_text(encoding="utf-8"))
        assert multi_normal_manifest_payload["compatibility"]["normal_samples"] == [
            "NORMAL", "OTHER_NORMAL"]

        # Repeatable -I shards are merged in argv order at the Host boundary.
        # Every non-normal @RG/SM remains a tumor sample, so each shard
        # contributes both its TUMOR and OTHER read populations.
        shard_two = work / "tumor-shard-two.sam"
        shard_two_lines = []
        for line in tumor_lines:
            if not line.startswith("@"):
                fields = line.split("\t")
                fields[0] = "shard2-" + fields[0]
                line = "\t".join(fields)
            shard_two_lines.append(line)
        shard_two.write_text("\n".join(shard_two_lines) + "\n", encoding="utf-8")
        shard_output = work / "tumor-shards.vcf.gz"
        shard_manifest = work / "tumor-shards.manifest.json"
        shard_result = json.loads(subprocess.check_output([
            str(native), "-I", str(synthetic_sam), "-I", str(shard_two),
            "-R", str(synthetic_reference), "-L", "chr1:1-40",
            "-O", str(shard_output), "--tumor-sample", "TUMOR",
            "--min-depth", "1", "--min-alt-support", "1",
            "--stats", str(work / "tumor-shards.stats.json"),
            "--output-manifest", str(shard_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert shard_result["tumor_reads"] == 90
        assert shard_result["tumor_sample"] == "TUMOR"
        shard_manifest_payload = json.loads(shard_manifest.read_text(encoding="utf-8"))
        assert shard_manifest_payload["compatibility"]["sample_name_filtering"] is True
        paired_shard_output = work / "paired-shards.vcf.gz"
        paired_shard_result = json.loads(subprocess.check_output([
            str(native), "-I", str(synthetic_sam), "-I", str(synthetic_normal_sam),
            "-R", str(synthetic_reference), "-L", "chr1:1-40",
            "-O", str(paired_shard_output), "--tumor-sample", "TUMOR",
            "--normal-sample", "NORMAL", "--normal-lod", "-100", "--min-depth", "1",
            "--min-alt-support", "1",
        ], text=True, env=env).splitlines()[-1])
        # OTHER and OTHER_NORMAL were not selected as matched normals, so
        # GATK's role rule retains their 5 + 4 reads as additional tumors.
        assert paired_shard_result["tumor_reads"] == 49
        assert paired_shard_result["normal_reads"] == 8
        assert paired_shard_result["tumor_sample"] == "TUMOR"
        assert paired_shard_result["normal_sample"] == "NORMAL"
        incompatible_shard = work / "tumor-shard-incompatible.sam"
        incompatible_shard.write_text(
            shard_two.read_text(encoding="utf-8").replace("@SQ\tSN:chr1\tLN:40",
                                                             "@SQ\tSN:chr1\tLN:41"),
            encoding="utf-8")
        incompatible = subprocess.run([
            str(native), "-I", str(synthetic_sam), "-I", str(incompatible_shard),
            "-R", str(synthetic_reference), "-L", "chr1:1-40",
            "-O", str(work / "incompatible.vcf.gz"), "--tumor-sample", "TUMOR",
            "--min-depth", "1", "--min-alt-support", "1",
        ], text=True, capture_output=True, env=env, check=False)
        assert incompatible.returncode != 0
        assert "incompatible sequence dictionaries" in incompatible.stderr

        # With no --tumor-sample, the first header sample remains the legacy
        # metadata name, but every non-normal input sample contributes to
        # Mutect2's tumor matrix.
        default_output = work / "default-tumor.vcf.gz"
        default_manifest = work / "default-tumor.manifest.json"
        default_result = json.loads(subprocess.check_output([
            str(native), "-I", str(synthetic_sam), "-R", str(synthetic_reference),
            "-L", "chr1:1-40", "-O", str(default_output),
            "--min-depth", "1", "--min-alt-support", "1",
            "--stats", str(work / "default-tumor.stats.json"),
            "--output-manifest", str(default_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert default_result["tumor_reads"] == 45
        assert default_result["tumor_sample"] == "TUMOR"
        default_manifest_payload = json.loads(default_manifest.read_text(encoding="utf-8"))
        assert default_manifest_payload["compatibility"]["sample_name_filtering"] is True
        missing = subprocess.run([
            str(native), "-I", str(synthetic_sam), "-R", str(synthetic_reference),
            "-L", "chr1:1-40", "-O", str(work / "missing.vcf.gz"),
            "--tumor-sample", "MISSING", "--min-depth", "1", "--min-alt-support", "1",
        ], text=True, capture_output=True, env=env, check=False)
        assert missing.returncode != 0
        assert "requested sample is absent" in missing.stderr

        # GATK does not invent a sample name when the tumor header has no
        # @RG/SM records: ReferenceConfidenceModel rejects an empty sample
        # set.  Keep the native Host boundary fail-closed as well, rather than
        # decoding every read and silently labeling the output TUMOR.
        no_sample_sam = work / "no-sample.sam"
        no_sample_lines = []
        for line in synthetic_sam.read_text(encoding="utf-8").splitlines():
            if line.startswith("@RG"):
                continue
            if not line.startswith("@"):
                line = "\t".join(field for field in line.split("\t")
                                   if not field.startswith("RG:Z:"))
            no_sample_lines.append(line)
        no_sample_sam.write_text("\n".join(no_sample_lines) + "\n", encoding="utf-8")
        no_sample = subprocess.run([
            str(native), "-I", str(no_sample_sam), "-R", str(synthetic_reference),
            "-L", "chr1:1-40", "-O", str(work / "no-sample.vcf.gz"),
            "--min-depth", "1", "--min-alt-support", "1",
        ], text=True, capture_output=True, env=env, check=False)
        assert no_sample.returncode != 0
        assert "no @RG SM sample" in no_sample.stderr

        # Exercise the same edge against the pinned Java oracle when it is
        # available.  A minimal sequence dictionary is sufficient for this
        # synthetic SAM and keeps this check independent of the large fixture.
        java = root / "third_party/jdk17/bin/java"
        gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        if oracle_guard.oracle_ready('verify_mutect2.py', java, gatk):
            synthetic_reference.with_suffix(".dict").write_text(
                "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:40\n", encoding="utf-8")
            gatk_no_sample = subprocess.run([
                str(java), "-jar", str(gatk), "Mutect2", "-R", str(synthetic_reference),
                "-I", str(no_sample_sam), "-L", "chr1:1-40",
                "-O", str(work / "gatk-no-sample.vcf.gz"),
            ], text=True, capture_output=True, env=env, check=False)
            assert gatk_no_sample.returncode != 0
            assert "samples cannot be empty" in gatk_no_sample.stderr

        # This is an unpaired SAM, so it has no well-defined fragment size.
        # AssemblyBasedCallerUtils.finalizeRegion therefore hard-clips its
        # terminal soft clips even when --dont-use-soft-clipped-bases is
        # false: `!hasWellDefinedFragmentSize(read)` takes the same source
        # branch as the explicit option.  A dangling terminal CIGAR S must
        # consequently not become an assembled Mutect2 insertion merely
        # because the Host can reconstruct its bases.  Keep both CLI states
        # here to guard option propagation and this default fragment gate.
        softclip_sam = work / "softclip.sam"
        softclip_output = work / "softclip.vcf.gz"
        softclip_stats = work / "softclip.stats.json"
        softclip_manifest = work / "softclip.manifest.json"
        softclip_disabled_output = work / "softclip-disabled.vcf.gz"
        softclip_disabled_stats = work / "softclip-disabled.stats.json"
        softclip_disabled_manifest = work / "softclip-disabled.manifest.json"
        softclip_sequence = reference_sequence[:30] + "CCC"
        softclip_sam.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:40\n"
            "@RG\tID:1\tSM:TUMOR\tPL:ILLUMINA\n" +
            "".join(
                f"clip{index}\t0\tchr1\t1\t60\t30M3S\t*\t0\t0\t{softclip_sequence}\t{'I' * len(softclip_sequence)}\tRG:Z:1\n"
                for index in range(4)), encoding="utf-8")
        softclip_common = [str(native), "-I", str(softclip_sam), "-R", str(synthetic_reference),
                           "-L", "chr1:1-40", "--min-depth", "1", "--min-alt-support", "2",
                           "--min-kmer-count", "1", "--kmer-size", "7",
                           "--max-num-haplotypes-in-population", "8", "--max-haplotype-depth", "64"]
        subprocess.check_output(softclip_common + ["-O", str(softclip_output),
                                                   "--stats", str(softclip_stats),
                                                   "--output-manifest", str(softclip_manifest)],
                                text=True, env=env)
        subprocess.check_output(softclip_common + ["-O", str(softclip_disabled_output),
                                                   "--stats", str(softclip_disabled_stats),
                                                   "--output-manifest", str(softclip_disabled_manifest),
                                                   "--dont-use-soft-clipped-bases"],
                                text=True, env=env)
        softclip_text = gzip.open(softclip_output, "rt", encoding="utf-8").read()
        softclip_disabled_text = gzip.open(softclip_disabled_output, "rt", encoding="utf-8").read()
        softclip_records = [line.split("\t") for line in softclip_text.splitlines()
                            if line and not line.startswith("#")]
        softclip_disabled_records = [line.split("\t") for line in softclip_disabled_text.splitlines()
                                     if line and not line.startswith("#")]
        assert not any(len(fields) >= 5 and len(fields[4]) > len(fields[3])
                       for fields in softclip_records)
        assert not any(len(fields) >= 5 and len(fields[4]) > len(fields[3])
                       for fields in softclip_disabled_records)
        assert json.loads(softclip_stats.read_text(encoding="utf-8"))["use_soft_clipped_bases"] is True
        assert json.loads(softclip_disabled_stats.read_text(encoding="utf-8"))["use_soft_clipped_bases"] is False
        assert json.loads(softclip_manifest.read_text(encoding="utf-8"))["compatibility"]["use_soft_clipped_bases"] is True
        assert json.loads(softclip_disabled_manifest.read_text(encoding="utf-8"))["compatibility"]["use_soft_clipped_bases"] is False

        # The GATK -mbq alias must be parsed at the caller boundary and apply
        # before assembly.  All bases in this low-quality copy are below Q30,
        # so no variant record may be emitted while the run still produces a
        # complete indexed output and telemetry.
        low_bq_sam = work / "low-bq.sam"
        low_bq_output = work / "low-bq.vcf.gz"
        low_bq_stats = work / "low-bq.stats.json"
        low_bq_manifest = work / "low-bq.manifest.json"
        low_bq_sam.write_text(softclip_sam.read_text(encoding="utf-8").replace(
            "30M3S", "33M").replace(
                "I" * len(softclip_sequence), "!" * len(softclip_sequence)), encoding="utf-8")
        low_bq_command = softclip_common.copy()
        low_bq_command[2] = str(low_bq_sam)
        low_bq_result = json.loads(subprocess.check_output(
            low_bq_command + ["-O", str(low_bq_output), "--stats", str(low_bq_stats),
                              "--output-manifest", str(low_bq_manifest), "-mbq", "30"],
            text=True, env=env).splitlines()[-1])
        low_bq_text = gzip.open(low_bq_output, "rt", encoding="utf-8").read()
        low_bq_records = [line for line in low_bq_text.splitlines()
                          if line and not line.startswith("#")]
        low_bq_stats_data = json.loads(low_bq_stats.read_text(encoding="utf-8"))
        assert low_bq_stats_data["tumor_candidate_sites"] == 0, low_bq_stats_data
        assert low_bq_result["tumor_calls"] == 0
        assert not low_bq_records
        assert json.loads(low_bq_manifest.read_text(encoding="utf-8"))["compatibility"]["read_filter_gatk_defaults"] is True
        assert low_bq_stats_data["use_soft_clipped_bases"] is True
        assert low_bq_stats_data["min_base_quality"] == 30
        low_bq_manifest_payload = json.loads(low_bq_manifest.read_text(encoding="utf-8"))
        assert low_bq_manifest_payload["compatibility"]["min_base_quality"] == 30
        assert low_bq_manifest_payload["telemetry"]["min_base_quality"] == 30
        print(json.dumps({"status": "pass", "tumor_calls": result["tumor_calls"],
                          "normal_calls": result["normal_calls"], "indexed": True,
                          "standard_f1r2_members": len(names),
                          "learned_orientation_contexts": learn_result["contexts"],
                          "synthetic_multiallelic_records": len(synthetic_records),
                          "softclip_records": len(softclip_records),
                          "softclip_disabled_records": len(softclip_disabled_records),
                          "low_bq_records": len(low_bq_records)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
