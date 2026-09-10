#!/usr/bin/env python3
"""Verify the native GVCF -> VCF materialization contract."""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    hc = Path(os.environ.get("FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    genotype = Path(os.environ.get("FASTGATK_GENOTYPE_BINARY", root / "fastgatk-native/build/fastgatk-genotype-gvcf"))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    java = root / "third_party/jdk17/bin/java"
    gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    with tempfile.TemporaryDirectory(prefix="fastgatk-genotype-gvcf-") as temp:
        work = Path(temp)
        gvcf = work / "input.g.vcf.gz"
        output = work / "output.vcf.gz"
        manifest = work / "output.manifest.json"
        subprocess.run([
            str(hc), "-I", str(bam), "-R", str(reference), "-L", "17:69000-69100",
            "-O", str(gvcf), "-ERC", "GVCF", "--min-depth", "1", "--min-alt-support", "1",
        ], check=True, stdout=subprocess.DEVNULL)
        summary = json.loads(subprocess.check_output([
            str(genotype), "-V", str(gvcf), "-O", str(output),
            "--output-manifest", str(manifest),
        ], text=True).splitlines()[-1])
        assert summary["status"] == "contract-compatible"
        assert summary["output_records"] > 0
        assert output.exists() and Path(f"{output}.tbi").exists()
        text = gzip.open(output, "rt", encoding="utf-8").read()
        records = [line for line in text.splitlines() if line and not line.startswith("#")]
        assert records and all("<NON_REF>" not in line for line in records)
        assert all("AD" in line.split("\t")[8] for line in records)
        metadata = json.loads(manifest.read_text())
        assert metadata["compatibility"]["reference_block_materialization"] is True
        assert metadata["compatibility"]["vcf_index"] is True
        assert metadata["compatibility"]["format_compaction"] is True
        assert metadata["compatibility"]["multi_allelic_materialization"] is True
        assert metadata["compatibility"]["allele_union"] is True
        assert metadata["compatibility"]["kokkos_genotype_pl"] is True
        assert metadata["compatibility"]["arbitrary_ploidy_gt_gq"] is True
        assert metadata["compatibility"]["kokkos_allele_counts"] is True
        assert metadata["telemetry"]["genotype_kernel_calls"] > 0
        assert metadata["telemetry"]["genotype_kernel_seconds"] > 0.0
        assert metadata["telemetry"]["genotype_kernel_execution_space"]
        assert metadata["telemetry"]["kokkos_default_concurrency"] >= 1
        assert metadata["telemetry"]["pl_remap_kernel_calls"] > 0
        assert metadata["telemetry"]["pl_remap_kernel_seconds"] > 0.0
        assert metadata["telemetry"]["pl_remap_kernel_execution_space"]
        assert metadata["telemetry"]["allele_field_remap_kernel_calls"] > 0
        assert metadata["telemetry"]["allele_field_remap_kernel_execution_space"]
        assert metadata["telemetry"]["allele_count_kernel_calls"] > 0
        assert metadata["telemetry"]["allele_count_kernel_seconds"] > 0.0
        assert metadata["telemetry"]["allele_count_kernel_execution_space"]
        assert metadata["compatibility"]["genotype_prior_calculator_assumingHW"] is True
        assert metadata["compatibility"]["posterior_qual_opt_in"] is True
        assert metadata["compatibility"]["cohort_af_calculator"] is True
        assert metadata["compatibility"]["allele_frequency_calculator_em"] is True
        assert metadata["compatibility"]["spanning_deletion_nonvariant_set"] is True
        assert metadata["compatibility"]["excess_het_exact_diploid"] is True
        assert metadata["compatibility"]["cross_sample_reference_confidence"] is True
        assert metadata["telemetry"]["cohort_af_kernel_calls"] > 0
        assert metadata["telemetry"]["excess_het_calls"] > 0
        assert metadata["telemetry"]["cohort_af_samples"] > 0
        assert metadata["telemetry"]["cohort_af_converged"] > 0
        assert metadata["telemetry"]["cross_sample_reference_kernel_calls"] > 0
        assert metadata["telemetry"]["cross_sample_reference_kernel_seconds"] > 0.0
        assert metadata["telemetry"]["cross_sample_reference_execution_space"]
        assert metadata["telemetry"]["posterior_kernel_calls"] == 0
        assert metadata["compatibility"]["three_stage_pipeline"] is True
        pipeline = metadata["telemetry"]
        assert pipeline["pipeline_lifecycle"] == (
            "Host decode->bounded queue->Kokkos compute->encode->sink"
        )
        assert pipeline["pipeline_decoded_items"] == pipeline["pipeline_computed_items"]
        assert pipeline["pipeline_computed_items"] == pipeline["pipeline_encoded_items"] == summary["output_records"]
        assert pipeline["pipeline_peak_decoded_bytes"] > 0
        assert pipeline["pipeline_peak_computed_bytes"] > 0
        assert pipeline["pipeline_peak_encoded_bytes"] > 0
        assert all("RCQ=" in record[7] and "RCP=" in record[7] for record in (
            line.split("\t") for line in records))

        # Plain VCF output uses the same shared Tribble LinearIndex v3 writer
        # as the other native tools.  GATK must be able to open that sidecar
        # and answer an interval query without a conversion step.
        plain_output = work / "output-plain.vcf"
        plain_manifest = work / "output-plain.manifest.json"
        plain_summary = json.loads(subprocess.check_output([
            str(genotype), "-V", str(gvcf), "-O", str(plain_output),
            "--output-manifest", str(plain_manifest),
        ], text=True).splitlines()[-1])
        assert plain_summary["output_records"] == summary["output_records"]
        assert plain_output.exists() and Path(f"{plain_output}.idx").exists()
        assert not Path(f"{plain_output}.tbi").exists()
        plain_query = work / "output-plain-query.vcf"
        plain_query_result = subprocess.run([
            str(java), "-jar", str(gatk_jar), "SelectVariants", "-V", str(plain_output),
            "-L", "17:69000-69100", "-O", str(plain_query),
        ], text=True, capture_output=True, check=False)
        assert plain_query_result.returncode == 0, plain_query_result.stderr
        assert sum(1 for line in plain_query.read_text(encoding="utf-8").splitlines()
                   if line and not line.startswith("#")) == summary["output_records"]
        plain_metadata = json.loads(plain_manifest.read_text(encoding="utf-8"))
        assert plain_metadata["compatibility"]["vcf_index"] is True

        # Optional-boolean index control must reject malformed values and
        # suppress only the sidecar, leaving the materialized records intact.
        no_index_output = work / "no-index-output.vcf.gz"
        no_index_manifest = work / "no-index-output.manifest.json"
        no_index_summary = json.loads(subprocess.check_output([
            str(genotype), "-V", str(gvcf), "-O", str(no_index_output),
            "--create-output-variant-index=false",
            "--output-manifest", str(no_index_manifest),
        ], text=True).splitlines()[-1])
        assert no_index_summary["output_records"] == summary["output_records"]
        assert no_index_output.exists() and not Path(f"{no_index_output}.tbi").exists()
        no_index_metadata = json.loads(no_index_manifest.read_text(encoding="utf-8"))
        assert no_index_metadata["compatibility"]["vcf_index"] is False
        assert gzip.open(no_index_output, "rt", encoding="utf-8").read() == text
        invalid_index = subprocess.run([
            str(genotype), "-V", str(gvcf), "-O", str(work / "invalid-index.vcf.gz"),
            "--create-output-variant-index=maybe",
        ], text=True, capture_output=True)
        assert invalid_index.returncode != 0 and "invalid boolean" in invalid_index.stderr

        # Indexed interval traversal must avoid the sequential whole-file
        # scan when the BGZF gVCF carries a TBI/CSI sidecar.  The output
        # records still go through the identical Host staging and Kokkos
        # genotype path; only the decode/traversal boundary changes.
        indexed_output = work / "indexed-output.vcf.gz"
        indexed_manifest = work / "indexed-output.manifest.json"
        indexed_summary = json.loads(subprocess.check_output([
            str(genotype), "-V", str(gvcf), "-L", "17:69000-69100",
            "-O", str(indexed_output), "--output-manifest", str(indexed_manifest),
        ], text=True).splitlines()[-1])
        assert indexed_summary["output_records"] > 0
        indexed_metadata = json.loads(indexed_manifest.read_text(encoding="utf-8"))
        assert indexed_metadata["compatibility"]["indexed_interval_traversal"] is True
        assert indexed_metadata["telemetry"]["indexed_inputs"] == 1
        assert indexed_metadata["telemetry"]["indexed_interval_queries"] >= 1
        assert indexed_metadata["telemetry"]["interval_skipped"] == 0
        assert indexed_output.exists() and Path(f"{indexed_output}.tbi").exists()

        # Repeatable -L selectors must honor GATK's set rule.  The default is
        # UNION; INTERSECTION keeps only the overlap of selector groups and
        # records the chosen rule in both the summary and OutputManifest.
        intersection_output = work / "intersection-output.vcf.gz"
        intersection_manifest = work / "intersection-output.manifest.json"
        intersection_summary = json.loads(subprocess.check_output([
            str(genotype), "-V", str(gvcf),
            "-L", "17:69000-69080", "-L", "17:69040-69100",
            "--interval-set-rule", "INTERSECTION",
            "-O", str(intersection_output),
            "--output-manifest", str(intersection_manifest),
        ], text=True).splitlines()[-1])
        assert intersection_summary["interval_set_rule"] == "INTERSECTION"
        intersection_metadata = json.loads(intersection_manifest.read_text(encoding="utf-8"))
        assert intersection_metadata["compatibility"]["interval_set_rule"] == "INTERSECTION"
        assert intersection_metadata["telemetry"]["interval_set_rule"] == "INTERSECTION"
        intersection_records = [line.split("\t") for line in gzip.open(
            intersection_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert intersection_records
        assert all(69040 <= int(record[1]) < 69080 for record in intersection_records)

        # Explicit stream-by-locus mode must produce the same VCF bytes as the
        # aggregate joint-materialization path while keeping only one locus
        # per source plus bounded pipeline queues in memory.
        stream_output = work / "stream-output.vcf.gz"
        stream_manifest = work / "stream-output.manifest.json"
        stream_summary = json.loads(subprocess.check_output([
            str(genotype), "-V", str(gvcf), "-O", str(stream_output),
            "--stream-by-locus", "--output-manifest", str(stream_manifest),
        ], text=True).splitlines()[-1])
        assert stream_summary["stream_by_locus"] is True
        assert stream_summary["streamed_loci"] == stream_summary["output_records"]
        stream_metadata = json.loads(stream_manifest.read_text(encoding="utf-8"))
        assert stream_metadata["compatibility"]["stream_by_locus"] is True
        assert stream_metadata["compatibility"]["bounded_k_way_merge"] is True
        assert stream_metadata["compatibility"]["three_stage_pipeline"] is True
        assert stream_metadata["telemetry"]["streamed_loci"] > 0
        assert stream_metadata["telemetry"]["streamed_peak_host_bytes"] > 0
        assert stream_metadata["telemetry"]["stream_max_inflight_records"] >= 1
        assert stream_metadata["telemetry"]["pipeline_decoded_items"] == stream_summary["output_records"]
        assert gzip.open(stream_output, "rt", encoding="utf-8").read() == text

        # The bounded path must retain its indexed traversal contract when an
        # interval is supplied, not silently fall back to a whole-file scan.
        stream_indexed_output = work / "stream-indexed-output.vcf.gz"
        stream_indexed_manifest = work / "stream-indexed-output.manifest.json"
        stream_indexed_summary = json.loads(subprocess.check_output([
            str(genotype), "-V", str(gvcf), "-L", "17:69000-69100",
            "-O", str(stream_indexed_output), "--stream-by-locus",
            "--output-manifest", str(stream_indexed_manifest),
        ], text=True).splitlines()[-1])
        assert stream_indexed_summary["stream_by_locus"] is True
        assert stream_indexed_summary["output_records"] == indexed_summary["output_records"]
        stream_indexed_metadata = json.loads(stream_indexed_manifest.read_text(encoding="utf-8"))
        assert stream_indexed_metadata["compatibility"]["indexed_interval_traversal"] is True
        assert stream_indexed_metadata["telemetry"]["indexed_inputs"] == 1
        assert stream_indexed_metadata["telemetry"]["indexed_interval_queries"] >= 1
        assert stream_indexed_metadata["telemetry"]["interval_skipped"] == 0
        assert gzip.open(stream_indexed_output, "rt", encoding="utf-8").read() == gzip.open(
            indexed_output, "rt", encoding="utf-8").read()

        # The explicit GATK-compatible annotation profile keeps the standard
        # finalized INFO/FORMAT set while removing native diagnostics.  This
        # is a field-set contract; the separate GATK oracle owns raw text
        # formatting and INFO ordering.
        strict_output = work / "strict-output.vcf.gz"
        strict_manifest = work / "strict-output.manifest.json"
        strict_summary = json.loads(subprocess.check_output([
            str(genotype), "-V", str(gvcf), "-O", str(strict_output),
            "--gatk-compatible-annotations",
            "--output-manifest", str(strict_manifest),
        ], text=True).splitlines()[-1])
        assert strict_summary["output_records"] == summary["output_records"]
        strict_metadata = json.loads(strict_manifest.read_text(encoding="utf-8"))
        assert strict_metadata["compatibility"]["gatk_annotation_compatibility"] is True
        strict_records = [line.split("\t") for line in gzip.open(
            strict_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert strict_records
        for strict_record in strict_records:
            assert all(tag not in strict_record[7] for tag in
                       ("RAW_MQandDP=", "RCQ=", "RCP="))
            assert "SB" not in strict_record[8].split(":")
        strict_variant = next((record for record in strict_records
                               if record[4] != "."), None)
        # The native-HC fixture does not carry reducible RAW_MQandDP/SB input;
        # the GATK-produced gVCF oracle covers the derived FS/MQ/QD/SOR path.
        has_reducible_annotations = any(
            "RAW_MQandDP=" in record[7] or "SB" in record[8]
            for record in records)
        if strict_variant is not None and has_reducible_annotations:
            assert all(tag in strict_variant[7] for tag in ("FS=", "MQ=", "QD=", "SOR="))

        # The adapter-produced workspace index is a portable bridge for the
        # native path. Opaque GenomicsDB workspaces remain fail-closed and are
        # routed to the configured GATK backend by the dispatcher.
        gendb_workspace = work / "gendb-workspace"
        gendb_workspace.mkdir()
        (gendb_workspace / "fastgatk-inputs.tsv").write_text(
            "# fastgatk GenomicsDB adapter input index v1\n" + str(gvcf.resolve()) + "\n",
            encoding="utf-8",
        )
        gendb_output = work / "gendb-output.vcf.gz"
        gendb_manifest = work / "gendb-output.manifest.json"
        gendb_summary = json.loads(subprocess.check_output([
            str(genotype), "-V", f"gendb://{gendb_workspace}", "-O", str(gendb_output),
            "--output-manifest", str(gendb_manifest),
        ], text=True).splitlines()[-1])
        assert gendb_summary["output_records"] == summary["output_records"]
        assert gendb_output.exists() and Path(f"{gendb_output}.tbi").exists()
        gendb_metadata = json.loads(gendb_manifest.read_text(encoding="utf-8"))
        assert gendb_metadata["telemetry"]["input_files"] == 1
        assert gendb_metadata["compatibility"]["genomicsdb_workspace"] is True
        assert gendb_metadata["telemetry"]["genomicsdb_workspace_input_index"] is True

        native_gendb_workspace = work / "native-gendb-workspace"
        native_gendb_workspace.mkdir()
        native_input = str(gvcf.resolve())
        (native_gendb_workspace / "fastgatk-inputs.tsv").write_text(
            "# fastgatk GenomicsDB adapter input index v1\n" + native_input + "\n",
            encoding="utf-8")
        (native_gendb_workspace / "fastgatk-workspace.json").write_text(
            json.dumps({"schema_version": 1, "backend": "fastgatk-sparse-index",
                        "input_count": 1, "batch_size": 1, "reader_threads": 1,
                        "inputs": [native_input]}) + "\n", encoding="utf-8")
        native_gendb_output = work / "native-gendb-output.vcf.gz"
        native_gendb_summary = json.loads(subprocess.check_output([
            str(genotype), "-V", f"gendb://{native_gendb_workspace}",
            "-O", str(native_gendb_output),
        ], text=True).splitlines()[-1])
        assert native_gendb_summary["output_records"] == summary["output_records"]
        assert native_gendb_output.exists() and Path(f"{native_gendb_output}.tbi").exists()

        # A published native sparse workspace carries a sidecar contract.  A
        # malformed backend tag must fail closed instead of being interpreted
        # as a plain VCF path; index-only legacy workspaces remain accepted by
        # the previous check above.
        malformed_workspace = work / "malformed-native-workspace"
        malformed_workspace.mkdir()
        (malformed_workspace / "fastgatk-inputs.tsv").write_text(
            str(gvcf.resolve()) + "\n", encoding="utf-8")
        (malformed_workspace / "fastgatk-workspace.json").write_text(
            '{"schema_version":1,"backend":"other-backend","input_count":1}\n',
            encoding="utf-8")
        malformed_output = work / "malformed-output.vcf.gz"
        malformed = subprocess.run([
            str(genotype), "-V", f"gendb://{malformed_workspace}", "-O", str(malformed_output),
        ], text=True, capture_output=True)
        assert malformed.returncode != 0
        assert "native GenomicsDB workspace" in malformed.stderr

        # GATK's --gp-qual opt-in is conditional: with the default
        # PL-assignment path and no GP/PP posterior field in the input, the
        # Java engine keeps the ordinary PL-derived QUAL and does not run a
        # posterior-quality update.
        posterior_output = work / "posterior.vcf.gz"
        posterior_manifest = work / "posterior.manifest.json"
        posterior_result = subprocess.run([
            str(genotype), "-V", str(gvcf), "-O", str(posterior_output),
            "--gp-qual", "--heterozygosity", "1e-3",
            "--indel-heterozygosity", "0.000125",
            "--output-manifest", str(posterior_manifest),
        ], text=True, capture_output=True)
        assert posterior_result.returncode == 0, posterior_result.stderr
        posterior_metadata = json.loads(posterior_manifest.read_text(encoding="utf-8"))
        assert posterior_metadata["compatibility"]["genotype_prior_calculator_assumingHW"] is True
        assert posterior_metadata["telemetry"]["posterior_kernel_calls"] == 0
        assert posterior_metadata["telemetry"]["posterior_samples"] == 0
        posterior_records = [line.split("\t") for line in gzip.open(
            posterior_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert posterior_records and all(float(record[5]) >= 0.0 for record in posterior_records)

        multi_input = work / "multi.g.vcf.gz"
        multi_output = work / "multi.vcf.gz"
        multi_manifest = work / "multi.manifest.json"
        multi_header = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##INFO=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
"""
        with gzip.open(multi_input, "wt", encoding="utf-8") as handle:
            handle.write(multi_header)
            # Deliberately make GT disagree with PL; joint materialization
            # must derive the final diploid call and GQ from the likelihoods.
            handle.write("chr1\t1\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=22\tGT:DP:AD:PL\t0/0:10:5,5,0:50,0,50,99,99,99\t0/0:12:0,12,0:100,80,0,99,99,99\n")
        multi_summary = json.loads(subprocess.check_output([
            str(genotype), "-V", str(multi_input), "-O", str(multi_output),
            "--output-manifest", str(multi_manifest),
        ], text=True).splitlines()[-1])
        assert multi_summary["sample_count"] == 2
        multi_text = gzip.open(multi_output, "rt", encoding="utf-8").read()
        assert "\tFORMAT\tS1\tS2" in multi_text
        multi_records = [line.split("\t") for line in multi_text.splitlines()
                         if line and not line.startswith("#")]
        assert len(multi_records) == 1 and ",<NON_REF>" not in multi_records[0][4]
        format_names = multi_records[0][8].split(":")
        assert "GQ" in format_names
        gt_index = format_names.index("GT")
        gq_index = format_names.index("GQ")
        assert multi_records[0][9].split(":")[gt_index] == "0/1"
        assert multi_records[0][10].split(":")[gt_index] == "1/1"
        assert multi_records[0][9].split(":")[gq_index] == "50"
        assert multi_records[0][10].split(":")[gq_index] == "80"
        assert len(multi_records[0][9].split(":")[-2].split(",")) == 3
        assert len(multi_records[0][10].split(":")[-2].split(",")) == 3
        info = {item.split("=", 1)[0]: item.split("=", 1)[1]
                for item in multi_records[0][7].split(";") if "=" in item}
        assert info["AC"] == "3" and info["AN"] == "4"
        assert abs(float(info["AF"]) - 0.75) < 1e-6

        # A spanning-deletion ALT is part of GATK's non-variant genotype set
        # for cohort AF/QUAL, but remains a concrete ALT in the output.  Keep
        # it adjacent to a regular ALT and <NON_REF> to exercise both union
        # remapping and the shared Kokkos spanning-deletion reduction.
        star_input = work / "spanning-deletion.g.vcf.gz"
        star_output = work / "spanning-deletion.vcf.gz"
        star_header = multi_header.replace("\tS1\tS2\n", "\tSTAR\n")
        with gzip.open(star_input, "wt", encoding="utf-8") as handle:
            handle.write(star_header)
            handle.write(
                "chr1\t2\t.\tA\t*,G,<NON_REF>\t.\tPASS\tDP=20\t"
                "GT:DP:AD:PL\t0/1:20:12,8,0,0:0,0,100,100,100,100,100,100,100,100\n"
            )
        star_manifest = work / "spanning-deletion.manifest.json"
        star_result = subprocess.run([
            str(genotype), "-V", str(star_input), "-O", str(star_output),
            "--output-manifest", str(star_manifest),
        ], text=True, capture_output=True)
        assert star_result.returncode == 0, star_result.stderr
        star_records = [line.split("\t") for line in gzip.open(
            star_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(star_records) == 1 and "*" in star_records[0][4]
        star_metadata = json.loads(star_manifest.read_text(encoding="utf-8"))
        assert star_metadata["compatibility"]["spanning_deletion_nonvariant_set"] is True
        star_posterior_output = work / "spanning-deletion-posterior.vcf.gz"
        star_posterior_manifest = work / "spanning-deletion-posterior.manifest.json"
        star_posterior_result = subprocess.run([
            str(genotype), "-V", str(star_input), "-O", str(star_posterior_output),
            "--gp-qual", "--output-manifest", str(star_posterior_manifest),
        ], text=True, capture_output=True)
        assert star_posterior_result.returncode == 0, star_posterior_result.stderr
        star_posterior_records = [line.split("\t") for line in gzip.open(
            star_posterior_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(star_posterior_records) == 1
        assert float(star_posterior_records[0][5]) <= float(star_records[0][5])
        star_posterior_metadata = json.loads(star_posterior_manifest.read_text(encoding="utf-8"))
        assert star_posterior_metadata["compatibility"]["spanning_deletion_posterior_qual"] is True
        assert star_posterior_metadata["telemetry"]["posterior_kernel_calls"] == 0
        assert star_posterior_metadata["telemetry"]["posterior_samples"] == 0

        # ExcessHet uses GATK's diploid, biallelic Wigginton right-tail test.
        # Four all-heterozygous samples have an exact p-value of 0.2285714,
        # hence a Phred score of 6.4098.  This exercises the cohort annotation
        # rather than the single-sample degenerate 0.0000 case.
        excess_input = work / "excess-het.g.vcf.gz"
        excess_output = work / "excess-het.vcf.gz"
        excess_header = multi_header.replace("\tS1\tS2\n", "\tS1\tS2\tS3\tS4\n")
        with gzip.open(excess_input, "wt", encoding="utf-8") as handle:
            handle.write(excess_header)
            sample = "0/0:10:5,5,0:99,0,99,99,99,99"
            handle.write("chr1\t1\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=40\t"
                         "GT:DP:AD:PL\t" + "\t".join([sample] * 4) + "\n")
        excess_result = subprocess.run([
            str(genotype), "-V", str(excess_input), "-O", str(excess_output),
            "--gatk-compatible-annotations",
        ], text=True, capture_output=True)
        assert excess_result.returncode == 0, excess_result.stderr
        excess_text = gzip.open(excess_output, "rt", encoding="utf-8").read()
        assert "##INFO=<ID=ExcessHet," in excess_text
        excess_records = [line.split("\t") for line in excess_text.splitlines()
                          if line and not line.startswith("#")]
        assert len(excess_records) == 1
        excess_info = {item.split("=", 1)[0]: item.split("=", 1)[1]
                       for item in excess_records[0][7].split(";") if "=" in item}
        assert abs(float(excess_info["ExcessHet"]) - 6.4098) < 1.0e-3

        # The default remains GATK's PREFER_PLS assignment.  The explicit
        # USE_POSTERIOR_PROBABILITIES mode uses the final cohort AF posterior
        # and can therefore select a different call for an intentionally
        # ambiguous PL row; its Kokkos assignment telemetry must be visible.
        posterior_assignment_input = work / "posterior-assignment.g.vcf.gz"
        posterior_assignment_output = work / "posterior-assignment.vcf.gz"
        posterior_assignment_manifest = work / "posterior-assignment.manifest.json"
        posterior_assignment_header = multi_header.replace("\tS1\tS2\n", "\tPA\n")
        with gzip.open(posterior_assignment_input, "wt", encoding="utf-8") as handle:
            handle.write(posterior_assignment_header)
            handle.write("chr1\t1\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=10\t"
                         "GT:DP:AD:PL\t0/1:10:5,5,0:5,0,5,99,99,99\n")
        posterior_assignment_result = subprocess.run([
            str(genotype), "-V", str(posterior_assignment_input),
            "-O", str(posterior_assignment_output),
            "--genotype-assignment-method", "USE_POSTERIOR_PROBABILITIES",
            "--gp-qual",
            "--output-manifest", str(posterior_assignment_manifest),
        ], text=True, capture_output=True)
        assert posterior_assignment_result.returncode == 0, posterior_assignment_result.stderr
        posterior_assignment_metadata = json.loads(
            posterior_assignment_manifest.read_text(encoding="utf-8"))
        assert posterior_assignment_metadata["compatibility"]["posterior_genotype_assignment"] is True
        assert posterior_assignment_metadata["telemetry"]["posterior_assignment_kernel_calls"] > 0
        assert posterior_assignment_metadata["telemetry"]["posterior_assignment_samples"] > 0
        assert posterior_assignment_metadata["telemetry"]["posterior_kernel_calls"] > 0
        posterior_assignment_records = [line.split("\t") for line in gzip.open(
            posterior_assignment_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(posterior_assignment_records) == 1
        posterior_assignment_format = posterior_assignment_records[0][8].split(":")
        posterior_assignment_sample = posterior_assignment_records[0][9].split(":")
        assert posterior_assignment_sample[posterior_assignment_format.index("GT")] == "0/0"
        assert "GP" in posterior_assignment_format and "PG" in posterior_assignment_format
        gp_values = posterior_assignment_sample[posterior_assignment_format.index("GP")].split(",")
        pg_values = posterior_assignment_sample[posterior_assignment_format.index("PG")].split(",")
        assert len(gp_values) == 3 and len(pg_values) == 3
        assert all(float(value) >= 0.0 for value in gp_values + pg_values)

        # USE_POSTERIORS_ANNOTATION consumes the existing GATK FORMAT/PP
        # vector after the same Number=G allele projection.  PP deliberately
        # disagrees with PL here, so the selected GT must follow PP and the
        # annotation must remain in the output record.
        posterior_annotation_input = work / "posterior-annotation.g.vcf.gz"
        posterior_annotation_output = work / "posterior-annotation.vcf.gz"
        posterior_annotation_manifest = work / "posterior-annotation.manifest.json"
        posterior_annotation_header = multi_header.replace(
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>\n",
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>\n"
            "##FORMAT=<ID=PP,Number=G,Type=Integer,Description=Phred-scaled Posterior Genotype Probabilities>\n")
        posterior_annotation_header = posterior_annotation_header.replace("\tS1\tS2\n", "\tPA\n")
        with gzip.open(posterior_annotation_input, "wt", encoding="utf-8") as handle:
            handle.write(posterior_annotation_header)
            handle.write("chr1\t1\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=10\t"
                         "GT:DP:AD:PL:PP\t0/0:10:5,5,0:0,50,100,99,99,99:80,0,20,99,99,99\n")
        posterior_annotation_result = subprocess.run([
            str(genotype), "-V", str(posterior_annotation_input),
            "-O", str(posterior_annotation_output),
            "--genotype-assignment-method", "USE_POSTERIORS_ANNOTATION",
            "--output-manifest", str(posterior_annotation_manifest),
        ], text=True, capture_output=True)
        assert posterior_annotation_result.returncode == 0, posterior_annotation_result.stderr
        posterior_annotation_metadata = json.loads(
            posterior_annotation_manifest.read_text(encoding="utf-8"))
        assert posterior_annotation_metadata["compatibility"]["use_posteriors_annotation"] is True
        assert posterior_annotation_metadata["telemetry"][
            "posterior_annotation_assignment_kernel_calls"] > 0
        posterior_annotation_records = [line.split("\t") for line in gzip.open(
            posterior_annotation_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(posterior_annotation_records) == 1
        posterior_annotation_format = posterior_annotation_records[0][8].split(":")
        posterior_annotation_sample = posterior_annotation_records[0][9].split(":")
        assert posterior_annotation_sample[posterior_annotation_format.index("GT")] == "0/1"
        assert posterior_annotation_sample[posterior_annotation_format.index("PP")] == "80,0,20"
        assert posterior_annotation_sample[posterior_annotation_format.index("GQ")] == "20"

        # SET_TO_NO_CALL_NO_ANNOTATIONS follows the Java GenotypeBuilder
        # contract: GT becomes no-call and every genotype-level annotation is
        # removed, including FORMAT fields not materialized by Record.
        no_annotation_output = work / "no-annotation.vcf.gz"
        no_annotation_input = work / "no-annotation.g.vcf.gz"
        no_annotation_manifest = work / "no-annotation.manifest.json"
        no_annotation_header = multi_header.replace(
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>\n",
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>\n"
            "##FORMAT=<ID=PP,Number=G,Type=Integer,Description=Posterior>\n"
            "##FORMAT=<ID=GP,Number=G,Type=Float,Description=Posterior>\n"
            "##FORMAT=<ID=PG,Number=G,Type=Float,Description=Prior>\n"
            "##FORMAT=<ID=FT,Number=1,Type=String,Description=Genotype filter>\n"
            "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=Minimum depth>\n")
        no_annotation_header = no_annotation_header.replace("\tS1\tS2\n", "\tNA\n")
        with gzip.open(no_annotation_input, "wt", encoding="utf-8") as handle:
            handle.write(no_annotation_header)
            handle.write("chr1\t1\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=10\t"
                         "GT:DP:AD:PL:PP:GP:PG:FT:MIN_DP\t"
                         "0/1:10:5,5,0:0,50,100,99,99,99:80,0,20,99,99,99:"
                         "80,0,20:80,0,20:LowQual:3\n")
        no_annotation_result = subprocess.run([
            str(genotype), "-V", str(no_annotation_input),
            "-O", str(no_annotation_output),
            "--genotype-assignment-method", "SET_TO_NO_CALL_NO_ANNOTATIONS",
            "--output-manifest", str(no_annotation_manifest),
        ], text=True, capture_output=True)
        assert no_annotation_result.returncode == 0, no_annotation_result.stderr
        no_annotation_records = [line.split("\t") for line in gzip.open(
            no_annotation_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(no_annotation_records) == 1
        no_annotation_format = no_annotation_records[0][8].split(":")
        no_annotation_sample = no_annotation_records[0][9].split(":")
        assert no_annotation_sample[no_annotation_format.index("GT")] == "./."
        assert no_annotation_format == ["GT"]
        no_annotation_metadata = json.loads(no_annotation_manifest.read_text(encoding="utf-8"))
        assert no_annotation_metadata["compatibility"]["all_format_annotation_cleanup"] is True
        assert no_annotation_metadata["telemetry"]["no_annotation_format_fields"] >= 8

        # BEST_MATCH_TO_ORIGINAL keeps the input GT even when PL prefers a
        # different genotype; GQ=0 plus PL[0]=0 follows GATK's no-call rule.
        best_match_input = work / "best-match.g.vcf.gz"
        best_match_output = work / "best-match.vcf.gz"
        best_match_header = multi_header.replace("\tS1\tS2\n", "\tBM\n")
        with gzip.open(best_match_input, "wt", encoding="utf-8") as handle:
            handle.write(best_match_header)
            handle.write("chr1\t1\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=10\t"
                         "GT:DP:AD:PL\t1/1:10:0,0,10:0,40,80,99,99,99\n")
        best_match_result = subprocess.run([
            str(genotype), "-V", str(best_match_input), "-O", str(best_match_output),
            "--genotype-assignment-method", "BEST_MATCH_TO_ORIGINAL",
        ], text=True, capture_output=True)
        assert best_match_result.returncode == 0, best_match_result.stderr
        best_match_records = [line.split("\t") for line in gzip.open(
            best_match_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(best_match_records) == 1
        best_match_format = best_match_records[0][8].split(":")
        best_match_sample = best_match_records[0][9].split(":")
        assert best_match_sample[best_match_format.index("GT")] == "1/1"
        best_match_zero_input = work / "best-match-zero.g.vcf.gz"
        best_match_zero_output = work / "best-match-zero.vcf.gz"
        best_match_zero_header = best_match_header.replace(
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n",
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n"
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>\n")
        with gzip.open(best_match_zero_input, "wt", encoding="utf-8") as handle:
            handle.write(best_match_zero_header)
            handle.write("chr1\t1\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=10\t"
                         "GT:GQ:DP:AD:PL\t1/1:0:10:0,0,10:0,40,80,99,99,99\n")
        best_match_zero_result = subprocess.run([
            str(genotype), "-V", str(best_match_zero_input), "-O", str(best_match_zero_output),
            "--genotype-assignment-method", "BEST_MATCH_TO_ORIGINAL",
        ], text=True, capture_output=True)
        assert best_match_zero_result.returncode == 0, best_match_zero_result.stderr
        best_match_zero_records = [line.split("\t") for line in gzip.open(
            best_match_zero_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        zero_format = best_match_zero_records[0][8].split(":")
        zero_sample = best_match_zero_records[0][9].split(":")
        assert zero_sample[zero_format.index("GT")] == "./."

        # The same triangular/combinatorial kernel must handle non-diploid
        # VCF PL ordering (triploid biallelic order: 000,001,011,111).
        triploid_input = work / "triploid.g.vcf.gz"
        triploid_output = work / "triploid.vcf.gz"
        triploid_header = multi_header.replace("\tS1\tS2\n", "\tT1\n")
        with gzip.open(triploid_input, "wt", encoding="utf-8") as handle:
            handle.write(triploid_header)
            handle.write("chr1\t1\t.\tA\tG\t.\tPASS\tDP=10\tGT:DP:AD:PL\t"
                         "0/0/0:10:5,5:50,0,80,100\n")
        subprocess.run([
            str(genotype), "-V", str(triploid_input), "-O", str(triploid_output),
        ], check=True, stdout=subprocess.DEVNULL)
        triploid_records = [line.split("\t") for line in gzip.open(
            triploid_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(triploid_records) == 1
        triploid_format = triploid_records[0][8].split(":")
        triploid_sample = triploid_records[0][9].split(":")
        assert triploid_sample[triploid_format.index("GT")] == "0/0/1"
        assert triploid_sample[triploid_format.index("GQ")] == "50"
        triploid_info = {item.split("=", 1)[0]: item.split("=", 1)[1]
                         for item in triploid_records[0][7].split(";") if "=" in item}
        assert triploid_info["AC"] == "1" and triploid_info["AN"] == "3"

        # GATK's AlleleFrequencyCalculator also accepts a diploid hom-ref
        # genotype with GQ but no PL (common in reblocked/reference-confidence
        # inputs) and synthesizes an approximate PL row.  The native cohort
        # API must keep that compatibility path on the shared Kokkos backend.
        gq_only_input = work / "gq-only.g.vcf.gz"
        gq_only_output = work / "gq-only.vcf.gz"
        gq_only_manifest = work / "gq-only.manifest.json"
        gq_only_header = multi_header.replace(
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>\n", "")
        gq_only_header = gq_only_header.replace(
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n",
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n"
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>\n")
        gq_only_header = gq_only_header.replace("\tS1\tS2\n", "\tGQ_ONLY\n")
        with gzip.open(gq_only_input, "wt", encoding="utf-8") as handle:
            handle.write(gq_only_header)
            handle.write("chr1\t1\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=10\t"
                         "GT:GQ:AD\t0/0:30:10,0,0\n")
        gq_only_summary = json.loads(subprocess.check_output([
            str(genotype), "-V", str(gq_only_input), "-O", str(gq_only_output),
            "--output-manifest", str(gq_only_manifest),
        ], text=True).splitlines()[-1])
        assert gq_only_summary["output_records"] == 1
        gq_only_metadata = json.loads(gq_only_manifest.read_text(encoding="utf-8"))
        assert gq_only_metadata["compatibility"]["gq_only_homref_af_approximation"] is True
        assert gq_only_metadata["telemetry"]["cohort_af_approximate_gq_samples"] == 1
        assert gq_only_metadata["telemetry"]["cohort_af_samples"] == 1
        gq_only_records = [line.split("\t") for line in gzip.open(
            gq_only_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(gq_only_records) == 1
        gq_only_info = {item.split("=", 1)[0]: item.split("=", 1)[1]
                        for item in gq_only_records[0][7].split(";") if "=" in item}
        assert gq_only_info["AN"] == "2" and gq_only_info["AC"] == "0"

        # Preserve all concrete ALTs from a multi-allelic gVCF and remap the
        # triangular PL vector using the VCF genotype-index formula.
        multi_alt_input = work / "multi-alt.g.vcf.gz"
        multi_alt_output = work / "multi-alt.vcf.gz"
        with gzip.open(multi_alt_input, "wt", encoding="utf-8") as handle:
            handle.write(multi_header)
            handle.write("chr1\t1\t.\tA\tG,T,<NON_REF>\t.\tPASS\tDP=22\tGT:DP:AD:PL\t"
                         "0/2:10:5,2,3,0:50,20,40,30,10,0\t"
                         "1/2:12:0,4,8,0:90,80,70,60,50,0\n")
        subprocess.run([
            str(genotype), "-V", str(multi_alt_input), "-O", str(multi_alt_output),
        ], check=True, stdout=subprocess.DEVNULL)
        multi_alt_records = [line.split("\t") for line in gzip.open(
            multi_alt_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(multi_alt_records) == 1
        assert multi_alt_records[0][4] == "G,T"
        multi_alt_format = multi_alt_records[0][8].split(":")
        assert len(multi_alt_records[0][9].split(":")[multi_alt_format.index("AD")].split(",")) == 3
        assert len(multi_alt_records[0][9].split(":")[multi_alt_format.index("PL")].split(",")) == 6
        assert len(multi_alt_records[0][10].split(":")[multi_alt_format.index("AD")].split(",")) == 3
        assert len(multi_alt_records[0][10].split(":")[multi_alt_format.index("PL")].split(",")) == 6

        # Joint shards can carry disjoint ALT subsets at one locus.  The
        # merge must form the union and remap each shard's triangular PL into
        # the joint allele order before writing one record.
        shard_g = work / "shard-g.g.vcf.gz"
        shard_t = work / "shard-t.g.vcf.gz"
        header_g = multi_header.replace("\tS1\tS2\n", "\tSG\n")
        header_t = multi_header.replace("\tS1\tS2\n", "\tST\n")
        with gzip.open(shard_g, "wt", encoding="utf-8") as handle:
            handle.write(header_g)
            handle.write("chr1\t1\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=10\tGT:DP:AD:PL\t0/1:10:5,5,0:50,0,50,99,99,99\n")
        with gzip.open(shard_t, "wt", encoding="utf-8") as handle:
            handle.write(header_t)
            handle.write("chr1\t1\t.\tA\tT,<NON_REF>\t.\tPASS\tDP=11\tGT:DP:AD:PL\t1/1:11:0,0,11:100,99,0,99,99,99\n")
        union_output = work / "union.vcf.gz"
        union_result = subprocess.run([
            str(genotype), "-V", str(shard_g), "-V", str(shard_t),
            "-O", str(union_output),
        ], text=True, capture_output=True)
        assert union_result.returncode == 0, union_result.stderr
        union_records = [line.split("\t") for line in gzip.open(
            union_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(union_records) == 1 and union_records[0][4] == "G,T"
        union_format = union_records[0][8].split(":")
        assert len(union_records[0][9].split(":")[union_format.index("AD")].split(",")) == 3
        assert len(union_records[0][9].split(":")[union_format.index("PL")].split(",")) == 6
        assert len(union_records[0][10].split(":")[union_format.index("AD")].split(",")) == 3
        assert len(union_records[0][10].split(":")[union_format.index("PL")].split(",")) == 6
        union_info = {item.split("=", 1)[0]: item.split("=", 1)[1]
                      for item in union_records[0][7].split(";") if "=" in item}
        assert union_info["AC"] == "1,2" and union_info["AN"] == "4"
        assert [float(value) for value in union_info["AF"].split(",")] == [0.25, 0.5]

        # A shard may present the same ALT union in a different order.  The
        # Kokkos projection sorts the mapped source tuple before ranking, so
        # this remains valid instead of reading the wrong PL row.
        reverse_input = work / "reverse-order.g.vcf.gz"
        reverse_header = multi_header.replace("\tS1\tS2\n", "\tSR\n")
        with gzip.open(reverse_input, "wt", encoding="utf-8") as handle:
            handle.write(reverse_header)
            handle.write("chr1\t1\t.\tA\tT,G,<NON_REF>\t.\tPASS\tDP=11\t"
                         "GT:DP:AD:PL\t1/2:11:1,4,6,0:90,80,70,60,50,0\n")
        reverse_output = work / "reverse-order.vcf.gz"
        reverse_result = subprocess.run([
            str(genotype), "-V", str(shard_g), "-V", str(reverse_input),
            "-O", str(reverse_output),
        ], text=True, capture_output=True)
        assert reverse_result.returncode == 0, reverse_result.stderr
        reverse_records = [line.split("\t") for line in gzip.open(
            reverse_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(reverse_records) == 1 and reverse_records[0][4] == "G,T"
        reverse_format = reverse_records[0][8].split(":")
        assert len(reverse_records[0][10].split(":")[reverse_format.index("PL")].split(",")) == 6

        # A reference-only GVCF block from one sample must participate in a
        # concrete site emitted by another sample.  Its REF/<NON_REF> PL is
        # projected onto REF/concrete-ALT genotypes, yielding a finite 0/0
        # call rather than silently dropping the sample from the cohort.
        variant_output = work / "variant-plus-reference-block.vcf.gz"
        variant_header = multi_header.replace("\tS1\tS2\n", "\tVA\n")
        reference_header = multi_header.replace("\tS1\tS2\n", "\tRB\n")
        variant_part = work / "variant-part.g.vcf.gz"
        reference_part = work / "reference-part.g.vcf.gz"
        with gzip.open(variant_part, "wt", encoding="utf-8") as handle:
            handle.write(variant_header)
            handle.write("chr1\t1\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=10\tGT:DP:AD:PL\t0/1:10:5,5,0:50,0,50,99,99,99\n")
        with gzip.open(reference_part, "wt", encoding="utf-8") as handle:
            handle.write(reference_header)
            handle.write("chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tDP=12\tGT:DP:AD:PL\t0/0:12:12,0:0,60,99\n")
        reference_result = subprocess.run([
            str(genotype), "-V", str(variant_part), "-V", str(reference_part),
            "-O", str(variant_output),
        ], text=True, capture_output=True)
        assert reference_result.returncode == 0, reference_result.stderr
        reference_records = [line.split("\t") for line in gzip.open(
            variant_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(reference_records) == 1 and reference_records[0][4] == "G"
        reference_format = reference_records[0][8].split(":")
        assert "<NON_REF>" not in "\t".join(reference_records[0])
        assert reference_records[0][9].split(":")[reference_format.index("GT")] == "0/1"
        assert reference_records[0][10].split(":")[reference_format.index("GT")] == "0/0"
        assert reference_records[0][10].split(":")[reference_format.index("GQ")] == "60"
        assert len(reference_records[0][10].split(":")[reference_format.index("PL")].split(",")) == 3

        # A block spanning a concrete site from another sample is split using
        # the indexed reference, and the point segment is merged into the
        # variant locus instead of dropping the reference sample.
        span_reference = work / "chr1.fa"
        span_reference.write_text(">chr1\n" + "A" * 100 + "\n", encoding="utf-8")
        span_reference.with_suffix(span_reference.suffix + ".fai").write_text(
            "chr1\t100\t6\t100\t101\n", encoding="utf-8")
        span_header = multi_header.replace(
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>\n",
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=End position>\n")
        span_variant = work / "span-variant.g.vcf.gz"
        span_block = work / "span-block.g.vcf.gz"
        span_variant_header = span_header.replace("\tS1\tS2\n", "\tSV\n")
        span_block_header = span_header.replace("\tS1\tS2\n", "\tSB\n")
        with gzip.open(span_variant, "wt", encoding="utf-8") as handle:
            handle.write(span_variant_header)
            handle.write("chr1\t2\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=10\t"
                         "GT:DP:AD:PL\t0/1:10:5,5,0:50,0,50,99,99,99\n")
        with gzip.open(span_block, "wt", encoding="utf-8") as handle:
            handle.write(span_block_header)
            handle.write("chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=4;DP=12\t"
                         "GT:DP:AD:PL\t0/0:12:12,0:0,60,99\n")
        span_output = work / "span-joint.vcf.gz"
        span_result = subprocess.run([
            str(genotype), "-R", str(span_reference), "-V", str(span_variant),
            "-V", str(span_block), "-O", str(span_output),
        ], text=True, capture_output=True)
        assert span_result.returncode == 0, span_result.stderr
        span_records = [line.split("\t") for line in gzip.open(
            span_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(span_records) == 1 and span_records[0][1] == "2"
        span_format = span_records[0][8].split(":")
        assert span_records[0][10].split(":")[span_format.index("GT")] == "0/0"

        # GATK's optional --include-non-variant-sites path retains pure
        # reference-confidence blocks as REF-only records instead of dropping
        # the locus group.  The native contract compacts each single-sample
        # PL to the sole REF genotype and preserves sample DP/AD/GQ.
        ref_only_input = work / "reference-only.g.vcf.gz"
        ref_only_output = work / "reference-only.vcf.gz"
        ref_only_manifest = work / "reference-only.manifest.json"
        ref_only_header = multi_header.replace("\tS1\tS2\n", "\tR1\n")
        with gzip.open(ref_only_input, "wt", encoding="utf-8") as handle:
            handle.write(ref_only_header)
            handle.write("chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tDP=10\tGT:DP:AD:PL\t0/0:10:10,0:0,60,99\n")
            handle.write("chr1\t2\t.\tC\t<NON_REF>\t.\tPASS\tDP=8\tGT:DP:AD:PL\t0/0:8:8,0:0,50,99\n")
        ref_only_result = subprocess.run([
            str(genotype), "-V", str(ref_only_input), "-O", str(ref_only_output),
            "--include-non-variant-sites", "--output-manifest", str(ref_only_manifest),
        ], text=True, capture_output=True)
        assert ref_only_result.returncode == 0, ref_only_result.stderr
        ref_only_records = [line.split("\t") for line in gzip.open(
            ref_only_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(ref_only_records) == 2
        for record in ref_only_records:
            assert record[4] == "."
            format_names = record[8].split(":")
            sample = record[9].split(":")
            assert sample[format_names.index("GT")] == "0/0"
            assert len(sample[format_names.index("PL")].split(",")) == 1
        ref_only_metadata = json.loads(ref_only_manifest.read_text(encoding="utf-8"))
        assert ref_only_metadata["compatibility"]["include_non_variant_sites"] is True
        assert ref_only_metadata["telemetry"]["include_non_variant_sites"] is True

        # Dense expansion is also available in the bounded stream path.  A
        # multi-base reference block must be emitted one coordinate at a time
        # and remain byte-identical to aggregate materialization.
        dense_aggregate = work / "dense-aggregate.vcf.gz"
        dense_stream = work / "dense-stream.vcf.gz"
        dense_stream_manifest = work / "dense-stream.manifest.json"
        dense_aggregate_result = subprocess.run([
            str(genotype), "-R", str(span_reference), "-V", str(span_block),
            "-O", str(dense_aggregate), "--include-non-variant-sites",
        ], text=True, capture_output=True)
        assert dense_aggregate_result.returncode == 0, dense_aggregate_result.stderr
        dense_stream_result = subprocess.run([
            str(genotype), "-R", str(span_reference), "-V", str(span_block),
            "-O", str(dense_stream), "--include-non-variant-sites",
            "--stream-by-locus", "--output-manifest", str(dense_stream_manifest),
        ], text=True, capture_output=True)
        assert dense_stream_result.returncode == 0, dense_stream_result.stderr
        assert gzip.open(dense_stream, "rt", encoding="utf-8").read() == gzip.open(
            dense_aggregate, "rt", encoding="utf-8").read()
        dense_metadata = json.loads(dense_stream_manifest.read_text(encoding="utf-8"))
        assert dense_metadata["compatibility"]["include_non_variant_sites"] is True
        assert dense_metadata["telemetry"]["include_non_variant_sites"] is True
        assert dense_metadata["telemetry"]["streamed_loci"] == 4

        # Shards of one already-merged cohort remain deterministic, while
        # disjoint sample headers are merged into one joint record with
        # missing fields preserved for samples absent from a shard.
        shard_output = work / "shards.vcf.gz"
        shard_result = subprocess.run([
            str(genotype), "-V", str(multi_input), "-V", str(multi_input),
            "-O", str(shard_output),
        ], text=True, capture_output=True)
        assert shard_result.returncode == 0, shard_result.stderr

        different_input = work / "different-samples.g.vcf.gz"
        different_header = multi_header.replace("\tS1\tS2\n", "\tS3\tS4\n")
        with gzip.open(different_input, "wt", encoding="utf-8") as handle:
            handle.write(different_header)
            handle.write("chr1\t1\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=22\tGT:DP:AD:PL\t0/1:10:5,5,0:50,0,50,99,99,99\t1/1:12:0,12,0:100,80,0,99,99,99\n")
        different_output = work / "different.vcf.gz"
        different_result = subprocess.run([
            str(genotype), "-V", str(multi_input), "-V", str(different_input),
            "-O", str(different_output),
        ], text=True, capture_output=True)
        assert different_result.returncode == 0, different_result.stderr
        different_text = gzip.open(different_output, "rt", encoding="utf-8").read()
        assert "\tFORMAT\tS1\tS2\tS3\tS4" in different_text
        different_records = [line.split("\t") for line in different_text.splitlines()
                             if line and not line.startswith("#")]
        assert len(different_records) == 1 and len(different_records[0]) == 13

        # Repeated interval selectors use union semantics and are applied
        # before locus grouping/materialization.
        interval_input = work / "intervals.g.vcf.gz"
        interval_output = work / "intervals.vcf.gz"
        interval_manifest = work / "intervals.manifest.json"
        with gzip.open(interval_input, "wt", encoding="utf-8") as handle:
            handle.write(multi_header)
            for position in (1, 2, 3):
                handle.write(
                    f"chr1\t{position}\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=10\t"
                    "GT:DP:AD:PL\t0/1:10:5,5,0:50,0,50,99,99,99\t"
                    "0/0:10:10,0,0:0,40,80,99,99,99\n"
                )
        interval_result = subprocess.run([
            str(genotype), "-V", str(interval_input), "-L", "chr1:1-1",
            "--intervals", "chr1:2-2", "-O", str(interval_output),
            "--output-manifest", str(interval_manifest),
        ], text=True, capture_output=True)
        assert interval_result.returncode == 0, interval_result.stderr
        interval_summary = json.loads(interval_result.stdout.splitlines()[-1])
        assert interval_summary["interval_skipped"] == 1
        interval_records = [line.split("\t") for line in gzip.open(
            interval_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [record[1] for record in interval_records] == ["1", "2"]
        interval_metadata = json.loads(interval_manifest.read_text(encoding="utf-8"))
        assert interval_metadata["compatibility"]["interval_subset"] is True
        assert interval_metadata["telemetry"]["intervals"] == 2
        interval_list = work / "genotype.interval_list"
        interval_list.write_text("@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:100\n"
                                 "chr1\t1\t1\t+\tfirst\nchr1\t2\t2\t+\tsecond\n",
                                 encoding="utf-8")
        interval_file_output = work / "interval-file.vcf.gz"
        interval_file_manifest = work / "interval-file.manifest.json"
        interval_file_result = subprocess.run([
            str(genotype), "-V", str(interval_input), "-L", str(interval_list),
            "-O", str(interval_file_output), "--output-manifest", str(interval_file_manifest),
        ], text=True, capture_output=True)
        assert interval_file_result.returncode == 0, interval_file_result.stderr
        interval_file_metadata = json.loads(interval_file_manifest.read_text(encoding="utf-8"))
        assert interval_file_metadata["telemetry"]["interval_list_inputs"] == 1
        assert interval_file_metadata["telemetry"]["interval_list_records"] == 2

        # gVCF interval selection is span-overlap based: a reference block at
        # POS=1 with END=5 is retained by a chr1:4-4 selector.
        block_input = work / "block.g.vcf.gz"
        block_output = work / "block.vcf.gz"
        block_header = (
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n"
            "##ALT=<ID=NON_REF,Description=Any alternate>\n"
            "##INFO=<ID=END,Number=1,Type=Integer,Description=End position>\n"
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depth>\n"
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        with gzip.open(block_input, "wt", encoding="utf-8") as handle:
            handle.write(block_header +
                         "chr1\t1\t.\tA\tG,<NON_REF>\t.\tPASS\tEND=5;DP=10\t"
                         "GT:DP:AD:PL\t0/1:10:5,5,0:50,0,50,99,99,99\n")
        block_result = subprocess.run([
            str(genotype), "-V", str(block_input), "-L", "chr1:4-4", "-O", str(block_output),
        ], text=True, capture_output=True)
        assert block_result.returncode == 0, block_result.stderr
        block_records = [line for line in gzip.open(
            block_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(block_records) == 1 and block_records[0].split("\t")[1] == "1"

        # A pure REF/<NON_REF> block is compacted to a REF-only output row in
        # --include-non-variant-sites mode.  The writer-facing PL then has
        # width one, but the cross-sample reference posterior must still use
        # the original width-three rows from every contributing sample.  This
        # fixture verifies that RCQ/RCP are produced after the sample merge and
        # that the Kokkos telemetry counts both samples rather than only the
        # first source row.
        reference_only_input = work / "reference-only-multisample.g.vcf.gz"
        reference_only_output = work / "reference-only-multisample.vcf.gz"
        reference_only_manifest = work / "reference-only-multisample.manifest.json"
        reference_only_header = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##ALT=<ID=NON_REF,Description=Any alternate>
##INFO=<ID=END,Number=1,Type=Integer,Description=End position>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depth>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
"""
        with gzip.open(reference_only_input, "wt", encoding="utf-8") as handle:
            handle.write(reference_only_header)
            handle.write(
                "chr1\t10\t.\tA\t<NON_REF>\t.\tPASS\tEND=20;DP=20\t"
                "GT:DP:AD:GQ:PL\t0/0:10:10,0:30:0,30,90\t"
                "0/0:10:10,0:20:0,20,60\n"
            )
        reference_only_result = subprocess.run([
            str(genotype), "-V", str(reference_only_input), "-O", str(reference_only_output),
            "--include-non-variant-sites", "--output-manifest", str(reference_only_manifest),
        ], text=True, capture_output=True)
        assert reference_only_result.returncode == 0, reference_only_result.stderr
        reference_only_text = gzip.open(
            reference_only_output, "rt", encoding="utf-8").read()
        reference_only_records = [line.split("\t") for line in reference_only_text.splitlines()
                                  if line and not line.startswith("#")]
        # Dense non-variant mode expands the source block into one record per
        # reference base; inspect the first materialized point while checking
        # that every emitted point carries the cross-sample diagnostics.
        assert len(reference_only_records) == 11
        assert all("RCQ=" in record[7] and "RCP=" in record[7]
                   for record in reference_only_records)
        reference_only_record = reference_only_records[0]
        assert reference_only_record[3:5] == ["A", "."]
        reference_only_format = reference_only_record[8].split(":")
        assert reference_only_format == ["GT", "DP", "AD", "GQ", "PL"]
        assert reference_only_record[9].split(":")[reference_only_format.index("GT")] == "0/0"
        assert reference_only_record[10].split(":")[reference_only_format.index("GT")] == "0/0"
        reference_only_metadata = json.loads(reference_only_manifest.read_text(encoding="utf-8"))
        assert reference_only_metadata["telemetry"]["cross_sample_reference_kernel_calls"] > 0
        assert reference_only_metadata["telemetry"]["cross_sample_reference_samples"] == 22
    print(json.dumps({"status": "pass", "output_records": summary["output_records"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
