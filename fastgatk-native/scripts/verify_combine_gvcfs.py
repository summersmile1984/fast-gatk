#!/usr/bin/env python3
"""Verify deterministic shard merge and reference-block coalescing."""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=END,Number=1,Type=Integer,Description=End position of reference block>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Read depth>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def write_vcf(path: Path, body: str, header: str = HEADER) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as output:
        output.write(header)
        output.write(body)


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build")) / "fastgatk-combine-gvcfs"
    java = root / "third_party/jdk17/bin/java"
    gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    with tempfile.TemporaryDirectory(prefix="fastgatk-combine-gvcfs-") as directory:
        work = Path(directory)
        shard1 = work / "shard1.g.vcf.gz"
        shard2 = work / "shard2.g.vcf.gz"
        sample1 = work / "sample1.g.vcf.gz"
        sample2 = work / "sample2.g.vcf.gz"
        output = work / "combined.g.vcf.gz"
        manifest = work / "combined.g.vcf.gz.manifest.json"
        write_vcf(shard1, "chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=3;DP=10\n"
                  "chr1\t8\t.\tA\tG\t.\tPASS\tDP=10\n")
        # Distinct INFO DP values must not prevent adjacent reference-block
        # coalescing; FORMAT/sample calls are handled separately below.
        write_vcf(shard2, "chr1\t4\t.\tA\t<NON_REF>\t.\tPASS\tEND=7;DP=12\n"
                  "chr1\t8\t.\tA\tG\t.\tPASS\tDP=10\n")
        sample_header = HEADER.replace(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>\n"
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n"
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n")
        write_vcf(sample1, "chr1\t8\t.\tA\tG\t.\tPASS\tDP=10\tGT:DP:AD:PL\t0/1:10:5,5:50,0,50\n",
                  sample_header)
        write_vcf(sample2, "chr1\t8\t.\tA\tG\t.\tPASS\tDP=12\tGT:DP:AD:PL\t1/1:12:0,12:100,80,0\n",
                  sample_header.replace("SAMPLE1", "SAMPLE2"))
        command = [str(binary), "-V", str(shard1), "-V", str(shard2),
                   "-O", str(output), "--output-manifest", str(manifest)]
        env = os.environ.copy()
        result = json.loads(subprocess.check_output(command, text=True, env=env).splitlines()[-1])
        assert result["tool"] == "CombineGVCFs" and result["status"] == "prototype"
        assert result["input_files"] == 2 and result["input_records"] == 4
        assert result["output_records"] == 2
        assert output.exists() and Path(f"{output}.tbi").exists()
        text = gzip.open(output, "rt", encoding="utf-8").read()
        records = [line for line in text.splitlines() if line and not line.startswith("#")]
        assert len(records) == 2
        assert "END=7" in records[0] and "DP=10" in records[0]
        assert "\t8\t.\tA\tG\t" in records[1]
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["execution_space"] == "Host"
        assert metadata["determinism"] == "strict"
        assert metadata["compatibility"]["site_level_merge"] is True
        assert metadata["compatibility"]["allele_union"] is True
        assert metadata["compatibility"]["gq_rematerialization"] is True
        assert metadata["compatibility"]["reference_block_merge"] is True
        assert metadata["compatibility"]["kokkos_genotype_assignment"] is False
        assert metadata["compatibility"]["native_genotype_materialization"] is False
        assert metadata["telemetry"]["genotype_kernel_calls"] == 0

        # The shared GATK writer option must remove FORMAT/sample payload only
        # at the publication boundary; merging still consumed both samples.
        sites_only = work / "combined-sites-only.g.vcf.gz"
        sites_only_manifest = work / "combined-sites-only.manifest.json"
        sites_only_run = subprocess.run([
            str(binary), "-V", str(sample1), "-V", str(sample2),
            "-O", str(sites_only), "--sites-only-vcf-output=true",
            "--output-manifest", str(sites_only_manifest),
        ], text=True, capture_output=True, check=False, env=env)
        assert sites_only_run.returncode == 0, sites_only_run.stderr
        site_lines = [line for line in gzip.open(sites_only, "rt", encoding="utf-8")
                      if line and not line.startswith("#")]
        header_lines = [line for line in gzip.open(sites_only, "rt", encoding="utf-8")
                        if line.startswith("#CHROM")]
        assert header_lines and header_lines[0].rstrip("\n").split("\t") == [
            "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
        ]
        assert site_lines and all(len(line.rstrip("\n").split("\t")) == 8 for line in site_lines)
        sites_only_metadata = json.loads(sites_only_manifest.read_text(encoding="utf-8"))
        assert sites_only_metadata["compatibility"]["sites_only_vcf_output"] is True
        assert sites_only_metadata["telemetry"]["sites_only_vcf_output"] is True

        # Plain VCF/GVCF output publishes a Tribble `.idx` sidecar and is
        # directly consumable by a bundled GATK interval traversal.
        plain_output = work / "combined-plain.g.vcf"
        plain_manifest = work / "combined-plain.manifest.json"
        plain_result = subprocess.run([
            str(binary), "-V", str(shard1), "-V", str(shard2), "-O", str(plain_output),
            "--output-manifest", str(plain_manifest),
        ], text=True, capture_output=True, check=False, env=env)
        assert plain_result.returncode == 0, plain_result.stderr
        plain_summary = json.loads(plain_result.stdout.splitlines()[-1])
        assert plain_summary["output_records"] == result["output_records"]
        assert plain_output.exists() and Path(f"{plain_output}.idx").exists()
        assert not Path(f"{plain_output}.tbi").exists()
        plain_query = work / "combined-plain-query.vcf"
        plain_query_result = subprocess.run([
            str(java), "-jar", str(gatk_jar), "SelectVariants", "-V", str(plain_output),
            "-L", "chr1:8-8", "-O", str(plain_query),
        ], text=True, capture_output=True, check=False)
        assert plain_query_result.returncode == 0, plain_query_result.stderr
        assert sum(1 for line in plain_query.read_text(encoding="utf-8").splitlines()
                   if line and not line.startswith("#")) == 1
        plain_metadata = json.loads(plain_manifest.read_text(encoding="utf-8"))
        assert plain_metadata["compatibility"]["vcf_index"] is True
        # The bounded k-way merge keeps one decoded record per input and must
        # preserve the aggregate merge bytes exactly (including block coalescing).
        stream_output = work / "combined-stream.g.vcf.gz"
        stream_manifest = work / "combined-stream.g.vcf.gz.manifest.json"
        stream_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(shard1), "-V", str(shard2),
            "-O", str(stream_output), "--stream-merge",
            "--output-manifest", str(stream_manifest),
        ], text=True, env=env).splitlines()[-1])
        assert stream_result["stream_merge"] is True
        assert stream_result["input_records"] == 4 and stream_result["output_records"] == 2
        assert gzip.open(output, "rb").read() == gzip.open(stream_output, "rb").read()
        stream_metadata = json.loads(stream_manifest.read_text(encoding="utf-8"))
        assert stream_metadata["compatibility"]["stream_merge"] is True
        assert stream_metadata["compatibility"]["k_way_merge"] is True
        assert stream_metadata["telemetry"]["max_inflight_records"] <= 3
        assert stream_metadata["telemetry"]["peak_queued_bytes"] > 0

        # An indexed VCF input must switch the stream cursor to HTSlib's
        # interval iterator (rather than scanning the complete file).  The
        # aggregate and bounded stream paths remain byte-identical for the
        # selected interval and expose the traversal in the manifest.
        indexed_aggregate_output = work / "indexed-aggregate.g.vcf.gz"
        indexed_stream_output = work / "indexed-stream.g.vcf.gz"
        indexed_stream_manifest = work / "indexed-stream.manifest.json"
        subprocess.run([
            str(binary), "-V", str(output), "-L", "chr1:8-8",
            "-O", str(indexed_aggregate_output),
        ], check=True, env=env, stdout=subprocess.DEVNULL)
        indexed_stream_result = subprocess.run([
            str(binary), "-V", str(output), "-L", "chr1:8-8",
            "-O", str(indexed_stream_output), "--stream-merge",
            "--output-manifest", str(indexed_stream_manifest),
        ], text=True, capture_output=True, env=env)
        assert indexed_stream_result.returncode == 0, indexed_stream_result.stderr
        assert gzip.open(indexed_aggregate_output, "rb").read() == gzip.open(
            indexed_stream_output, "rb").read()
        indexed_summary = json.loads(indexed_stream_result.stdout.splitlines()[-1])
        assert indexed_summary["output_records"] == 1
        indexed_metadata = json.loads(indexed_stream_manifest.read_text(encoding="utf-8"))
        assert indexed_metadata["compatibility"]["indexed_interval_traversal"] is True
        assert indexed_metadata["telemetry"]["indexed_inputs"] == 1
        assert indexed_metadata["telemetry"]["indexed_interval_queries"] == 1
        assert indexed_metadata["telemetry"]["interval_skipped"] == 0

        # Repeated/overlapping selectors have union semantics in indexed
        # stream mode too; normalization must prevent the same BGZF records
        # from being emitted once per query.
        repeated_output = work / "indexed-repeated-stream.g.vcf.gz"
        repeated_manifest = work / "indexed-repeated-stream.manifest.json"
        repeated_result = subprocess.run([
            str(binary), "-V", str(output), "-L", "chr1:8-8",
            "-L", "chr1:8-8", "-O", str(repeated_output), "--stream-merge",
            "--output-manifest", str(repeated_manifest),
        ], text=True, capture_output=True, env=env)
        assert repeated_result.returncode == 0, repeated_result.stderr
        repeated_rows = [line for line in gzip.open(
            repeated_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(repeated_rows) == 1
        repeated_metadata = json.loads(repeated_manifest.read_text(encoding="utf-8"))
        assert repeated_metadata["telemetry"]["indexed_interval_queries"] == 1

        # GATK's INTERSECTION rule applies a set operation across repeated
        # selectors.  The overlap of 1..8 and 8..8 contains only the variant
        # at coordinate 8; both aggregate and bounded stream paths must agree.
        intersection_output = work / "indexed-intersection.g.vcf.gz"
        intersection_manifest = work / "indexed-intersection.manifest.json"
        intersection_result = subprocess.run([
            str(binary), "-V", str(output), "-L", "chr1:1-8",
            "-L", "chr1:8-8", "--interval-set-rule", "INTERSECTION",
            "-O", str(intersection_output), "--output-manifest", str(intersection_manifest),
        ], text=True, capture_output=True, env=env)
        assert intersection_result.returncode == 0, intersection_result.stderr
        intersection_summary = json.loads(intersection_result.stdout.splitlines()[-1])
        assert intersection_summary["interval_set_rule"] == "INTERSECTION"
        intersection_rows = [line for line in gzip.open(
            intersection_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(intersection_rows) == 1 and "\t8\t.\tA\tG\t" in intersection_rows[0]
        intersection_metadata = json.loads(intersection_manifest.read_text(encoding="utf-8"))
        assert intersection_metadata["compatibility"]["interval_set_rule"] == "INTERSECTION"
        assert intersection_metadata["telemetry"]["interval_set_rule"] == "INTERSECTION"
        intersection_stream_output = work / "indexed-intersection-stream.g.vcf.gz"
        intersection_stream_manifest = work / "indexed-intersection-stream.manifest.json"
        intersection_stream_result = subprocess.run([
            str(binary), "-V", str(output), "-L", "chr1:1-8",
            "-L", "chr1:8-8", "--interval-set-rule=INTERSECTION",
            "-O", str(intersection_stream_output), "--stream-merge",
            "--output-manifest", str(intersection_stream_manifest),
        ], text=True, capture_output=True, env=env)
        assert intersection_stream_result.returncode == 0, intersection_stream_result.stderr
        assert gzip.open(intersection_output, "rb").read() == gzip.open(
            intersection_stream_output, "rb").read()
        intersection_stream_summary = json.loads(
            intersection_stream_result.stdout.splitlines()[-1])
        assert intersection_stream_summary["interval_set_rule"] == "INTERSECTION"

        split_aggregate_output = work / "split-aggregate.g.vcf.gz"
        split_stream_output = work / "split-stream.g.vcf.gz"
        split_stream_manifest = work / "split-stream.manifest.json"
        split_aggregate = subprocess.run([
            str(binary), "-V", str(shard1), "-O", str(split_aggregate_output),
            "--convert-to-base-pair-resolution",
        ], text=True, capture_output=True, env=env)
        assert split_aggregate.returncode == 0, split_aggregate.stderr
        split_stream = subprocess.run([
            str(binary), "-V", str(shard1), "-O", str(split_stream_output),
            "--stream-merge", "--convert-to-base-pair-resolution",
            "--output-manifest", str(split_stream_manifest),
        ], text=True, capture_output=True, env=env)
        assert split_stream.returncode == 0, split_stream.stderr
        assert gzip.open(split_stream_output, "rb").read() == gzip.open(
            split_aggregate_output, "rb").read()
        split_stream_metadata = json.loads(split_stream_manifest.read_text(encoding="utf-8"))
        assert split_stream_metadata["compatibility"]["lazy_reference_band_split"] is True
        assert split_stream_metadata["telemetry"]["lazy_reference_blocks"] == 1
        assert split_stream_metadata["telemetry"]["lazy_reference_segments"] == 3
        sample_output = work / "samples.g.vcf.gz"
        sample_manifest = work / "samples.g.vcf.gz.manifest.json"
        sample_result = json.loads(subprocess.check_output([
            str(binary), "-V", str(sample1), "-V", str(sample2),
            "-O", str(sample_output), "--output-manifest", str(sample_manifest),
            "--fastgatk-materialize-genotypes",
        ], text=True, env=env).splitlines()[-1])
        assert sample_result["sample_count"] == 2
        sample_text = gzip.open(sample_output, "rt", encoding="utf-8").read()
        assert "\tFORMAT\tSAMPLE1\tSAMPLE2" in sample_text
        sample_records = [line.split("\t") for line in sample_text.splitlines()
                          if line and not line.startswith("#")]
        assert len(sample_records) == 1
        assert sample_records[0][9].startswith("0/1:") and sample_records[0][10].startswith("1/1:")
        sample_metadata = json.loads(sample_manifest.read_text(encoding="utf-8"))
        assert sample_metadata["compatibility"]["native_genotype_materialization"] is True
        assert sample_metadata["compatibility"]["gatk_format_semantics"] is False
        assert sample_metadata["compatibility"]["format_sample_merge"] is True
        assert sample_metadata["telemetry"]["genotype_kernel_calls"] > 0
        assert sample_metadata["telemetry"]["genotype_kernel_seconds"] > 0.0
        assert sample_metadata["telemetry"]["genotype_kernel_execution_space"]
        assert sample_metadata["telemetry"]["pl_remap_kernel_calls"] == 0
        sample_stream_output = work / "samples-stream.g.vcf.gz"
        sample_stream_manifest = work / "samples-stream.manifest.json"
        sample_stream_result = subprocess.run([
            str(binary), "-V", str(sample1), "-V", str(sample2),
            "-O", str(sample_stream_output), "--fastgatk-materialize-genotypes",
            "--stream-merge", "--output-manifest", str(sample_stream_manifest),
        ], text=True, capture_output=True, env=env)
        assert sample_stream_result.returncode == 0, sample_stream_result.stderr
        assert gzip.open(sample_output, "rb").read() == gzip.open(sample_stream_output, "rb").read()
        sample_stream_metadata = json.loads(sample_stream_manifest.read_text(encoding="utf-8"))
        assert sample_stream_metadata["compatibility"]["stream_merge"] is True

        # The standard GATK spelling must be a drop-in alias for the native
        # materialization switch.  Also verify deterministic reference-band
        # breakpoints: a block covering 1..12 split at multiples of five is
        # emitted as 1..4, 5..9, and 10..12 (1-based coordinates).
        alias_output = work / "call-genotypes-alias.g.vcf.gz"
        alias_manifest = work / "call-genotypes-alias.manifest.json"
        alias_result = subprocess.run([
            str(binary), "-V", sample1, "-O", alias_output,
            "--call-genotypes", "--break-bands-at-multiples-of", "5",
            "--output-manifest", alias_manifest,
        ], text=True, capture_output=True, env=env)
        assert alias_result.returncode == 0, alias_result.stderr
        alias_records = [line.split("\t") for line in gzip.open(
            alias_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(alias_records) == 1  # ref_sample1 is a one-locus variant
        alias_metadata = json.loads(alias_manifest.read_text(encoding="utf-8"))
        assert alias_metadata["compatibility"]["call_genotypes"] is True
        assert alias_metadata["telemetry"]["call_genotypes"] is True
        false_call_output = work / "call-genotypes-false.g.vcf.gz"
        false_call_manifest = work / "call-genotypes-false.manifest.json"
        false_call_result = subprocess.run([
            str(binary), "-V", sample1, "-O", false_call_output,
            "--call-genotypes", "false", "--output-manifest", false_call_manifest,
        ], text=True, capture_output=True, env=env)
        assert false_call_result.returncode == 0, false_call_result.stderr
        false_call_summary = json.loads(false_call_result.stdout.splitlines()[-1])
        assert false_call_summary["call_genotypes"] is False
        false_call_metadata = json.loads(false_call_manifest.read_text(encoding="utf-8"))
        assert false_call_metadata["compatibility"]["call_genotypes"] is False
        # The explicit false value must not be confused with the bare switch;
        # native experimental materialization remains disabled in this path.
        assert false_call_metadata["compatibility"]["native_genotype_materialization"] is False
        block_for_split = work / "split-block.g.vcf.gz"
        split_output = work / "split-block-out.g.vcf.gz"
        write_vcf(block_for_split,
                  "chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=12;DP=10\n")
        split_result = subprocess.run([
            str(binary), "-V", str(block_for_split), "-O", str(split_output),
            "--break-bands-at-multiples-of", "5",
        ], text=True, capture_output=True, env=env)
        assert split_result.returncode == 0, split_result.stderr
        split_records = [line.split("\t") for line in gzip.open(
            split_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [int(row[1]) for row in split_records] == [1, 5, 10]
        split_info = [dict(item.split("=", 1) for item in row[7].split(";") if "=" in item)
                      for row in split_records]
        assert [info.get("END") for info in split_info] == ["4", "9", "12"]
        bp_output = work / "bp-block-out.g.vcf.gz"
        bp_result = subprocess.run([
            str(binary), "-V", str(block_for_split), "-O", str(bp_output),
            "--convert-to-base-pair-resolution",
        ], text=True, capture_output=True, env=env)
        assert bp_result.returncode == 0, bp_result.stderr
        bp_records = [line for line in gzip.open(
            bp_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(bp_records) == 12
        assert all("END=" not in row.split("\t")[7] for row in bp_records)

        # CombineGVCFs now uses the shared arbitrary-ploidy Kokkos genotype
        # kernel, not a diploid-only Host triangular loop.
        triploid = work / "triploid.g.vcf.gz"
        triploid_output = work / "triploid-combined.g.vcf.gz"
        triploid_header = sample_header.replace("SAMPLE1", "TRIPLOID")
        write_vcf(triploid, "chr1\t8\t.\tA\tG\t.\tPASS\tDP=10\t"
                  "GT:DP:AD:PL\t0/0/0:10:5,5:50,0,80,100\n", triploid_header)
        subprocess.run([
            str(binary), "-V", str(triploid), "-O", str(triploid_output),
            "--fastgatk-materialize-genotypes",
        ], check=True, env=env, stdout=subprocess.DEVNULL)
        triploid_records = [line.split("\t") for line in gzip.open(
            triploid_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(triploid_records) == 1
        triploid_format = triploid_records[0][8].split(":")
        triploid_sample = triploid_records[0][9].split(":")
        assert triploid_sample[triploid_format.index("GT")] == "0/0/1"
        assert triploid_sample[triploid_format.index("GQ")] == "50"

        # Reference blocks from distinct samples must coalesce while retaining
        # each sample's FORMAT call, even when their INFO/FORMAT depths differ.
        ref_sample1 = work / "ref-sample1.g.vcf.gz"
        ref_sample2 = work / "ref-sample2.g.vcf.gz"
        ref_header = sample_header.replace("SAMPLE1", "REF1")
        write_vcf(ref_sample1,
                  "chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=3;DP=10\t"
                  "GT:DP:AD:PL\t0/0:10:10,0:0,99,99\n", ref_header)
        write_vcf(ref_sample2,
                  "chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=3;DP=12\t"
                  "GT:DP:AD:PL\t0/0:12:12,0:0,99,99\n",
                  ref_header.replace("REF1", "REF2"))
        ref_output = work / "ref-blocks.g.vcf.gz"
        subprocess.run([
            str(binary), "-V", str(ref_sample1), "-V", str(ref_sample2),
            "-O", str(ref_output), "--fastgatk-materialize-genotypes",
        ], check=True, env=env, stdout=subprocess.DEVNULL)
        ref_lines = [line.split("\t") for line in gzip.open(
            ref_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(ref_lines) == 1 and "END=3" in ref_lines[0][7]
        assert "\tFORMAT\tREF1\tREF2" in gzip.open(
            ref_output, "rt", encoding="utf-8").read()
        ref_format = ref_lines[0][8].split(":")
        ref_sample1 = ref_lines[0][9].split(":")
        ref_sample2 = ref_lines[0][10].split(":")
        assert ref_sample1[ref_format.index("GT")] == "0/0"
        assert ref_sample1[ref_format.index("DP")] == "10"
        assert ref_sample2[ref_format.index("GT")] == "0/0"
        assert ref_sample2[ref_format.index("DP")] == "12"

        # A one-base hole must not be bridged by reference-block merging.
        gap1 = work / "gap1.g.vcf.gz"
        gap2 = work / "gap2.g.vcf.gz"
        write_vcf(gap1, "chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=3;DP=10\n")
        write_vcf(gap2, "chr1\t5\t.\tA\t<NON_REF>\t.\tPASS\tEND=7;DP=10\n")
        gap_output = work / "gap.g.vcf.gz"
        subprocess.run([
            str(binary), "-V", str(gap1), "-V", str(gap2), "-O", str(gap_output),
        ], check=True, env=env, stdout=subprocess.DEVNULL)
        gap_lines = [line for line in gzip.open(
            gap_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(gap_lines) == 2

        # Repeated interval selectors use union semantics before site merge;
        # the position-1 block is fully contained by 1..3, while the
        # position-4 block is excluded and both position-8 shards are still
        # coalesced/deduplicated.
        interval_output = work / "interval.g.vcf.gz"
        interval_manifest = work / "interval.manifest.json"
        interval_result = subprocess.run([
            str(binary), "-V", str(shard1), "-V", str(shard2),
            "-L", "chr1:1-3", "--intervals", "chr1:8-8",
            "-O", str(interval_output), "--output-manifest", str(interval_manifest),
        ], text=True, capture_output=True, env=env)
        assert interval_result.returncode == 0, interval_result.stderr
        interval_summary = json.loads(interval_result.stdout.splitlines()[-1])
        assert interval_summary["interval_skipped"] == 1
        interval_lines = [line.split("\t") for line in gzip.open(
            interval_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert [line[1] for line in interval_lines] == ["1", "8"]
        interval_metadata = json.loads(interval_manifest.read_text(encoding="utf-8"))
        assert interval_metadata["compatibility"]["interval_subset"] is True
        assert interval_metadata["telemetry"]["intervals"] == 2
        interval_list = work / "combine.interval_list"
        interval_list.write_text("@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:100\n"
                                 "chr1\t1\t3\t+\tfirst\nchr1\t8\t8\t+\tsecond\n",
                                 encoding="utf-8")
        interval_file_output = work / "interval-file.g.vcf.gz"
        interval_file_manifest = work / "interval-file.manifest.json"
        interval_file_result = subprocess.run([
            str(binary), "-V", str(shard1), "-V", str(shard2), "-L", str(interval_list),
            "-O", str(interval_file_output), "--output-manifest", str(interval_file_manifest),
        ], text=True, capture_output=True, env=env)
        assert interval_file_result.returncode == 0, interval_file_result.stderr
        interval_file_metadata = json.loads(interval_file_manifest.read_text(encoding="utf-8"))
        assert interval_file_metadata["telemetry"]["interval_list_inputs"] == 1
        assert interval_file_metadata["telemetry"]["interval_list_records"] == 2

        # A reference block is retained when its END is contained by the
        # interval; Java drops a block when the interval cuts its END.
        block_input = work / "block.g.vcf.gz"
        block_output = work / "block-out.g.vcf.gz"
        write_vcf(block_input, "chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=5;DP=10\n")
        block_result = subprocess.run([
            str(binary), "-V", str(block_input), "-L", "chr1:4-5", "-O", str(block_output),
        ], text=True, capture_output=True, env=env)
        assert block_result.returncode == 0, block_result.stderr
        block_records = [line for line in gzip.open(
            block_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(block_records) == 1 and block_records[0].split("\t")[1] == "1"

        # Combining sample shards with disjoint concrete ALTs must produce a
        # deterministic ALT union and triangular PL layout.
        alt_g = work / "alt-g.g.vcf.gz"
        alt_t = work / "alt-t.g.vcf.gz"
        write_vcf(alt_g, "chr1\t9\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=10\tGT:DP:AD:PL\t0/1:10:5,5,0:50,0,50,99,99,99\n",
                  sample_header.replace("SAMPLE1", "ALT_G"))
        write_vcf(alt_t, "chr1\t9\t.\tA\tT,<NON_REF>\t.\tPASS\tDP=11\tGT:DP:AD:PL\t1/1:11:0,0,11:100,99,0,99,99,99\n",
                  sample_header.replace("SAMPLE1", "ALT_T"))
        union_output = work / "union.g.vcf.gz"
        union_manifest = work / "union.manifest.json"
        subprocess.run([
            str(binary), "-V", str(alt_g), "-V", str(alt_t), "-O", str(union_output),
            "--fastgatk-materialize-genotypes",
            "--output-manifest", str(union_manifest),
        ], check=True, env=env, stdout=subprocess.DEVNULL)
        union_records = [line.split("\t") for line in gzip.open(
            union_output, "rt", encoding="utf-8").read().splitlines()
            if line and not line.startswith("#")]
        assert len(union_records) == 1 and union_records[0][4] == "G,T,<NON_REF>"
        union_format = union_records[0][8].split(":")
        assert len(union_records[0][9].split(":")[union_format.index("AD")].split(",")) == 4
        assert len(union_records[0][9].split(":")[union_format.index("PL")].split(",")) == 10
        assert len(union_records[0][10].split(":")[union_format.index("AD")].split(",")) == 4
        assert len(union_records[0][10].split(":")[union_format.index("PL")].split(",")) == 10
        union_metadata = json.loads(union_manifest.read_text(encoding="utf-8"))
        assert union_metadata["telemetry"]["pl_remap_kernel_calls"] > 0
        assert union_metadata["telemetry"]["allele_field_remap_kernel_calls"] > 0
        assert union_metadata["telemetry"]["pl_remap_kernel_execution_space"]
        duplicate = subprocess.run([
            str(binary), "-V", str(sample1), "-V", str(sample1),
            "-O", str(work / "duplicate.g.vcf.gz")
        ], text=True, capture_output=True, env=env)
        assert duplicate.returncode != 0 and "duplicate" in duplicate.stderr
        print(json.dumps({"status": "pass", "records": len(records), "indexed": True}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
