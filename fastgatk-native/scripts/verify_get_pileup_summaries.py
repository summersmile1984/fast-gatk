#!/usr/bin/env python3
"""Verify the native GetPileupSummaries table and Kokkos count contract."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def run(binary: Path, args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run([str(binary), *args], text=True, capture_output=True, check=False)


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_GET_PILEUP_BINARY",
        str(root / "fastgatk-native/build/fastgatk-get-pileup-summaries"),
    ))
    assert binary.is_file() and os.access(binary, os.X_OK)
    with tempfile.TemporaryDirectory(prefix="fastgatk-get-pileup-summaries-") as directory:
        work = Path(directory)
        reads = work / "reads.sam"
        reads.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:100\n"
            "@RG\tID:rg1\tSM:TUMOR\n"
            "r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAA\tIIIIIIIIII\tRG:Z:rg1\n"
            "r2\t0\tchr1\t1\t60\t10M\t*\t0\t0\tAGGTTCGTAA\tIIIIIIIIII\tRG:Z:rg1\n"
            "lowmq\t0\tchr1\t1\t10\t10M\t*\t0\t0\tAGGTTCGTAA\tIIIIIIIIII\tRG:Z:rg1\n"
            "duplicate\t1024\tchr1\t1\t60\t10M\t*\t0\t0\tAGGTTCGTAA\tIIIIIIIIII\tRG:Z:rg1\n",
            encoding="utf-8",
        )
        variants = work / "sites.vcf"
        variants.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=100>\n"
            "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t2\trs2\tC\tG\t50\tPASS\tAF=0.10\n"
            "chr1\t5\trs5\tA\tT\t50\tPASS\tAF=0.15\n"
            "chr1\t8\trs8\tT\tC\t50\tPASS\tAF=0.50\n",
            encoding="utf-8",
        )
        output = work / "pileups.table"
        manifest = work / "pileups.manifest.json"
        result = run(binary, [
            "-I", str(reads), "-V", str(variants), "-L", str(variants),
            "-O", str(output), "--output-manifest", str(manifest), "--batch-records", "1",
        ])
        assert result.returncode == 0, result.stderr
        lines = output.read_text(encoding="utf-8").splitlines()
        assert lines[0] == "#<METADATA>SAMPLE=TUMOR"
        assert lines[1] == "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency"
        assert lines[2].startswith("chr1\t2\t1\t1\t0\t0.1")
        assert lines[3].startswith("chr1\t5\t1\t1\t0\t0.15")
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["status"] == "prototype"
        assert metadata["compatibility"]["kokkos_count_kernel"] is True
        assert metadata["telemetry"]["read_filters"] == []
        assert metadata["telemetry"]["reads_seen"] == 4
        assert metadata["telemetry"]["reads_used"] == 2
        assert metadata["telemetry"]["bases_used"] == 4
        assert metadata["telemetry"]["persistent_buffer_allocations"] >= 2
        assert metadata["telemetry"]["persistent_buffer_reuses"] > 0
        assert metadata["telemetry"]["persistent_buffer_capacity_records"] >= 1
        telemetry = metadata["telemetry"]
        assert telemetry["pipeline_lifecycle"] == "decode->compute->encode->sink"
        assert telemetry["pipeline_decoded_items"] == telemetry["pipeline_computed_items"]
        assert telemetry["pipeline_computed_items"] == telemetry["pipeline_encoded_items"]
        assert telemetry["pipeline_decoded_items"] >= 1
        assert telemetry["pipeline_peak_decoded_bytes"] > 0
        assert telemetry["pipeline_peak_computed_bytes"] > 0
        assert telemetry["pipeline_peak_encoded_bytes"] > 0
        expected_execution_space = os.environ.get("FASTGATK_EXPECTED_EXECUTION_SPACE")
        if expected_execution_space:
            assert metadata["telemetry"]["count_kernel_execution_space"] == expected_execution_space

        # GoodCigarReadFilter is applied per read.  A malformed CIGAR in one
        # batch must not make otherwise valid records disappear from the
        # pileup.  The batch size of two forces both records through one
        # decode/count boundary.
        malformed_reads = work / "malformed-cigar.sam"
        malformed_reads.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:100\n"
            "@RG\tID:rg1\tSM:TUMOR\n"
            "valid\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAA\tIIIIIIIIII\tRG:Z:rg1\n"
            # HTSlib accepts a zero-length CIGAR element and retains the
            # record, while GoodCigarReadFilter rejects it.  This gives the
            # end-to-end test a malformed record in the same decoded batch
            # without making sam_read1 stop at the record boundary.
            "bad-cigar\t0\tchr1\t1\t60\t0M\t*\t0\t0\t*\t*\tRG:Z:rg1\n",
            encoding="utf-8",
        )
        malformed_output = work / "malformed-cigar.table"
        malformed_manifest = work / "malformed-cigar.manifest.json"
        malformed_result = run(binary, [
            "-I", str(malformed_reads), "-V", str(variants), "-L", str(variants),
            "-O", str(malformed_output), "--batch-records", "2",
            "--output-manifest", str(malformed_manifest),
        ])
        assert malformed_result.returncode == 0, malformed_result.stderr
        malformed_lines = malformed_output.read_text(encoding="utf-8").splitlines()
        assert malformed_lines[2].startswith("chr1\t2\t1\t0\t0\t0.1"), \
            "a malformed CIGAR must not discard a valid read in the same batch"
        malformed_metadata = json.loads(malformed_manifest.read_text(encoding="utf-8"))
        assert malformed_metadata["telemetry"]["reads_seen"] == 2
        assert malformed_metadata["telemetry"]["reads_used"] == 1

        # Repeatable -L selectors support GATK UNION (default) or
        # INTERSECTION.  The intersection of chr1:1-5 and chr1:2-2 keeps only
        # the population site at POS=2 and records the rule in the manifest.
        intersection_output = work / "intersection.table"
        intersection_manifest = work / "intersection.manifest.json"
        intersection_result = run(binary, [
            "-I", str(reads), "-V", str(variants),
            "-L", "chr1:1-5", "-L", "chr1:2-2",
            "--interval-set-rule", "INTERSECTION",
            "-O", str(intersection_output), "--output-manifest", str(intersection_manifest),
        ])
        assert intersection_result.returncode == 0, intersection_result.stderr
        intersection_lines = intersection_output.read_text(encoding="utf-8").splitlines()
        assert len(intersection_lines) == 3 and intersection_lines[2].startswith("chr1\t2\t1\t1\t0\t0.1")
        intersection_metadata = json.loads(intersection_manifest.read_text(encoding="utf-8"))
        assert intersection_metadata["compatibility"]["interval_set_rule"] == "INTERSECTION"
        assert intersection_metadata["telemetry"]["interval_set_rule"] == "INTERSECTION"

        # Exercise a real tabix-indexed VCF resource.  The first contig with a
        # SNP in this checked-in GATK fixture is chr20:61098; the native path
        # must use the TBI sequence-name map (rather than assuming VCF header
        # rid values line up with tabix tids) and report one interval query.
        indexed_variants = root / (
            "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/"
            "walkers/fasta/FastaAlternateReferenceMaker/"
            "NA12878.WGS.b37.chr20.firstMB.vcf.gz"
        )
        assert indexed_variants.is_file()
        indexed_reads = work / "indexed-reads.sam"
        indexed_reads.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:20\tLN:63025520\n"
            "@RG\tID:rg1\tSM:TUMOR\n"
            "r1\t0\t20\t61090\t60\t20M\t*\t0\t0\tAAAAAAAACAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1\n"
            "r2\t0\t20\t61090\t60\t20M\t*\t0\t0\tAAAAAAAATAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII\tRG:Z:rg1\n",
            encoding="utf-8",
        )
        indexed_output = work / "indexed.table"
        indexed_manifest = work / "indexed.manifest.json"
        result = run(binary, [
            "-I", str(indexed_reads), "-V", str(indexed_variants),
            "-L", "20:61098-61098", "--max-af", "0.6",
            "-O", str(indexed_output), "--output-manifest", str(indexed_manifest),
        ])
        assert result.returncode == 0, result.stderr
        indexed_lines = indexed_output.read_text(encoding="utf-8").splitlines()
        assert indexed_lines[2].startswith("20\t61098\t1\t1\t0\t0.5")
        indexed_metadata = json.loads(indexed_manifest.read_text(encoding="utf-8"))
        assert indexed_metadata["telemetry"]["indexed_variant_traversal"] is True
        assert indexed_metadata["telemetry"]["variant_interval_queries"] == 1

        indexed_selector = work / "indexed-selector.vcf"
        indexed_selector.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=20,length=63025520>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "20\t61098\trs1\tC\tT\t50\tPASS\t.\n",
            encoding="utf-8",
        )
        selector_output = work / "indexed-selector.table"
        selector_manifest = work / "indexed-selector.manifest.json"
        result = run(binary, [
            "-I", str(indexed_reads), "-V", str(indexed_variants),
            "-L", str(indexed_selector), "--max-af", "0.6",
            "-O", str(selector_output), "--output-manifest", str(selector_manifest),
        ])
        assert result.returncode == 0, result.stderr
        selector_lines = selector_output.read_text(encoding="utf-8").splitlines()
        assert selector_lines[2].startswith("20\t61098\t1\t1\t0\t0.5")
        selector_metadata = json.loads(selector_manifest.read_text(encoding="utf-8"))
        assert selector_metadata["telemetry"]["indexed_variant_traversal"] is True
        assert selector_metadata["telemetry"]["variant_interval_queries"] == 1

        include_duplicates = work / "include-duplicates.table"
        result = run(binary, [
            "-I", str(reads), "-V", str(variants), "-L", str(variants),
            "-O", str(include_duplicates), "--include-duplicates",
        ])
        assert result.returncode == 0, result.stderr
        duplicate_lines = include_duplicates.read_text(encoding="utf-8").splitlines()
        assert duplicate_lines[2].startswith("chr1\t2\t1\t2\t0\t0.1")

        custom_filter = work / "custom-filter.table"
        result = run(binary, [
            "-I", str(reads), "-V", str(variants), "-L", str(variants),
            "-O", str(custom_filter), "--disable-tool-default-read-filters",
            "-RF", "MappingQualityReadFilter",
        ])
        assert result.returncode == 0, result.stderr
        custom_lines = custom_filter.read_text(encoding="utf-8").splitlines()
        assert custom_lines[2].startswith("chr1\t2\t1\t2\t0\t0.1")

        # The parameterized ReadLengthReadFilter follows the same typed
        # bounds used by HC/Mutect2/CountReads/FlagStat.  Use a separate
        # two-read fixture so the default oracle above remains unchanged.
        length_reads = work / "length-reads.sam"
        length_reads.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:100\n"
            "@RG\tID:rg1\tSM:TUMOR\n"
            "long\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAA\tIIIIIIIIII\tRG:Z:rg1\n"
            "short\t0\tchr1\t1\t60\t5M\t*\t0\t0\tAGGTA\tIIIII\tRG:Z:rg1\n",
            encoding="utf-8",
        )
        length_filter = work / "length-filter.table"
        length_manifest = work / "length-filter.manifest.json"
        result = run(binary, [
            "-I", str(length_reads), "-V", str(variants), "-L", str(variants),
            "-O", str(length_filter), "--output-manifest", str(length_manifest),
            "--disable-tool-default-read-filters", "--read-filter", "ReadLengthReadFilter",
            "--min-read-length", "10", "--max-read-length", "10",
        ])
        assert result.returncode == 0, result.stderr
        length_lines = length_filter.read_text(encoding="utf-8").splitlines()
        assert length_lines[2].startswith("chr1\t2\t1\t0\t0\t0.1")
        length_metadata = json.loads(length_manifest.read_text(encoding="utf-8"))
        assert length_metadata["telemetry"]["read_filter_require_read_length"] is True
        assert length_metadata["telemetry"]["read_filter_min_read_length"] == 10
        assert length_metadata["telemetry"]["read_filter_max_read_length"] == 10
        assert length_metadata["telemetry"]["reads_seen"] == 2
        assert length_metadata["telemetry"]["reads_used"] == 1

        disabled_filter = work / "disabled-filter.table"
        disabled_manifest = work / "disabled-filter.manifest.json"
        result = run(binary, [
            "-I", str(reads), "-V", str(variants), "-L", str(variants),
            "-O", str(disabled_filter), "--output-manifest", str(disabled_manifest),
            "-DF", "NotDuplicateReadFilter",
        ])
        assert result.returncode == 0, result.stderr
        disabled_lines = disabled_filter.read_text(encoding="utf-8").splitlines()
        assert disabled_lines[2].startswith("chr1\t2\t1\t2\t0\t0.1")
        disabled_metadata = json.loads(disabled_manifest.read_text(encoding="utf-8"))
        assert disabled_metadata["telemetry"]["disabled_read_filters"] == ["NotDuplicateReadFilter"]

        no_default_filters = work / "no-default-filters.table"
        result = run(binary, [
            "-I", str(reads), "-V", str(variants), "-L", str(variants),
            "-O", str(no_default_filters), "--disable-tool-default-read-filters",
        ])
        assert result.returncode == 0, result.stderr
        no_filter_lines = no_default_filters.read_text(encoding="utf-8").splitlines()
        assert no_filter_lines[2].startswith("chr1\t2\t1\t3\t0\t0.1")

        filtered = work / "filtered.table"
        interval = work / "sites.interval_list"
        interval.write_text("@HD\tVN:1.6\nchr1\t2\t2\t+\tchr1\n", encoding="utf-8")
        result = run(binary, [
            "-I", str(reads), "-V", str(variants), "-L", str(interval),
            "-O", str(filtered), "--min-af", "0.05", "--max-af", "0.2",
        ])
        assert result.returncode == 0, result.stderr
        assert len(filtered.read_text(encoding="utf-8").splitlines()) == 3

        unsupported = run(binary, [
            "-I", str(reads), "-V", str(variants), "-L", str(variants),
            "-O", str(work / "unsupported.table"), "--read-filter", "NoSuchReadFilter",
        ])
        assert unsupported.returncode != 0
        assert "UNSUPPORTED_PARAMETER" in unsupported.stderr

        no_interval = run(binary, [
            "-I", str(reads), "-V", str(variants), "-O", str(work / "no-interval.table"),
        ])
        assert no_interval.returncode != 0
        assert "-L/--intervals is required" in no_interval.stderr

        mismatched = work / "mismatched-dict.vcf"
        mismatched.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=101>\n"
            "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t2\trs2\tC\tG\t50\tPASS\tAF=0.10\n",
            encoding="utf-8",
        )
        mismatch = run(binary, [
            "-I", str(reads), "-V", str(mismatched), "-L", str(mismatched),
            "-O", str(work / "mismatch.table"),
        ])
        assert mismatch.returncode != 0
        assert "sequence dictionary length mismatch" in mismatch.stderr
        disabled_dictionary = run(binary, [
            "-I", str(reads), "-V", str(mismatched), "-L", str(mismatched),
            "-O", str(work / "mismatch-disabled.table"),
            "--disable-sequence-dictionary-validation",
        ])
        assert disabled_dictionary.returncode == 0, disabled_dictionary.stderr

        missing_af = work / "missing-af.vcf"
        missing_af.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=100>\n"
            "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>\n"
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t2\trs2\tC\tG\t50\tPASS\tDP=10\n",
            encoding="utf-8",
        )
        missing = run(binary, [
            "-I", str(reads), "-V", str(missing_af), "-L", str(missing_af),
            "-O", str(work / "missing-af.table"),
        ])
        assert missing.returncode != 0
        assert "no variants in population VCF had an AF field" in missing.stderr

        # GATK succeeds when AF-bearing records exist but all are outside the
        # requested population-AF range: the output contains metadata and the
        # six-column header with no data rows.
        out_of_range = work / "out-of-range.vcf"
        out_of_range.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=100>\n"
            "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t2\trs2\tC\tG\t50\tPASS\tAF=0.90\n",
            encoding="utf-8",
        )
        empty_output = work / "out-of-range.table"
        empty_manifest = work / "out-of-range.manifest.json"
        empty = run(binary, [
            "-I", str(reads), "-V", str(out_of_range), "-L", str(out_of_range),
            "-O", str(empty_output), "--output-manifest", str(empty_manifest),
        ])
        assert empty.returncode == 0, empty.stderr
        assert empty_output.read_text(encoding="utf-8") == (
            "#<METADATA>SAMPLE=TUMOR\n"
            "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n"
        )
        empty_metadata = json.loads(empty_manifest.read_text(encoding="utf-8"))
        assert empty_metadata["telemetry"]["input_sites"] == 1
        assert empty_metadata["telemetry"]["selected_sites"] == 0

        # GATK obtains SAMPLE metadata from the first read-group SM tag and
        # fails when the BAM has no sample at all; UNKNOWN must not be emitted
        # as a silently consumable contamination table.
        no_sample_reads = work / "no-sample.sam"
        no_sample_reads.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:100\n"
            "r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAA\tIIIIIIIIII\n",
            encoding="utf-8",
        )
        no_sample = run(binary, [
            "-I", str(no_sample_reads), "-V", str(variants), "-L", str(variants),
            "-O", str(work / "no-sample.table"),
        ])
        assert no_sample.returncode != 0
        assert "no sample" in no_sample.stderr.lower()

    print(json.dumps({"status": "pass", "sites": 2, "reads_used": 2, "bases_used": 4}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
