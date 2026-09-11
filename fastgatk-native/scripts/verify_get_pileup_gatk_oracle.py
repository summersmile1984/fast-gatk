#!/usr/bin/env python3
"""Compare native GetPileupSummaries with the pinned GATK table writer."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    binary = Path(os.environ.get(
        "FASTGATK_GET_PILEUP_BINARY",
        str(root / "fastgatk-native/build/fastgatk-get-pileup-summaries"),
    ))
    required = (java, gatk, binary)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_get_pileup_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK GetPileupSummaries oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-get-pileup-gatk-oracle-") as directory:
        work = Path(directory)
        sam = work / "reads.sam"
        sam.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:100\n"
            "@RG\tID:rg1\tSM:TUMOR\n"
            "r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAA\tIIIIIIIIII\tRG:Z:rg1\n"
            "r2\t0\tchr1\t1\t60\t10M\t*\t0\t0\tAGGTTCGTAA\tIIIIIIIIII\tRG:Z:rg1\n"
            "lowmq\t0\tchr1\t1\t10\t10M\t*\t0\t0\tAGGTTCGTAA\tIIIIIIIIII\tRG:Z:rg1\n"
            "duplicate\t1024\tchr1\t1\t60\t10M\t*\t0\t0\tAGGTTCGTAA\tIIIIIIIIII\tRG:Z:rg1\n",
            encoding="utf-8",
        )
        reads = work / "reads.bam"
        sites = work / "sites.vcf"
        sites.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=100>\n"
            "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t2\trs2\tC\tG\t50\tPASS\tAF=0.123456789\n"
            "chr1\t5\trs5\tA\tT\t50\tPASS\tAF=1.2e-2\n"
            "chr1\t8\trs8\tT\tC\t50\tPASS\tAF=0.10000000149\n",
            encoding="utf-8",
        )
        sort_run = run([
            str(java), "-jar", str(gatk), "SortSam", "-I", str(sam), "-O", str(reads),
            "-SO", "coordinate", "--CREATE_INDEX", "true",
        ])
        assert sort_run.returncode == 0, sort_run.stderr
        index_run = run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(sites)])
        assert index_run.returncode == 0, index_run.stderr

        gatk_output = work / "gatk.table"
        native_output = work / "native.table"
        gatk_run = run([
            str(java), "-jar", str(gatk), "GetPileupSummaries", "-I", str(reads),
            "-V", str(sites), "-L", str(sites), "-O", str(gatk_output),
        ])
        assert gatk_run.returncode == 0, gatk_run.stderr
        native_run = run([
            str(binary), "-I", str(reads), "-V", str(sites), "-L", str(sites),
            "-O", str(native_output),
        ])
        assert native_run.returncode == 0, native_run.stderr

        expected = gatk_output.read_text(encoding="utf-8")
        actual = native_output.read_text(encoding="utf-8")
        assert actual == expected, {"expected": expected, "actual": actual}

        # No AF-qualified sites is a successful GATK traversal, not a bad
        # input: both writers must preserve metadata/header and emit no rows.
        gatk_empty = work / "gatk-empty.table"
        native_empty = work / "native-empty.table"
        gatk_empty_run = run([
            str(java), "-jar", str(gatk), "GetPileupSummaries", "-I", str(reads),
            "-V", str(sites), "-L", str(sites), "--min-af", "0", "--max-af", "0.001",
            "-O", str(gatk_empty),
        ])
        assert gatk_empty_run.returncode == 0, gatk_empty_run.stderr
        native_empty_run = run([
            str(binary), "-I", str(reads), "-V", str(sites), "-L", str(sites),
            "--min-af", "0", "--max-af", "0.001", "-O", str(native_empty),
        ])
        assert native_empty_run.returncode == 0, native_empty_run.stderr
        assert native_empty.read_text(encoding="utf-8") == gatk_empty.read_text(encoding="utf-8")

        # GATK applies -XL after -L and honors interval padding on the
        # half-open traversal coordinates.  Keep both a normal exclusion and
        # the padded single-locus boundary pinned to the Java writer.
        gatk_excluded = work / "gatk-excluded.table"
        native_excluded = work / "native-excluded.table"
        exclusion_args = ["-L", "chr1:1-8", "-XL", "chr1:2-5"]
        gatk_exclusion_run = run([
            str(java), "-jar", str(gatk), "GetPileupSummaries", "-I", str(reads),
            "-V", str(sites), *exclusion_args, "-O", str(gatk_excluded),
        ])
        assert gatk_exclusion_run.returncode == 0, gatk_exclusion_run.stderr
        native_exclusion_run = run([
            str(binary), "-I", str(reads), "-V", str(sites), *exclusion_args,
            "-O", str(native_excluded),
        ])
        assert native_exclusion_run.returncode == 0, native_exclusion_run.stderr
        assert native_excluded.read_text(encoding="utf-8") == gatk_excluded.read_text(encoding="utf-8")

        gatk_padded = work / "gatk-padded-exclusion.table"
        native_padded = work / "native-padded-exclusion.table"
        padded_args = ["-L", "chr1:1-8", "-XL", "chr1:7-7", "-ixp", "1"]
        gatk_padded_run = run([
            str(java), "-jar", str(gatk), "GetPileupSummaries", "-I", str(reads),
            "-V", str(sites), *padded_args, "-O", str(gatk_padded),
        ])
        assert gatk_padded_run.returncode == 0, gatk_padded_run.stderr
        native_padded_run = run([
            str(binary), "-I", str(reads), "-V", str(sites), *padded_args,
            "-O", str(native_padded),
        ])
        assert native_padded_run.returncode == 0, native_padded_run.stderr
        assert native_padded.read_text(encoding="utf-8") == gatk_padded.read_text(encoding="utf-8")

        # Boolean arguments in Nextflow command lines are commonly emitted as
        # a separate token; accept that form exactly like GATK's Barclay
        # parser, not only the historical bare-flag spelling.
        gatk_no_filters = work / "gatk-no-filters.table"
        native_no_filters = work / "native-no-filters.table"
        no_filter_args = ["--disable-tool-default-read-filters", "true"]
        gatk_no_filter_run = run([
            str(java), "-jar", str(gatk), "GetPileupSummaries", "-I", str(reads),
            "-V", str(sites), "-L", str(sites), *no_filter_args, "-O", str(gatk_no_filters),
        ])
        assert gatk_no_filter_run.returncode == 0, gatk_no_filter_run.stderr
        native_no_filter_run = run([
            str(binary), "-I", str(reads), "-V", str(sites), "-L", str(sites),
            *no_filter_args, "-O", str(native_no_filters),
        ])
        assert native_no_filter_run.returncode == 0, native_no_filter_run.stderr
        assert native_no_filters.read_text(encoding="utf-8") == gatk_no_filters.read_text(encoding="utf-8")

        # GATK accepts the short -XRF spelling.  Start from its empty filter
        # set, then retain only the low-MAPQ read through the inverse.  (GATK
        # correctly rejects an attempt to both -DF and -XRF the same named
        # default filter.)  This is a visible table-writer boundary, not
        # merely parser acceptance.
        gatk_inverted = work / "gatk-inverted-mapq.table"
        native_inverted = work / "native-inverted-mapq.table"
        native_inverted_manifest = work / "native-inverted-mapq.manifest.json"
        inverted_args = ["--disable-tool-default-read-filters", "true",
                         "-XRF", "MappingQualityReadFilter", "--minimum-mapping-quality", "50"]
        gatk_inverted_run = run([
            str(java), "-jar", str(gatk), "GetPileupSummaries", "-I", str(reads),
            "-V", str(sites), "-L", str(sites), *inverted_args, "-O", str(gatk_inverted),
        ])
        assert gatk_inverted_run.returncode == 0, gatk_inverted_run.stderr
        native_inverted_run = run([
            str(binary), "-I", str(reads), "-V", str(sites), "-L", str(sites),
            *inverted_args, "-O", str(native_inverted), "--output-manifest", str(native_inverted_manifest),
        ])
        assert native_inverted_run.returncode == 0, native_inverted_run.stderr
        assert native_inverted.read_text(encoding="utf-8") == gatk_inverted.read_text(encoding="utf-8")
        inverted_manifest = json.loads(native_inverted_manifest.read_text(encoding="utf-8"))
        assert inverted_manifest["compatibility"]["inverted_read_filter"] is True
        assert inverted_manifest["telemetry"]["inverted_read_filter_count"] == 1

        # GATK obtains SAMPLE from the first read-group SM tag.  A header with
        # no sample must fail in both implementations instead of producing a
        # silently consumable SAMPLE=UNKNOWN table.
        no_sample_sam = work / "no-sample.sam"
        no_sample_sam.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:100\n"
            "r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAA\tIIIIIIIIII\n",
            encoding="utf-8",
        )
        no_sample_reads = work / "no-sample.bam"
        no_sample_sort = run([
            str(java), "-jar", str(gatk), "SortSam", "-I", str(no_sample_sam),
            "-O", str(no_sample_reads), "-SO", "coordinate", "--CREATE_INDEX", "true",
        ])
        assert no_sample_sort.returncode == 0, no_sample_sort.stderr
        gatk_no_sample = run([
            str(java), "-jar", str(gatk), "GetPileupSummaries", "-I", str(no_sample_reads),
            "-V", str(sites), "-L", str(sites), "-O", str(work / "gatk-no-sample.table"),
        ])
        native_no_sample = run([
            str(binary), "-I", str(no_sample_reads), "-V", str(sites), "-L", str(sites),
            "-O", str(work / "native-no-sample.table"),
        ])
        assert gatk_no_sample.returncode != 0, gatk_no_sample.stdout
        assert native_no_sample.returncode != 0, native_no_sample.stdout
        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "records": 3,
            "table_bytes_exact": True,
            "empty_site_table_exact": True,
            "exclude_intervals_exact": True,
            "exclude_interval_padding_exact": True,
            "boolean_read_filter_argument_exact": True,
            "inverted_read_filter_short_alias_exact": True,
            "read_filter_counts_exact": True,
            "af_format_exact": True,
            "no_sample_fail_closed": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
