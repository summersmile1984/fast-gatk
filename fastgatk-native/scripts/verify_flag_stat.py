#!/usr/bin/env python3
"""Contract/oracle test for the native FlagStat path."""

from __future__ import annotations

import json
import gzip
import os
import re
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BAM = ROOT / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
BINARY = Path(os.environ.get("FASTGATK_FLAG_STAT_BINARY", ROOT / "fastgatk-native/build/fastgatk-flag-stat"))
EXPECTED = [
    "493 in total",
    "53 QC failure",
    "59 duplicates",
    "487 mapped (98.78%)",
    "493 paired in sequencing",
    "247 read1",
    "246 read2",
    "477 properly paired (96.75%)",
    "477 with itself and mate mapped",
    "10 singletons (2.03%)",
    "0 with mate mapped to a different chr",
    "0 with mate mapped to a different chr (mapQ>=5)",
]


def main() -> int:
    if not BINARY.is_file():
        print(json.dumps({"status": "skip", "suite": "flag-stat", "reason": f"missing binary: {BINARY}"}))
        return 0
    with tempfile.TemporaryDirectory(prefix="fastgatk-flag-stat-") as directory:
        work = Path(directory)
        output = work / "flagstat.txt"
        manifest = work / "flagstat.manifest.json"
        result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-O", str(output), "--batch-records", "41",
             "--threads", "2", "--output-manifest", str(manifest)],
            text=True, capture_output=True,
        )
        assert result.returncode == 0, (result.stdout, result.stderr)
        actual = result.stdout.splitlines()
        assert actual == EXPECTED, (actual, EXPECTED)
        assert output.read_text().splitlines() == EXPECTED
        payload = json.loads(manifest.read_text())
        assert payload["tool"] == "FlagStat"
        # Keep the executable's manifest/telemetry status in lockstep with
        # the dispatcher registry's contract-compatible promotion.
        assert payload["status"] == "contract-compatible"
        assert payload["input"] == str(BAM)
        assert payload["input_count"] == 1
        assert payload["inputs"] == [str(BAM)]
        assert payload["records"] == 493
        assert payload["default_read_filter_count"] == 1
        assert payload["default_read_filters"] == ["WellformedReadFilter"]
        assert payload["telemetry"]["persistent_buffer_allocations"] == 1
        assert payload["telemetry"]["persistent_buffer_reuses"] > 0
        assert payload["telemetry"]["persistent_buffer_capacity_records"] >= 41
        assert payload["telemetry"]["kernel_execution_policy"] == "RangePolicy+atomic"
        assert payload["telemetry"]["pipeline_lifecycle"] == "decode->compute->encode->sink"
        assert payload["telemetry"]["pipeline_decoded_items"] == payload["batches"]
        assert payload["telemetry"]["pipeline_computed_items"] == payload["batches"]
        assert payload["telemetry"]["pipeline_encoded_items"] == payload["batches"]
        assert payload["telemetry"]["pipeline_peak_decoded_bytes"] > 0
        assert payload["telemetry"]["pipeline_peak_computed_bytes"] > 0
        assert payload["telemetry"]["pipeline_peak_encoded_bytes"] > 0
        assert json.loads(result.stderr.splitlines()[-1])["status"] == "contract-compatible"

        # Repeated -I is the shard aggregation form used by GATK and by the
        # surrounding SLURM/Nextflow workflow.  Every FlagStat counter and
        # the percentages must aggregate over both inputs.
        multi_manifest = work / "multi.manifest.json"
        multi_result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "--input", str(BAM),
             "--batch-records", "41", "--output-manifest", str(multi_manifest)],
            text=True, capture_output=True,
        )
        assert multi_result.returncode == 0, (multi_result.stdout, multi_result.stderr)
        assert multi_result.stdout.splitlines() == [
            "986 in total", "106 QC failure", "118 duplicates",
            "974 mapped (98.78%)", "986 paired in sequencing", "494 read1",
            "492 read2", "954 properly paired (96.75%)",
            "954 with itself and mate mapped", "20 singletons (2.03%)",
            "0 with mate mapped to a different chr",
            "0 with mate mapped to a different chr (mapQ>=5)",
        ]
        multi_payload = json.loads(multi_manifest.read_text())
        assert multi_payload["records"] == 986
        assert multi_payload["input_count"] == 2
        assert multi_payload["inputs"] == [str(BAM), str(BAM)]
        assert multi_payload["telemetry"]["pipeline_decoded_items"] == \
            multi_payload["telemetry"]["pipeline_computed_items"]
        assert multi_payload["telemetry"]["pipeline_computed_items"] == \
            multi_payload["telemetry"]["pipeline_encoded_items"]

        interval_result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", "17:69000-69500", "--batch-records", "29"],
            text=True, capture_output=True,
        )
        assert interval_result.returncode == 0, (interval_result.stdout, interval_result.stderr)
        assert interval_result.stdout.splitlines() == [
            "338 in total", "34 QC failure", "44 duplicates", "334 mapped (98.82%)",
            "338 paired in sequencing", "168 read1", "170 read2", "328 properly paired (97.04%)",
            "328 with itself and mate mapped", "6 singletons (1.78%)",
            "0 with mate mapped to a different chr", "0 with mate mapped to a different chr (mapQ>=5)",
        ]

        intersection_result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", "17:69000-69500", "-L", "17:69200-69300",
             "-isr", "INTERSECTION", "--output-manifest", str(work / "intersection.manifest.json")],
            text=True, capture_output=True,
        )
        assert intersection_result.returncode == 0 and intersection_result.stdout.splitlines() == [
            "11 in total", "0 QC failure", "1 duplicates", "11 mapped (100.00%)",
            "11 paired in sequencing", "7 read1", "4 read2", "11 properly paired (100.00%)",
            "11 with itself and mate mapped", "0 singletons (0.00%)",
            "0 with mate mapped to a different chr", "0 with mate mapped to a different chr (mapQ>=5)",
        ], (intersection_result.stdout, intersection_result.stderr)
        intersection_payload = json.loads((work / "intersection.manifest.json").read_text())
        assert intersection_payload["interval_set_rule"] == "INTERSECTION"

        invalid_rule = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", "17:69000-69500",
             "--interval-set-rule", "XOR"],
            text=True, capture_output=True,
        )
        assert invalid_rule.returncode != 0
        assert "interval-set-rule must be UNION or INTERSECTION" in invalid_rule.stderr

        compressed_interval = work / "regions.interval_list.gz"
        with gzip.open(compressed_interval, "wt", encoding="utf-8") as stream:
            stream.write("17:69000-69500\n")
        compressed_result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", str(compressed_interval), "--batch-records", "29"],
            text=True, capture_output=True,
        )
        assert compressed_result.returncode == 0, (compressed_result.stdout, compressed_result.stderr)
        assert compressed_result.stdout.splitlines() == interval_result.stdout.splitlines()

        compressed_bed = work / "regions.bed.gz"
        with gzip.open(compressed_bed, "wt", encoding="utf-8") as stream:
            stream.write("track name=targets\n17\t68999\t69500\n")
        compressed_bed_result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", str(compressed_bed), "--batch-records", "29"],
            text=True, capture_output=True,
        )
        assert compressed_bed_result.returncode == 0, (
            compressed_bed_result.stdout, compressed_bed_result.stderr)
        assert compressed_bed_result.stdout.splitlines() == interval_result.stdout.splitlines()

        malformed_compressed = work / "malformed.interval_list.gz"
        with gzip.open(malformed_compressed, "wt", encoding="utf-8") as stream:
            stream.write("17\t69000\n")
        malformed_result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", str(malformed_compressed)],
            text=True, capture_output=True,
        )
        assert malformed_result.returncode != 0 and "malformed interval file" in malformed_result.stderr

        filtered_result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-RF", "NotDuplicateReadFilter", "--batch-records", "37"],
            text=True, capture_output=True,
        )
        assert filtered_result.returncode == 0, (filtered_result.stdout, filtered_result.stderr)
        assert filtered_result.stdout.splitlines() == [
            "434 in total", "49 QC failure", "0 duplicates", "428 mapped (98.62%)",
            "434 paired in sequencing", "217 read1", "217 read2", "419 properly paired (96.54%)",
            "419 with itself and mate mapped", "9 singletons (2.07%)",
            "0 with mate mapped to a different chr", "0 with mate mapped to a different chr (mapQ>=5)",
        ]

        mapping_quality = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-RF", "MappingQualityReadFilter",
             "--minimum-mapping-quality", "30", "--maximum-mapping-quality", "50",
             "--batch-records", "37"], text=True, capture_output=True,
        )
        assert mapping_quality.returncode == 0, (mapping_quality.stdout, mapping_quality.stderr)
        assert mapping_quality.stdout.splitlines() == [
            "18 in total", "2 QC failure", "2 duplicates", "18 mapped (100.00%)",
            "18 paired in sequencing", "11 read1", "7 read2", "9 properly paired (50.00%)",
            "9 with itself and mate mapped", "9 singletons (50.00%)",
            "0 with mate mapped to a different chr",
            "0 with mate mapped to a different chr (mapQ>=5)",
        ]
        # -XRF is GATK's InvertedReadFilter surface.  The fixture contains
        # reads on both sides of this parameterized MAPQ predicate.
        inverted_mapping_quality = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-XRF", "MappingQualityReadFilter",
             "--minimum-mapping-quality", "30", "--maximum-mapping-quality", "50",
             "--batch-records", "37"], text=True, capture_output=True,
        )
        assert inverted_mapping_quality.returncode == 0, (
            inverted_mapping_quality.stdout, inverted_mapping_quality.stderr)
        assert inverted_mapping_quality.stdout.splitlines() == [
            "475 in total", "51 QC failure", "57 duplicates", "469 mapped (98.74%)",
            "475 paired in sequencing", "236 read1", "239 read2",
            "468 properly paired (98.53%)", "468 with itself and mate mapped",
            "1 singletons (0.21%)", "0 with mate mapped to a different chr",
            "0 with mate mapped to a different chr (mapQ>=5)",
        ]

        padded = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", "17:69990-69990", "-ip", "10",
             "--output-manifest", str(work / "padded.manifest.json")],
            text=True, capture_output=True,
        )
        assert padded.returncode == 0 and padded.stdout.splitlines() == [
            "17 in total", "1 QC failure", "0 duplicates", "17 mapped (100.00%)",
            "17 paired in sequencing", "9 read1", "8 read2", "17 properly paired (100.00%)",
            "17 with itself and mate mapped", "0 singletons (0.00%)",
            "0 with mate mapped to a different chr", "0 with mate mapped to a different chr (mapQ>=5)",
        ], (padded.stdout, padded.stderr)
        padded_payload = json.loads((work / "padded.manifest.json").read_text())
        assert padded_payload["interval_padding"] == 10
        assert padded_payload["interval_exclusion_padding"] == 0

        excluded_padding = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", "17:69980-70000",
             "-XL", "17:69990-69990", "-ixp", "2"],
            text=True, capture_output=True,
        )
        assert excluded_padding.returncode == 0 and excluded_padding.stdout.splitlines() == [
            "17 in total", "1 QC failure", "0 duplicates", "17 mapped (100.00%)",
            "17 paired in sequencing", "9 read1", "8 read2", "17 properly paired (100.00%)",
            "17 with itself and mate mapped", "0 singletons (0.00%)",
            "0 with mate mapped to a different chr", "0 with mate mapped to a different chr (mapQ>=5)",
        ], (excluded_padding.stdout, excluded_padding.stderr)
        exclude_only = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-XL", "17:69500-69500", "-ixp", "50"],
            text=True, capture_output=True,
        )
        assert exclude_only.returncode == 0 and exclude_only.stdout.splitlines() == [
            "415 in total", "42 QC failure", "48 duplicates", "411 mapped (99.04%)",
            "415 paired in sequencing", "214 read1", "201 read2", "404 properly paired (97.35%)",
            "404 with itself and mate mapped", "7 singletons (1.69%)",
            "0 with mate mapped to a different chr", "0 with mate mapped to a different chr (mapQ>=5)",
        ], (exclude_only.stdout, exclude_only.stderr)
        merged = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", "17:69990-69990",
             "-L", "17:69991-69991", "--interval-merging-rule=OVERLAPPING_ONLY",
             "--output-manifest", str(work / "overlap-only.manifest.json")],
            text=True, capture_output=True,
        )
        assert merged.returncode == 0 and merged.stdout.splitlines() == [
            "11 in total", "0 QC failure", "0 duplicates", "11 mapped (100.00%)",
            "11 paired in sequencing", "7 read1", "4 read2", "11 properly paired (100.00%)",
            "11 with itself and mate mapped", "0 singletons (0.00%)",
            "0 with mate mapped to a different chr", "0 with mate mapped to a different chr (mapQ>=5)",
        ], (merged.stdout, merged.stderr)
        merged_payload = json.loads((work / "overlap-only.manifest.json").read_text())
        assert merged_payload["interval_merging_rule"] == "OVERLAPPING_ONLY"
        assert merged_payload["intervals"] == 2
        for option, message in (("-ip", "interval-padding must be non-negative"),
                                ("-ixp", "interval-exclusion-padding must be non-negative")):
            invalid_padding = subprocess.run(
                [str(BINARY), "-I", str(BAM), option, "-1"],
                text=True, capture_output=True,
            )
            assert invalid_padding.returncode != 0 and message in invalid_padding.stderr, (
                option, invalid_padding.stdout, invalid_padding.stderr)
        invalid_merge = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-imr", "TOUCHING"],
            text=True, capture_output=True,
        )
        assert invalid_merge.returncode != 0 and "interval-merging-rule must be ALL or OVERLAPPING_ONLY" in invalid_merge.stderr
        all_excluded = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", "17:69000-69005",
             "-XL", "17:69000-69005"], text=True, capture_output=True,
        )
        assert all_excluded.returncode != 0 and "removed all territory specified by -L" in all_excluded.stderr, (
            all_excluded.stdout, all_excluded.stderr)
        mapping_quality_without_filter = subprocess.run(
            [str(BINARY), "-I", str(BAM), "--maximum-mapping-quality", "50"],
            text=True, capture_output=True,
        )
        assert mapping_quality_without_filter.returncode != 0
        assert "requires MappingQualityReadFilter" in mapping_quality_without_filter.stderr
        read_name = "809R9ABXX101220:5:21:7638:11612"
        read_name_result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-RF", "ReadNameReadFilter",
             "--read-name", read_name,
             "--output-manifest", str(work / "read-name.manifest.json")],
            text=True, capture_output=True,
        )
        assert read_name_result.returncode == 0 and read_name_result.stdout.splitlines() == [
            "2 in total", "0 QC failure", "0 duplicates", "2 mapped (100.00%)",
            "2 paired in sequencing", "1 read1", "1 read2", "2 properly paired (100.00%)",
            "2 with itself and mate mapped", "0 singletons (0.00%)",
            "0 with mate mapped to a different chr", "0 with mate mapped to a different chr (mapQ>=5)",
        ], (read_name_result.stdout, read_name_result.stderr)
        read_name_payload = json.loads((work / "read-name.manifest.json").read_text())
        assert read_name_payload["read_name_count"] == 1
        read_name_without_filter = subprocess.run(
            [str(BINARY), "-I", str(BAM), "--read-name", read_name],
            text=True, capture_output=True,
        )
        assert read_name_without_filter.returncode != 0 and "requires ReadNameReadFilter" in read_name_without_filter.stderr
        read_group_result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-RF", "ReadGroupReadFilter",
             "--keep-read-group", "809R9.5"], text=True, capture_output=True,
        )
        assert read_group_result.returncode == 0 and read_group_result.stdout.splitlines() == EXPECTED, (
            read_group_result.stdout, read_group_result.stderr)
        read_group_without_filter = subprocess.run(
            [str(BINARY), "-I", str(BAM), "--keep-read-group", "809R9.5"],
            text=True, capture_output=True,
        )
        assert read_group_without_filter.returncode != 0 and "requires ReadGroupReadFilter" in read_group_without_filter.stderr
        read_tag_result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-RF", "ReadTagValueFilter",
             "--read-filter-tag", "NM", "--read-filter-tag-comp", "1",
             "--read-filter-tag-op", "GREATER_OR_EQUAL",
             "--output-manifest", str(work / "read-tag.manifest.json")],
            text=True, capture_output=True,
        )
        assert read_tag_result.returncode == 0 and read_tag_result.stdout.splitlines() == [
            "189 in total", "27 QC failure", "24 duplicates", "189 mapped (100.00%)",
            "189 paired in sequencing", "84 read1", "105 read2", "185 properly paired (97.88%)",
            "185 with itself and mate mapped", "4 singletons (2.12%)",
            "0 with mate mapped to a different chr", "0 with mate mapped to a different chr (mapQ>=5)",
        ], (read_tag_result.stdout, read_tag_result.stderr)
        read_tag_payload = json.loads((work / "read-tag.manifest.json").read_text())
        assert read_tag_payload["read_filter_tag"] == "NM"
        assert read_tag_payload["read_filter_tag_comp"] == 1
        assert read_tag_payload["read_filter_tag_op"] == "GREATER_OR_EQUAL"
        for args, message in ((["--read-filter-tag", "NM"], "requires ReadTagValueFilter"),
                              (["-RF", "ReadTagValueFilter"], "requires a two-character"),
                              (["-RF", "ReadTagValueFilter", "--read-filter-tag", "N", "--read-filter-tag-comp", "1"], "requires a two-character"),
                              (["-RF", "ReadTagValueFilter", "--read-filter-tag", "NM", "--read-filter-tag-op", "BOGUS"], "read-filter-tag-op must be")):
            invalid_tag = subprocess.run([str(BINARY), "-I", str(BAM), *args], text=True, capture_output=True)
            assert invalid_tag.returncode != 0 and message in invalid_tag.stderr, (
                args, invalid_tag.stdout, invalid_tag.stderr)
        read_group_blacklist = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-RF", "ReadGroupBlackListReadFilter",
             "--read-group-black-list", "RG:809R9.5"], text=True, capture_output=True,
        )
        assert read_group_blacklist.returncode == 0 and read_group_blacklist.stdout.splitlines() == [
            "0 in total", "0 QC failure", "0 duplicates", "0 mapped (NaN%)",
            "0 paired in sequencing", "0 read1", "0 read2", "0 properly paired (NaN%)",
            "0 with itself and mate mapped", "0 singletons (NaN%)",
            "0 with mate mapped to a different chr", "0 with mate mapped to a different chr (mapQ>=5)",
        ], (read_group_blacklist.stdout, read_group_blacklist.stderr)
        blacklist_without_filter = subprocess.run(
            [str(BINARY), "-I", str(BAM), "--read-group-black-list", "RG:809R9.5"],
            text=True, capture_output=True,
        )
        assert blacklist_without_filter.returncode != 0 and "requires ReadGroupBlackListReadFilter" in blacklist_without_filter.stderr
        nonzero_fragment = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-RF", "NonZeroFragmentLengthReadFilter",
             "--batch-records", "37"], text=True, capture_output=True,
        )
        assert nonzero_fragment.returncode == 0
        assert nonzero_fragment.stdout.splitlines() == [
            "477 in total", "51 QC failure", "58 duplicates", "477 mapped (100.00%)",
            "477 paired in sequencing", "238 read1", "239 read2", "477 properly paired (100.00%)",
            "477 with itself and mate mapped", "0 singletons (0.00%)",
            "0 with mate mapped to a different chr",
            "0 with mate mapped to a different chr (mapQ>=5)",
        ], (nonzero_fragment.stdout, nonzero_fragment.stderr)
        mate_same = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-RF", "MateOnSameContigOrNoMappedMateReadFilter"],
            text=True, capture_output=True,
        )
        assert mate_same.returncode == 0 and mate_same.stdout.splitlines()[0] == "487 in total", (
            mate_same.stdout, mate_same.stderr)
        fragment_filter = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-RF", "FragmentLengthReadFilter",
             "--min-fragment-length", "100", "--max-fragment-length", "300"],
            text=True, capture_output=True,
        )
        assert fragment_filter.returncode == 0
        assert fragment_filter.stdout.splitlines() == [
            "413 in total", "45 QC failure", "56 duplicates", "413 mapped (100.00%)",
            "413 paired in sequencing", "206 read1", "207 read2", "413 properly paired (100.00%)",
            "413 with itself and mate mapped", "0 singletons (0.00%)",
            "0 with mate mapped to a different chr",
            "0 with mate mapped to a different chr (mapQ>=5)",
        ], (fragment_filter.stdout, fragment_filter.stderr)
        fragment_without_filter = subprocess.run(
            [str(BINARY), "-I", str(BAM), "--min-fragment-length", "100"],
            text=True, capture_output=True,
        )
        assert fragment_without_filter.returncode != 0
        assert "requires FragmentLengthReadFilter" in fragment_without_filter.stderr

        length_filtered = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-RF", "ReadLengthReadFilter",
             "--min-read-length", "80", "--max-read-length", "100", "--batch-records", "37"],
            text=True, capture_output=True,
        )
        assert length_filtered.returncode == 0, (length_filtered.stdout, length_filtered.stderr)
        assert length_filtered.stdout.splitlines() == [
            "0 in total", "0 QC failure", "0 duplicates", "0 mapped (NaN%)",
            "0 paired in sequencing", "0 read1", "0 read2", "0 properly paired (NaN%)",
            "0 with itself and mate mapped", "0 singletons (NaN%)",
            "0 with mate mapped to a different chr", "0 with mate mapped to a different chr (mapQ>=5)",
        ]
        no_defaults = subprocess.run(
            [str(BINARY), "-I", str(BAM), "--disable-tool-default-read-filters",
             "--output-manifest", str(work / "no-default-filters.manifest.json")],
            text=True, capture_output=True,
        )
        assert no_defaults.returncode == 0 and no_defaults.stdout.splitlines() == EXPECTED, (
            no_defaults.stdout, no_defaults.stderr)
        no_defaults_payload = json.loads(
            (work / "no-default-filters.manifest.json").read_text())
        assert no_defaults_payload["disable_tool_default_read_filters"] is True
        assert no_defaults_payload["default_read_filter_count"] == 0
        assert no_defaults_payload["default_read_filters"] == []

        malformed = work / "wellformed-boundary.sam"
        malformed.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:100\n"
            "@RG\tID:rg1\tSM:TUMOR\n"
            "good\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n"
            "no_rg\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
            "skip\t0\tchr1\t20\t60\t2M2N2M\t*\t0\t0\tACGT\tIIII\tRG:Z:rg1\n"
            "empty\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n",
            encoding="utf-8",
        )
        wellformed = subprocess.run(
            [str(BINARY), "-I", str(malformed)], text=True, capture_output=True,
        )
        assert wellformed.returncode == 0 and wellformed.stdout.splitlines() == [
            "1 in total", "0 QC failure", "0 duplicates", "1 mapped (100.00%)",
            "0 paired in sequencing", "0 read1", "0 read2", "0 properly paired (0.00%)",
            "0 with itself and mate mapped", "0 singletons (0.00%)",
            "0 with mate mapped to a different chr", "0 with mate mapped to a different chr (mapQ>=5)",
        ], (wellformed.stdout, wellformed.stderr)
        raw_records = subprocess.run(
            [str(BINARY), "-I", str(malformed), "--disable-tool-default-read-filters"],
            text=True, capture_output=True,
        )
        assert raw_records.returncode == 0 and raw_records.stdout.splitlines()[0] == "4 in total", (
            raw_records.stdout, raw_records.stderr)
        raw_filter_disable = subprocess.run(
            [str(BINARY), "-I", str(malformed), "-DF", "WellformedReadFilter"],
            text=True, capture_output=True,
        )
        assert raw_filter_disable.returncode == 0 and raw_filter_disable.stdout.splitlines()[0] == "4 in total", (
            raw_filter_disable.stdout, raw_filter_disable.stderr)

        java = ROOT / "third_party/jdk17/bin/java"
        jar = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        java_oracle = False
        if java.is_file() and jar.is_file():
            oracle = subprocess.run([str(java), "-jar", str(jar), "FlagStat", "-I", str(BAM)],
                                    text=True, capture_output=True)
            if oracle.returncode == 0:
                java_lines = [line for line in oracle.stdout.splitlines() if re.match(r"^\d+ .*", line)]
                assert java_lines == EXPECTED, (java_lines, EXPECTED, oracle.stderr)
                multi_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "FlagStat", "-I", str(BAM),
                     "-I", str(BAM)], text=True, capture_output=True)
                assert multi_oracle.returncode == 0, (
                    multi_oracle.stdout, multi_oracle.stderr)
                assert [line for line in multi_oracle.stdout.splitlines()
                        if re.match(r"^\d+ .*", line)] == [
                            "986 in total", "106 QC failure", "118 duplicates",
                            "974 mapped (98.78%)", "986 paired in sequencing",
                            "494 read1", "492 read2",
                            "954 properly paired (96.75%)",
                            "954 with itself and mate mapped",
                            "20 singletons (2.03%)",
                            "0 with mate mapped to a different chr",
                            "0 with mate mapped to a different chr (mapQ>=5)",
                        ], (multi_oracle.stdout, multi_oracle.stderr)
                interval_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "FlagStat", "-I", str(BAM),
                     "-L", "17:69000-69500"], text=True, capture_output=True)
                assert interval_oracle.returncode == 0
                intersection_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "FlagStat", "-I", str(BAM),
                     "-L", "17:69000-69500", "-L", "17:69200-69300",
                     "--interval-set-rule", "INTERSECTION"], text=True, capture_output=True)
                assert intersection_oracle.returncode == 0
                intersection_lines = [line for line in intersection_oracle.stdout.splitlines()
                                      if re.match(r"^\d+ .*", line)]
                assert intersection_lines == intersection_result.stdout.splitlines(), (
                    intersection_lines, intersection_result.stdout, intersection_oracle.stderr)
                padded_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "FlagStat", "-I", str(BAM),
                     "-L", "17:69990-69990", "-ip", "10"],
                    text=True, capture_output=True)
                assert padded_oracle.returncode == 0
                padded_lines = [line for line in padded_oracle.stdout.splitlines()
                                if re.match(r"^\d+ .*", line)]
                assert padded_lines == padded.stdout.splitlines(), (
                    padded_lines, padded.stdout, padded_oracle.stderr)
                excluded_padding_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "FlagStat", "-I", str(BAM),
                     "-L", "17:69980-70000", "-XL", "17:69990-69990", "-ixp", "2"],
                    text=True, capture_output=True)
                assert excluded_padding_oracle.returncode == 0
                excluded_padding_lines = [line for line in excluded_padding_oracle.stdout.splitlines()
                                          if re.match(r"^\d+ .*", line)]
                assert excluded_padding_lines == excluded_padding.stdout.splitlines(), (
                    excluded_padding_lines, excluded_padding.stdout, excluded_padding_oracle.stderr)
                exclude_only_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "FlagStat", "-I", str(BAM),
                     "-XL", "17:69500-69500", "-ixp", "50"],
                    text=True, capture_output=True)
                assert exclude_only_oracle.returncode == 0
                exclude_only_lines = [line for line in exclude_only_oracle.stdout.splitlines()
                                      if re.match(r"^\d+ .*", line)]
                assert exclude_only_lines == exclude_only.stdout.splitlines(), (
                    exclude_only_lines, exclude_only.stdout, exclude_only_oracle.stderr)
                overlapping_only_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "FlagStat", "-I", str(BAM),
                     "-L", "17:69990-69990", "-L", "17:69991-69991",
                     "-imr", "OVERLAPPING_ONLY"], text=True, capture_output=True)
                assert overlapping_only_oracle.returncode == 0
                overlapping_only_lines = [line for line in overlapping_only_oracle.stdout.splitlines()
                                          if re.match(r"^\d+ .*", line)]
                assert overlapping_only_lines == merged.stdout.splitlines(), (
                    overlapping_only_lines, merged.stdout, overlapping_only_oracle.stderr)
                all_excluded_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "FlagStat", "-I", str(BAM),
                     "-L", "17:69000-69005", "-XL", "17:69000-69005"],
                    text=True, capture_output=True)
                assert all_excluded_oracle.returncode != 0 and "removed all territory" in (
                    all_excluded_oracle.stderr + all_excluded_oracle.stdout)
                length_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "FlagStat", "-I", str(BAM),
                     "-RF", "ReadLengthReadFilter", "--min-read-length", "80",
                     "--max-read-length", "100"], text=True, capture_output=True)
                assert length_oracle.returncode == 0
                assert [line for line in length_oracle.stdout.splitlines()
                        if re.match(r"^\d+ .*", line)] == [
                            "0 in total", "0 QC failure", "0 duplicates",
                            "0 mapped (NaN%)", "0 paired in sequencing", "0 read1",
                            "0 read2", "0 properly paired (NaN%)",
                            "0 with itself and mate mapped", "0 singletons (NaN%)",
                            "0 with mate mapped to a different chr",
                            "0 with mate mapped to a different chr (mapQ>=5)",
                        ], (length_oracle.stdout, length_oracle.stderr)
                malformed_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "FlagStat", "-I", str(malformed)],
                    text=True, capture_output=True,
                )
                assert malformed_oracle.returncode == 0, (
                    malformed_oracle.stdout, malformed_oracle.stderr)
                malformed_lines = [line for line in malformed_oracle.stdout.splitlines()
                                   if re.match(r"^\d+ .*", line)]
                assert malformed_lines and malformed_lines[0] == "1 in total", (
                    malformed_oracle.stdout, malformed_oracle.stderr)
                malformed_raw_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "FlagStat", "-I", str(malformed),
                     "--disable-tool-default-read-filters"],
                    text=True, capture_output=True,
                )
                assert malformed_raw_oracle.returncode == 0, (
                    malformed_raw_oracle.stdout, malformed_raw_oracle.stderr)
                malformed_raw_lines = [line for line in malformed_raw_oracle.stdout.splitlines()
                                      if re.match(r"^\d+ .*", line)]
                assert malformed_raw_lines and malformed_raw_lines[0] == "4 in total", (
                    malformed_raw_oracle.stdout, malformed_raw_oracle.stderr)
                interval_lines = [line for line in interval_oracle.stdout.splitlines()
                                  if re.match(r"^\d+ .*", line)]
                assert interval_lines == interval_result.stdout.splitlines(), (
                    interval_lines, interval_result.stdout, interval_oracle.stderr)
                mapping_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "FlagStat", "-I", str(BAM),
                     "-RF", "MappingQualityReadFilter", "--minimum-mapping-quality", "30",
                     "--maximum-mapping-quality", "50"], text=True, capture_output=True)
                assert mapping_oracle.returncode == 0, (mapping_oracle.stdout, mapping_oracle.stderr)
                mapping_lines = [line for line in mapping_oracle.stdout.splitlines()
                                 if re.match(r"^\d+ .*", line)]
                assert mapping_lines == [
                    "18 in total", "2 QC failure", "2 duplicates", "18 mapped (100.00%)",
                    "18 paired in sequencing", "11 read1", "7 read2", "9 properly paired (50.00%)",
                    "9 with itself and mate mapped", "9 singletons (50.00%)",
                    "0 with mate mapped to a different chr",
                    "0 with mate mapped to a different chr (mapQ>=5)",
                ], (mapping_lines, mapping_oracle.stdout, mapping_oracle.stderr)
                inverted_mapping_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "FlagStat", "-I", str(BAM),
                     "-XRF", "MappingQualityReadFilter", "--minimum-mapping-quality", "30",
                     "--maximum-mapping-quality", "50"], text=True, capture_output=True)
                assert inverted_mapping_oracle.returncode == 0, (
                    inverted_mapping_oracle.stdout, inverted_mapping_oracle.stderr)
                inverted_mapping_lines = [line for line in inverted_mapping_oracle.stdout.splitlines()
                                          if re.match(r"^\d+ .*", line)]
                assert inverted_mapping_lines == inverted_mapping_quality.stdout.splitlines(), (
                    inverted_mapping_lines, inverted_mapping_quality.stdout,
                    inverted_mapping_oracle.stderr)
                read_name_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "FlagStat", "-I", str(BAM),
                     "-RF", "ReadNameReadFilter", "--read-name", read_name],
                    text=True, capture_output=True)
                assert read_name_oracle.returncode == 0
                read_name_lines = [line for line in read_name_oracle.stdout.splitlines()
                                   if re.match(r"^\d+ .*", line)]
                assert read_name_lines == read_name_result.stdout.splitlines(), (
                    read_name_lines, read_name_result.stdout, read_name_oracle.stderr)
                read_group_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "FlagStat", "-I", str(BAM),
                     "-RF", "ReadGroupReadFilter", "--keep-read-group", "809R9.5"],
                    text=True, capture_output=True)
                assert read_group_oracle.returncode == 0
                read_group_lines = [line for line in read_group_oracle.stdout.splitlines()
                                    if re.match(r"^\d+ .*", line)]
                assert read_group_lines == EXPECTED, (read_group_lines, read_group_oracle.stderr)
                read_tag_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "FlagStat", "-I", str(BAM),
                     "-RF", "ReadTagValueFilter", "--read-filter-tag", "NM",
                     "--read-filter-tag-comp", "1", "--read-filter-tag-op", "GREATER_OR_EQUAL"],
                    text=True, capture_output=True)
                assert read_tag_oracle.returncode == 0
                read_tag_lines = [line for line in read_tag_oracle.stdout.splitlines()
                                  if re.match(r"^\d+ .*", line)]
                assert read_tag_lines == read_tag_result.stdout.splitlines(), (
                    read_tag_lines, read_tag_result.stdout, read_tag_oracle.stderr)
                blacklist_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "FlagStat", "-I", str(BAM),
                     "-RF", "ReadGroupBlackListReadFilter", "--read-group-black-list", "RG:809R9.5"],
                    text=True, capture_output=True)
                assert blacklist_oracle.returncode == 0
                blacklist_lines = [line for line in blacklist_oracle.stdout.splitlines()
                                   if re.match(r"^\d+ .*", line)]
                assert blacklist_lines == read_group_blacklist.stdout.splitlines(), (
                    blacklist_lines, read_group_blacklist.stdout, blacklist_oracle.stderr)
                nonzero_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "FlagStat", "-I", str(BAM),
                     "-RF", "NonZeroFragmentLengthReadFilter"], text=True, capture_output=True)
                assert nonzero_oracle.returncode == 0, (
                    nonzero_oracle.stdout, nonzero_oracle.stderr)
                nonzero_lines = [line for line in nonzero_oracle.stdout.splitlines()
                                 if re.match(r"^\d+ .*", line)]
                assert nonzero_lines == [
                    "477 in total", "51 QC failure", "58 duplicates", "477 mapped (100.00%)",
                    "477 paired in sequencing", "238 read1", "239 read2", "477 properly paired (100.00%)",
                    "477 with itself and mate mapped", "0 singletons (0.00%)",
                    "0 with mate mapped to a different chr",
                    "0 with mate mapped to a different chr (mapQ>=5)",
                ], (nonzero_lines, nonzero_oracle.stdout, nonzero_oracle.stderr)
                fragment_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "FlagStat", "-I", str(BAM),
                     "-RF", "FragmentLengthReadFilter", "--min-fragment-length", "100",
                     "--max-fragment-length", "300"], text=True, capture_output=True)
                assert fragment_oracle.returncode == 0, (
                    fragment_oracle.stdout, fragment_oracle.stderr)
                fragment_lines = [line for line in fragment_oracle.stdout.splitlines()
                                  if re.match(r"^\d+ .*", line)]
                assert fragment_lines == [
                    "413 in total", "45 QC failure", "56 duplicates", "413 mapped (100.00%)",
                    "413 paired in sequencing", "206 read1", "207 read2", "413 properly paired (100.00%)",
                    "413 with itself and mate mapped", "0 singletons (0.00%)",
                    "0 with mate mapped to a different chr",
                    "0 with mate mapped to a different chr (mapQ>=5)",
                ], (fragment_lines, fragment_oracle.stdout, fragment_oracle.stderr)
                java_oracle = True
    print(json.dumps({"status": "pass", "suite": "flag-stat", "records": 493,
                      "java_oracle": java_oracle}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
