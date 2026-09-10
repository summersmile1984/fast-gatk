#!/usr/bin/env python3
"""Contract/oracle test for the native CountReads path."""

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
BINARY = Path(os.environ.get("FASTGATK_COUNT_READS_BINARY", ROOT / "fastgatk-native/build/fastgatk-count-reads"))


def main() -> int:
    if not BINARY.is_file():
        print(json.dumps({"status": "skip", "suite": "count-reads", "reason": f"missing binary: {BINARY}"}))
        return 0
    with tempfile.TemporaryDirectory(prefix="fastgatk-count-reads-") as directory:
        work = Path(directory)
        output = work / "count.txt"
        manifest = work / "count.manifest.json"
        result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-O", str(output), "--batch-records", "37",
             "--threads", "2", "--output-manifest", str(manifest)],
            text=True, capture_output=True,
        )
        assert result.returncode == 0, (result.stdout, result.stderr)
        assert result.stdout.strip() == "493"
        assert output.read_text().strip() == "493"
        payload = json.loads(manifest.read_text())
        assert payload["tool"] == "CountReads"
        # Registry and binary-side OutputManifest must expose the same
        # compatibility level; otherwise a stale binary could be selected as
        # contract-compatible by the dispatcher without proving its output
        # contract.
        assert payload["status"] == "contract-compatible"
        assert payload["input"] == str(BAM)
        assert payload["input_count"] == 1
        assert payload["inputs"] == [str(BAM)]
        assert payload["records"] == 493
        assert payload["batches"] > 1
        assert payload["default_read_filter_count"] == 1
        assert payload["default_read_filters"] == ["WellformedReadFilter"]
        assert payload["telemetry"]["persistent_buffer_allocations"] == 1
        assert payload["telemetry"]["persistent_buffer_reuses"] > 0
        assert payload["telemetry"]["persistent_buffer_capacity_records"] >= 37
        assert payload["telemetry"]["kernel_execution_policy"] == "RangePolicy+reduction"
        assert payload["telemetry"]["kernel_lifecycle"] == (
            "HostBatch->KernelPlan.prepare->Kokkos Views->execute->collect"
        )
        assert payload["telemetry"]["pipeline_lifecycle"] == "decode->compute->encode->sink"
        assert payload["telemetry"]["pipeline_decoded_items"] == payload["batches"]
        assert payload["telemetry"]["pipeline_computed_items"] == payload["batches"]
        assert payload["telemetry"]["pipeline_encoded_items"] == payload["batches"]
        assert payload["telemetry"]["pipeline_peak_decoded_bytes"] > 0
        assert payload["telemetry"]["pipeline_peak_computed_bytes"] > 0
        assert payload["telemetry"]["pipeline_peak_encoded_bytes"] > 0
        assert json.loads(result.stderr.splitlines()[-1])["status"] == "contract-compatible"

        # GATK ReadWalker accepts repeated -I values and aggregates each
        # shard in input order.  This is the normal SLURM/Nextflow gather
        # shape, so exercise both the count and manifest contract.
        multi_manifest = work / "multi.manifest.json"
        multi_result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "--input", str(BAM),
             "--batch-records", "37", "--output-manifest", str(multi_manifest)],
            text=True, capture_output=True,
        )
        assert multi_result.returncode == 0, (multi_result.stdout, multi_result.stderr)
        assert multi_result.stdout.strip() == "986"
        multi_payload = json.loads(multi_manifest.read_text())
        assert multi_payload["records"] == 986
        assert multi_payload["input_count"] == 2
        assert multi_payload["inputs"] == [str(BAM), str(BAM)]
        assert multi_payload["telemetry"]["pipeline_decoded_items"] == \
            multi_payload["telemetry"]["pipeline_computed_items"]
        assert multi_payload["telemetry"]["pipeline_computed_items"] == \
            multi_payload["telemetry"]["pipeline_encoded_items"]

        # Exercise the CRAM-compatible reference injection path on a BAM
        # fixture as well; HTSlib accepts the configured .fai even though this
        # record-count walk does not need to fetch bases from it.
        reference = ROOT / "gatk-source/src/test/resources/Homo_sapiens_assembly38_chrM_only.fasta"
        reference_result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-R", str(reference), "--batch-records", "128"],
            text=True, capture_output=True,
        )
        assert reference_result.returncode == 0, (reference_result.stdout, reference_result.stderr)
        assert reference_result.stdout.strip() == "493"

        interval_result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", "17:69000-69500", "--batch-records", "29",
             "--output-manifest", str(work / "interval.manifest.json")],
            text=True, capture_output=True,
        )
        assert interval_result.returncode == 0, (interval_result.stdout, interval_result.stderr)
        assert interval_result.stdout.strip() == "338"
        interval_payload = json.loads((work / "interval.manifest.json").read_text())
        assert interval_payload["intervals"] == 1
        assert interval_payload["indexed_intervals"] is True

        padded = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", "17:69990-69990", "-ip", "10",
             "--output-manifest", str(work / "padded.manifest.json")],
            text=True, capture_output=True,
        )
        assert padded.returncode == 0 and padded.stdout.strip() == "17", (
            padded.stdout, padded.stderr)
        padded_payload = json.loads((work / "padded.manifest.json").read_text())
        assert padded_payload["interval_padding"] == 10
        assert padded_payload["interval_exclusion_padding"] == 0

        excluded_padding = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", "17:69980-70000",
             "-XL", "17:69990-69990", "-ixp", "2"],
            text=True, capture_output=True,
        )
        assert excluded_padding.returncode == 0 and excluded_padding.stdout.strip() == "17", (
            excluded_padding.stdout, excluded_padding.stderr)
        exclude_only = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-XL", "17:69500-69500", "-ixp", "50"],
            text=True, capture_output=True,
        )
        assert exclude_only.returncode == 0 and exclude_only.stdout.strip() == "415", (
            exclude_only.stdout, exclude_only.stderr)
        merged = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", "17:69990-69990",
             "-L", "17:69991-69991", "--interval-merging-rule=OVERLAPPING_ONLY",
             "--output-manifest", str(work / "overlap-only.manifest.json")],
            text=True, capture_output=True,
        )
        assert merged.returncode == 0 and merged.stdout.strip() == "11", (
            merged.stdout, merged.stderr)
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

        # Repeatable -L selectors are intersected as selector sets when the
        # GATK interval-set-rule is INTERSECTION (an interval file remains one
        # selector, even when it contains multiple records).
        intersection = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", "17:69000-69500", "-L", "17:69200-69300",
             "-isr", "INTERSECTION", "--output-manifest", str(work / "intersection.manifest.json")],
            text=True, capture_output=True,
        )
        assert intersection.returncode == 0 and intersection.stdout.strip() == "11", (
            intersection.stdout, intersection.stderr)
        intersection_payload = json.loads((work / "intersection.manifest.json").read_text())
        assert intersection_payload["interval_set_rule"] == "INTERSECTION"

        invalid_rule = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", "17:69000-69500",
             "--interval-set-rule", "XOR"],
            text=True, capture_output=True,
        )
        assert invalid_rule.returncode != 0
        assert "interval-set-rule must be UNION or INTERSECTION" in invalid_rule.stderr

        interval_file = work / "regions.list"
        interval_file.write_text("17:69000-69500\n", encoding="utf-8")
        file_result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", str(interval_file), "-XL", "17:69000-69050",
             "--batch-records", "31"], text=True, capture_output=True,
        )
        assert file_result.returncode == 0, (file_result.stdout, file_result.stderr)
        assert file_result.stdout.strip() == "334"

        compressed_interval = work / "regions.interval_list.gz"
        with gzip.open(compressed_interval, "wt", encoding="utf-8") as stream:
            stream.write("17:69000-69500\n")
        compressed_result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", str(compressed_interval), "--batch-records", "31"],
            text=True, capture_output=True,
        )
        assert compressed_result.returncode == 0, (compressed_result.stdout, compressed_result.stderr)
        assert compressed_result.stdout.strip() == "338"

        compressed_bed = work / "regions.bed.gz"
        with gzip.open(compressed_bed, "wt", encoding="utf-8") as stream:
            stream.write("track name=targets\n17\t68999\t69500\n")
        compressed_bed_result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", str(compressed_bed), "--batch-records", "31"],
            text=True, capture_output=True,
        )
        assert compressed_bed_result.returncode == 0, (
            compressed_bed_result.stdout, compressed_bed_result.stderr)
        assert compressed_bed_result.stdout.strip() == "338"

        malformed_compressed = work / "malformed.interval_list.gz"
        with gzip.open(malformed_compressed, "wt", encoding="utf-8") as stream:
            stream.write("17\t69000\n")
        malformed_result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-L", str(malformed_compressed)],
            text=True, capture_output=True,
        )
        assert malformed_result.returncode != 0 and "malformed interval file" in malformed_result.stderr

        filter_result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-RF", "NotDuplicateReadFilter", "-RF", "MappedReadFilter",
             "--batch-records", "43", "--output-manifest", str(work / "filter.manifest.json")],
            text=True, capture_output=True,
        )
        assert filter_result.returncode == 0, (filter_result.stdout, filter_result.stderr)
        assert filter_result.stdout.strip() == "428"
        filter_payload = json.loads((work / "filter.manifest.json").read_text())
        assert filter_payload["read_filter_count"] == 2
        mapping_quality = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-RF", "MappingQualityReadFilter",
             "--minimum-mapping-quality", "30", "--maximum-mapping-quality", "50",
             "--output-manifest", str(work / "mapping-quality.manifest.json")],
            text=True, capture_output=True,
        )
        assert mapping_quality.returncode == 0, (mapping_quality.stdout, mapping_quality.stderr)
        assert mapping_quality.stdout.strip() == "18"
        mapping_payload = json.loads((work / "mapping-quality.manifest.json").read_text())
        assert mapping_payload["read_filter_min_mapping_quality"] == 30
        assert mapping_payload["read_filter_max_mapping_quality"] == 50
        mapping_quality_default = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-RF", "MappingQualityReadFilter",
             "--minimum-mapping-quality", "30"], text=True, capture_output=True,
        )
        assert mapping_quality_default.returncode == 0
        assert mapping_quality_default.stdout.strip() == "410"
        mapping_quality_without_filter = subprocess.run(
            [str(BINARY), "-I", str(BAM), "--minimum-mapping-quality", "30"],
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
        assert read_name_result.returncode == 0 and read_name_result.stdout.strip() == "2", (
            read_name_result.stdout, read_name_result.stderr)
        read_name_payload = json.loads((work / "read-name.manifest.json").read_text())
        assert read_name_payload["read_name_count"] == 1
        read_name_without_filter = subprocess.run(
            [str(BINARY), "-I", str(BAM), "--read-name", read_name],
            text=True, capture_output=True,
        )
        assert read_name_without_filter.returncode != 0 and "requires ReadNameReadFilter" in read_name_without_filter.stderr
        read_group_result = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-RF", "ReadGroupReadFilter",
             "--keep-read-group", "809R9.5",
             "--output-manifest", str(work / "read-group.manifest.json")],
            text=True, capture_output=True,
        )
        assert read_group_result.returncode == 0 and read_group_result.stdout.strip() == "493", (
            read_group_result.stdout, read_group_result.stderr)
        read_group_payload = json.loads((work / "read-group.manifest.json").read_text())
        assert read_group_payload["keep_read_group_count"] == 1
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
        assert read_tag_result.returncode == 0 and read_tag_result.stdout.strip() == "189", (
            read_tag_result.stdout, read_tag_result.stderr)
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
             "--read-group-black-list", "RG:809R9.5",
             "--output-manifest", str(work / "read-group-blacklist.manifest.json")],
            text=True, capture_output=True,
        )
        assert read_group_blacklist.returncode == 0 and read_group_blacklist.stdout.strip() == "0", (
            read_group_blacklist.stdout, read_group_blacklist.stderr)
        blacklist_payload = json.loads((work / "read-group-blacklist.manifest.json").read_text())
        assert blacklist_payload["read_group_blacklist_count"] == 1
        blacklist_without_filter = subprocess.run(
            [str(BINARY), "-I", str(BAM), "--read-group-black-list", "RG:809R9.5"],
            text=True, capture_output=True,
        )
        assert blacklist_without_filter.returncode != 0 and "requires ReadGroupBlackListReadFilter" in blacklist_without_filter.stderr
        simple_filter_expectations = {
            "NonZeroFragmentLengthReadFilter": "477",
            "MateOnSameContigOrNoMappedMateReadFilter": "487",
            "CigarContainsNoNOperator": "493",
            "GoodCigarReadFilter": "493",
            "HasReadGroupReadFilter": "493",
            "MateDifferentStrandReadFilter": "477",
            "NotProperlyPairedReadFilter": "16",
        }
        for filter_name, expected_count in simple_filter_expectations.items():
            filtered = subprocess.run(
                [str(BINARY), "-I", str(BAM), "-RF", filter_name],
                text=True, capture_output=True,
            )
            assert filtered.returncode == 0 and filtered.stdout.strip() == expected_count, (
                filter_name, filtered.stdout, filtered.stderr)
        fragment_filter = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-RF", "FragmentLengthReadFilter",
             "--min-fragment-length", "100", "--max-fragment-length", "300",
             "--output-manifest", str(work / "fragment-length.manifest.json")],
            text=True, capture_output=True,
        )
        assert fragment_filter.returncode == 0 and fragment_filter.stdout.strip() == "413", (
            fragment_filter.stdout, fragment_filter.stderr)
        fragment_payload = json.loads((work / "fragment-length.manifest.json").read_text())
        assert fragment_payload["read_filter_min_fragment_length"] == 100
        assert fragment_payload["read_filter_max_fragment_length"] == 300
        fragment_without_filter = subprocess.run(
            [str(BINARY), "-I", str(BAM), "--max-fragment-length", "300"],
            text=True, capture_output=True,
        )
        assert fragment_without_filter.returncode != 0
        assert "requires FragmentLengthReadFilter" in fragment_without_filter.stderr
        length_filter = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-RF", "ReadLengthReadFilter",
             "--min-read-length", "80", "--max-read-length", "100",
             "--output-manifest", str(work / "length-filter.manifest.json")],
            text=True, capture_output=True,
        )
        assert length_filter.returncode == 0, (length_filter.stdout, length_filter.stderr)
        assert length_filter.stdout.strip() == "0"
        length_payload = json.loads((work / "length-filter.manifest.json").read_text())
        assert length_payload["read_filter_require_read_length"] is True
        assert length_payload["read_filter_min_read_length"] == 80
        assert length_payload["read_filter_max_read_length"] == 100
        no_defaults = subprocess.run(
            [str(BINARY), "-I", str(BAM), "--disable-tool-default-read-filters",
             "--output-manifest", str(work / "no-default-filters.manifest.json")],
            text=True, capture_output=True,
        )
        assert no_defaults.returncode == 0 and no_defaults.stdout.strip() == "493", (
            no_defaults.stdout, no_defaults.stderr)
        no_defaults_payload = json.loads(
            (work / "no-default-filters.manifest.json").read_text())
        assert no_defaults_payload["disable_tool_default_read_filters"] is True
        assert no_defaults_payload["default_read_filter_count"] == 0
        assert no_defaults_payload["default_read_filters"] == []

        # CountReads inherits ReadWalker.WellformedReadFilter by default.  A
        # valid SAM containing a missing RG, a skipped-reference CIGAR, and an
        # empty sequence must therefore count only the well-formed record;
        # disabling tool defaults restores the raw record count.
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
        assert wellformed.returncode == 0 and wellformed.stdout.strip() == "1", (
            wellformed.stdout, wellformed.stderr)
        raw_records = subprocess.run(
            [str(BINARY), "-I", str(malformed), "--disable-tool-default-read-filters"],
            text=True, capture_output=True,
        )
        assert raw_records.returncode == 0 and raw_records.stdout.strip() == "4", (
            raw_records.stdout, raw_records.stderr)
        raw_filter_disable = subprocess.run(
            [str(BINARY), "-I", str(malformed), "-DF", "WellformedReadFilter"],
            text=True, capture_output=True,
        )
        assert raw_filter_disable.returncode == 0 and raw_filter_disable.stdout.strip() == "4", (
            raw_filter_disable.stdout, raw_filter_disable.stderr)

        unknown_filter = subprocess.run(
            [str(BINARY), "-I", str(BAM), "-RF", "UnsupportedReadFilter"],
            text=True, capture_output=True,
        )
        assert unknown_filter.returncode != 0 and "UNSUPPORTED_PARAMETER" in unknown_filter.stderr

        java = ROOT / "third_party/jdk17/bin/java"
        jar = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        java_oracle = False
        if java.is_file() and jar.is_file():
            oracle = subprocess.run([str(java), "-jar", str(jar), "CountReads", "-I", str(BAM)],
                                    text=True, capture_output=True)
            if oracle.returncode == 0:
                matches = re.findall(r"(?m)^\s*(\d+)\s*$", oracle.stdout)
                assert matches and matches[-1] == "493", (oracle.stdout, oracle.stderr)
                multi_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "CountReads", "-I", str(BAM),
                     "-I", str(BAM)], text=True, capture_output=True)
                assert multi_oracle.returncode == 0, (
                    multi_oracle.stdout, multi_oracle.stderr)
                multi_matches = re.findall(r"(?m)^\s*(\d+)\s*$", multi_oracle.stdout)
                assert multi_matches and multi_matches[-1] == "986", (
                    multi_oracle.stdout, multi_oracle.stderr)
                interval_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "CountReads", "-I", str(BAM),
                     "-L", "17:69000-69500"], text=True, capture_output=True)
                assert interval_oracle.returncode == 0
                intersection_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "CountReads", "-I", str(BAM),
                     "-L", "17:69000-69500", "-L", "17:69200-69300",
                     "--interval-set-rule", "INTERSECTION"], text=True, capture_output=True)
                assert intersection_oracle.returncode == 0
                intersection_matches = re.findall(r"(?m)^\s*(\d+)\s*$", intersection_oracle.stdout)
                assert intersection_matches and intersection_matches[-1] == "11", (
                    intersection_oracle.stdout, intersection_oracle.stderr)
                padded_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "CountReads", "-I", str(BAM),
                     "-L", "17:69990-69990", "-ip", "10"],
                    text=True, capture_output=True)
                assert padded_oracle.returncode == 0
                padded_matches = re.findall(r"(?m)^\s*(\d+)\s*$", padded_oracle.stdout)
                assert padded_matches and padded_matches[-1] == "17", (
                    padded_oracle.stdout, padded_oracle.stderr)
                excluded_padding_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "CountReads", "-I", str(BAM),
                     "-L", "17:69980-70000", "-XL", "17:69990-69990", "-ixp", "2"],
                    text=True, capture_output=True)
                assert excluded_padding_oracle.returncode == 0
                excluded_padding_matches = re.findall(
                    r"(?m)^\s*(\d+)\s*$", excluded_padding_oracle.stdout)
                assert excluded_padding_matches and excluded_padding_matches[-1] == "17", (
                    excluded_padding_oracle.stdout, excluded_padding_oracle.stderr)
                exclude_only_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "CountReads", "-I", str(BAM),
                     "-XL", "17:69500-69500", "-ixp", "50"],
                    text=True, capture_output=True)
                assert exclude_only_oracle.returncode == 0
                exclude_only_matches = re.findall(r"(?m)^\s*(\d+)\s*$", exclude_only_oracle.stdout)
                assert exclude_only_matches and exclude_only_matches[-1] == "415", (
                    exclude_only_oracle.stdout, exclude_only_oracle.stderr)
                overlapping_only_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "CountReads", "-I", str(BAM),
                     "-L", "17:69990-69990", "-L", "17:69991-69991",
                     "-imr", "OVERLAPPING_ONLY"], text=True, capture_output=True)
                assert overlapping_only_oracle.returncode == 0
                overlapping_only_matches = re.findall(
                    r"(?m)^\s*(\d+)\s*$", overlapping_only_oracle.stdout)
                assert overlapping_only_matches and overlapping_only_matches[-1] == "11", (
                    overlapping_only_oracle.stdout, overlapping_only_oracle.stderr)
                all_excluded_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "CountReads", "-I", str(BAM),
                     "-L", "17:69000-69005", "-XL", "17:69000-69005"],
                    text=True, capture_output=True)
                assert all_excluded_oracle.returncode != 0 and "removed all territory" in (
                    all_excluded_oracle.stderr + all_excluded_oracle.stdout)
                length_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "CountReads", "-I", str(BAM),
                     "-RF", "ReadLengthReadFilter", "--min-read-length", "80",
                     "--max-read-length", "100"], text=True, capture_output=True)
                assert length_oracle.returncode == 0
                length_matches = re.findall(r"(?m)^\s*(\d+)\s*$", length_oracle.stdout)
                assert length_matches and length_matches[-1] == "0", (
                    length_oracle.stdout, length_oracle.stderr)
                malformed_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "CountReads", "-I", str(malformed)],
                    text=True, capture_output=True,
                )
                assert malformed_oracle.returncode == 0, (
                    malformed_oracle.stdout, malformed_oracle.stderr)
                malformed_matches = re.findall(r"(?m)^\s*(\d+)\s*$", malformed_oracle.stdout)
                assert malformed_matches and malformed_matches[-1] == "1", (
                    malformed_oracle.stdout, malformed_oracle.stderr)
                malformed_raw_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "CountReads", "-I", str(malformed),
                     "--disable-tool-default-read-filters"],
                    text=True, capture_output=True,
                )
                assert malformed_raw_oracle.returncode == 0, (
                    malformed_raw_oracle.stdout, malformed_raw_oracle.stderr)
                malformed_raw_matches = re.findall(
                    r"(?m)^\s*(\d+)\s*$", malformed_raw_oracle.stdout)
                assert malformed_raw_matches and malformed_raw_matches[-1] == "4", (
                    malformed_raw_oracle.stdout, malformed_raw_oracle.stderr)
                interval_matches = re.findall(r"(?m)^\s*(\d+)\s*$", interval_oracle.stdout)
                assert interval_matches and interval_matches[-1] == "338", (
                    interval_oracle.stdout, interval_oracle.stderr)
                exclusion_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "CountReads", "-I", str(BAM),
                     "-L", "17:69000-69500", "-XL", "17:69000-69050"],
                    text=True, capture_output=True)
                assert exclusion_oracle.returncode == 0
                exclusion_matches = re.findall(r"(?m)^\s*(\d+)\s*$", exclusion_oracle.stdout)
                assert exclusion_matches and exclusion_matches[-1] == "334", (
                    exclusion_oracle.stdout, exclusion_oracle.stderr)
                mapping_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "CountReads", "-I", str(BAM),
                     "-RF", "MappingQualityReadFilter", "--minimum-mapping-quality", "30",
                     "--maximum-mapping-quality", "50"], text=True, capture_output=True)
                assert mapping_oracle.returncode == 0, (mapping_oracle.stdout, mapping_oracle.stderr)
                mapping_matches = re.findall(r"(?m)^\s*(\d+)\s*$", mapping_oracle.stdout)
                assert mapping_matches and mapping_matches[-1] == "18", (
                    mapping_oracle.stdout, mapping_oracle.stderr)
                read_name_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "CountReads", "-I", str(BAM),
                     "-RF", "ReadNameReadFilter", "--read-name", read_name],
                    text=True, capture_output=True)
                assert read_name_oracle.returncode == 0
                read_name_matches = re.findall(r"(?m)^\s*(\d+)\s*$", read_name_oracle.stdout)
                assert read_name_matches and read_name_matches[-1] == "2", (
                    read_name_oracle.stdout, read_name_oracle.stderr)
                read_group_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "CountReads", "-I", str(BAM),
                     "-RF", "ReadGroupReadFilter", "--keep-read-group", "809R9.5"],
                    text=True, capture_output=True)
                assert read_group_oracle.returncode == 0
                read_group_matches = re.findall(r"(?m)^\s*(\d+)\s*$", read_group_oracle.stdout)
                assert read_group_matches and read_group_matches[-1] == "493", (
                    read_group_oracle.stdout, read_group_oracle.stderr)
                read_tag_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "CountReads", "-I", str(BAM),
                     "-RF", "ReadTagValueFilter", "--read-filter-tag", "NM",
                     "--read-filter-tag-comp", "1", "--read-filter-tag-op", "GREATER_OR_EQUAL"],
                    text=True, capture_output=True)
                assert read_tag_oracle.returncode == 0
                read_tag_matches = re.findall(r"(?m)^\s*(\d+)\s*$", read_tag_oracle.stdout)
                assert read_tag_matches and read_tag_matches[-1] == "189", (
                    read_tag_oracle.stdout, read_tag_oracle.stderr)
                blacklist_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "CountReads", "-I", str(BAM),
                     "-RF", "ReadGroupBlackListReadFilter", "--read-group-black-list", "RG:809R9.5"],
                    text=True, capture_output=True)
                assert blacklist_oracle.returncode == 0
                blacklist_matches = re.findall(r"(?m)^\s*(\d+)\s*$", blacklist_oracle.stdout)
                assert blacklist_matches and blacklist_matches[-1] == "0", (
                    blacklist_oracle.stdout, blacklist_oracle.stderr)
                for filter_name, expected_count in simple_filter_expectations.items():
                    filter_oracle = subprocess.run(
                        [str(java), "-jar", str(jar), "CountReads", "-I", str(BAM),
                         "-RF", filter_name], text=True, capture_output=True)
                    assert filter_oracle.returncode == 0, (
                        filter_name, filter_oracle.stdout, filter_oracle.stderr)
                    filter_matches = re.findall(r"(?m)^\s*(\d+)\s*$", filter_oracle.stdout)
                    assert filter_matches and filter_matches[-1] == expected_count, (
                        filter_name, filter_oracle.stdout, filter_oracle.stderr)
                fragment_oracle = subprocess.run(
                    [str(java), "-jar", str(jar), "CountReads", "-I", str(BAM),
                     "-RF", "FragmentLengthReadFilter", "--min-fragment-length", "100",
                     "--max-fragment-length", "300"], text=True, capture_output=True)
                assert fragment_oracle.returncode == 0, (
                    fragment_oracle.stdout, fragment_oracle.stderr)
                fragment_matches = re.findall(r"(?m)^\s*(\d+)\s*$", fragment_oracle.stdout)
                assert fragment_matches and fragment_matches[-1] == "413", (
                    fragment_oracle.stdout, fragment_oracle.stderr)
                java_oracle = True
    print(json.dumps({"status": "pass", "suite": "count-reads", "records": 493,
                      "java_oracle": java_oracle}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
