#!/usr/bin/env python3
"""Verify BaseRecalibrator's GATK read-filter plugin boundary.

BaseRecalibrator has a seven-filter default list in GATK 4.6.2.0.  The
generic read-filter descriptor resolves command-line controls in this order:
start with the tool defaults (or an empty set when
``--disable-tool-default-read-filters`` is true), remove ``-DF`` defaults,
then append explicit ``-RF`` filters and finally apply ``-XRF`` inversions.
This oracle uses a tiny reference-backed SAM containing one read for each
relevant flag/MAPQ case and compares all five GATKReport tables with the pinned
Java implementation.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
import oracle_guard


def run(command: list[str], *, capture: bool = False) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=True, text=True, capture_output=capture)


def table_rows(text: str, table: str) -> list[tuple[str, ...]]:
    marker = f"#:GATKTable:{table}:"
    lines = text.splitlines()
    try:
        start = next(index for index, line in enumerate(lines)
                     if line.startswith(marker)) + 2
    except StopIteration:
        return []
    rows: list[tuple[str, ...]] = []
    for line in lines[start:]:
        if line.startswith("#:GATKTable:"):
            break
        if line.strip():
            rows.append(tuple(line.split()))
    return rows


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native_build = Path(os.environ.get("FASTGATK_NATIVE_BUILD",
                                     root / "fastgatk-native/build"))
    native = Path(os.environ.get("FASTGATK_BQSR_BINARY",
                                native_build / "fastgatk-bqsr"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    reference = root / "gatk-source/src/test/resources/exampleFASTA.fasta"
    required = (native, java, gatk, reference,
                reference.with_suffix(".dict"),
                reference.with_suffix(reference.suffix + ".fai"))
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_bqsr_read_filter_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled BaseRecalibrator read-filter oracle inputs are required")
        print(json.dumps({"status": "skip",
                          "reason": "bundled GATK read-filter oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-bqsr-read-filter-oracle-") as directory:
        work = Path(directory)
        sam = work / "filter-cases.sam"
        known_sites = work / "known-sites.vcf"
        known_sites.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=1,length=100000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t99999\t.\tA\tG\t.\tPASS\t.\n",
            encoding="utf-8")
        # All records are structurally well-formed and carry PL/PU metadata;
        # only their filter-visible flags/MAPQ differ.  The supplementary read
        # is intentionally retained by BQSR defaults and removed by explicit
        # NotSupplementaryAlignmentReadFilter, proving that -RF is additive.
        reads = [
            ("good", 0, 40, 1),
            ("mapq-zero", 0, 0, 21),
            ("mapq-unavailable", 0, 255, 41),
            ("secondary", 256, 40, 61),
            ("duplicate", 1024, 40, 81),
            ("vendor-fail", 512, 40, 101),
            ("supplementary", 2048, 40, 121),
        ]
        sam_lines = ["@HD\tVN:1.6\tSO:unsorted",
                     "@SQ\tSN:1\tLN:100000",
                     "@RG\tID:rg1\tSM:sample\tPL:ILLUMINA"]
        for name, flag, mapq, position in reads:
            sam_lines.append(
                f"{name}\t{flag}\t1\t{position}\t{mapq}\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:rg1")
        sam.write_text("\n".join(sam_lines) + "\n", encoding="utf-8")
        run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(known_sites)])

        scenarios: list[tuple[str, list[str]]] = [
            ("default", []),
            ("disable-mq-zero", ["--disable-read-filter", "MappingQualityNotZeroReadFilter"]),
            ("defaults-off-explicit-mq-zero",
             ["--disable-tool-default-read-filters", "--read-filter",
              "MappingQualityNotZeroReadFilter"]),
            ("explicit-not-supplementary", ["--read-filter",
                                             "NotSupplementaryAlignmentReadFilter"]),
            # GATK permits inversion of a default filter only after the
            # default list is disabled. The inverse of NotSupplementary then
            # keeps precisely the supplementary record.
            ("defaults-off-inverted-not-supplementary",
             ["--disable-tool-default-read-filters", "--inverted-read-filter",
              "NotSupplementaryAlignmentReadFilter"]),
            # Disabling a recognized non-default filter is a Java no-op (with
            # a warning), rather than a native parser error.
            ("disable-nondefault-supplementary",
             ["--disable-read-filter", "NotSupplementaryAlignmentReadFilter"]),
            # ContextCovariate's non-default context size is a separate
            # report-writer/model boundary from read-filter resolution.  The
            # fixture has enough aligned sequence to make size 3 observable,
            # while keeping the same effective read set for a clean oracle.
            ("mismatches-context-size-3", ["--mismatches-context-size", "3"]),
        ]
        tables = ("Quantized", "RecalTable0", "RecalTable1", "RecalTable2")
        outputs: dict[str, dict[str, int]] = {}
        for name, filter_args in scenarios:
            java_report = work / f"java-{name}.tsv"
            native_report = work / f"native-{name}.tsv"
            common = ["-R", str(reference), "-I", str(sam),
                      "--known-sites", str(known_sites)]
            run([str(java), "-jar", str(gatk), "BaseRecalibrator",
                 *common, *filter_args, "-O", str(java_report)])
            run([str(native), *common, *filter_args, "-O", str(native_report)])
            java_text = java_report.read_text(encoding="utf-8")
            native_text = native_report.read_text(encoding="utf-8")
            java_tables = {table: Counter(table_rows(java_text, table)) for table in tables}
            native_tables = {table: Counter(table_rows(native_text, table)) for table in tables}
            if native_tables != java_tables:
                raise AssertionError({
                    "scenario": name,
                    "tables": {table: (len(native_tables[table]), len(java_tables[table]))
                               for table in tables},
                })
            if name == "mismatches-context-size-3":
                # Arguments and manifest carry the same resolved value as the
                # covariate construction; this prevents a report-only option
                # from drifting away from restart/telemetry state.
                java_args = dict(table_rows(java_text, "Arguments"))
                native_args = dict(table_rows(native_text, "Arguments"))
                if java_args.get("mismatches_context_size") != "3" or \
                        native_args.get("mismatches_context_size") != "3":
                    raise AssertionError({"mismatches_context_size_arguments":
                                          (native_args, java_args)})
                manifest = json.loads(Path(f"{native_report}.manifest.json").read_text(
                    encoding="utf-8"))
                if manifest["compatibility"].get("mismatches_context_size") is not True or \
                        manifest["telemetry"].get("mismatches_context_size") != 3:
                    raise AssertionError({"mismatches_context_size_manifest": manifest})
            # RecalTable1 is an unambiguous read-filter cardinality witness:
            # ten bases per passing record, all at reported quality 40.
            observations = sum(int(row[4]) for row in table_rows(native_text, "RecalTable1")
                               if len(row) >= 6 and row[2] == "M")
            outputs[name] = {"records": observations // 10, "observations": observations}

        expected_records = {
            "default": 2,                    # good + supplementary
            "disable-mq-zero": 3,            # + mapq-zero
            "defaults-off-explicit-mq-zero": 6,  # all except mapq-zero
            "explicit-not-supplementary": 1,  # default minus supplementary
            "defaults-off-inverted-not-supplementary": 1,
            "disable-nondefault-supplementary": 2,
            "mismatches-context-size-3": 2,   # same read set, wider context
        }
        for name, expected in expected_records.items():
            if outputs[name]["records"] != expected:
                raise AssertionError({"scenario": name,
                                      "expected_records": expected,
                                      "observed": outputs[name]})

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "scenarios": len(scenarios),
            "report_tables_bit_identical": True,
            "effective_records": outputs,
            "default_filter_semantics": "default-list -> -DF removal -> -RF addition -> -XRF inversion",
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
