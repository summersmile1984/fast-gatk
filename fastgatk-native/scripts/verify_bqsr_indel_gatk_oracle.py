#!/usr/bin/env python3
"""Pin BaseRecalibrator's BI/BD indel-quality event model to GATK 4.6.2.0.

The ordinary BQSR oracle uses the shipped chr17 fixture, which has no BI/BD
tags and therefore only exercises the Q45 fallback.  This oracle materializes
the same reads as SAM, adds deterministic per-base BI/BD qualities, and
compares all four GATKReport tables.  A truncated tag is also required to
fail closed in the native Host reader instead of silently shifting offsets.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from collections import Counter
from pathlib import Path


def run(command: list[str], *, capture: bool = False,
        check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=check, text=True, capture_output=capture)


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
    build = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
    bqsr = Path(os.environ.get("FASTGATK_BQSR_BINARY", build / "fastgatk-bqsr"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (bqsr, java, gatk, bam, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            missing = [str(path) for path in required if not path.is_file()]
            raise SystemExit(f"bundled BQSR indel oracle inputs are required: {missing}")
        print(json.dumps({"status": "skip",
                          "reason": "bundled GATK BQSR indel oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-bqsr-indel-oracle-") as directory:
        work = Path(directory)
        source_sam = work / "source.sam"
        tagged_sam = work / "tagged.sam"
        malformed_sam = work / "malformed.sam"
        known_sites = work / "known-sites.vcf"
        gatk_report = work / "gatk.recal.tsv"
        native_report = work / "native.recal.tsv"
        malformed_report = work / "malformed.recal.tsv"

        # PrintReads gives us a valid SAM header and the exact pinned read
        # records without requiring pysam/samtools.  BI/BD are then replaced
        # (rather than appended) so the input contains one authoritative tag.
        run([str(java), "-jar", str(gatk), "PrintReads", "-I", str(bam),
             "--create-output-bam-index", "false", "-O", str(source_sam)])
        tagged_lines: list[str] = []
        malformed_lines: list[str] = []
        for line in source_sam.read_text(encoding="utf-8").splitlines():
            if not line or line.startswith("@"):
                tagged_lines.append(line)
                malformed_lines.append(line)
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                raise AssertionError(f"malformed source SAM record: {line[:120]}")
            sequence = fields[9]
            fields[11:] = [field for field in fields[11:]
                           if not field.startswith(("BI:", "BD:"))]
            # Four distinct qualities per tag make the test sensitive to
            # per-offset lookup instead of merely checking one aggregate.
            insertion = "".join(chr(33 + 10 + (index % 4) * 10)
                                 for index in range(len(sequence)))
            deletion = "".join(chr(33 + 11 + (index % 4) * 10)
                               for index in range(len(sequence)))
            tagged = fields + [f"BI:Z:{insertion}", f"BD:Z:{deletion}"]
            tagged_lines.append("\t".join(tagged))
            malformed = fields + [f"BI:Z:{insertion[:-1]}", f"BD:Z:{deletion}"]
            malformed_lines.append("\t".join(malformed))
        tagged_sam.write_text("\n".join(tagged_lines) + "\n", encoding="utf-8")
        malformed_sam.write_text("\n".join(malformed_lines) + "\n", encoding="utf-8")
        # Keep the generated input shape visible in a failing CI log.
        first_tagged = next(line for line in tagged_lines if line and not line.startswith("@"))
        if "BI:Z:" not in first_tagged or "BD:Z:" not in first_tagged:
            raise AssertionError("generated SAM did not carry BI/BD tags")

        known_sites.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=17,length=1000000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "17\t69067\t.\tT\tG\t.\tPASS\t.\n",
            encoding="utf-8")
        run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(known_sites)])

        run([str(java), "-jar", str(gatk), "BaseRecalibrator", "-R", str(reference),
             "-I", str(tagged_sam), "--known-sites", str(known_sites),
             "--indels", "-O", str(gatk_report)])
        native = run([str(bqsr), "-I", str(tagged_sam), "-R", str(reference),
                      "--known-sites", str(known_sites),
                      "--compute-indel-bqsr-tables", "-O", str(native_report)],
                     capture=True)
        native_summary = json.loads(native.stdout.strip().splitlines()[-1])

        expected_tables = ("Quantized", "RecalTable0", "RecalTable1", "RecalTable2")
        gatk_tables = {table: Counter(table_rows(gatk_report.read_text(encoding="utf-8"), table))
                       for table in expected_tables}
        native_tables = {table: Counter(table_rows(native_report.read_text(encoding="utf-8"), table))
                         for table in expected_tables}
        if native_tables != gatk_tables:
            raise AssertionError({
                "report_tables": {
                    table: (len(native_tables[table]), len(gatk_tables[table]))
                    for table in expected_tables
                },
                "missing": {
                    table: list((gatk_tables[table] - native_tables[table]).elements())[:8]
                    for table in expected_tables if gatk_tables[table] != native_tables[table]
                },
                "extra": {
                    table: list((native_tables[table] - gatk_tables[table]).elements())[:8]
                    for table in expected_tables if gatk_tables[table] != native_tables[table]
                },
                "recal_table1_missing": list((gatk_tables["RecalTable1"] - native_tables["RecalTable1"]).elements())[:20],
                "recal_table1_extra": list((native_tables["RecalTable1"] - gatk_tables["RecalTable1"]).elements())[:20],
            })

        # Confirm that both I and D tables contain the four non-default
        # qualities.  This catches a regression to the old all-Q45 path even
        # if a future oracle accidentally compares only row counts.
        indel_rows = table_rows(native_report.read_text(encoding="utf-8"), "RecalTable1")
        for event, expected in (("I", {10, 20, 30, 40}), ("D", {11, 21, 31, 41})):
            observed = {int(row[1]) for row in indel_rows if row[2] == event}
            if observed != expected:
                raise AssertionError({"event": event, "observed_qualities": sorted(observed)})

        # The public RecalibrationArgumentCollection exposes -ics as the
        # short spelling of --indels-context-size.  Direct replacement must
        # preserve that spelling while retaining the same report rows.
        gatk_alias_report = work / "gatk-alias.recal.tsv"
        native_alias_report = work / "native-alias.recal.tsv"
        run([str(java), "-jar", str(gatk), "BaseRecalibrator", "-R", str(reference),
             "-I", str(tagged_sam), "--known-sites", str(known_sites),
             "-ics", "4", "--indels", "-O", str(gatk_alias_report)])
        run([str(bqsr), "-I", str(tagged_sam), "-R", str(reference),
             "--known-sites", str(known_sites), "-ics", "4",
             "--compute-indel-bqsr-tables", "-O", str(native_alias_report)], capture=True)
        for table in expected_tables:
            assert Counter(table_rows(gatk_alias_report.read_text(encoding="utf-8"), table)) == \
                   Counter(table_rows(native_alias_report.read_text(encoding="utf-8"), table))

        # A present but truncated BI/BD tag is malformed for a GATKRead and
        # must not be accepted as a shorter flat span by the native reader.
        malformed_run = run([str(bqsr), "-I", str(malformed_sam), "-R", str(reference),
                             "--known-sites", str(known_sites),
                             "--compute-indel-bqsr-tables", "-O", str(malformed_report)],
                            capture=True, check=False)
        if malformed_run.returncode == 0 or "indel quality length" not in malformed_run.stderr:
            raise AssertionError("truncated BI/BD tag did not fail closed")

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "traversal": "whole tagged SAM",
            "report_tables_bit_identical": True,
            "report_tables": {table: len(native_tables[table]) for table in expected_tables},
            "indel_qualities": {"I": [10, 20, 30, 40], "D": [11, 21, 31, 41]},
            "truncated_tag_rejected": True,
            "short_indel_context_alias": True,
            "native_records": native_summary["records"],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
