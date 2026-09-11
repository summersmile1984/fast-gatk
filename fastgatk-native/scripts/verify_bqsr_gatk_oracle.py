#!/usr/bin/env python3
"""Run a focused BaseRecalibrator/ApplyBQSR Java-vs-native oracle.

The broad ``verify_bqsr.py`` contract intentionally exercises many native
options.  This smaller oracle is cheap enough for every configured backend and
keeps one release-pinned semantic slice explicit: native BaseRecalibrator must
produce the same GATKReport covariate tables, and native ApplyBQSR must apply a
GATK-produced report without changing the Java-decoded SAM records.  The
ApplyBQSR comparison includes the OQ path, where the input qualities are read
from the original-quality tag before recalibration.
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
    return subprocess.run(command, check=True, text=True,
                          capture_output=capture)


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


def sam_records(path: Path) -> list[str]:
    return [line for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("@")]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native_build = Path(os.environ.get("FASTGATK_NATIVE_BUILD",
                                     root / "fastgatk-native/build"))
    bqsr = Path(os.environ.get("FASTGATK_BQSR_BINARY",
                               native_build / "fastgatk-bqsr"))
    apply = Path(os.environ.get("FASTGATK_APPLY_BQSR_BINARY",
                               native_build / "fastgatk-apply-bqsr"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (bqsr, apply, java, gatk, bam, reference)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_bqsr_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled BaseRecalibrator/ApplyBQSR oracle inputs are required")
        print(json.dumps({"status": "skip",
                          "reason": "bundled GATK BQSR oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-bqsr-gatk-oracle-") as directory:
        work = Path(directory)
        region = "17:69000-69100"
        known_sites = work / "known-sites.vcf"
        gatk_report = work / "gatk.recal.tsv"
        native_report = work / "native.recal.tsv"
        gatk_bam = work / "gatk.recal.bam"
        native_bam = work / "native.recal.bam"
        gatk_oq_bam = work / "gatk.oq.bam"
        native_oq_bam = work / "native.oq.bam"
        gatk_sam = work / "gatk.sam"
        native_sam = work / "native.sam"
        gatk_oq_sam = work / "gatk.oq.sam"
        native_oq_sam = work / "native.oq.sam"

        known_sites.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=17,length=1000000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "17\t69067\t.\tT\tG\t.\tPASS\t.\n",
            encoding="utf-8")
        run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(known_sites)])

        run([str(java), "-jar", str(gatk), "BaseRecalibrator", "-R", str(reference),
             "-I", str(bam), "-L", region, "--known-sites", str(known_sites),
             "-O", str(gatk_report)])
        native_result = run([str(bqsr), "-I", str(bam), "-R", str(reference),
                             "-L", region, "--known-sites", str(known_sites),
                             "-O", str(native_report)], capture=True)
        native_summary = json.loads(native_result.stdout.strip().splitlines()[-1])
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
                }
            })

        # First compare ordinary recalibration.  PrintReads uses the same
        # HTSJDK decoder for both outputs, avoiding BAM block/compression
        # differences while retaining every SAM record and quality character.
        run([str(java), "-jar", str(gatk), "ApplyBQSR", "-R", str(reference),
             "-I", str(bam), "--bqsr-recal-file", str(gatk_report),
             "--create-output-bam-index", "false", "-O", str(gatk_bam)])
        native_apply = run([str(apply), "-I", str(bam), "-R", str(reference),
                            "--bqsr-recal-file", str(gatk_report),
                            "--create-output-bam-index=false", "-O", str(native_bam)],
                           capture=True)
        native_apply_summary = json.loads(native_apply.stdout.strip().splitlines()[-1])
        for source, output in ((gatk_bam, gatk_sam), (native_bam, native_sam)):
            run([str(java), "-jar", str(gatk), "PrintReads", "-I", str(source),
                 "--create-output-bam-index", "false", "-O", str(output)])
        ordinary_gatk = sam_records(gatk_sam)
        ordinary_native = sam_records(native_sam)
        if ordinary_gatk != ordinary_native:
            raise AssertionError("ApplyBQSR ordinary Java/native SAM records differ")

        # Repeat through GATK's --use-original-qualities path.  This verifies
        # that OQ is retained and used as the transformer input, not merely
        # copied after ordinary recalibration.
        run([str(java), "-jar", str(gatk), "ApplyBQSR", "-R", str(reference),
             "-I", str(bam), "--bqsr-recal-file", str(gatk_report),
             "--use-original-qualities", "--create-output-bam-index", "false",
             "-O", str(gatk_oq_bam)])
        native_oq = run([str(apply), "-I", str(bam), "-R", str(reference),
                         "--bqsr-recal-file", str(gatk_report),
                         "--use-original-qualities", "--create-output-bam-index=false",
                         "-O", str(native_oq_bam)], capture=True)
        native_oq_summary = json.loads(native_oq.stdout.strip().splitlines()[-1])
        for source, output in ((gatk_oq_bam, gatk_oq_sam), (native_oq_bam, native_oq_sam)):
            run([str(java), "-jar", str(gatk), "PrintReads", "-I", str(source),
                 "--create-output-bam-index", "false", "-O", str(output)])
        oq_gatk = sam_records(gatk_oq_sam)
        oq_native = sam_records(native_oq_sam)
        if oq_gatk != oq_native:
            raise AssertionError("ApplyBQSR --use-original-qualities Java/native records differ")
        if len(oq_native) != len(ordinary_native) or not any("\tOQ:Z:" in line for line in oq_native):
            raise AssertionError("ApplyBQSR OQ output contract failed")

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "region": region,
            "report_tables_bit_identical": True,
            "report_tables": {table: len(native_tables[table]) for table in expected_tables},
            "apply_records_bit_identical": True,
            "apply_records_compared": len(ordinary_native),
            "apply_oq_records_bit_identical": True,
            "apply_oq_records_compared": len(oq_native),
            "apply_oq_present": True,
            "native_records": native_summary["records"],
            "native_apply_records": native_apply_summary["records"],
            "native_oq_records": native_oq_summary["records"],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
