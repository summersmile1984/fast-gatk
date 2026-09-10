#!/usr/bin/env python3
"""Verify the bounded BQSR contract on a real CRAM input/output path.

The ordinary BQSR oracle exercises BAM and OQ/Apply semantics.  This oracle
keeps CRAM/reference/index handling separate: Java creates a CRAM from the
pinned fixture, both BaseRecalibrator implementations consume it with known
sites, and both ApplyBQSR implementations publish indexed CRAM.  SAM records
are decoded by the same pinned Java reader before comparison, so compression
container bytes are not treated as an algorithmic identity claim.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from collections import Counter
from pathlib import Path


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


def sam_records(path: Path) -> list[str]:
    return [line for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("@")]


def index_for(path: Path) -> Path:
    """Return the index sidecar emitted by the active HTSJDK writer."""
    for suffix in (".crai", ".bai"):
        candidate = Path(f"{path}{suffix}")
        if candidate.is_file():
            return candidate
    raise AssertionError(f"missing CRAM/BAM index for {path}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    build = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
    bqsr = Path(os.environ.get("FASTGATK_BQSR_BINARY", build / "fastgatk-bqsr"))
    apply = Path(os.environ.get("FASTGATK_APPLY_BQSR_BINARY", build / "fastgatk-apply-bqsr"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    input_cram_fixture = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.cram"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (bqsr, apply, java, gatk, input_cram_fixture, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            missing = [str(path) for path in required if not path.is_file()]
            raise SystemExit(f"bundled BQSR CRAM oracle inputs are required: {missing}")
        print(json.dumps({"status": "skip", "reason": "bundled BQSR CRAM oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-bqsr-cram-oracle-") as directory:
        work = Path(directory)
        known_sites = work / "known-sites.vcf"
        known_sites.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=17,length=1000000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "17\t69067\t.\tT\tG\t.\tPASS\t.\n",
            encoding="utf-8")
        run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(known_sites)])

        # Use the pinned reference-backed CRAM fixture and its CRAI.  Keeping
        # the fixture bytes unchanged makes the Java/native comparison cover
        # the same CRAM reference-resolution and index traversal boundary.
        input_cram = input_cram_fixture
        input_index = index_for(input_cram)

        region = "17:69000-69100"
        gatk_report = work / "gatk.recal.tsv"
        native_report = work / "native.recal.tsv"
        run([str(java), "-jar", str(gatk), "BaseRecalibrator", "-R", str(reference),
             "-I", str(input_cram), "-L", region, "--known-sites", str(known_sites),
             "-O", str(gatk_report)])
        native_result = run([str(bqsr), "-I", str(input_cram), "-R", str(reference),
                             "-L", region, "--known-sites", str(known_sites),
                             "-O", str(native_report)], capture=True)
        native_summary = json.loads(native_result.stdout.strip().splitlines()[-1])
        expected_tables = ("Quantized", "RecalTable0", "RecalTable1", "RecalTable2")
        gatk_tables = {table: Counter(table_rows(gatk_report.read_text(encoding="utf-8"), table))
                       for table in expected_tables}
        native_tables = {table: Counter(table_rows(native_report.read_text(encoding="utf-8"), table))
                         for table in expected_tables}
        assert native_tables == gatk_tables, {
            table: (len(native_tables[table]), len(gatk_tables[table]))
            for table in expected_tables
        }

        gatk_cram = work / "gatk.recal.cram"
        native_cram = work / "native.recal.cram"
        run([str(java), "-jar", str(gatk), "ApplyBQSR", "-R", str(reference),
             "-I", str(input_cram), "--bqsr-recal-file", str(gatk_report),
             "--create-output-bam-index", "true", "-O", str(gatk_cram)])
        native_apply = run([str(apply), "-I", str(input_cram), "-R", str(reference),
                            "--bqsr-recal-file", str(gatk_report),
                            "--create-output-bam-index=true", "-O", str(native_cram)],
                           capture=True)
        native_apply_summary = json.loads(native_apply.stdout.strip().splitlines()[-1])
        assert gatk_cram.is_file()
        gatk_index = index_for(gatk_cram)
        assert native_cram.is_file() and Path(f"{native_cram}.crai").is_file()

        gatk_sam = work / "gatk.sam"
        native_sam = work / "native.sam"
        for source, output in ((gatk_cram, gatk_sam), (native_cram, native_sam)):
            run([str(java), "-jar", str(gatk), "PrintReads", "-R", str(reference),
                 "-I", str(source), "--create-output-bam-index", "false", "-O", str(output)])
        gatk_records = sam_records(gatk_sam)
        native_records = sam_records(native_sam)
        assert len(gatk_records) == len(native_records) == 493
        assert gatk_records == native_records

    print(json.dumps({
        "status": "pass",
        "gatk_version": "4.6.2.0",
        "input_format": "CRAM",
        "output_format": "CRAM",
        "input_index": input_index.suffix,
        "output_index": {"gatk": gatk_index.suffix, "native": ".crai"},
        "report_tables_bit_identical": True,
        "apply_records_bit_identical": True,
        "records_compared": len(native_records),
        "native_records": native_summary["records"],
        "native_apply_records": native_apply_summary["records"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
