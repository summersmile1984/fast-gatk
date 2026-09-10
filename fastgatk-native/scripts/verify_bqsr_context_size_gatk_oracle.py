#!/usr/bin/env python3
"""Pin report-owned ContextCovariate width across BaseRecalibrator/ApplyBQSR."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from collections import Counter
from pathlib import Path


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=True, text=True, capture_output=True)


def table_rows(path: Path, table: str) -> Counter[tuple[str, ...]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    marker = f"#:GATKTable:{table}:"
    start = next(index for index, line in enumerate(lines)
                 if line.startswith(marker)) + 2
    rows: Counter[tuple[str, ...]] = Counter()
    for line in lines[start:]:
        if line.startswith("#:GATKTable:"):
            break
        if line.strip():
            rows[tuple(line.split())] += 1
    return rows


def sam_records(path: Path) -> list[str]:
    return [line for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("@")]


def make_fixture(reference: Path, sam: Path, known: Path) -> None:
    genome = "".join(line.strip() for line in reference.read_text(encoding="utf-8").splitlines()
                     if not line.startswith(">"))
    start = 69000
    template = genome[start - 1:start - 1 + 50]
    substitutions = {"A": "C", "C": "A", "G": "T", "T": "G"}
    records: list[str] = []
    for index in range(200):
        bases = list(template)
        if index < 100:
            bases[20] = substitutions[bases[20]]
        records.append("\t".join((f"context-{index:03d}", "0", "17", str(start),
                                  "60", "50M", "*", "0", "0", "".join(bases),
                                  "?" * 50, "RG:Z:RG1")))
    sam.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:17\tLN:1000000\n"
        "@RG\tID:RG1\tSM:S1\tPL:ILLUMINA\n" + "\n".join(records) + "\n",
        encoding="utf-8")
    known.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=17,length=1000000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        encoding="utf-8")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    build = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
    base = Path(os.environ.get("FASTGATK_BQSR_BINARY", build / "fastgatk-bqsr"))
    apply = Path(os.environ.get("FASTGATK_APPLY_BQSR_BINARY",
                                build / "fastgatk-apply-bqsr"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not all(path.is_file() for path in (base, apply, java, gatk, reference)):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled BQSR context-size oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "GATK context-size oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-bqsr-context-size-") as directory:
        work = Path(directory)
        reads, known = work / "reads.sam", work / "known.vcf"
        java_report, native_report = work / "java.table", work / "native.table"
        make_fixture(reference, reads, known)
        run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(known)])
        common_java = [str(java), "-jar", str(gatk), "BaseRecalibrator",
                       "-R", str(reference), "-I", str(reads), "--known-sites", str(known),
                       "--mismatches-context-size", "3", "--verbosity", "ERROR"]
        common_native = [str(base), "-R", str(reference), "-I", str(reads),
                         "--known-sites", str(known), "--mismatches-context-size", "3"]
        run(common_java + ["-O", str(java_report)])
        native_base = run(common_native + ["-O", str(native_report)])
        tables = ("Arguments", "Quantized", "RecalTable0", "RecalTable1", "RecalTable2")
        java_tables = {table: table_rows(java_report, table) for table in tables}
        native_tables = {table: table_rows(native_report, table) for table in tables}
        if java_tables != native_tables:
            raise AssertionError({table: {
                "java_only": sum((java_tables[table] - native_tables[table]).values()),
                "native_only": sum((native_tables[table] - java_tables[table]).values())}
                for table in tables})

        java_bam, native_bam = work / "java.bam", work / "native.bam"
        run([str(java), "-jar", str(gatk), "ApplyBQSR", "-R", str(reference),
             "-I", str(reads), "--bqsr-recal-file", str(java_report),
             "--create-output-bam-index", "false", "-O", str(java_bam),
             "--verbosity", "ERROR"])
        native_apply = run([str(apply), "-R", str(reference), "-I", str(reads),
                            "--bqsr-recal-file", str(java_report),
                            "--create-output-bam-index=false", "-O", str(native_bam)])
        apply_manifest = json.loads(Path(str(native_bam) + ".manifest.json").read_text(
            encoding="utf-8"))
        manifest_context_size = apply_manifest.get("telemetry", {}).get(
            "mismatches_context_size")
        if manifest_context_size != 3:
            raise AssertionError("ApplyBQSR manifest lost report-owned context size")
        decoded: list[list[str]] = []
        for source, output in ((java_bam, work / "java.sam"),
                               (native_bam, work / "native.sam")):
            run([str(java), "-jar", str(gatk), "PrintReads", "-I", str(source),
                 "--create-output-bam-index", "false", "-O", str(output),
                 "--verbosity", "ERROR"])
            decoded.append(sam_records(output))
        if decoded[0] != decoded[1]:
            raise AssertionError("report-owned size-3 ContextCovariate changed ApplyBQSR QUAL")
        if len(decoded[0]) != 200 or not any(row.split("\t")[10] != "?" * 50
                                             for row in decoded[0]):
            raise AssertionError("context fixture did not exercise a non-trivial recalibrated quality")

        print(json.dumps({
            "status": "pass", "gatk_version": "4.6.2.0",
            "mismatches_context_size": 3,
            "manifest_context_size": manifest_context_size,
            "report_tables_bit_identical": True,
            "apply_records_bit_identical": True,
            "apply_records_compared": len(decoded[0]),
            "native_base_records": json.loads(native_base.stdout.strip().splitlines()[-1])["records"],
            "native_apply_records": json.loads(native_apply.stdout.strip().splitlines()[-1])["records"],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
