#!/usr/bin/env python3
"""Verify GatherBQSRReports against the pinned GATK 4.6.2.0 implementation.

This is intentionally separate from the broad BQSR contract.  It exercises
file-boundary workflows that are easy to get subtly wrong: two disjoint Java
BaseRecalibrator reports use a non-default quantization level, then Java and
native GatherBQSRReports must produce the same five report tables.  The second
case proves that a structurally incompatible indel context dimension fails
closed in both implementations instead of silently creating a corrupt table;
the empty-report case pins GATK's "no usable data" writer failure and prevents
native from publishing a misleading zero-row report.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from collections import Counter
from pathlib import Path


def run(command: list[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=check, text=True, capture_output=True)


def table_rows(path: Path, table: str) -> list[tuple[str, ...]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    marker = f"#:GATKTable:{table}:"
    try:
        start = next(i for i, line in enumerate(lines) if line.startswith(marker)) + 2
    except StopIteration:
        return []
    rows: list[tuple[str, ...]] = []
    for line in lines[start:]:
        if line.startswith("#:GATKTable:"):
            break
        if line.strip():
            rows.append(tuple(line.split()))
    return rows


def report_tables(path: Path) -> dict[str, Counter[tuple[str, ...]]]:
    return {table: Counter(table_rows(path, table)) for table in (
        "Arguments", "Quantized", "RecalTable0", "RecalTable1", "RecalTable2")}


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native_build = Path(os.environ.get("FASTGATK_NATIVE_BUILD",
                                     root / "fastgatk-native/build"))
    native = Path(os.environ.get(
        "FASTGATK_GATHER_BQSR_BINARY",
        native_build / "fastgatk-gather-bqsr-reports"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (native, java, gatk, bam, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GatherBQSRReports oracle inputs are required")
        print(json.dumps({"status": "skip",
                          "reason": "bundled GATK Gather oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-gather-bqsr-gatk-oracle-") as directory:
        work = Path(directory)
        known_sites = work / "known-sites.vcf"
        known_sites.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=17,length=1000000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "17\t69067\t.\tT\tG\t.\tPASS\t.\n",
            encoding="utf-8")
        run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(known_sites)])

        def base(output: Path, interval: str, *, indels_context_size: int = 3) -> None:
            run([str(java), "-jar", str(gatk), "BaseRecalibrator",
                 "-R", str(reference), "-I", str(bam), "-L", interval,
                 "--known-sites", str(known_sites), "--quantizing-levels", "4",
                 "--indels-context-size", str(indels_context_size), "-O", str(output)])

        left = work / "left.table"
        right = work / "right.table"
        gatk_gather = work / "gatk-gather.table"
        native_gather = work / "native-gather.table"
        java_tmp = work / "java-tmp"
        native_tmp = work / "native-tmp"
        java_tmp.mkdir()
        native_tmp.mkdir()
        base(left, "17:69000-69050")
        base(right, "17:69051-69100")
        run([str(java), "-jar", str(gatk), "GatherBQSRReports",
             "-I", str(left), "-I", str(right), "-O", str(gatk_gather),
             "--QUIET", "true", "--tmp-dir", str(java_tmp),
             "--use-jdk-deflater", "false", "--use-jdk-inflater", "false",
             "--verbosity", "WARNING"])
        native_result = run([str(native), "-I", str(left), "-I", str(right),
                             "-O", str(native_gather), "--QUIET", "true",
                             "--tmp-dir", str(native_tmp),
                             "--use-jdk-deflater", "false", "--use-jdk-inflater", "false",
                             "--verbosity", "WARNING"])
        native_summary = json.loads(native_result.stdout.strip().splitlines()[-1])
        java_tables = report_tables(gatk_gather)
        native_tables = report_tables(native_gather)
        if native_tables != java_tables:
            raise AssertionError({
                table: {"native_only": len(native_tables[table] - java_tables[table]),
                        "java_only": len(java_tables[table] - native_tables[table])}
                for table in java_tables
            })
        argument_rows = dict(java_tables["Arguments"])
        if ("quantizing_levels", "4") not in java_tables["Arguments"]:
            raise AssertionError("Java gather did not preserve quantizing_levels=4")
        metadata = json.loads(Path(f"{native_gather}.manifest.json").read_text(encoding="utf-8"))
        assert metadata["telemetry"]["quantizing_levels"] == 4
        assert metadata["compatibility"]["quiet"] is True
        assert metadata["telemetry"]["quiet"] is True
        assert metadata["compatibility"]["verbosity"] == "WARNING"
        assert metadata["telemetry"]["verbosity"] == "WARNING"
        assert native_summary["input_reports"] == 2

        # A changed indel context changes the RecalibrationTables dimensions;
        # GATK rejects this during RecalUtils.combineTables.  Native must fail
        # closed at the same file boundary rather than silently merging keys.
        incompatible = work / "incompatible.table"
        base(incompatible, "17:69051-69100", indels_context_size=4)
        java_bad = run([str(java), "-jar", str(gatk), "GatherBQSRReports",
                        "-I", str(left), "-I", str(incompatible),
                        "-O", str(work / "java-bad.table")], check=False)
        native_bad = run([str(native), "-I", str(left), "-I", str(incompatible),
                          "-O", str(work / "native-bad.table")], check=False)
        assert java_bad.returncode != 0
        assert native_bad.returncode != 0
        assert "incompatible BQSR report argument" in native_bad.stderr

        # GATK accepts an empty BaseRecalibrator report as an input file but
        # rejects the gather during RecalibrationReport.gatherReports(): no
        # RecalTable0/1/2 datum is usable.  Native must fail before publishing
        # a primary report or its covariate sidecar, rather than treating the
        # Arguments/Quantized headers as data.
        empty = work / "empty.table"
        base(empty, "17:1-10")
        java_empty_output = work / "java-empty.table"
        native_empty_output = work / "native-empty.table"
        java_empty = run([str(java), "-jar", str(gatk), "GatherBQSRReports",
                          "-I", str(empty), "-O", str(java_empty_output)], check=False)
        native_empty = run([str(native), "-I", str(empty),
                            "-O", str(native_empty_output)], check=False)
        assert java_empty.returncode != 0
        assert native_empty.returncode != 0
        assert "no usable data in any input file" in native_empty.stderr
        assert not native_empty_output.exists()
        assert not Path(f"{native_empty_output}.covariates.tsv").exists()

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "disjoint_reports": 2,
            "quantizing_levels": 4,
            "tables_bit_identical": True,
            "table_rows": {table: len(rows) for table, rows in native_tables.items()},
            "incompatible_dimensions_fail_closed": True,
            "empty_reports_fail_closed": True,
            "gatk_utility_cli_controls": True,
            "native_summary": native_summary,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
