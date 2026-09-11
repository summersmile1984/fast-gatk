#!/usr/bin/env python3
"""Pin BaseRecalibrator's long-read CycleCovariate boundary.

GATK's CycleCovariate is not an unbounded integer map: the default maximum
absolute cycle is 500 and a longer read fails, while an explicit
``--maximum-cycle-value`` permits the larger key domain.  The native path
previously had no way to select that domain and also omitted inserted read
bases from the substitution quality table.  This oracle uses a small SAM
fixture containing a 1,500-base read and a CIGAR insertion, checks the default
fail-closed boundary, then compares all materialized GATK report tables for a
larger configured maximum on Java and native implementations.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
import oracle_guard


def run(command: list[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=check, text=True, capture_output=True)


def table_rows(path: Path, table: str) -> Counter[tuple[str, ...]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    marker = f"#:GATKTable:{table}:"
    try:
        start = next(index for index, line in enumerate(lines)
                     if line.startswith(marker)) + 2
    except StopIteration:
        return Counter()
    rows: Counter[tuple[str, ...]] = Counter()
    for line in lines[start:]:
        if line.startswith("#:GATKTable:"):
            break
        if line.strip():
            rows[tuple(line.split())] += 1
    return rows


def make_fixture(reference: Path, sam: Path, known_sites: Path) -> None:
    sequence = "".join(line.strip() for line in reference.read_text(encoding="utf-8").splitlines()
                        if not line.startswith(">"))

    def record(name: str, pos: int, cigar: str, bases: str) -> str:
        return "\t".join((name, "0", "17", str(pos), "60", cigar, "*", "0", "0",
                           bases, "I" * len(bases), "RG:Z:LONG"))

    long_start = 69000
    long_bases = sequence[long_start - 1:long_start - 1 + 1500]
    insertion_start = 70600
    insertion_ref = sequence[insertion_start - 1:insertion_start - 1 + 98]
    insertion_bases = insertion_ref[:50] + "GG" + insertion_ref[50:]
    header = ("@HD\tVN:1.6\tSO:coordinate\n"
              "@SQ\tSN:17\tLN:1000000\n"
              "@RG\tID:LONG\tSM:S1\tPL:ILLUMINA\n")
    sam.write_text(header + record("long1500", long_start, "1500M", long_bases) + "\n" +
                   record("insertion", insertion_start, "50M2I48M", insertion_bases) + "\n",
                   encoding="utf-8")
    known_sites.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=17,length=1000000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "17\t69010\t.\tA\tG\t.\tPASS\t.\n", encoding="utf-8")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native_build = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
    native = Path(os.environ.get("FASTGATK_BQSR_BINARY", native_build / "fastgatk-bqsr"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (native, java, gatk, reference)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_bqsr_long_read_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("long-read BQSR oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK long-read oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-bqsr-long-gatk-oracle-") as directory:
        work = Path(directory)
        sam = work / "long.sam"
        known = work / "known.vcf"
        make_fixture(reference, sam, known)
        run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(known)])

        common_java = [str(java), "-jar", str(gatk), "BaseRecalibrator", "-R", str(reference),
                       "-I", str(sam), "--known-sites", str(known), "--verbosity", "ERROR"]
        common_native = [str(native), "-R", str(reference), "-I", str(sam),
                         "--known-sites", str(known)]

        java_default = run(common_java + ["-O", str(work / "java-default.table")], check=False)
        native_default = run(common_native + ["-O", str(work / "native-default.table")], check=False)
        if java_default.returncode == 0 or native_default.returncode == 0:
            raise AssertionError("default maximum cycle must reject a 1,500-base read")

        java_report = work / "java-long.table"
        native_report = work / "native-long.table"
        java_long = run(common_java + ["--maximum-cycle-value", "1500", "-O", str(java_report)])
        native_long = run(common_native + ["--maximum-cycle-value", "1500", "-O", str(native_report)])
        native_summary = json.loads(native_long.stdout.strip().splitlines()[-1])
        expected_tables = ("Quantized", "RecalTable0", "RecalTable1", "RecalTable2")
        java_tables = {table: table_rows(java_report, table) for table in expected_tables}
        native_tables = {table: table_rows(native_report, table) for table in expected_tables}
        if native_tables != java_tables:
            raise AssertionError({
                table: {"java_only": sum((java_tables[table] - native_tables[table]).values()),
                        "native_only": sum((native_tables[table] - java_tables[table]).values())}
                for table in expected_tables
            })

        # The short alias selects the same larger CycleCovariate domain.  It
        # is intentionally compared through report tables, not just parser
        # success, so the alias cannot silently revert to the default model.
        native_alias_report = work / "native-alias.table"
        native_alias = run(common_native + ["--max-cycle", "1500", "-O", str(native_alias_report)])
        alias_tables = {table: table_rows(native_alias_report, table) for table in expected_tables}
        if alias_tables != native_tables:
            raise AssertionError("--max-cycle did not reproduce --maximum-cycle-value report")

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "long_read_length": 1500,
            "insertion_read_length": 100,
            "default_maximum_cycle_rejected_java": True,
            "default_maximum_cycle_rejected_native": True,
            "maximum_cycle_value": 1500,
            "report_tables_bit_identical": True,
            "report_table_rows": {table: sum(rows.values()) for table, rows in native_tables.items()},
            "short_alias_report_bit_identical": True,
            "native_summary": native_summary,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
