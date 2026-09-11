#!/usr/bin/env python3
"""Pin the no-tag BI/BD terminal quality to GATK 4.6.2.0.

GATK's missing BI/BD fallback starts at Q45.  Under the default conservative
PCR model only positions before the terminal read base are lowered; the Java
loop deliberately leaves the final insertion/deletion quality at Q45.  This
oracle checks the serialized Java PairHMM input and the native telemetry for
both the effective Q40 adjustment and the raw terminal Q45 value.  Explicit
BI/BD decoding is covered separately by the existing BQSR indel oracle and is
not changed by this bounded default-only fix.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def run(command: list[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=check)


def pairhmm_rows(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 7:
            raise AssertionError(f"unexpected PairHMM row width: {len(fields)}")
        rows.append(fields)
    return rows


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    build = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
    binary = Path(os.environ.get("FASTGATK_MUTECT2_BINARY", build / "fastgatk-mutect2"))
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (binary, java, gatk, bam, reference)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_pairhmm_default_indel_quality_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            missing = [str(path) for path in required if not path.is_file()]
            raise SystemExit(f"PairHMM default-quality oracle inputs are required: {missing}")
        print(json.dumps({"status": "skip", "reason": "bundled GATK PairHMM oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-pairhmm-default-quality-") as directory:
        work = Path(directory)
        java_pairhmm = work / "gatk.pairhmm.txt"
        java_vcf = work / "gatk.vcf.gz"
        native_vcf = work / "native.vcf.gz"
        native_stats = work / "native.stats.json"
        run([
            str(java), "-jar", str(gatk), "Mutect2", "-R", str(reference), "-I", str(bam),
            "--tumor-sample", "NA12878", "-L", "17:69000-70000", "-O", str(java_vcf),
            "--pair-hmm-results-file", str(java_pairhmm),
        ])
        rows = pairhmm_rows(java_pairhmm)
        if not rows:
            raise AssertionError("GATK emitted no PairHMM rows")
        # Columns are hap-bases, read-bases, read-qual, BI, BD, GCP, result.
        # This fixture has no BI/BD optional tags.  GATK fills the missing
        # arrays with raw Q45, then applies the conservative repeat model to
        # original-read coordinates before AssemblyRegion clipping.  A clipped
        # terminal may therefore retain a lower quality (Q20 here) because it
        # was an interior repeat position in the source read.
        terminal_insertions = {ord(row[3][-1]) - 33 for row in rows if row[3]}
        terminal_deletions = {ord(row[4][-1]) - 33 for row in rows if row[4]}
        expected_terminal = {20, 45}
        if terminal_insertions != expected_terminal or terminal_deletions != expected_terminal:
            raise AssertionError({
                "terminal_insertion_qualities": sorted(terminal_insertions),
                "terminal_deletion_qualities": sorted(terminal_deletions),
                "expected_terminal_qualities": sorted(expected_terminal),
            })
        if not any((ord(row[3][0]) - 33) < 45 for row in rows if row[3]):
            raise AssertionError("GATK fixture did not exercise conservative PCR adjustment")

        native_run = run([
            str(binary), "-R", str(reference), "-I", str(bam), "--tumor-sample", "NA12878",
            "-L", "17:69000-70000", "-O", str(native_vcf), "--stats", str(native_stats),
            "--min-depth", "1", "--min-alt-support", "1",
            "--phred-scaled-global-read-mismapping-rate", "45",
        ])
        stats = json.loads(native_stats.read_text(encoding="utf-8"))
        if stats["pairhmm_default_indel_quality"] != 40:
            raise AssertionError(stats["pairhmm_default_indel_quality"])
        if stats["pairhmm_default_indel_terminal_quality"] != 45:
            raise AssertionError(stats["pairhmm_default_indel_terminal_quality"])
        if stats["pairhmm_pairs"] <= 0 or stats["pairhmm_pcr_adjusted_positions"] <= 0:
            raise AssertionError({
                "pairhmm_pairs": stats["pairhmm_pairs"],
                "pcr_adjusted_positions": stats["pairhmm_pcr_adjusted_positions"],
            })

        # NONE is the explicit no-PCR path: both the effective and terminal
        # default are Q45.  This also ensures the new raw default is not hidden
        # behind the Conservative telemetry compatibility scalar.
        none_stats_path = work / "native-none.stats.json"
        run([
            str(binary), "-R", str(reference), "-I", str(bam), "--tumor-sample", "NA12878",
            "-L", "17:69000-70000", "-O", str(work / "native-none.vcf.gz"),
            "--stats", str(none_stats_path), "--pcr-indel-model", "NONE",
            "--min-depth", "1", "--min-alt-support", "1",
            "--phred-scaled-global-read-mismapping-rate", "45",
        ])
        none_stats = json.loads(none_stats_path.read_text(encoding="utf-8"))
        if (none_stats["pairhmm_default_indel_quality"],
                none_stats["pairhmm_default_indel_terminal_quality"]) != (45, 45):
            raise AssertionError({
                "default": none_stats["pairhmm_default_indel_quality"],
                "terminal": none_stats["pairhmm_default_indel_terminal_quality"],
            })

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "java_pairhmm_rows": len(rows),
            "java_terminal_indel_quality": sorted(expected_terminal),
            "native_conservative_effective_quality": 40,
            "native_conservative_terminal_quality": 45,
            "native_none_effective_quality": 45,
            "native_none_terminal_quality": 45,
            "native_pairhmm_pairs": stats["pairhmm_pairs"],
            "explicit_bi_bd_path_unchanged": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
