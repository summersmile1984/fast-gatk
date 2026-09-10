#!/usr/bin/env python3
"""Pinned GATK oracle for HCC1143 reads overhanging the truncated chr20 reference.

The shipped HCC1143 BAMs retain reads whose CIGAR ends past the artificial
1,000,000-base chr20 reference.  GATK's AssemblyRegionWalker traverses its
reference dictionary (not every CIGAR-projected coordinate), so this bounded
tail interval must complete and emit the same empty VCF on Java and native.
"""

from __future__ import annotations

import gzip
import os
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
JAVA = ROOT / "third_party" / "jdk17" / "bin" / "java"
GATK = ROOT / "third_party" / "gatk-package" / "gatk-4.6.2.0" / "gatk-package-4.6.2.0-local.jar"
NATIVE = Path(os.environ.get(
    "FASTGATK_MUTECT2_BINARY", str(ROOT / "fastgatk-native" / "build" / "fastgatk-mutect2")))
REFERENCE = ROOT / "testdata" / "real" / "cnv_somatic" / "human_g1k_v37.chr-20.truncated.fasta"
TUMOR = ROOT / "testdata" / "real" / "cnv_somatic" / "chr20" / "HCC1143_tumor.bam"
NORMAL = ROOT / "testdata" / "real" / "cnv_somatic" / "chr20" / "HCC1143_normal.bam"
TUMOR_SAMPLE = "HCC1143"
NORMAL_SAMPLE = "HCC1143 BL"
# This is one-based and inclusive at the GATK CLI boundary.  Its right edge
# is the truncated reference end, while selected reads extend beyond it.
INTERVAL = "20:999900-1000000"


def rows(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle if line and not line.startswith("#")]


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def main() -> int:
    required = (JAVA, GATK, NATIVE, REFERENCE, TUMOR, NORMAL)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise AssertionError(f"missing HCC1143 reference-boundary inputs: {missing}")
        print('{"status":"skip","reason":"HCC1143 reference-boundary inputs unavailable"}')
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-hcc1143-reference-boundary-") as directory:
        work = Path(directory)
        gatk_vcf = work / "gatk.vcf.gz"
        native_vcf = work / "native.vcf.gz"
        gatk = run([
            str(JAVA), "-jar", str(GATK), "Mutect2",
            "-R", str(REFERENCE), "-I", str(TUMOR), "-I", str(NORMAL),
            "-tumor", TUMOR_SAMPLE, "-normal", NORMAL_SAMPLE,
            "-L", INTERVAL, "-O", str(gatk_vcf),
        ])
        if gatk.returncode != 0:
            raise AssertionError(f"GATK Mutect2 failed:\n{gatk.stderr[-4000:]}")
        native = run([
            str(NATIVE), "-R", str(REFERENCE), "-I", str(TUMOR), "-I", str(NORMAL),
            "--tumor-sample", TUMOR_SAMPLE, "--normal-sample", NORMAL_SAMPLE,
            "-L", INTERVAL, "--native-pair-hmm-threads", "1", "-O", str(native_vcf),
        ])
        if native.returncode != 0:
            raise AssertionError(f"native Mutect2 failed:\n{native.stderr[-4000:]}")
        gatk_rows = rows(gatk_vcf)
        native_rows = rows(native_vcf)
        if gatk_rows != native_rows or gatk_rows:
            raise AssertionError({"gatk_rows": gatk_rows, "native_rows": native_rows})

    print('{"status":"pass","fixture":"hcc1143_reference_end","vcf_rows":0}')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
