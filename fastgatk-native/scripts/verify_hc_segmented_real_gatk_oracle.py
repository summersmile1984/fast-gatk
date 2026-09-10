#!/usr/bin/env python3
"""Strict GATK oracle across adjacent real-BAM HaplotypeCaller intervals.

The regular chr17 fixture verifies a complete 1 kb traversal.  This oracle
splits that same span at three independent -L boundaries to make interval
selection, AssemblyRegion padding and active-core ownership observable.  It
compares the complete VCF data rows, intentionally excluding only producer
headers and execution provenance.
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
JAVA = ROOT / "third_party/jdk17/bin/java"
GATK = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
NATIVE = Path(os.environ.get(
    "FASTGATK_HC_BINARY", ROOT / "fastgatk-native/build/fastgatk-hc-call"))
BAM = ROOT / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
REFERENCE = ROOT / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"

# The fourth interval deliberately has no call.  Pinning both populated and
# empty spans catches a call that leaks across a requested output boundary.
SEGMENTS = (
    ("17:69000-69250", 1),
    ("17:69251-69500", 1),
    ("17:69501-69750", 1),
    ("17:69751-70000", 0),
)


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False,
                          env={**os.environ, "OMP_NUM_THREADS": "1"})


def rows(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle
                if line and not line.startswith("#")]


def main() -> int:
    required = (JAVA, GATK, NATIVE, BAM, REFERENCE)
    assert all(path.is_file() for path in required), required
    assert os.access(NATIVE, os.X_OK), NATIVE
    observed: dict[str, int] = {}
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-segmented-oracle-") as directory:
        work = Path(directory)
        for ordinal, (interval, expected_rows) in enumerate(SEGMENTS):
            gatk_vcf = work / f"gatk-{ordinal}.vcf.gz"
            native_vcf = work / f"native-{ordinal}.vcf.gz"
            gatk = run([
                str(JAVA), "-jar", str(GATK), "HaplotypeCaller",
                "-R", str(REFERENCE), "-I", str(BAM), "-L", interval,
                "-O", str(gatk_vcf), "--add-output-vcf-command-line", "false",
            ])
            assert gatk.returncode == 0, gatk.stderr
            native = run([
                str(NATIVE), "-R", str(REFERENCE), "-I", str(BAM), "-L", interval,
                "-O", str(native_vcf), "--add-output-vcf-command-line", "false",
                "--native-pair-hmm-threads", "1",
            ])
            assert native.returncode == 0, native.stderr
            gatk_rows = rows(gatk_vcf)
            native_rows = rows(native_vcf)
            assert len(gatk_rows) == expected_rows, (interval, gatk_rows)
            assert native_rows == gatk_rows, (interval, gatk_rows, native_rows)
            observed[interval] = len(gatk_rows)
    print(json.dumps({"status": "pass", "segments": observed,
                      "vcf_data_rows_exact": True}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
