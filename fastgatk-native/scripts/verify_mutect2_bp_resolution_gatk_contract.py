#!/usr/bin/env python3
"""Pin tumor-only Mutect2 ``-ERC BP_RESOLUTION`` to GATK 4.6.2.0.

Unlike GVCF band compression, BP resolution emits one <NON_REF> record per
reference position.  This is a separate writer path: it must retain the INFO
column, GATK's FORMAT-header order, and the per-base AD/DP/TLOD values.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def without_provenance(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle
                if not line.startswith("##GATKCommandLine=")]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party" / "jdk17" / "bin" / "java"
    gatk = root / "third_party" / "gatk-package" / "gatk-4.6.2.0" / \
        "gatk-package-4.6.2.0-local.jar"
    binary = Path(os.environ.get(
        "FASTGATK_MUTECT2_BINARY", str(root / "fastgatk-native/build/fastgatk-mutect2")))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (java, gatk, binary, bam, reference)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_mutect2_bp_resolution_gatk_contract.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK Mutect2 BP-resolution inputs are required")
        print(json.dumps({"status": "skip", "reason": "oracle inputs unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-bp-resolution-") as directory:
        work = Path(directory)
        gatk_output = work / "gatk.bp.vcf.gz"
        native_output = work / "native.bp.vcf.gz"
        common = ["-R", str(reference), "-I", str(bam), "--tumor-sample", "NA12878",
                  "-L", "17:69000-69066", "-ERC", "BP_RESOLUTION",
                  "--create-output-variant-index", "false",
                  "--add-output-vcf-command-line", "false"]
        gatk_run = run([str(java), "-Xmx2g", "-jar", str(gatk), "Mutect2",
                        *common, "-O", str(gatk_output)])
        assert gatk_run.returncode == 0, gatk_run.stderr[-4000:]
        native_run = run([str(binary), *common, "-O", str(native_output)])
        assert native_run.returncode == 0, native_run.stderr[-4000:]
        expected = without_provenance(gatk_output)
        observed = without_provenance(native_output)
        assert observed == expected
        data = [line for line in observed if not line.startswith("#")]
        assert len(data) == 67, data
        assert all(line.split("\t")[4] == "<NON_REF>" for line in data), data
        assert [int(line.split("\t")[1]) for line in data] == list(range(69000, 69067))

    print(json.dumps({
        "status": "pass",
        "mode": "Mutect2 -ERC BP_RESOLUTION",
        "oracle": "GATK-4.6.2.0",
        "records": 67,
        "vcf_exact_without_execution_provenance": True,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
