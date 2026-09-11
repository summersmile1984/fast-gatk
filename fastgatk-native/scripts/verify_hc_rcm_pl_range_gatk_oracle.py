#!/usr/bin/env python3
"""Verify that HC RCM PL values retain the full GATK integer range.

The PairHMM candidate writer has a bounded PL envelope, but GATK's
ReferenceConfidenceModel does not cap PLs at 999.  This small BP_RESOLUTION
window is intentionally high-depth so the distinction is observable while
keeping the oracle independent of broader HC assembly differences.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def records(path: Path) -> list[tuple[str, ...]]:
    result: list[tuple[str, ...]] = []
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            values = dict(zip(fields[8].split(":"), fields[9].split(":")))
            result.append((
                fields[0], fields[1], fields[3], fields[4],
                values.get("GT", ""), values.get("AD", ""),
                values.get("DP", ""), values.get("GQ", ""),
                values.get("PL", ""),
            ))
    return result


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", str(root / "fastgatk-native/build/fastgatk-hc-call")))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (java, gatk, native, bam, reference)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_hc_rcm_pl_range_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled HC RCM PL-range oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    region = "17:69646-69647"
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-rcm-pl-range-") as directory:
        work = Path(directory)
        gatk_output = work / "gatk.bp.g.vcf.gz"
        native_output = work / "native.bp.g.vcf.gz"
        gatk_run = subprocess.run([
            str(java), "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", region,
            "--emit-ref-confidence", "BP_RESOLUTION", "-O", str(gatk_output),
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--seconds-between-progress-updates", "1",
            "--indel-size-to-eliminate-in-ref-model", "1",
        ], text=True, capture_output=True, check=False)
        assert gatk_run.returncode == 0, gatk_run.stderr
        native_run = subprocess.run([
            str(native), "-R", str(reference), "-I", str(bam), "-L", region,
            "-ERC", "BP_RESOLUTION", "-O", str(native_output),
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--min-depth", "1", "--min-alt-support", "1", "--threads", "2",
            "--indel-size-to-eliminate-in-ref-model", "1",
        ], text=True, capture_output=True, check=False)
        assert native_run.returncode == 0, native_run.stderr
        gatk_records = records(gatk_output)
        native_records = records(native_output)
        assert gatk_records == native_records, {
            "gatk_records": gatk_records, "native_records": native_records,
        }
        assert [int(row[1]) for row in native_records] == [69646, 69647]
        pls = [int(value) for row in native_records for value in row[8].split(",")]
        assert max(pls) > 999, pls
        print(json.dumps({
            "status": "pass", "gatk_version": "4.6.2.0",
            "mode": "BP_RESOLUTION", "indel_size_to_eliminate_in_ref_model": 1,
            "records": len(native_records), "sample_fields_exact": True,
            "pl_range_unclipped": True, "maximum_pl": max(pls),
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
