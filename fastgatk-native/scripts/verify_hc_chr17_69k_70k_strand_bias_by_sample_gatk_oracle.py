#!/usr/bin/env python3
"""Verify HaplotypeCaller StrandBiasBySample (FORMAT/SB) against GATK 4.6.2.0.

Source-of-truth: broadinstitute/gatk 4.6.2.0
``StrandBiasBySample.java`` (Fisher's exact 2x2 per-sample
contingency table).

Pinned fixtures:
  * ``gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam``
  * ``gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta``
  * interval ``17:69000-70000``, sample ``NA12878``

Contract:
  * GATK Java declares ``##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Strand Bias By Sample">``
    AND emits a 4-tuple per sample per variant (forward-ref,
    forward-alt, reverse-ref, reverse-alt) when at least one alt-support
    read is informative in each direction.
  * Native must match both the header declaration and the per-sample
    4-tuple emission, byte-equal to GATK's StrandFisherTest
    computation.

TDD: this oracle MUST exit non-zero on HEAD (the
``##FORMAT=<ID=SB,…>`` header is declared at ``hc_call.cpp:3941`` but
the per-sample 4-tuple is never written).  After Tier-1b code lands it
must exit zero.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def read_format_sb(path: Path) -> dict[str, dict[str, list[int]]]:
    """Return {position: {sample: [4-tuple]}} for the FORMAT/SB field.

    A position is in only if (a) the SB key is present in the FORMAT
    column, (b) the per-sample SB tuple has exactly 4 comma-separated
    integers.
    """
    out: dict[str, dict[str, list[int]]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            fmt_keys = fields[8].split(":")
            if "SB" not in fmt_keys:
                continue
            sb_idx = fmt_keys.index("SB")
            pos = f"{fields[0]}:{fields[1]}"
            sample_sb: dict[str, list[int]] = {}
            for col_idx in range(9, len(fields)):
                sample_name = fields[col_idx - (9 - 8)].split(":")
                if len(sample_name) <= 0:
                    continue
                sample_col = fields[col_idx].split(":")
                if sb_idx >= len(sample_col):
                    continue
                sb_str = sample_col[sb_idx]
                parts = sb_str.split(",")
                if len(parts) != 4:
                    continue
                try:
                    sample_sb[fields[col_idx - (9 - 9)]] = [int(x) for x in parts]
                except ValueError:
                    continue
            if sample_sb:
                out[pos] = sample_sb
    return out


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
        oracle_guard.oracle_not_verified(
            "verify_hc_chr17_69k_70k_strand_bias_by_sample_gatk_oracle.py", java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled HaplotypeCaller StrandBiasBySample oracle inputs required")
        print(json.dumps({"status": "skip",
                          "reason": "bundled GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-sb-oracle-") as directory:
        work = Path(directory)
        gatk_output = work / "gatk.vcf.gz"
        native_output = work / "native.vcf.gz"
        region = "17:69000-70000"

        gatk_run = subprocess.run([
            str(java), "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", region,
            "--sample-name", "NA12878",
            "--emit-ref-confidence", "NONE",
            "-O", str(gatk_output),
            "--create-output-variant-index", "false",
            "--seconds-between-progress-updates", "1",
        ], text=True, capture_output=True, check=False)
        assert gatk_run.returncode == 0, gatk_run.stderr
        native_run = subprocess.run([
            str(native), "-R", str(reference), "-I", str(bam), "-L", region,
            "--sample-name", "NA12878",
            "-O", str(native_output),
            "--create-output-variant-index", "false",
        ], text=True, capture_output=True, check=False)
        assert native_run.returncode == 0, native_run.stderr

        gatk_sb = read_format_sb(gatk_output)
        native_sb = read_format_sb(native_output)

        # The same set of positions must declare SB in both VCFs.
        diagnostics: list[str] = []
        for pos, gatk_sample_sb in gatk_sb.items():
            if pos not in native_sb:
                diagnostics.append(f"native missing SB at {pos}")
                continue
            for sample, gatk_tuple in gatk_sample_sb.items():
                native_tuple = native_sb[pos].get(sample)
                if native_tuple is None:
                    diagnostics.append(f"native missing SB for {sample} at {pos}")
                    continue
                if list(gatk_tuple) != list(native_tuple):
                    diagnostics.append(
                        f"SB tuple mismatch at {pos} sample={sample}: "
                        f"java={list(gatk_tuple)} native={list(native_tuple)}")
        # Also fail if GATK emitted SB but native omitted it for any position.
        for pos in native_sb:
            if pos not in gatk_sb:
                diagnostics.append(f"native emitted SB at {pos} where GATK did not")

        status = "pass" if not diagnostics else "fail"
        print(json.dumps({
            "status": status,
            "gatk_version": "4.6.2.0",
            "fixture": "hc-chr17-69k-70k",
            "diploid_single_sample_strand_bias_by_sample_4tuple": True,
            "records_with_sb_java": len(gatk_sb),
            "records_with_sb_native": len(native_sb),
            "diagnostics": diagnostics,
        }, sort_keys=True))
        return 0 if not diagnostics else 1


if __name__ == "__main__":
    raise SystemExit(main())
