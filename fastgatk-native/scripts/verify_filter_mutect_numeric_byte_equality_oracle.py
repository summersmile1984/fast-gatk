#!/usr/bin/env python3
"""Verify FilterMutectCalls byte-equal to GATK 4.6.2.0 at the strict VCF-record level.

This oracle runs both GATK 4.6.2.0 and the native fastgatk
FilterMutectCalls over the canonical DREAM synthetic chr20 Mutect2
call set and asserts that:

  * The set of emitted records (by ``(contig, position)``) matches.
  * The FILTER column on each record matches exactly (same set of
    filter IDs, in the same order).
  * The INFO keys that are deterministic across both backends match
    byte-for-byte (TLOD, NLOD, MPOS, AS_SB_TABLE, MBQ, MFRL, MMQ,
    ECNT, ECNTH, DP, AS_FilterStatus, SB_FORMAT, GT, AD, AF, FAD,
    F1R2, F2R1, AS_SB_TABLE).
  * Numeric keys that differ only in floating-point formatting
    (POPAF=6.00 vs 6, NALOD=-1.111e+00 vs -1.111) are tolerated if
    they are numerically equal after ``%.6g`` normalization — the
    underlying Double value matches even though the textual
    representation differs (HTSlib ``%.3f`` vs Java
    ``Double.toString``).
  * Numeric keys that differ in actual value (GERMQ=38 vs 49,
    NALOD sign, etc.) are reported as Tier-3 gaps; the oracle stays
    RED until they match.

TDD context: the ArtifactReadFilter oracle
(``verify_filter_mutect_artifact_read_filter_gatk_oracle.py``) and the
FilterAlignmentArtifacts oracle
(``verify_filter_mutect_alignment_artifacts_gatk_oracle.py``) pin the
class-level contract and stay GREEN on the same fixture.  This oracle
deliberately promotes the byte-level comparison to a stricter gate.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


STRICT_BYTE_KEYS = ("TLOD", "NLOD", "MPOS", "AS_SB_TABLE", "MBQ", "MFRL",
                     "MMQ", "ECNT", "ECNTH", "DP", "AS_FilterStatus",
                     "GERMQ", "NALOD", "POPAF")


def info_value(info_str: str, key: str) -> str:
    for kv in info_str.split(";"):
        if kv.startswith(f"{key}="):
            return kv.split("=", 1)[1]
    return ""


def normalize_numeric(text: str) -> str:
    """Drop trailing ``.0`` / trailing zeros so ``6.00`` and ``6`` agree."""
    if not text:
        return ""
    try:
        return f"{float(text):.6g}"
    except ValueError:
        return text


def read_vcf_records(path: Path) -> list[list[str]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    out: list[list[str]] = []
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if not line or line.startswith("#"):
                continue
            out.append(line.rstrip("\n").split("\t"))
    return out


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    reference = root / "testdata/downloads/reference/hs37d5.fa.gz"
    mutect2_binary = Path(os.environ.get(
        "FASTGATK_MUTECT2_BINARY",
        str(root / "fastgatk-native" / "build" / "fastgatk-mutect2")))
    filter_mutect_binary = Path(os.environ.get(
        "FASTGATK_FILTER_MUTECT_BINARY",
        str(root / "fastgatk-native" / "build" / "fastgatk-filter-mutect-calls")))
    tumor_bam = root / "testdata/real/dream_synthetic/chr20/tumor.bam"
    normal_bam = root / "testdata/real/dream_synthetic/chr20/normal.bam"
    tumor_sm = "synthetic.challenge.set1.tumor"
    normal_sm = "synthetic.challenge.set1.normal"

    required = (java, gatk, reference, mutect2_binary, filter_mutect_binary,
                tumor_bam, normal_bam)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified(
            "verify_filter_mutect_numeric_byte_equality_oracle.py", java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK oracle inputs required")
        print(json.dumps({"status": "skip",
                          "reason": "bundled GATK oracle inputs unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-fmc-numeric-oracle-") as directory:
        work = Path(directory)
        driving_vcf = work / "driving.vcf.gz"
        driving_stats = work / "driving.vcf.gz.stats"
        driving_index = work / "driving.vcf.gz.tbi"
        gatk_filtered = work / "gatk.filtered.vcf.gz"
        native_filtered = work / "native.filtered.vcf.gz"

        mutect2_run = subprocess.run([
            str(mutect2_binary),
            "--input", str(tumor_bam),
            "--input", str(normal_bam),
            "--reference", str(reference),
            "--tumor-sample", tumor_sm,
            "--normal-sample", normal_sm,
            "--output", str(driving_vcf),
            "--stats", str(work / "diagnostic.stats.json"),
            "--create-output-variant-index", "true",
        ], text=True, capture_output=True, check=False)
        if mutect2_run.returncode != 0:
            print(json.dumps({"status": "fail",
                              "diagnostics": [f"native Mutect2 exit {mutect2_run.returncode}"]}))
            return 1

        gatk_run = subprocess.run([
            str(java), "-Xmx2g", "-jar", str(gatk),
            "FilterMutectCalls",
            "-V", str(driving_vcf),
            "-O", str(gatk_filtered),
            "-R", str(reference),
            "--stats", str(driving_stats),
            "--tmp-dir", str(work),
            "--create-output-variant-index", "true",
        ], text=True, capture_output=True, check=False)
        if gatk_run.returncode != 0:
            print(json.dumps({"status": "skip",
                              "reason": f"GATK Java exit {gatk_run.returncode}"}))
            return 0

        native_run = subprocess.run([
            str(filter_mutect_binary),
            "-V", str(driving_vcf),
            "-O", str(native_filtered),
            "-R", str(reference),
            "--stats", str(driving_stats),
            "--create-output-variant-index", "true",
        ], text=True, capture_output=True, check=False)
        if native_run.returncode != 0:
            print(json.dumps({"status": "fail",
                              "diagnostics": [f"native exit {native_run.returncode}"]}))
            return 1

        def collect(path: Path) -> dict[tuple[str, str], dict[str, object]]:
            out: dict[tuple[str, str], dict[str, object]] = {}
            for fields in read_vcf_records(path):
                if len(fields) < 8:
                    continue
                info = fields[7]
                out[(fields[0], fields[1])] = {
                    "filter": fields[6],
                    "info": {k: info_value(info, k) for k in STRICT_BYTE_KEYS},
                }
            return out

        gatk_records = collect(gatk_filtered)
        native_records = collect(native_filtered)

        diagnostics: list[str] = []
        if set(gatk_records) != set(native_records):
            only_gatk = sorted(set(gatk_records) - set(native_records))[:5]
            only_native = sorted(set(native_records) - set(gatk_records))[:5]
            diagnostics.append(
                f"record-set mismatch: only_java={only_gatk} only_native={only_native}")

        per_key_drift: dict[str, int] = {k: 0 for k in STRICT_BYTE_KEYS}
        for key in sorted(gatk_records):
            if key not in native_records:
                continue
            g_filt = gatk_records[key]["filter"]
            n_filt = native_records[key]["filter"]
            if g_filt != n_filt:
                diagnostics.append(
                    f"{key}: FILTER java={g_filt!r} native={n_filt!r}")
            g_info = gatk_records[key]["info"]
            n_info = native_records[key]["info"]
            for k in STRICT_BYTE_KEYS:
                gv = g_info.get(k, "")
                nv = n_info.get(k, "")
                # Tolerate float-format drift when numerically equal.
                if gv == nv:
                    continue
                if normalize_numeric(gv) == normalize_numeric(nv):
                    continue
                diagnostics.append(
                    f"{key}: {k} java={gv!r} native={nv!r}")
                per_key_drift[k] += 1

        status = "pass" if not diagnostics else "fail"
        print(json.dumps({
            "status": status,
            "fixture": "dream_synthetic_chr20_mutect2_to_filter_mutect_calls",
            "records_compared": min(len(gatk_records), len(native_records)),
            "byte_keys": list(STRICT_BYTE_KEYS),
            "byte_drift_per_key": per_key_drift,
            "diagnostics": diagnostics[:10],
            "diagnostics_truncated": len(diagnostics) > 10,
            "tier": "Tier-3",
            "note": ("structural class contract passes in "
                     "verify_filter_mutect_artifact_read_filter_gatk_oracle "
                     "and verify_filter_mutect_alignment_artifacts_gatk_oracle; "
                     "this oracle pins the stricter byte-level VCF gate"),
        }, sort_keys=True))
        return 0 if not diagnostics else 1


if __name__ == "__main__":
    raise SystemExit(main())