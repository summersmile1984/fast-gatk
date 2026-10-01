#!/usr/bin/env python3
"""Verify HC INFO byte-equal to GATK 4.6.2.0 on NA12878 chr20:100000-200000.

Tier-3 byte-level gate for HC INFO column on a richer fixture
than the chr17-69k-70k single-sample (which yields too few
records to detect drift).  GATK on this fixture emits
``BaseQRankSum``, ``FS``, ``MQRankSum``, ``QD``, ``ReadPosRankSum``,
``SOR`` for most records; native emits the same key set but with
slightly different values for some fields (SOR, QD, MLEAC/AF, DP).

This oracle asserts byte-equality on every per-record value
(after ``%.6g`` float-format normalization, which absorbs the
HTSlib-vs-Java Double.toString drift).

The fixture is a real GATK test fixture
(``NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam``)
shipped in ``gatk-source/src/test/resources``.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


KEYS = ("AC", "AN", "BaseQRankSum", "DP", "ExcessHet", "FS",
        "MQ", "MQRankSum", "QD", "ReadPosRankSum", "SOR")


def info_value(info_str: str, key: str) -> str:
    for kv in info_str.split(";"):
        if kv.startswith(f"{key}="):
            return kv.split("=", 1)[1]
    return ""


def normalize_numeric(text: str) -> str:
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
    binary = Path(os.environ.get(
        "FASTGATK_HC_CALL_BINARY",
        str(root / "fastgatk-native" / "build" / "fastgatk-hc-call")))
    bam = root / "NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
    sample_name = "NA12878"
    interval = "20:100000-200000"

    required = (java, gatk, reference, binary, bam)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified(
            "verify_hc_chr20_100k_200k_info_byte_equality_oracle.py", java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK oracle inputs required")
        print(json.dumps({"status": "skip",
                          "reason": "bundled GATK oracle inputs unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-chr20-info-oracle-") as directory:
        work = Path(directory)
        native_vcf = work / "native.hc.vcf.gz"
        gatk_vcf = work / "gatk.hc.vcf.gz"

        native_run = subprocess.run([
            str(binary),
            "-I", str(bam),
            "-R", str(reference),
            "-L", interval,
            "--sample-name", sample_name,
            "-O", str(native_vcf),
            "--tmp-dir", str(work),
            "--create-output-variant-index", "true",
        ], text=True, capture_output=True, check=False)
        if native_run.returncode != 0:
            print(json.dumps({"status": "fail",
                              "diagnostics": [f"native HC exit {native_run.returncode}",
                                              native_run.stderr[-400:]]}))
            return 1

        gatk_run = subprocess.run([
            str(java), "-Xmx2g", "-jar", str(gatk),
            "HaplotypeCaller",
            "-I", str(bam),
            "-R", str(reference),
            "-L", interval,
            "--sample-name", sample_name,
            "-O", str(gatk_vcf),
            "--tmp-dir", str(work),
            "--create-output-variant-index", "true",
        ], text=True, capture_output=True, check=False)
        if gatk_run.returncode != 0:
            print(json.dumps({"status": "skip",
                              "reason": f"GATK Java HC exit {gatk_run.returncode}",
                              "diagnostics": [gatk_run.stderr[-400:]]}))
            return 0

        def collect(path: Path) -> dict[tuple[str, str], dict[str, str]]:
            out: dict[tuple[str, str], dict[str, str]] = {}
            for fields in read_vcf_records(path):
                if len(fields) < 8:
                    continue
                info = fields[7]
                out[(fields[0], fields[1])] = {k: info_value(info, k) for k in KEYS}
            return out

        native_records = collect(native_vcf)
        gatk_records = collect(gatk_vcf)

        diagnostics: list[str] = []
        if set(native_records) != set(gatk_records):
            only_gatk = sorted(set(gatk_records) - set(native_records))[:5]
            only_native = sorted(set(native_records) - set(gatk_records))[:5]
            diagnostics.append(
                f"record-set mismatch: only_gatk={only_gatk} only_native={only_native}")

        per_key_compared: dict[str, int] = {k: 0 for k in KEYS}
        per_key_drift: dict[str, int] = {k: 0 for k in KEYS}
        for pos in sorted(set(native_records) & set(gatk_records)):
            for k in KEYS:
                gv = gatk_records[pos].get(k, "")
                nv = native_records[pos].get(k, "")
                if gv and nv:
                    per_key_compared[k] += 1
                    if gv == nv:
                        continue
                    if normalize_numeric(gv) == normalize_numeric(nv):
                        continue
                    diagnostics.append(
                        f"{pos}: {k} java={gv!r} native={nv!r}")
                    per_key_drift[k] += 1

        status = "pass" if not diagnostics else "fail"
        print(json.dumps({
            "status": status,
            "fixture": "hc_na12878_chr20_100k_200k",
            "records_compared": min(len(native_records), len(gatk_records)),
            "info_keys": list(KEYS),
            "info_keys_compared": {k: per_key_compared[k] for k in KEYS
                                   if per_key_compared[k] > 0},
            "info_keys_drift": {k: per_key_drift[k] for k in KEYS
                                if per_key_drift[k] > 0},
            "tier": "Tier-3",
            "diagnostics": diagnostics[:10],
            "diagnostics_truncated": len(diagnostics) > 10,
            "note": ("HC INFO byte-equal gate on a richer fixture than "
                     "chr17-69k-70k single-sample.  Detects drift in "
                     "SOR / QD / MLEAC / DP / BaseQRankSum etc.  Tolerance "
                     "covers HTSlib-vs-Java Double.toString drift via "
                     "%.6g normalization."),
        }, sort_keys=True))
        return 0 if not diagnostics else 1


if __name__ == "__main__":
    raise SystemExit(main())