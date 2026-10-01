#!/usr/bin/env python3
"""Verify Mutect2 AS_* family byte-equal to GATK 4.6.2.0 at the value level.

The structural set-equality oracle
(``verify_mutect2_as_annotation_oracle.py``) only checks which AS_*
keys each backend emits.  This value-level oracle goes further:

For every AS_* key that both backends emit on a record, the value
must be byte-identical.  Native currently emits `AS_SB_TABLE`,
`AS_MBQ`, `AS_MFRL`; GATK on the dream_synthetic chr20 fixture
emits only `AS_SB_TABLE`.  Value-level parity therefore only
constrains `AS_SB_TABLE` today, but the gate is structured so
that any additional AS_* key GATK emits on a richer fixture
(e.g. hcc1143 chr20, NA12878 trio) immediately becomes part of
the byte-equal contract.

Per-record AS_SB_TABLE format (GATK 4.6.2.0):
    ``AS_SB_TABLE=ref_f,ref_r|alt_f,alt_r|...``
One ALT-allele pair per comma-joined segment; commas separate the
two ALT strands within a segment; pipes separate REF from ALT.

Tolerance:
    Float-format drift (e.g. ``0.00`` vs ``0.000``) is recorded but
    not pinned: HTSlib uses ``%.3f`` for AS_* values while GATK uses
    Java's ``Double.toString``.  Tightening this is a separate
    Tier-3 byte-level formatting task.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


AS_KEYS = (
    "AS_BaseQRankSum",
    "AS_FS",
    "AS_FisherStrand",
    "AS_InbreedingCoeff",
    "AS_MBQ",
    "AS_MFRL",
    "AS_MPOS",
    "AS_MQRankSum",
    "AS_QualByDepth",
    "AS_ReadPosRankSum",
    "AS_SB_TABLE",
    "AS_SOR",
    "AS_StrandBiasBySample",
    "AS_StrandOddsRatio",
)


def info_value(info_str: str, key: str) -> str:
    for kv in info_str.split(";"):
        if kv.startswith(f"{key}="):
            return kv.split("=", 1)[1]
    return ""


def read_vcf_records(path: Path) -> list[list[str]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    out: list[list[str]] = []
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if not line or line.startswith("#"):
                continue
            out.append(line.rstrip("\n").split("\t"))
    return out


def germline_resource_path(root: Path) -> str:
    candidates = [
        root / "testdata/real/dream_synthetic/dream3-chr20.vcf",
        root / "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/mutect/dream/vcfs/dream3-chr20.vcf",
    ]
    for c in candidates:
        if c.is_file():
            return str(c)
    return ""


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    reference = root / "testdata/downloads/reference/hs37d5.fa.gz"
    mutect2_binary = Path(os.environ.get(
        "FASTGATK_MUTECT2_BINARY",
        str(root / "fastgatk-native" / "build" / "fastgatk-mutect2")))
    tumor_bam = root / "testdata/real/dream_synthetic/chr20/tumor.bam"
    normal_bam = root / "testdata/real/dream_synthetic/chr20/normal.bam"
    tumor_sm = "synthetic.challenge.set1.tumor"
    normal_sm = "synthetic.challenge.set1.normal"

    required = (java, gatk, reference, mutect2_binary, tumor_bam, normal_bam)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified(
            "verify_mutect2_as_value_byte_equality_oracle.py", java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK oracle inputs required")
        print(json.dumps({"status": "skip",
                          "reason": "bundled GATK oracle inputs unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-as-value-oracle-") as directory:
        work = Path(directory)
        native_vcf = work / "native.mutect2.vcf.gz"
        gatk_vcf = work / "gatk.mutect2.vcf.gz"

        native_run = subprocess.run([
            str(mutect2_binary),
            "--input", str(tumor_bam),
            "--input", str(normal_bam),
            "--reference", str(reference),
            "--tumor-sample", tumor_sm,
            "--normal-sample", normal_sm,
            "--output", str(native_vcf),
            "--stats", str(work / "native.diag.json"),
            "--create-output-variant-index", "true",
        ], text=True, capture_output=True, check=False)
        if native_run.returncode != 0:
            print(json.dumps({"status": "fail",
                              "diagnostics": [f"native Mutect2 exit {native_run.returncode}",
                                              native_run.stderr[-400:]]}))
            return 1

        germline_arg = []
        germline_res = germline_resource_path(root)
        if germline_res:
            germline_arg = ["--germline-resource", germline_res]

        gatk_run = subprocess.run([
            str(java), "-Xmx2g", "-jar", str(gatk),
            "Mutect2",
            "-I", str(tumor_bam),
            "-I", str(normal_bam),
            "-R", str(reference),
            "-tumor", tumor_sm,
            "-normal", normal_sm,
            "-O", str(gatk_vcf),
            *germline_arg,
            "--tmp-dir", str(work),
            "--create-output-variant-index", "true",
        ], text=True, capture_output=True, check=False)
        if gatk_run.returncode != 0:
            print(json.dumps({"status": "skip",
                              "reason": f"GATK Java Mutect2 exit {gatk_run.returncode}",
                              "diagnostics": [gatk_run.stderr[-400:]]}))
            return 0

        def collect(path: Path) -> dict[tuple[str, str], dict[str, str]]:
            out: dict[tuple[str, str], dict[str, str]] = {}
            for fields in read_vcf_records(path):
                if len(fields) < 8:
                    continue
                info = fields[7]
                out[(fields[0], fields[1])] = {
                    k: info_value(info, k) for k in AS_KEYS
                }
            return out

        native_records = collect(native_vcf)
        gatk_records = collect(gatk_vcf)

        diagnostics: list[str] = []
        # Per-key, per-record value byte-equal where both backends emit.
        per_key_compared: dict[str, int] = {k: 0 for k in AS_KEYS}
        per_key_drift: dict[str, int] = {k: 0 for k in AS_KEYS}
        for pos in sorted(set(native_records) & set(gatk_records)):
            for k in AS_KEYS:
                gv = gatk_records[pos].get(k, "")
                nv = native_records[pos].get(k, "")
                if gv and nv:
                    per_key_compared[k] += 1
                    if gv != nv:
                        diagnostics.append(
                            f"{pos}: {k} java={gv!r} native={nv!r}")
                        per_key_drift[k] += 1

        status = "pass" if not diagnostics else "fail"
        print(json.dumps({
            "status": status,
            "fixture": "dream_synthetic_chr20_mutect2",
            "records_compared": min(len(native_records), len(gatk_records)),
            "as_keys_compared": {k: per_key_compared[k] for k in AS_KEYS
                                if per_key_compared[k] > 0},
            "as_keys_drift": {k: per_key_drift[k] for k in AS_KEYS
                              if per_key_drift[k] > 0},
            "tier": "Tier-1a",
            "diagnostics": diagnostics[:10],
            "diagnostics_truncated": len(diagnostics) > 10,
            "note": ("Value-level byte-equal gate.  Compares the value of "
                     "each AS_* key that both backends emit on the same "
                     "record.  On dream_synthetic chr20 only AS_SB_TABLE "
                     "is shared, so this oracle pins that field's value; "
                     "richer fixtures (hcc1143, NA12878 trio) would extend "
                     "the comparison to the remaining AS_* keys."),
        }, sort_keys=True))
        return 0 if not diagnostics else 1


if __name__ == "__main__":
    raise SystemExit(main())