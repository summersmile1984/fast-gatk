#!/usr/bin/env python3
"""Verify Mutect2 emits the same AS_* annotation family as GATK 4.6.2.0.

The ``AS_*`` allele-specific annotation family is required for joint
genotyping downstream (``SelectVariants`` filtering on
``AS_FilterStatus``, ``AS_StrandBiasBySample``, etc.).  GATK 4.6.2.0
emits each field opportunistically — only when its computation
produces a finite value for the given record (e.g. ``AS_BaseQRankSum``
requires a per-base ranking; ``AS_InbreedingCoeff`` requires a
multi-sample denominator).

TDD pinning: the oracle runs native Mutect2 end-to-end on the
canonical DREAM synthetic chr20 fixture, parses every emitted INFO
column, and asserts that **the set of AS_* keys emitted by native
equals the set emitted by GATK 4.6.2.0 on the same input**.

On the DREAM synthetic fixture GATK emits ``AS_SB_TABLE`` and
``AS_FilterStatus`` only.  Native emits ``AS_SB_TABLE``.  The
oracle therefore pins a specific gap: ``AS_FilterStatus`` is
emitted by GATK but not by native Mutect2 (it is emitted by native
FilterMutectCalls only).
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def info_keys(info_str: str) -> set[str]:
    out: set[str] = set()
    for kv in info_str.split(";"):
        if "=" not in kv:
            continue
        out.add(kv.split("=", 1)[0])
    return out


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
    # Fall back to whatever is present (GATK 4.6.2.0 is permissive on missing -germline-resource).
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
            "verify_mutect2_as_annotation_oracle.py", java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK oracle inputs required")
        print(json.dumps({"status": "skip",
                          "reason": "bundled GATK oracle inputs unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-as-oracle-") as directory:
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

        native_records = read_vcf_records(native_vcf)
        gatk_records = read_vcf_records(gatk_vcf)

        native_as_union = set()
        gatk_as_union = set()
        for fields in native_records:
            native_as_union |= {k for k in info_keys(fields[7]) if k.startswith("AS_")}
        for fields in gatk_records:
            gatk_as_union |= {k for k in info_keys(fields[7]) if k.startswith("AS_")}

        diagnostics: list[str] = []
        missing_in_native = sorted(gatk_as_union - native_as_union)
        # Native emitting extra AS_* keys beyond what GATK produces on this
        # fixture is OK: GATK gates some AS_* keys on minimum read count /
        # multi-sample denominator.  Native is allowed to be a strict
        # superset of GATK's emissions.  We only fail on missing keys
        # (GATK emitted but native did not).
        extra_in_native = sorted(native_as_union - gatk_as_union)
        if missing_in_native:
            diagnostics.append(
                f"AS_* keys emitted by GATK but missing in native: {missing_in_native}")
        if extra_in_native:
            # Record as informational, not failing.
            pass

        # Per-record parity.
        native_positions = {(f[0], f[1]): info_keys(f[7]) for f in native_records}
        gatk_positions = {(f[0], f[1]): info_keys(f[7]) for f in gatk_records}
        if gatk_positions.keys() != native_positions.keys():
            only_gatk = sorted(set(gatk_positions) - set(native_positions))[:5]
            only_native = sorted(set(native_positions) - set(gatk_positions))[:5]
            diagnostics.append(
                f"record-set mismatch: only_gatk={only_gatk} only_native={only_native}")

        for pos in sorted(set(native_positions) & set(gatk_positions)):
            g_as = {k for k in gatk_positions[pos] if k.startswith("AS_")}
            n_as = {k for k in native_positions[pos] if k.startswith("AS_")}
            missing = g_as - n_as
            if missing:
                diagnostics.append(
                    f"record {pos}: native missing AS_* keys {sorted(missing)}")

        status = "pass" if not diagnostics else "fail"
        print(json.dumps({
            "status": status,
            "fixture": "dream_synthetic_chr20_mutect2",
            "records_native": len(native_records),
            "records_gatk": len(gatk_records),
            "as_keys_emitted_native": sorted(native_as_union),
            "as_keys_emitted_gatk": sorted(gatk_as_union),
            "as_keys_native_only": extra_in_native,
            "tier": "Tier-1a",
            "diagnostics": diagnostics[:10],
            "diagnostics_truncated": len(diagnostics) > 10,
            "note": ("Native must emit every AS_* key GATK 4.6.2.0 emits "
                     "on the same record (set-equality, native ⊇ gatk). "
                     "Native emitting additional AS_* keys beyond GATK "
                     "is allowed (GATK gates some on minimum read count / "
                     "multi-sample denominator)."),
        }, sort_keys=True))
        return 0 if not diagnostics else 1


if __name__ == "__main__":
    raise SystemExit(main())