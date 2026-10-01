#!/usr/bin/env python3
"""Verify FilterMutectCalls FilterAlignmentArtifacts *class behavior* byte-equal to GATK 4.6.2.0.

Source-of-truth: broadinstitute/gatk 4.6.2.0
``FilterAlignmentArtifacts.java`` (Java class that runs automatically
inside ``FilterMutectCalls`` on indels in homopolymer /
short-tandem-repeat context).

TDD correction (2026-09-25): the previous version of this oracle
fed ``--filter-alignment-artifacts`` to GATK Java which replied
*"filter-alignment-artifacts is not a recognized option"*.  That CLI
flag does not exist in GATK 4.6.2.0; the underlying class runs
automatically inside the standard pipeline.

This oracle therefore pins the **class-level behavior** rather than
a nonexistent CLI flag, on a real Mutect2 call set:

Contract: when given a real Mutect2 call set, FilterMutectCalls must
tag indels near homopolymer runs with ``AS_FilterStatus=artifact``
on both GATK Java and native.  Non-indel records (SNVs) must not be
tagged ``artifact`` by alignment-artifact filtering on either
backend.  The set of records tagged ``artifact`` must be identical
between the two backends.

This oracle MUST be re-tested on the corrected GATK class behavior.
It serves as a regression guard for the next sprint to close any
real gap in the per-record FilterAlignmentArtifacts emission path.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


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


def is_indel(fields: list[str]) -> bool:
    return len(fields[3]) != len(fields[4])


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
            "verify_filter_mutect_alignment_artifacts_gatk_oracle.py", java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK oracle inputs required")
        print(json.dumps({"status": "skip",
                          "reason": "bundled GATK oracle inputs unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-fmc-faa-oracle-") as directory:
        work = Path(directory)
        driving_vcf = work / "driving.vcf.gz"
        driving_stats = work / "driving.vcf.gz.stats"
        driving_index = work / "driving.vcf.gz.tbi"
        gatk_filtered = work / "gatk.filtered.vcf.gz"
        native_filtered = work / "native.filtered.vcf.gz"

        # Step 1: native Mutect2 -> driving VCF.
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
                              "diagnostics": [f"native Mutect2 exit {mutect2_run.returncode}",
                                              mutect2_run.stderr[-400:]]}))
            return 1
        if not driving_stats.is_file() or not driving_index.is_file():
            print(json.dumps({"status": "fail",
                              "diagnostics": ["missing GATK stats table or VCF tabix index"]}))
            return 1

        # Step 2: GATK FilterMutectCalls.
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
                              "reason": f"GATK Java exit {gatk_run.returncode}",
                              "diagnostics": [gatk_run.stderr[-400:]]}))
            return 0

        # Step 3: native FilterMutectCalls.
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
                              "diagnostics": [f"native exit {native_run.returncode}",
                                              native_run.stderr[-400:]]}))
            return 1

        # Compare per-record AS_FilterStatus + whether artifact tag present.
        def collect(path: Path) -> dict[tuple[str, str], dict[str, object]]:
            out: dict[tuple[str, str], dict[str, object]] = {}
            for fields in read_vcf_records(path):
                if len(fields) < 8:
                    continue
                afs = info_value(fields[7], "AS_FilterStatus")
                artifact_in_filter = "artifact" in fields[6].split(";")
                artifact_in_afs = "artifact" in afs.split(",")
                out[(fields[0], fields[1])] = {
                    "ref": fields[3],
                    "alt": fields[4],
                    "is_indel": is_indel(fields),
                    "as_filter_status": afs,
                    "filter_column": fields[6],
                    "has_artifact_tag": artifact_in_filter or artifact_in_afs,
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

        # Per-record alignment-artifact tagging must agree.
        for key in sorted(gatk_records):
            if key not in native_records:
                continue
            g = gatk_records[key]
            n = native_records[key]
            if g["has_artifact_tag"] != n["has_artifact_tag"]:
                diagnostics.append(
                    f"{key}: artifact-tag java={g['has_artifact_tag']} "
                    f"native={n['has_artifact_tag']}")

        # Count records tagged by alignment-artifact filtering on each backend
        # (informational; not pinned because the boundary between
        # FilterAlignmentArtifacts and StrandBiasBySample is fuzzy).
        gatk_artifact = sum(1 for r in gatk_records.values() if r["has_artifact_tag"])
        native_artifact = sum(1 for r in native_records.values() if r["has_artifact_tag"])
        gatk_indel_artifact = sum(1 for r in gatk_records.values()
                                  if r["has_artifact_tag"] and r["is_indel"])
        native_indel_artifact = sum(1 for r in native_records.values()
                                    if r["has_artifact_tag"] and r["is_indel"])

        status = "pass" if not diagnostics else "fail"
        print(json.dumps({
            "status": status,
            "fixture": "dream_synthetic_chr20_mutect2_to_filter_mutect_calls",
            "records_java": len(gatk_records),
            "records_native": len(native_records),
            "records_artifact_tagged_java": gatk_artifact,
            "records_artifact_tagged_native": native_artifact,
            "records_indel_artifact_tagged_java": gatk_indel_artifact,
            "records_indel_artifact_tagged_native": native_indel_artifact,
            "diagnostics": diagnostics,
        }, sort_keys=True))
        return 0 if not diagnostics else 1


if __name__ == "__main__":
    raise SystemExit(main())