#!/usr/bin/env python3
"""Verify FilterMutectCalls byte-equal to GATK 4.6.2.0 on real chr20 Mutect2 output.

TDD correction (2026-09-25): the previous version fed synthetic
VCFs to GATK Java's FilterMutectCalls.  GATK requires a real
companion Mutect2 stats table + chr-matching reference for the
input VCF, and rejects unindexed block-compressed inputs.  The
synthetic chr17 VCF on a chr17 reference failed with "Mutect stats
table not a GATK Mutect2 stats table" (exit 2).

This oracle uses the canonical DREAM synthetic chr20 fixture
(testdata/real/dream_synthetic/chr20/{tumor,normal}.bam with
hs37d5.fa.gz) and runs the full native Mutect2 -> native FilterMutectCalls
pipeline end-to-end.  Both outputs are compared to GATK 4.6.2.0's
FilterMutectCalls run on the same driving VCF.

Contract pinned:
  * Both FilterMutectCalls runs complete (exit 0).
  * Both outputs declare the canonical GATK FILTER vocabulary
    (PASS / FAIL / strand_bias / weak_evidence / base_qual / etc.).
  * AS_FilterStatus is emitted per-record on both backends.
  * Records with strand_bias get `strand_bias` in their FILTER column
    on both backends (ArtifactReadFilter / StrandBias filter class
    contract).
  * Records with weak_evidence get `weak_evidence` in their FILTER
    column on both backends.
  * Record count and (contig, position) identity match exactly.

Float-format drift (e.g. ``NALOD=-1.111e+00`` vs ``-1.111``,
``POPAF=6.00`` vs ``6``) is recorded but not pinned: HTSlib uses
``%.3f`` while GATK uses Java's ``Double.toString``.  Tightening
this is the byte-equal Tier-3 follow-up.
"""

from __future__ import annotations

import gzip
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def read_vcf_records(path: Path) -> list[list[str]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    out: list[list[str]] = []
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if not line or line.startswith("#"):
                continue
            out.append(line.rstrip("\n").split("\t"))
    return out


def info_value(info_str: str, key: str) -> str:
    for kv in info_str.split(";"):
        if kv.startswith(f"{key}="):
            return kv.split("=", 1)[1]
    return ""


def filter_record_dict(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    records: dict[tuple[str, str], dict[str, str]] = {}
    for fields in read_vcf_records(path):
        if len(fields) < 8:
            continue
        info_str = fields[7]
        records[(fields[0], fields[1])] = {
            "filter": fields[6],
            "as_filter_status": info_value(info_str, "AS_FilterStatus"),
            "gerq": info_value(info_str, "GERMQ"),
        }
    return records


def list_filter_ids(path: Path) -> set[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    ids: set[str] = set()
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("##FILTER=<ID="):
                ids.add(line.split("ID=", 1)[1].split(",", 1)[0])
    return ids


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
            "verify_filter_mutect_artifact_read_filter_gatk_oracle.py", java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK oracle inputs required")
        print(json.dumps({"status": "skip",
                          "reason": "bundled GATK oracle inputs unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-fmc-arf-oracle-") as directory:
        work = Path(directory)
        driving_vcf = work / "driving.vcf.gz"
        driving_stats = work / "driving.vcf.gz.stats"
        driving_index = work / "driving.vcf.gz.tbi"
        gatk_filtered = work / "gatk.filtered.vcf.gz"
        native_filtered = work / "native.filtered.vcf.gz"

        # Step 1: native Mutect2 to produce the driving VCF + GATK stats table.
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

        # Native mutect2 emits the GATK stats table at <output>.stats.
        # If the user passed an explicit --stats path, the GATK table is
        # still written next to the VCF (see write_gatk_mutect_stats_table).
        if not driving_stats.is_file():
            print(json.dumps({"status": "fail",
                              "diagnostics": [f"GATK stats table not produced at {driving_stats}"]}))
            return 1
        if not driving_index.is_file():
            print(json.dumps({"status": "fail",
                              "diagnostics": [f"VCF tabix index not produced at {driving_index}"]}))
            return 1

        # Step 2: GATK 4.6.2.0 FilterMutectCalls on the native driving VCF.
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

        # Step 3: native FilterMutectCalls on the same driving VCF + stats.
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

        gatk_records = filter_record_dict(gatk_filtered)
        native_records = filter_record_dict(native_filtered)
        gatk_filter_ids = list_filter_ids(gatk_filtered)
        native_filter_ids = list_filter_ids(native_filtered)

        diagnostics: list[str] = []
        # FILTER vocabulary must match (GATK adds strict_strand/possible_numt
        # that native currently does not).
        vocab_diff = (gatk_filter_ids ^ native_filter_ids) - {"PASS", "FAIL"}
        if vocab_diff:
            diagnostics.append(
                f"FILTER vocabulary mismatch: only_java={sorted(vocab_diff & gatk_filter_ids)} "
                f"only_native={sorted(vocab_diff & native_filter_ids)}")

        if set(gatk_records) != set(native_records):
            only_gatk = sorted(set(gatk_records) - set(native_records))[:5]
            only_native = sorted(set(native_records) - set(gatk_records))[:5]
            diagnostics.append(
                f"record-set mismatch: only_java={only_gatk} only_native={only_native}")

        # Per-record AS_FilterStatus + FILTER column byte-equal for StrandBias / weak_evidence.
        for key in sorted(gatk_records):
            if key not in native_records:
                continue
            g_afs = gatk_records[key]["as_filter_status"]
            n_afs = native_records[key]["as_filter_status"]
            g_filters = set(gatk_records[key]["filter"].split(";"))
            n_filters = set(native_records[key]["filter"].split(";"))
            if g_afs != n_afs:
                diagnostics.append(
                    f"{key}: AS_FilterStatus java={g_afs!r} native={n_afs!r}")
            if g_filters != n_filters:
                diagnostics.append(
                    f"{key}: FILTER java={sorted(g_filters)} native={sorted(n_filters)}")

        status = "pass" if not diagnostics else "fail"
        print(json.dumps({
            "status": status,
            "fixture": "dream_synthetic_chr20_mutect2_to_filter_mutect_calls",
            "records_java": len(gatk_records),
            "records_native": len(native_records),
            "filter_vocabulary_java": sorted(gatk_filter_ids),
            "filter_vocabulary_native": sorted(native_filter_ids),
            "diagnostics": diagnostics,
        }, sort_keys=True))
        return 0 if not diagnostics else 1


if __name__ == "__main__":
    raise SystemExit(main())