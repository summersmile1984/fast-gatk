#!/usr/bin/env python3
"""Verify HaplotypeCaller InbreedingCoeff (INFO) against pinned GATK 4.6.2.0.

Source-of-truth: broadinstitute/gatk 4.6.2.0
``HaplotypeCaller.java`` + ``InbreedingCoeff.java`` (per-sample
heterozygosity estimate from `1 - sum_het/(n*(n-1)/2)`).

Pinned fixtures:
  * ``gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam``
  * ``gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta``
  * interval ``17:69000-70000``, sample ``NA12878``

Contract: GATK Java 4.6.2.0 declares
``##INFO=<ID=InbreedingCoeff,Number=1,…>`` in the VCF header but does
NOT emit ``InbreedingCoeff=X`` in any record body for a single-sample
run (the multi-sample denominator ``n*(n-1)/2`` is zero so the value is
uninformative; GATK omits the per-record key).  Native must match this
exact behavior byte-for-byte: header declared + per-record key absent.

For a future multi-sample fixture (≥2 samples), GATK emits the key
with the formula ``1 - sum_het/(n*(n-1)/2)`` per ``InbreedingCoeff.java``.
This oracle does NOT pin the multi-sample formula; it's pinned by the
joint-genotyping ``verify_mutect2_gatk_oracle.py``.

TDD: this oracle MUST exit non-zero against the current native binary
(the per-record InbreedingCoeff key is **over-emitted** as a literal
``InbreedingCoeff=0.0000`` placeholder at ``hc_call.cpp:3214, 5166`` —
which doesn't match GATK Java's actual absent-key behavior).
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def read_records(path: Path) -> list[tuple[str, ...]]:
    """Return ``[(chrom, pos, ref, alt, filter, has_inbreeding_coeff)]``.

    ``has_inbreeding_coeff`` is True iff the INFO field contains a
    key literally named ``InbreedingCoeff=...``.
    """
    rows: list[tuple[str, ...]] = []
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info_str = fields[7] if len(fields) > 7 else ""
            keys = set()
            for kv in info_str.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    keys.add(k)
                else:
                    keys.add(kv)
            rows.append((fields[0], fields[1], fields[3], fields[4], fields[6],
                         "InbreedingCoeff" in keys))
    return rows


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
            "verify_hc_chr17_69k_70k_inbreeding_coeff_gatk_oracle.py", java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled HaplotypeCaller InbreedingCoeff oracle inputs are required")
        print(json.dumps({"status": "skip",
                          "reason": "bundled GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-ic-oracle-") as directory:
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

        gatk_records = read_records(gatk_output)
        native_records = read_records(native_output)

        # Both runs emit the same set of variants (REF/ALT).  Native
        # must NOT emit ``InbreedingCoeff=`` for single-sample runs.
        diagnostics: list[str] = []
        assert len(gatk_records) == len(native_records), \
            f"GATK {len(gatk_records)} records vs native {len(native_records)}"
        for gatk_rec, native_rec in zip(gatk_records, native_records):
            if gatk_rec[:5] != native_rec[:5]:
                diagnostics.append(f"record mismatch: java={gatk_rec[:5]} native={native_rec[:5]}")
            if native_rec[5]:  # native emits InbreedingCoeff=X
                diagnostics.append(
                    f"native emits InbreedingCoeff key at {native_rec[0]}:{native_rec[1]} "
                    f"where GATK Java omits it for single-sample run")
        # Sanity: GATK should not emit the key either (denominator zero).
        for gatk_rec in gatk_records:
            if gatk_rec[5]:
                diagnostics.append(
                    f"GATK unexpectedly emits InbreedingCoeff at "
                    f"{gatk_rec[0]}:{gatk_rec[1]} — this changes the contract")

        status = "pass" if not diagnostics else "fail"
        print(json.dumps({
            "status": status,
            "gatk_version": "4.6.2.0",
            "fixture": "hc-chr17-69k-70k",
            "records": len(native_records),
            "diploid_single_sample_inbreeding_coefficient_emit_policy": "absent",
            "diagnostics": diagnostics,
        }, sort_keys=True))
        return 0 if not diagnostics else 1


if __name__ == "__main__":
    raise SystemExit(main())
