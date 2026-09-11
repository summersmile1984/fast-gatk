#!/usr/bin/env python3
"""Byte oracle for chr20 min-pruning=1 HaplotypeCaller gVCF semantics.

This real-data boundary combines three GATK HaplotypeCaller behaviors that
must remain coupled: a numerically monomorphic EventMap allele folds into an
RCM block, a retained EventMap allele can be concrete with GT=0/0 and a
nonzero full REF/ALT/<NON_REF> QUAL, and those concrete hom-ref rows take part
in the physical-phasing called-haplotype set.  The last condition makes the
nearby 10020679--10020681 calls intentionally unphased in GATK 4.6.2.0.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def rows(path: Path) -> list[list[str]]:
    return [line.rstrip("\n").split("\t")
            for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")]


def row_at(records: list[list[str]], position: int) -> list[str]:
    for record in records:
        if int(record[1]) == position:
            return record
    raise AssertionError(f"missing gVCF row at {position}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    reference = root / "fixtures/chr20/ref20mnp.fasta"
    bam = root / "fixtures/chr20/mnp.bam"
    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    gatk = root / (
        "third_party/gatk-package/gatk-4.6.2.0/"
        "gatk-package-4.6.2.0-local.jar")
    region = "20:10019901-10020710"

    if not native.exists() or not reference.exists() or not bam.exists():
        raise SystemExit("missing native build or chr20 min-pruning fixture")
    if not gatk.exists():
        oracle_guard.oracle_not_verified('verify_hc_chr20_min_pruning_gvcf_gatk_oracle.py', gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing GATK oracle jar: {gatk}")
        print(json.dumps({"status": "skipped", "reason": "GATK jar not present"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-chr20-min-pruning-gvcf-") as directory:
        work = Path(directory)
        gatk_gvcf = work / "gatk.g.vcf"
        native_gvcf = work / "native.g.vcf"
        manifest = work / "native.json"
        subprocess.run([
            java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", region,
            "-ERC", "GVCF", "--min-pruning", "1", "-O", str(gatk_gvcf),
            "--native-pair-hmm-threads", "1",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-R", str(reference), "-I", str(bam), "-L", region,
            "--emit-ref-confidence", "GVCF", "--min-pruning", "1", "--threads", "1",
            "-O", str(native_gvcf), "--output-manifest", str(manifest),
            "--add-output-vcf-command-line", "false",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        assert gatk_gvcf.read_bytes() == native_gvcf.read_bytes(), (
            "chr20 min-pruning=1 gVCF is not byte-identical with provenance disabled")
        records = rows(native_gvcf)
        assert len(records) == 135, f"unexpected GATK-compatible gVCF row count: {len(records)}"

        # A minuscule positive floating-point QUAL must not turn this
        # monomorphic EventMap allele into a concrete REF/G/<NON_REF> record.
        monomorphic = row_at(records, 10020152)
        assert monomorphic[3:8] == ["T", "<NON_REF>", ".", ".", "END=10020152"]

        # EventMap retention is structural: this row has GT=0/0 and AD ALT=0,
        # but its complete gVCF posterior carries QUAL 0.01.
        hom_ref_event = row_at(records, 10020493)
        assert hom_ref_event[3:6] == ["T", "A,<NON_REF>", "0.01"]
        assert hom_ref_event[8:] == ["GT:AD:DP:GQ:PL:SB",
                                     "0/0:1,0,0:1:3:0,3,37,3,37,37:0,1,0,0"]

        # The earlier concrete hom-ref rows participate in calledHaplotypes,
        # making this otherwise adjacent cluster unphasable.
        for position in (10020679, 10020680, 10020681):
            phased = row_at(records, position)
            assert "PGT" not in phased[8] and "PID" not in phased[8] and "PS" not in phased[8]
            assert phased[9].startswith("0/1:"), (
                f"GATK chr20 min-pruning=1 phase/GT mismatch at {position}: {phased[9]}")

        telemetry = json.loads(manifest.read_text(encoding="utf-8"))["telemetry"]
        assert telemetry["graph_min_pruning"] == 1
        assert telemetry["gvcf_requested"] is True
        print(json.dumps({
            "status": "pass",
            "release": "GATK 4.6.2.0",
            "region": region,
            "min_pruning": 1,
            "gvcf_records": len(records),
            "gvcf_byte_identical_without_provenance": True,
            "monomorphic_event_folded_to_rcm_block": True,
            "hom_ref_event_full_posterior_exact": True,
            "called_haplotype_phasing_exact": True,
        }, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
