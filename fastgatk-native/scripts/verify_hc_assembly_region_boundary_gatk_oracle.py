#!/usr/bin/env python3
"""Pinned GATK oracle for HC AssemblyRegion boundary/default semantics.

GATK 4.6.2.0 uses a 100-base assembly halo and a 50-base activity
probability propagation distance by default.  This boundary is normally
unobservable on a broad region, but becomes observable when the active-region
cap is exercised: the ownership of reads at a shard/AssemblyRegion edge can
change the PairHMM AD/DP evidence.  The native invocation must therefore
match GATK on a deliberately small ``--max-assembly-region-size`` window.

This is a bounded oracle, not a claim of complete HaplotypeCaller parity.  It
locks the active-region/trimmed-read boundary and leaves graph assembly and
other release-specific differences to their dedicated checks.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    gatk_jar = root / (
        "third_party/gatk-package/gatk-4.6.2.0/"
        "gatk-package-4.6.2.0-local.jar")
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    region = "17:69000-70000"
    max_region_size = "50"
    if not native.exists() or not bam.exists() or not reference.exists():
        raise SystemExit("missing native build or AssemblyRegion boundary fixture")
    if not gatk_jar.exists():
        oracle_guard.oracle_not_verified('verify_hc_assembly_region_boundary_gatk_oracle.py', gatk_jar)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit(f"missing GATK oracle jar: {gatk_jar}")
        print(json.dumps({"status": "skipped", "reason": "GATK jar not present"}))
        return 0

    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-assembly-boundary-oracle-") as directory:
        work = Path(directory)
        gatk_vcf = work / "gatk.vcf"
        native_vcf = work / "native.vcf"
        manifest = work / "native.json"
        common_java = [
            java, "-Xmx1g", "-jar", str(gatk_jar), "HaplotypeCaller",
            "-R", str(reference), "-I", str(bam), "-L", region,
            "--max-assembly-region-size", max_region_size,
            "--native-pair-hmm-threads", "2",
            "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--seconds-between-progress-updates", "1",
        ]
        subprocess.run(common_java + ["-O", str(gatk_vcf)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            str(native), "-I", str(bam), "-R", str(reference), "-L", region,
            "-O", str(native_vcf), "--max-assembly-region-size", max_region_size,
            "--threads", "2", "--create-output-variant-index", "false",
            "--add-output-vcf-command-line", "false",
            "--output-manifest", str(manifest),
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # With command-line provenance disabled the bounded VCF is expected
        # to be byte-identical, including the schema/header and all PairHMM
        # fields.  This catches accidental drift in either writer or trim
        # boundary without weakening the oracle to call-set-only equality.
        assert gatk_vcf.read_bytes() == native_vcf.read_bytes(), (
            "max-size=50 HC VCF differs from pinned GATK output")
        records = [
            line.rstrip("\n").split("\t")
            for line in native_vcf.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")
        ]
        by_position = {int(record[1]): record for record in records}
        expected_positions = {69067, 69368, 69631}
        assert set(by_position) == expected_positions, (
            f"unexpected bounded HC call set: {sorted(by_position)}")
        assert by_position[69368][9].startswith("0/1:20,22:42:"), (
            "17:69368 must retain both ALT-bearing boundary reads")

        manifest_data = json.loads(manifest.read_text(encoding="utf-8"))
        telemetry = manifest_data["telemetry"]
        assert telemetry["max_assembly_region_size"] == 50
        assert telemetry["assembly_region_padding"] == 100
        assert telemetry["max_probability_propagation_distance"] == 50
        assert telemetry["assembly_region_partitioned"] is True
        assert telemetry["pairhmm_used"] is True

        print(json.dumps({
            "status": "pass",
            "release": "GATK 4.6.2.0",
            "region": region,
            "max_assembly_region_size": 50,
            "assembly_region_padding": telemetry["assembly_region_padding"],
            "max_probability_propagation_distance": telemetry[
                "max_probability_propagation_distance"],
            "call_set": sorted(expected_positions),
            "vcf_byte_identical_without_provenance": True,
            "boundary_alt_depth_exact": "17:69368 AD=20,22 DP=42",
            "scope": "AssemblyRegion active ownership, trim and PairHMM boundary",
            "remaining_scope": "full HC graph/assembly equivalence remains outside this bounded oracle",
        }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
