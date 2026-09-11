#!/usr/bin/env python3
"""Verify HaplotypeCaller repeatable k-mer options against GATK 4.6.2.0.

ReadThreadingAssembler's ``--kmer-size`` is a repeatable argument and its
release default is the ordered set ``[10,25]``.  The native graph API keeps a
legacy scalar k for source compatibility, so this oracle protects the CLI
boundary: default and explicitly repeated/reordered requests must be passed
to the same graph retry policy and produce the same VCF as the pinned Java
caller.  Graph telemetry additionally exposes the normalized request list.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def run(command: list[str]) -> None:
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stdout={result.stdout[-2000:]}\nstderr={result.stderr[-4000:]}"
        )


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / (
        "third_party/gatk-package/gatk-4.6.2.0/"
        "gatk-package-4.6.2.0-local.jar")
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not all(path.exists() for path in (native, java, gatk, bam, reference)):
        oracle_guard.oracle_not_verified('verify_hc_kmer_list_gatk_oracle.py', java, gatk)
        raise SystemExit("missing pinned HC k-mer-list oracle assets")

    common_java = [
        str(java), "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
        "-R", str(reference), "-I", str(bam), "-L", "17:69000-70000",
        "--native-pair-hmm-threads", "2",
        "--create-output-variant-index", "false",
        "--add-output-vcf-command-line", "false",
        "--seconds-between-progress-updates", "1",
    ]
    common_native = [
        str(native), "-R", str(reference), "-I", str(bam), "-L", "17:69000-70000",
        "--threads", "2", "--add-output-vcf-command-line", "false",
    ]
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-kmer-list-oracle-") as directory:
        work = Path(directory)
        # First check the release defaults.  GATK's ReadThreadingAssembler
        # defaults to [10,25], while native must expose that list in its
        # manifest even when the caller does not spell the option out.
        gatk_default = work / "gatk-default.vcf"
        native_default = work / "native-default.vcf"
        native_default_manifest = work / "native-default.json"
        run([*common_java, "-O", str(gatk_default)])
        run([*common_native, "-O", str(native_default),
             "--output-manifest", str(native_default_manifest)])
        if gatk_default.read_bytes() != native_default.read_bytes():
            raise AssertionError("default [10,25] k-mer VCF is not byte-identical")

        # Reordered repeated arguments must normalize like GATK's sorted Java
        # List<Integer>, while the historical scalar telemetry retains the
        # final explicit value for old consumers.
        repeated = ["--kmer-size", "25", "--kmer-size", "10"]
        gatk_repeated = work / "gatk-repeated.vcf"
        native_repeated = work / "native-repeated.vcf"
        native_repeated_manifest = work / "native-repeated.json"
        run([*common_java, *repeated, "-O", str(gatk_repeated)])
        run([*common_native, *repeated, "-O", str(native_repeated),
             "--output-manifest", str(native_repeated_manifest)])
        if gatk_repeated.read_bytes() != native_repeated.read_bytes():
            raise AssertionError("repeated/reordered [25,10] k-mer VCF differs from GATK")

        default_telemetry = json.loads(
            native_default_manifest.read_text(encoding="utf-8"))["telemetry"]
        repeated_telemetry = json.loads(
            native_repeated_manifest.read_text(encoding="utf-8"))["telemetry"]
        if default_telemetry.get("graph_kmer_sizes") != [10, 25]:
            raise AssertionError(
                f"native default k-mer list mismatch: {default_telemetry.get('graph_kmer_sizes')}")
        if repeated_telemetry.get("graph_kmer_sizes") != [10, 25]:
            raise AssertionError(
                f"native repeated k-mer list mismatch: {repeated_telemetry.get('graph_kmer_sizes')}")
        if repeated_telemetry.get("graph_kmer_size") != 10:
            raise AssertionError(
                "legacy scalar graph_kmer_size must retain the final explicit value")

        print(json.dumps({
            "status": "pass",
            "release": "GATK 4.6.2.0",
            "region": "17:69000-70000",
            "default_kmer_sizes": [10, 25],
            "reordered_explicit_kmer_sizes": [25, 10],
            "normalized_kmer_sizes": [10, 25],
            "default_vcf_byte_identical_without_provenance": True,
            "repeated_vcf_byte_identical_without_provenance": True,
            "legacy_scalar_after_reordered_request": 10,
            "selected_default_kmer": default_telemetry.get("graph_kmer_size_selected"),
            "selected_repeated_kmer": repeated_telemetry.get("graph_kmer_size_selected"),
        }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
