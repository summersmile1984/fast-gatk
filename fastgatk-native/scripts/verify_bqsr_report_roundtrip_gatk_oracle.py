#!/usr/bin/env python3
"""Verify that a native BaseRecalibrator report is consumable by GATK.

The table-value oracle compares the four recalibration tables, but it does not
exercise the GATKReport fixed-width reader.  This gate closes that direct
replacement boundary: Java ApplyBQSR and native ApplyBQSR both consume the
native report, then the pinned Java reader decodes both outputs for an exact
SAM-record comparison.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


def run(command: list[str], *, capture: bool = False) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=True, text=True, capture_output=capture)


def sam_records(path: Path) -> list[str]:
    return [line for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("@")]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    build = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
    bqsr = Path(os.environ.get("FASTGATK_BQSR_BINARY", build / "fastgatk-bqsr"))
    apply = Path(os.environ.get("FASTGATK_APPLY_BQSR_BINARY", build / "fastgatk-apply-bqsr"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (bqsr, apply, java, gatk, bam, reference)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            missing = [str(path) for path in required if not path.is_file()]
            raise SystemExit(f"bundled BQSR report oracle inputs are required: {missing}")
        print(json.dumps({"status": "skip", "reason": "bundled BQSR report oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-bqsr-report-oracle-") as directory:
        work = Path(directory)
        known_sites = work / "known-sites.vcf"
        known_sites.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=17,length=1000000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "17\t69067\t.\tT\tG\t.\tPASS\t.\n",
            encoding="utf-8")
        run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(known_sites)])

        region = "17:69000-69100"
        native_report = work / "native.recal.tsv"
        run([str(bqsr), "-I", str(bam), "-R", str(reference), "-L", region,
             "--known-sites", str(known_sites), "-O", str(native_report)])

        # The Java side intentionally consumes the native report.  A failure
        # here is a real direct-replacement break, not a numerical tolerance.
        java_output = work / "java-from-native-report.bam"
        native_output = work / "native-from-native-report.bam"
        run([str(java), "-jar", str(gatk), "ApplyBQSR", "-R", str(reference),
             "-I", str(bam), "--bqsr-recal-file", str(native_report),
             "--create-output-bam-index", "false", "-O", str(java_output)])
        run([str(apply), "-I", str(bam), "-R", str(reference),
             "--bqsr-recal-file", str(native_report),
             "--create-output-bam-index=false", "-O", str(native_output)])

        java_sam = work / "java.sam"
        native_sam = work / "native.sam"
        for source, output in ((java_output, java_sam), (native_output, native_sam)):
            run([str(java), "-jar", str(gatk), "PrintReads", "-I", str(source),
                 "--create-output-bam-index", "false", "-O", str(output)])
        java_records = sam_records(java_sam)
        native_records = sam_records(native_sam)
        assert len(java_records) == len(native_records) == 493
        assert java_records == native_records

    print(json.dumps({
        "status": "pass",
        "gatk_version": "4.6.2.0",
        "report_producer": "native BaseRecalibrator",
        "report_consumer": ["GATK ApplyBQSR", "native ApplyBQSR"],
        "java_native_report_records_bit_identical": True,
        "records_compared": len(native_records),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
