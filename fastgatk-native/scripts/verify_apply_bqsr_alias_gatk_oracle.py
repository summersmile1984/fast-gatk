#!/usr/bin/env python3
"""Verify ApplyBQSR's pinned ``-bqsr`` short alias.

GATK 4.6.2.0 exposes ``-bqsr`` as the short alias of
``--bqsr-recal-file``.  This is a direct-replacement CLI boundary: the
native tool must accept the alias and it must select the same report and
produce the same transformed records as the long spelling.  Keep this gate
separate from the broader BQSR numerical oracle so a parser regression cannot
be hidden by a report/table comparison.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def run(command: list[str], *, capture: bool = False) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=True, text=True, capture_output=capture)


def sam_records(path: Path) -> list[str]:
    return [line for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("@")]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    build = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
    apply = Path(os.environ.get("FASTGATK_APPLY_BQSR_BINARY",
                               build / "fastgatk-apply-bqsr"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (apply, java, gatk, bam, reference)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_apply_bqsr_alias_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            missing = [str(path) for path in required if not path.is_file()]
            raise SystemExit(f"bundled ApplyBQSR alias oracle inputs are required: {missing}")
        print(json.dumps({"status": "skip",
                          "reason": "bundled GATK ApplyBQSR alias oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-apply-bqsr-alias-oracle-") as directory:
        work = Path(directory)
        known_sites = work / "known-sites.vcf"
        recal = work / "gatk.recal.tsv"
        known_sites.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=17,length=1000000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "17\t69067\t.\tT\tG\t.\tPASS\t.\n",
            encoding="utf-8")
        run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(known_sites)])
        run([str(java), "-jar", str(gatk), "BaseRecalibrator", "-R", str(reference),
             "-I", str(bam), "-L", "17:69000-69100", "--known-sites", str(known_sites),
             "-O", str(recal)])

        java_long = work / "java-long.bam"
        java_alias = work / "java-alias.bam"
        native_long = work / "native-long.bam"
        native_alias = work / "native-alias.bam"
        common_java = [str(java), "-jar", str(gatk), "ApplyBQSR", "-R", str(reference),
                       "-I", str(bam), "--create-output-bam-index", "false"]
        run(common_java + ["--bqsr-recal-file", str(recal), "-O", str(java_long)])
        run(common_java + ["-bqsr", str(recal), "-O", str(java_alias)])
        common_native = [str(apply), "-I", str(bam), "-R", str(reference),
                         "--create-output-bam-index=false"]
        run(common_native + ["--bqsr-recal-file", str(recal), "-O", str(native_long)])
        run(common_native + ["-bqsr", str(recal), "-O", str(native_alias)])

        # HTSJDK decoding compares the semantic records while ignoring BAM
        # compression and the tool-generated @PG header text.
        java_long_sam = work / "java-long.sam"
        java_alias_sam = work / "java-alias.sam"
        native_long_sam = work / "native-long.sam"
        native_alias_sam = work / "native-alias.sam"
        for source, output in ((java_long, java_long_sam),
                               (java_alias, java_alias_sam),
                               (native_long, native_long_sam),
                               (native_alias, native_alias_sam)):
            run([str(java), "-jar", str(gatk), "PrintReads", "-I", str(source),
                 "--create-output-bam-index", "false", "-O", str(output)])
        records = {"java_long": sam_records(java_long_sam),
                   "java_alias": sam_records(java_alias_sam),
                   "native_long": sam_records(native_long_sam),
                   "native_alias": sam_records(native_alias_sam)}
        if not records["java_alias"]:
            raise AssertionError("ApplyBQSR alias oracle produced no records")
        if records["java_alias"] != records["java_long"]:
            raise AssertionError("GATK -bqsr and --bqsr-recal-file differ")
        if records["native_alias"] != records["native_long"]:
            raise AssertionError("native -bqsr and --bqsr-recal-file differ")
        if records["native_alias"] != records["java_alias"]:
            raise AssertionError("native -bqsr records differ from pinned GATK")

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "alias": "-bqsr",
            "long_option": "--bqsr-recal-file",
            "java_alias_long_records_bit_identical": True,
            "native_alias_long_records_bit_identical": True,
            "native_java_alias_records_bit_identical": True,
            "records_compared": len(records["native_alias"]),
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
