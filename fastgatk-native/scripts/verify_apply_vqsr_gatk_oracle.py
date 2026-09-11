#!/usr/bin/env python3
"""Pinned GATK 4.6.2.0 oracle for ApplyVQSR AS per-ALT recalibration.

The broad ApplyVQSR contract exercises native filters and tranche plumbing, but
its Java fixture is single-ALT.  GATK's real AS recalibration VCF is a set of
one-ALT records at a shared locus, each carrying scalar VQSLOD/culprit values.
This oracle locks that multi-ALT lookup and the Number=A output vectors.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=score>
##INFO=<ID=AS_VQSLOD,Number=A,Type=Float,Description=allele-score>
##INFO=<ID=AS_FilterStatus,Number=A,Type=String,Description=allele-status>
##INFO=<ID=AS_culprit,Number=A,Type=String,Description=allele-culprit>
##INFO=<ID=culprit,Number=1,Type=String,Description=culprit>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def semantic_records(path: Path) -> list[tuple[str, ...]]:
    records: list[tuple[str, ...]] = []
    with gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(
        "r", encoding="utf-8"
    ) as stream:
        for line in stream:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            info: dict[str, object] = {}
            if fields[7] != ".":
                for item in fields[7].split(";"):
                    key, separator, value = item.partition("=")
                    if not separator:
                        info[key] = True
                    elif key in {"AS_VQSLOD", "VQSLOD"}:
                        info[key] = tuple(float(token) for token in value.split(","))
                    else:
                        info[key] = value
            records.append(
                (fields[0], fields[1], fields[2], fields[3], fields[4], fields[5],
                 fields[6], json.dumps(info, sort_keys=True, separators=(",", ":")))
            )
    return records


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_APPLY_VQSR_BINARY",
        root / "fastgatk-native/build/fastgatk-apply-vqsr",
    ))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    required = (binary, java, gatk)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_apply_vqsr_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("ApplyVQSR Java oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-apply-vqsr-gatk-oracle-") as directory:
        work = Path(directory)
        input_vcf = work / "input.vcf"
        recal_vcf = work / "recal.vcf"
        input_vcf.write_text(
            HEADER
            + "chr1\t1\t.\tA\tG,T\t50\tPASS\t.\n"
            + "chr1\t2\t.\tC\tA,G\t50\tPASS\t.\n",
            encoding="utf-8",
        )
        # GATK AS recalibration records are one ALT per row.  Their scalar
        # culprit is copied into the corresponding AS_culprit output slot.
        recal_vcf.write_text(
            HEADER
            + "chr1\t1\t.\tA\tG\t.\tPASS\tVQSLOD=2.0;AS_VQSLOD=2.0;culprit=QD\n"
            + "chr1\t1\t.\tA\tT\t.\tPASS\tVQSLOD=-1.0;AS_VQSLOD=-1.0;culprit=MQ\n"
            + "chr1\t2\t.\tC\tA\t.\tPASS\tVQSLOD=-1.0;AS_VQSLOD=-1.0;culprit=MQ\n"
            + "chr1\t2\t.\tC\tG\t.\tPASS\tVQSLOD=1.0;AS_VQSLOD=1.0;culprit=QD\n",
            encoding="utf-8",
        )
        for path in (input_vcf, recal_vcf):
            indexed = run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(path)])
            if indexed.returncode != 0:
                raise RuntimeError(f"GATK IndexFeatureFile failed: {indexed.stderr[-2000:]}")

        common = [
            "--recal-file", str(recal_vcf),
            "--use-allele-specific-annotations", "true",
            "--lod-score-cutoff", "0", "--mode", "SNP",
        ]
        java_output = work / "java.vcf.gz"
        java_result = run([
            str(java), "-jar", str(gatk), "ApplyVQSR", "-V", str(input_vcf),
            *common, "-O", str(java_output), "--create-output-variant-index", "false",
        ])
        if java_result.returncode != 0:
            raise RuntimeError(f"GATK ApplyVQSR failed: {java_result.stderr[-3000:]}")

        native_output = work / "native.vcf.gz"
        native_result = run([
            str(binary), "-V", str(input_vcf), *common, "-O", str(native_output),
            "--create-output-variant-index", "false",
        ])
        if native_result.returncode != 0:
            raise RuntimeError(f"native ApplyVQSR failed: {native_result.stderr[-3000:]}")

        expected = semantic_records(java_output)
        actual = semantic_records(native_output)
        if actual != expected:
            raise AssertionError({"expected": expected, "actual": actual})
        if len(actual) != 2:
            raise AssertionError(f"expected two multi-ALT records, got {len(actual)}")
        if "AS_FilterStatus=PASS,LOW_VQSLOD" not in gzip.open(native_output, "rt", encoding="utf-8").read():
            raise AssertionError("first Number=A filter vector was not preserved")
        if "AS_culprit=QD,MQ" not in gzip.open(native_output, "rt", encoding="utf-8").read():
            raise AssertionError("per-ALT scalar culprit values were not projected")

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "multi_alt_records": len(actual),
            "number_a_vectors_exact": True,
            "per_alt_scalar_culprit_exact": True,
            "site_filter_most_lenient_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
