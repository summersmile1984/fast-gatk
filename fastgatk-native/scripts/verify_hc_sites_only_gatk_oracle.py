#!/usr/bin/env python3
"""Pinned GATK oracle for HaplotypeCaller sites-only VCF output.

The sites-only switch is a writer boundary: HC must still consume full
FORMAT/sample evidence for genotyping and reference-confidence computation,
then omit FORMAT/sample columns from the published VCF/GVCF.  This guard
compares site columns and validates bare/separated Boolean forms on the pinned
chr17 fixture.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


def records(path: Path) -> tuple[list[str], list[list[str]]]:
    if path.suffix == ".gz":
        stream_context = gzip.open(path, "rt", encoding="utf-8")
    else:
        stream_context = path.open("rt", encoding="utf-8")
    with stream_context as stream:
        lines = [line.rstrip("\n") for line in stream]
    header = next(line for line in lines if line.startswith("#CHROM"))
    return header.split("\t"), [
        line.split("\t") for line in lines if line and not line.startswith("#")
    ]


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not all(path.is_file() for path in (native, java, gatk, bam, reference)):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("HaplotypeCaller sites-only oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "native or pinned GATK assets unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-sites-only-oracle-") as directory:
        work = Path(directory)
        results: dict[str, dict[str, object]] = {}
        cases = {
            "vcf-true-separated": (False, ["--sites-only-vcf-output", "true"]),
            "vcf-true-bare": (False, ["--sites-only-vcf-output"]),
            "vcf-false-separated": (False, ["--sites-only-vcf-output", "false"]),
            "gvcf-true-separated": (True, ["--sites-only-vcf-output", "true"]),
            "gvcf-false-separated": (True, ["--sites-only-vcf-output", "false"]),
        }
        for label, (gvcf, boolean_args) in cases.items():
            suffix = ".g.vcf.gz" if gvcf else ".vcf"
            native_output = work / f"{label}.native{suffix}"
            gatk_output = work / f"{label}.gatk{suffix}"
            manifest = work / f"{label}.manifest.json"
            common = ["-R", str(reference), "-I", str(bam), "-L", "17:69000-70000"]
            if gvcf:
                common += ["-ERC", "GVCF"]
            native_run = run([
                str(native), *common, "-O", str(native_output),
                "--output-manifest", str(manifest),
                "--create-output-variant-index", "false",
                "--add-output-vcf-command-line", "false", *boolean_args,
            ])
            assert native_run.returncode == 0, f"{label}: native failed: {native_run.stderr[-4000:]}"
            gatk_run = run([
                str(java), "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller",
                *common, "-O", str(gatk_output),
                "--create-output-variant-index", "false",
                "--add-output-vcf-command-line", "false", *boolean_args,
            ])
            assert gatk_run.returncode == 0, f"{label}: GATK failed: {gatk_run.stderr[-4000:]}"
            native_header, native_records = records(native_output)
            gatk_header, gatk_records = records(gatk_output)
            expected_sites_only = boolean_args != ["--sites-only-vcf-output", "false"]
            expected_columns = 8 if expected_sites_only else 10
            assert len(native_header) == expected_columns, (label, native_header)
            assert len(gatk_header) == expected_columns, (label, gatk_header)
            # GVCF reference-block coalescing is intentionally outside this
            # writer-boundary gate (the native RCM uses a different bounded
            # partition).  Concrete candidate rows, however, must remain the
            # same; ordinary VCF rows are all concrete and are compared below.
            if gvcf:
                native_compare = [row for row in native_records if row[4] != "<NON_REF>"]
                gatk_compare = [row for row in gatk_records if row[4] != "<NON_REF>"]
            else:
                native_compare = native_records
                gatk_compare = gatk_records
            assert len(native_compare) == len(gatk_compare), (
                label, len(native_compare), len(gatk_compare)
            )
            for index, (native_record, gatk_record) in enumerate(zip(native_compare, gatk_compare)):
                assert len(native_record) == expected_columns, (label, index, native_record)
                assert len(gatk_record) == expected_columns, (label, index, gatk_record)
                assert native_record[:8] == gatk_record[:8], (
                    f"{label}/{index}: site columns differ: {native_record[:8]} != {gatk_record[:8]}"
                )
                if not expected_sites_only:
                    assert native_record[8:] == gatk_record[8:], (
                        f"{label}/{index}: FORMAT/sample differs: {native_record[8:]} != {gatk_record[8:]}"
                    )
            manifest_data = json.loads(manifest.read_text(encoding="utf-8"))
            assert manifest_data["compatibility"]["sites_only_vcf_output"] is expected_sites_only
            results[label] = {
                "columns": expected_columns,
                "records": len(native_records),
                "shape_exact": True,
            }

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "cases": results,
            "sites_only_shape_exact": True,
            "format_payload_preserved_when_false": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
