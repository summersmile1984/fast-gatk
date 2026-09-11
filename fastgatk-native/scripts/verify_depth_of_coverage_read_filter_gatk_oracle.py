#!/usr/bin/env python3
"""Pinned GATK oracle for DepthOfCoverage's default read-filter boundary.

DepthOfCoverage 4.6.2.0 resolves WellformedReadFilter, MappedReadFilter,
NotDuplicateReadFilter and NotSecondaryAlignmentReadFilter.  In particular,
QC-fail is not a default rejection, while WellformedReadFilter rejects a
whole read containing CIGAR ``N`` (CigarContainsNoNOperator).  This oracle
keeps both cases in one traversal so a native implementation cannot silently
count only the M islands of a rejected spliced read.
"""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile
import oracle_guard


ROOT = pathlib.Path(__file__).resolve().parents[2]
BUILD = pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", ROOT / "fastgatk-native/build"))
NATIVE = pathlib.Path(os.environ.get(
    "FASTGATK_DEPTH_OF_COVERAGE_BINARY", str(BUILD / "fastgatk-depth-of-coverage")
))
REFERENCE = ROOT / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
JAVA = ROOT / "third_party/jdk17/bin/java"
GATK = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


def main() -> None:
    if not NATIVE.exists() or not REFERENCE.exists() or not JAVA.exists() or not GATK.exists():
        oracle_guard.oracle_not_verified('verify_depth_of_coverage_read_filter_gatk_oracle.py', JAVA, GATK)
        raise SystemExit("missing DepthOfCoverage read-filter oracle assets")
    with tempfile.TemporaryDirectory(prefix="fastgatk-depth-of-coverage-filter-oracle-") as temporary:
        work = pathlib.Path(temporary)
        sam = work / "filter.sam"
        sam.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:17\tLN:1000000\n"
            "@RG\tID:rgB\tSM:FILTER\n"
            # QC-fail is retained by GATK's default DepthOfCoverage filters.
            "qcfail\t512\t17\t69000\t60\t4M\t*\t0\t0\tAAAA\tIIII\tRG:Z:rgB\n"
            # WellformedReadFilter rejects the entire read because of N.
            "spliced\t0\t17\t69004\t60\t2M2N2M\t*\t0\t0\tCCCC\tIIII\tRG:Z:rgB\n",
            encoding="ascii",
        )
        bam = work / "filter.bam"
        subprocess.run(
            [str(JAVA), "-jar", str(GATK), "SortSam", "-I", str(sam), "-O", str(bam),
             "-SO", "coordinate", "--CREATE_INDEX", "true"],
            check=True, text=True, capture_output=True,
        )
        native_output = work / "native.out"
        native_manifest = work / "native.manifest.json"
        subprocess.run(
            [str(NATIVE), "-R", str(REFERENCE), "-I", str(bam), "-L", "17:69000-69012",
             "-O", str(native_output), "--output-manifest", str(native_manifest),
             "--omit-interval-statistics", "true", "--omit-per-sample-statistics", "true"],
            check=True, text=True, capture_output=True,
        )
        java_output = work / "java.out"
        subprocess.run(
            [str(JAVA), "-jar", str(GATK), "DepthOfCoverage", "-R", str(REFERENCE),
             "-I", str(bam), "-L", "17:69000-69012", "-O", str(java_output),
             "--omit-locus-table", "false", "--omit-interval-statistics", "true",
             "--omit-per-sample-statistics", "true", "--omit-depth-output-at-each-base", "false"],
            check=True, text=True, capture_output=True,
        )
        assert native_output.read_bytes() == java_output.read_bytes(), (
            native_output.read_text(encoding="ascii"), java_output.read_text(encoding="ascii")
        )
        expected = [
            "17:69000,1,1.00,1", "17:69001,1,1.00,1", "17:69002,1,1.00,1",
            "17:69003,1,1.00,1", "17:69004,0,0.00,0", "17:69005,0,0.00,0",
            "17:69006,0,0.00,0", "17:69007,0,0.00,0", "17:69008,0,0.00,0",
            "17:69009,0,0.00,0", "17:69010,0,0.00,0", "17:69011,0,0.00,0",
            "17:69012,0,0.00,0",
        ]
        assert native_output.read_text(encoding="ascii").splitlines()[1:] == expected
        manifest = json.loads(native_manifest.read_text(encoding="ascii"))
        assert manifest["compatibility"]["count_type"] == "COUNT_READS"
        assert manifest["telemetry"]["reads_seen"] == 2
        assert manifest["telemetry"]["reads_used"] == 1
        assert manifest["telemetry"]["reads_filtered"] == 1
        print(json.dumps({
            "status": "pass",
            "tool": "DepthOfCoverage",
            "gatk_version": "4.6.2.0",
            "native_java_locus_table_exact": True,
            "qcfail_retained": True,
            "cigar_reference_skip_read_rejected": True,
            "reads_seen": manifest["telemetry"]["reads_seen"],
            "reads_used": manifest["telemetry"]["reads_used"],
            "reads_filtered": manifest["telemetry"]["reads_filtered"],
        }, sort_keys=True))


if __name__ == "__main__":
    main()
