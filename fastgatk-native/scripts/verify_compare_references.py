#!/usr/bin/env python3
"""GATK 4.6.2.0 contract/oracle checks for native CompareReferences."""

from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile


def main() -> int:
    root = pathlib.Path(__file__).resolve().parents[2]
    native = pathlib.Path(os.environ.get(
        "FASTGATK_COMPARE_REFERENCES_BINARY",
        str(pathlib.Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
            / "fastgatk-compare-references"),
    ))
    if not native.exists():
        raise SystemExit(f"missing native binary: {native}")
    base = root / "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/reference/CompareReferences"
    with tempfile.TemporaryDirectory(prefix="fastgatk-compare-references-") as directory:
        work = pathlib.Path(directory)
        table = work / "renamed.table"
        manifest = work / "renamed.manifest.json"
        command = [
            str(native), "-R", str(base / "hg19mini.fasta"),
            "-refcomp", str(base / "hg19mini_1renamed.fasta"), "-O", str(table),
            "--md5-calculation-mode", "ALWAYS_RECALCULATE", "--output-manifest", str(manifest),
        ]
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        assert '"status":"contract-compatible"' in result.stderr
        expected_table = (base / "expected.testCompareReferencesRenamedSequence.table").read_text(encoding="ascii")
        if table.read_text(encoding="ascii") != expected_table:
            raise AssertionError("renamed sequence table differs from GATK fixture")
        if "\tDIFFER_IN_SEQUENCE_NAMES\n" not in result.stdout:
            raise AssertionError(result.stdout)
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        assert metadata["status"] == "contract-compatible"
        assert metadata["rows"] == 4 and metadata["references"] == 2
        assert metadata["md5_calculation_mode"] == "ALWAYS_RECALCULATE"
        assert metadata["telemetry"]["host_reference_io"] is True

        snps_dir = work / "snps"
        snps_dir.mkdir()
        snps = subprocess.run([
            str(native), "-R", str(base / "hg19mini.fasta"),
            "-refcomp", str(base / "hg19mini_chr2multiplesnps.fasta"),
            "-O", str(work / "snps.table"), "--base-comparison", "FIND_SNPS_ONLY",
            "--base-comparison-output", str(snps_dir), "--threads", "2",
        ], check=True, text=True, capture_output=True)
        assert '"status":"contract-compatible"' in snps.stderr
        snps_output = snps_dir / "hg19mini.fasta_hg19mini_chr2multiplesnps.fasta_snps.tsv"
        expected_snps = (base / "expected.hg19mini.fasta_hg19mini_chr2multiplesnps.fasta_snps.tsv").read_text(encoding="ascii")
        if snps_output.read_text(encoding="ascii") != expected_snps:
            raise AssertionError("FIND_SNPS_ONLY output differs from GATK fixture")
        if '"mismatches":3' not in snps.stderr:
            raise AssertionError(snps.stderr)

        missing_md5 = subprocess.run([
            str(native), "-R", str(base / "hg19mini_missingmd5.fasta"),
            "-refcomp", str(base / "hg19mini.fasta"), "--md5-calculation-mode", "USE_DICT",
        ], text=True, capture_output=True)
        assert missing_md5.returncode != 0 and "MD5 is missing" in missing_md5.stderr

        java = root / "third_party/jdk17/bin/java"
        jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        java_oracle = False
        if java.exists() and jar.exists():
            java_table = work / "java.table"
            java_result = subprocess.run([
                str(java), "-jar", str(jar), "CompareReferences", "-R", str(base / "hg19mini.fasta"),
                "-refcomp", str(base / "hg19mini_1renamed.fasta"), "-O", str(java_table),
                "--md5-calculation-mode", "ALWAYS_RECALCULATE",
            ], check=True, text=True, capture_output=True)
            if java_table.read_text(encoding="ascii") != table.read_text(encoding="ascii"):
                raise AssertionError("native/Java comparison tables differ")
            if "DIFFER_IN_SEQUENCE_NAMES" not in java_result.stdout:
                raise AssertionError(java_result.stdout)
            java_oracle = True
        print(json.dumps({"status": "pass", "rows": 4, "snps": 3,
                          "java_oracle": java_oracle}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
