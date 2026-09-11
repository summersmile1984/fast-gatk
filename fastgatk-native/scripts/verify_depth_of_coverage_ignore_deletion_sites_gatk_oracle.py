#!/usr/bin/env python3
"""Pinned GATK oracle for DepthOfCoverage deletion-site suppression.

GATK exposes ``--ignore-deletion-sites`` as an advanced boolean.  Its effect
is coupled to ``--include-deletions``: when both are requested, the deletion
pileup is excluded from the locus depth and base-count output.  The native
count path previously accepted ``--include-deletions`` but silently ignored
this switch, so deletion-only loci remained at depth one.
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
JAVA = ROOT / "third_party/jdk17/bin/java"
GATK = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"


def write_reference(path: pathlib.Path) -> None:
    sequence = "ACGTACGTACGTACGTACGT"
    path.write_text(">chrD\n" + sequence + "\n", encoding="ascii")
    path.with_suffix(path.suffix + ".fai").write_text(
        f"chrD\t{len(sequence)}\t6\t{len(sequence)}\t{len(sequence) + 1}\n",
        encoding="ascii",
    )
    path.with_suffix(".dict").write_text(
        f"@HD\tVN:1.6\n@SQ\tSN:chrD\tLN:{len(sequence)}\n", encoding="ascii"
    )


def run_java(reference: pathlib.Path, bam: pathlib.Path, output: pathlib.Path, ignore: str) -> None:
    subprocess.run(
        [str(JAVA), "-jar", str(GATK), "DepthOfCoverage", "-R", str(reference),
         "-I", str(bam), "-L", "chrD:1-10", "--include-deletions", "true",
         "--ignore-deletion-sites", ignore, "-O", str(output),
         "--omit-locus-table", "false", "--omit-interval-statistics", "true",
         "--omit-per-sample-statistics", "true", "--omit-depth-output-at-each-base", "false"],
        check=True, text=True, capture_output=True,
    )


def run_native(reference: pathlib.Path, bam: pathlib.Path, output: pathlib.Path,
               manifest: pathlib.Path, ignore: str) -> None:
    subprocess.run(
        [str(NATIVE), "-R", str(reference), "-I", str(bam), "-L", "chrD:1-10",
         "--include-deletions", "--ignore-deletion-sites", ignore, "-O", str(output),
         "--output-manifest", str(manifest), "--omit-interval-statistics", "true",
         "--omit-per-sample-statistics", "true"],
        check=True, text=True, capture_output=True,
    )


def main() -> None:
    if not NATIVE.exists() or not JAVA.exists() or not GATK.exists():
        oracle_guard.oracle_not_verified('verify_depth_of_coverage_ignore_deletion_sites_gatk_oracle.py', JAVA, GATK)
        raise SystemExit("missing DepthOfCoverage deletion oracle assets")
    with tempfile.TemporaryDirectory(prefix="fastgatk-depth-of-coverage-deletion-oracle-") as temporary:
        work = pathlib.Path(temporary)
        reference = work / "reference.fasta"
        write_reference(reference)
        sam = work / "deletion.sam"
        sam.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chrD\tLN:20\n"
            "@RG\tID:rgD\tSM:DEL_SAMPLE\n"
            "del\t0\tchrD\t1\t60\t4M2D4M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:rgD\n",
            encoding="ascii",
        )
        bam = work / "deletion.bam"
        subprocess.run(
            [str(JAVA), "-jar", str(GATK), "SortSam", "-I", str(sam), "-O", str(bam),
             "-SO", "coordinate", "--CREATE_INDEX", "true"],
            check=True, text=True, capture_output=True,
        )

        results = {}
        for ignore in ("false", "true"):
            java_output = work / f"java-{ignore}.out"
            native_output = work / f"native-{ignore}.out"
            native_manifest = work / f"native-{ignore}.manifest.json"
            run_java(reference, bam, java_output, ignore)
            run_native(reference, bam, native_output, native_manifest, ignore)
            assert native_output.read_bytes() == java_output.read_bytes(), (
                ignore, native_output.read_text(encoding="ascii"),
                java_output.read_text(encoding="ascii"),
            )
            rows = native_output.read_text(encoding="ascii").splitlines()[1:]
            assert len(rows) == 10
            deletion_rows = rows[4:6]
            if ignore == "true":
                assert deletion_rows == ["chrD:5,0,0.00,0", "chrD:6,0,0.00,0"]
            else:
                assert deletion_rows == ["chrD:5,1,1.00,1", "chrD:6,1,1.00,1"]
            manifest = json.loads(native_manifest.read_text(encoding="ascii"))
            assert manifest["compatibility"]["include_deletions"] is True
            assert manifest["compatibility"]["ignore_deletion_sites"] == (ignore == "true")
            assert manifest["compatibility"]["effective_include_deletions"] == (ignore == "false")
            results[ignore] = {
                "java_native_locus_table_exact": True,
                "deletion_rows": deletion_rows,
                "effective_include_deletions": ignore == "false",
            }

        print(json.dumps({
            "status": "pass",
            "tool": "DepthOfCoverage",
            "gatk_version": "4.6.2.0",
            "ignore_deletion_sites": results["true"],
            "include_deletions_without_ignore": results["false"],
        }, sort_keys=True))


if __name__ == "__main__":
    main()
