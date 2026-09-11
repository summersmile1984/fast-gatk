#!/usr/bin/env python3
"""Pinned GATK 4.6.2.0 oracle for SelectVariants FILTER exclusion."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##FILTER=<ID=q10,Description=Low quality>
##FILTER=<ID=LowQual,Description=Low quality>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t1\tpass\tA\tG\t50\tPASS\t.\tGT\t0/1
chr1\t2\tfiltered\tC\tT\t40\tq10\t.\tGT\t0/1
chr1\t3\tdot\tG\tA\t30\t.\t.\tGT\t0/0
chr1\t4\tmulti\tT\tC\t20\tq10;LowQual\t.\tGT\t0/1
"""


def invoke(java: Path, jar: Path, args: list[str]) -> None:
    result = subprocess.run([str(java), "-jar", str(jar), "SelectVariants", *args],
                            text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"GATK SelectVariants failed: {result.stderr[-3000:]}")


def native(binary: Path, args: list[str]) -> dict:
    result = subprocess.run([str(binary), *args], text=True,
                            capture_output=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"native SelectVariants failed: {result.stderr[-3000:]}")
    return json.loads(result.stdout.splitlines()[-1])


def rows(path: Path) -> list[tuple[str, ...]]:
    result = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            # FILTER exclusion is a site-level operation.  Preserve the full
            # record payload to catch accidental FORMAT/INFO rewrites.
            result.append(tuple(fields))
    return result


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_SELECT_VARIANTS_BINARY",
        root / "fastgatk-native/build/fastgatk-select-variants"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    required = (binary, java, jar)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_select_variants_filtered_oracle.py', java, jar)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("SelectVariants FILTER oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-select-variants-filtered-") as directory:
        work = Path(directory)
        source = work / "input.vcf"
        source.write_text(VCF, encoding="utf-8")

        gatk_true = work / "gatk-true.vcf"
        invoke(java, jar, ["-V", str(source), "--exclude-filtered", "true",
                           "-O", str(gatk_true)])
        expected_true = rows(gatk_true)
        if [row[1] for row in expected_true] != ["1", "3"]:
            raise AssertionError(f"unexpected GATK filtered result: {expected_true}")

        true_results: dict[str, list[tuple[str, ...]]] = {}
        true_summaries: dict[str, dict] = {}
        true_filtered_counts: dict[str, int] = {}
        for label, option in (
                ("bare", ["--exclude-filtered"]),
                ("true", ["--exclude-filtered", "true"]),
                ("inline", ["--exclude-filtered=true"]),
                ("alias", ["--exclude-filtered-variants", "1"]),
        ):
            output = work / f"native-{label}.vcf"
            manifest = work / f"native-{label}.manifest.json"
            summary = native(binary, ["-V", str(source), *option,
                                     "--create-output-variant-index=false",
                                     "--output-manifest", str(manifest),
                                     "-O", str(output)])
            true_results[label] = rows(output)
            true_summaries[label] = summary
            if true_results[label] != expected_true:
                raise AssertionError(f"native {label} differs from GATK: {true_results[label]}")
            metadata = json.loads(manifest.read_text(encoding="utf-8"))
            if metadata["compatibility"].get("exclude_filtered") is not True:
                raise AssertionError(f"native {label} manifest lost FILTER choice: {metadata}")
            if metadata["telemetry"].get("filtered_records") != 2:
                raise AssertionError(f"native {label} filtered count mismatch: {metadata}")
            true_filtered_counts[label] = metadata["telemetry"]["filtered_records"]

        gatk_false = work / "gatk-false.vcf"
        invoke(java, jar, ["-V", str(source), "--exclude-filtered", "false",
                           "-O", str(gatk_false)])
        expected_false = rows(gatk_false)
        false_output = work / "native-false.vcf"
        false_summary = native(binary, ["-V", str(source), "--exclude-filtered=false",
                                        "--create-output-variant-index=false",
                                        "-O", str(false_output)])
        if rows(false_output) != expected_false or [row[1] for row in rows(false_output)] != ["1", "2", "3", "4"]:
            raise AssertionError(("explicit false FILTER selection differs from GATK",
                                  rows(false_output), expected_false))
        false_metadata = json.loads(Path(str(false_output) + ".manifest.json").read_text(encoding="utf-8"))
        if false_metadata["compatibility"].get("exclude_filtered") is not False:
            raise AssertionError(f"false form manifest mismatch: {false_metadata}")
        if false_metadata["telemetry"].get("filtered_records") != 0:
            raise AssertionError(f"false form unexpectedly filtered: {false_metadata}")

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "filtered_positions": [row[1] for row in expected_true],
            "retained_positions_false": [row[1] for row in expected_false],
            "optional_boolean_forms_exact": sorted(true_results),
            "filtered_records_true": true_filtered_counts["true"],
            "filtered_records_false": false_metadata["telemetry"].get("filtered_records"),
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
