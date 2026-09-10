#!/usr/bin/env python3
"""Focused GATK 4.6.2.0 SelectVariants type/non-variant oracle.

This is intentionally separate from the broad synthetic SelectVariants
contract.  It locks the common Nextflow-facing selection slice and the
Barclay optional-boolean surface for ``--exclude-non-variants`` while keeping
the Java startup cost bounded.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path


VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##ALT=<ID=NON_REF,Description=Any alternate allele>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t1\t.\tA\t.\t50\tPASS\t.\tGT\t0/0
chr1\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/0
chr1\t3\t.\tA\tG\t50\tPASS\t.\tGT\t./.
chr1\t4\t.\tA\tG\t50\tPASS\t.\tGT\t0/1
chr1\t5\t.\tA\t<NON_REF>\t50\tPASS\t.\tGT\t0/0
chr1\t6\t.\tA\t<NON_REF>\t50\tPASS\t.\tGT\t0/1
"""


def records(path: Path) -> list[tuple[str, ...]]:
    opener = path.open
    if path.suffix == ".gz":
        import gzip
        opener = lambda **kwargs: gzip.open(path, "rt", **kwargs)
    result: list[tuple[str, ...]] = []
    with opener(encoding="utf-8") as stream:
        for line in stream:
            if line and not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                # Normalize INFO/FORMAT ordering while retaining the genotype
                # payload that is relevant to this selection boundary.
                result.append((fields[0], fields[1], fields[3], fields[4],
                                fields[5], fields[6], fields[7],
                                fields[8], fields[9]))
    return result


def invoke(java: Path, jar: Path, tool: str, args: list[str]) -> None:
    result = subprocess.run([str(java), "-jar", str(jar), tool, *args],
                            text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"GATK {tool} failed: {result.stderr[-3000:]}")


def native(binary: Path, args: list[str]) -> dict:
    result = subprocess.run([str(binary), *args], text=True,
                            capture_output=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"native SelectVariants failed: {result.stderr[-3000:]}")
    return json.loads(result.stdout.splitlines()[-1])


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_SELECT_VARIANTS_BINARY",
        root / "fastgatk-native/build/fastgatk-select-variants"))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    required = (binary, java, jar)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("SelectVariants Java oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-select-variants-gatk-oracle-") as directory:
        work = Path(directory)
        source = work / "input.vcf"
        source.write_text(VCF, encoding="utf-8")
        common = ["-V", str(source)]

        # GATK requires an explicit boolean value.  Native supports that form
        # plus the bare/inline forms, all mapping to the same records.
        gatk_true = work / "gatk-true.vcf"
        invoke(java, jar, "SelectVariants", common + [
            "--exclude-non-variants", "true", "-O", str(gatk_true)])
        expected_true = records(gatk_true)
        if [row[1] for row in expected_true] != ["4", "6"]:
            raise AssertionError(f"unexpected GATK non-variant result: {expected_true}")

        native_results: dict[str, list[tuple[str, ...]]] = {}
        native_summaries: dict[str, dict] = {}
        for label, switch in (("bare", ["--exclude-non-variants"]),
                              ("true", ["--exclude-non-variants", "true"]),
                              ("inline", ["--exclude-non-variants=true"]),
                              ("alias", ["--exclude-non-variant-sites", "1"])):
            output = work / f"native-{label}.vcf"
            manifest = work / f"native-{label}.manifest.json"
            summary = native(binary, common + switch + [
                "--create-output-variant-index=false", "-O", str(output),
                "--output-manifest", str(manifest)])
            native_results[label] = records(output)
            native_summaries[label] = summary
            if native_results[label] != expected_true:
                raise AssertionError(f"native {label} differs from GATK: {native_results[label]}")
            metadata = json.loads(manifest.read_text(encoding="utf-8"))
            if metadata["compatibility"]["exclude_non_variants"] is not True:
                raise AssertionError(f"native {label} manifest lost the true option")

        # False must retain all records, and the SNP type selector is a
        # separate, composable filter.  Compare this path to Java as well.
        gatk_snp = work / "gatk-snp.vcf"
        invoke(java, jar, "SelectVariants", common + [
            "--exclude-non-variants", "false", "--select-type-to-include", "SNP",
            "-O", str(gatk_snp)])
        native_snp = work / "native-snp.vcf"
        native_snp_summary = native(binary, common + [
            "--exclude-non-variants=false", "--select-type-to-include", "SNP",
            "--create-output-variant-index=false", "-O", str(native_snp)])
        if records(native_snp) != records(gatk_snp):
            raise AssertionError("native SNP type selection differs from GATK")
        if [row[1] for row in records(native_snp)] != ["2", "3", "4"]:
            raise AssertionError("unexpected SNP selection positions")

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "exclude_non_variants_records": len(expected_true),
            "exclude_non_variants_positions": [row[1] for row in expected_true],
            "optional_boolean_forms_exact": sorted(native_results),
            "select_type_snp_records": len(records(native_snp)),
            "native_summaries": {key: value.get("output_records")
                                 for key, value in native_summaries.items()},
            "native_snp_summary": native_snp_summary.get("output_records"),
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
