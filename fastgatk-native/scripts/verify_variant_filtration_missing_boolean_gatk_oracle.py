#!/usr/bin/env python3
"""Pinned GATK oracle for VariantFiltration missing-value Boolean binding.

GATK exposes --missing-values-evaluate-as-failing as a Boolean-valued
argument, accepting a bare switch or a separated true/false token while
rejecting an embedded '=' spelling. This keeps the Java and native command
boundaries identical, while preserving the JEXL null check as a separate path.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path


VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=20>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1
chr1\t2\t.\tC\tT\t50\tPASS\tDP=5\tGT\t0/1
chr1\t3\t.\tG\tA\t50\tPASS\tDP=20\tGT\t0/1
"""


def records(path: Path) -> list[tuple[str, str]]:
    opener = gzip.open if path.name.endswith(".gz") else open
    result: list[tuple[str, str]] = []
    with opener(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            result.append((fields[1], fields[6]))
    return result


def invoke(command: list[str], label: str) -> None:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise AssertionError(f"{label} failed ({result.returncode}): {result.stderr[-3000:]}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = Path(os.environ.get(
        "FASTGATK_VARIANT_FILTRATION_BINARY",
        root / "fastgatk-native/build/fastgatk-variant-filtration",
    ))
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    if not all(path.is_file() for path in (native, java, gatk)):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("VariantFiltration missing-value oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-variant-filtration-missing-boolean-") as directory:
        work = Path(directory)
        source = work / "missing.vcf"
        source.write_text(VCF, encoding="utf-8")
        cases = {
            "bare": ["--missing-values-evaluate-as-failing"],
            "explicit_true": ["--missing-values-evaluate-as-failing", "true"],
            "explicit_false": ["--missing-values-evaluate-as-failing", "false"],
        }
        outputs: dict[str, list[tuple[str, str]]] = {}
        for label, flag in cases.items():
            native_output = work / f"native-{label}.vcf.gz"
            gatk_output = work / f"gatk-{label}.vcf"
            common = ["--filter-expression", "DP < 10", "--filter-name", "LowDP", *flag]
            invoke([str(native), "-V", str(source), "-O", str(native_output), *common],
                   f"native {label}")
            invoke([str(java), "-jar", str(gatk), "VariantFiltration", "-V", str(source),
                    "-O", str(gatk_output), *common], f"GATK {label}")
            native_records = records(native_output)
            gatk_records = records(gatk_output)
            if native_records != gatk_records:
                raise AssertionError({"case": label, "native": native_records, "gatk": gatk_records})
            outputs[label] = native_records

        expected_failing = [("1", "LowDP"), ("2", "LowDP"), ("3", "PASS")]
        expected_non_failing = [("1", "PASS"), ("2", "LowDP"), ("3", "PASS")]
        if outputs["bare"] != expected_failing or outputs["explicit_true"] != expected_failing:
            raise AssertionError({"expected_true": expected_failing, "observed": outputs})
        if outputs["explicit_false"] != expected_non_failing:
            raise AssertionError({"case": "explicit_false", "expected": expected_non_failing,
                                  "observed": outputs["explicit_false"]})

        # Unlike the shared native optional-Boolean forms used by some tools,
        # Barclay rejects an embedded '=' for this GATK argument. Preserve
        # that exact invalid-command boundary instead of accepting a superset.
        inline = ["--missing-values-evaluate-as-failing=true"]
        native_inline = subprocess.run(
            [str(native), "-V", str(source), "-O", str(work / "native-inline-invalid.vcf.gz"),
             "--filter-expression", "DP < 10", "--filter-name", "LowDP", *inline],
            text=True, capture_output=True, check=False)
        gatk_inline = subprocess.run(
            [str(java), "-jar", str(gatk), "VariantFiltration", "-V", str(source),
             "-O", str(work / "gatk-inline-invalid.vcf"),
             "--filter-expression", "DP < 10", "--filter-name", "LowDP", *inline],
            text=True, capture_output=True, check=False)
        if native_inline.returncode == 0 or gatk_inline.returncode == 0:
            raise AssertionError({"native_inline": native_inline.returncode,
                                  "gatk_inline": gatk_inline.returncode})

        # A null comparison is independent of the missing-values policy.  It
        # must select the absent INFO field in both explicit false and true
        # modes, rather than being treated as numeric zero.
        null_outputs = []
        for label, flag in (("null_false", ["--missing-values-evaluate-as-failing", "false"]),
                            ("null_true", ["--missing-values-evaluate-as-failing", "true"])):
            native_output = work / f"native-{label}.vcf.gz"
            gatk_output = work / f"gatk-{label}.vcf"
            common = ["--filter-expression", 'vc.getAttribute("DP") == null',
                      "--filter-name", "MissingDP", *flag]
            invoke([str(native), "-V", str(source), "-O", str(native_output), *common], f"native {label}")
            invoke([str(java), "-jar", str(gatk), "VariantFiltration", "-V", str(source),
                    "-O", str(gatk_output), *common], f"GATK {label}")
            native_records = records(native_output)
            gatk_records = records(gatk_output)
            if native_records != gatk_records:
                raise AssertionError({"case": label, "native": native_records, "gatk": gatk_records})
            null_outputs.append(native_records)
        if null_outputs[0] != null_outputs[1] or null_outputs[0][0] != ("1", "MissingDP"):
            raise AssertionError({"null_outputs": null_outputs})

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "separated_boolean_forms_exact": True,
            "inline_equal_rejected": True,
            "missing_true_records": outputs["explicit_true"],
            "missing_false_records": outputs["explicit_false"],
            "null_policy_independent": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
