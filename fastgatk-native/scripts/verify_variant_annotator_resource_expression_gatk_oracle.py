#!/usr/bin/env python3
"""VariantAnnotator oracle: --comp / --dbsnp / --resource + -E expressions.

Runs GATK 4.6.2.0's VariantAnnotator and the native fastgatk-variant-
annotator on identical self-contained fixtures and requires the output
VCF to be byte-identical (headers aside from execution provenance, data
rows verbatim) — the without_provenance convention shared by the other
GATK contract oracles.

Covered semantics (VariantAnnotatorEngine / VariantOverlapAnnotator /
VAExpression in gatk-source .../walkers/annotator/):
  * --comp[:NAME]     : NAME membership flag, matched by minimum-
    representation biallelic equality (same ref+alt) over unfiltered
    source records starting at the site.
  * --dbsnp           : DB membership flag + rsID merge into the ID field
    (verbatim getID() join; appended with ';' only when not contained).
  * --resource:NAME + -E NAME.FIELD: the first resource record starting at
    the site contributes FIELD (ID/ALT/FILTER special cases); A/R-counted
    INFO fields are mapped onto the input's minimum-representation
    biallelics with "0" for unmatched alleles and omitted entirely when
    nothing matched; other fields transfer the raw value.
  * header lines: htsjdk's canonical lines for GQ (and the other default-
    annotation-owned keys), sorted FORMAT/INFO groups, contigs last.

The fixtures are written by this script (the gatk-source "large" fixtures
are unmaterialized git-lfs pointers in this checkout).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import oracle_guard


ROOT = Path(__file__).resolve().parents[2]
JAVA = ROOT / "third_party/jdk17/bin/java"
GATK = ROOT / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
NATIVE = Path(os.environ.get(
    "FASTGATK_VARIANT_ANNOTATOR_BINARY",
    str(ROOT / "fastgatk-native/build/fastgatk-variant-annotator")))

INPUT_VCF = """\
##fileformat=VCFv4.1
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Read Depth (only filtered reads used for calling)">
##source=UnifiedGenotyper
##contig=<ID=20,length=63025520>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
20\t10002353\trs372143558\tG\tT\t.\t.\t.
20\t10002374\t.\tC\tT\t.\t.\t.
20\t10002443\trs1;rs2\tT\tC\t.\t.\t.
20\t10002458\t.\tG\tGTT,GTTT\t.\t.\t.
20\t10002460\trstestnot\tTC\tT,TTT\t.\t.\t.
20\t10002470\trstest\tC\tT\t.\t.\t.
20\t10002470\trstest\tC\tCTTT\t.\t.\t.
20\t10002478\t.\tA\tATT,ATTT\t.\t.\t.
"""

# Covers GATK's VariantAnnotatorIntegrationTest match classes: same-tag
# rsID, no-tag rsID, filtered source (flags suppressed, expressions kept),
# multi-allelic expression value mapping ("1,0"), and a source record whose
# alleles only partially match the input.
RESOURCE_VCF = """\
##fileformat=VCFv4.2
##contig=<ID=20,length=63025520>
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples">
##FILTER=<ID=LowQual,Description="Low quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
20\t10002353\trs372143558\tG\tT\t.\tPASS\tAC=1;NS=1
20\t10002443\trs188831105\tT\tC\t.\tLowQual\tAC=2,3;NS=2
20\t10002458\trs34527371\tG\tGTT\t.\tPASS\tAC=1;NS=1
20\t10002470\trs2327260\tC\tT\t.\tPASS\tAC=1;NS=3
20\t10002478\trs33961276\tA\tATT\t.\tLowQual\tAC=1;NS=1
"""

EXPRESSIONS = ["foo.FILTER", "foo.ID", "foo.AC", "foo.NS"]


def annotator_args(work: Path) -> list[str]:
    resource = str(work / "resource.vcf")
    return ["--resource:foo", resource,
            *[item for expression in EXPRESSIONS for item in ("-E", expression)],
            "--comp:bar", resource,
            "--dbsnp", resource]


def without_provenance(path: Path) -> list[str]:
    with open(path, encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle
                if not line.startswith(("##GATKCommandLine=", "##fileDate="))]


def run(command: list[str], timeout: int = 600) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, capture_output=True, text=True, check=False,
                          timeout=timeout)


def main() -> int:
    if not JAVA.is_file() or not GATK.is_file() or not NATIVE.is_file():
        oracle_guard.oracle_not_verified(
            'verify_variant_annotator_resource_expression_gatk_oracle.py', JAVA, GATK)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK VariantAnnotator inputs are required")
        print(json.dumps({"status": "skip",
                          "reason": "bundled GATK or native VariantAnnotator unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-variant-annotator-") as directory:
        work = Path(directory)
        (work / "input.vcf").write_text(INPUT_VCF, encoding="utf-8")
        (work / "resource.vcf").write_text(RESOURCE_VCF, encoding="utf-8")

        # GATK's FeatureInputs require random access; the bundled
        # IndexFeatureFile builds the tribble index.
        index = run([str(JAVA), "-Xmx512m", "-jar", str(GATK), "IndexFeatureFile",
                     "-I", str(work / "resource.vcf")])
        if index.returncode != 0:
            print(json.dumps({"status": "fail", "stage": "index",
                              "stderr": index.stderr[-2000:]}))
            return 1

        gatk_output = work / "gatk.vcf"
        native_output = work / "native.vcf"
        gatk = run([str(JAVA), "-Xmx512m", "-jar", str(GATK), "VariantAnnotator",
                    "-V", str(work / "input.vcf"), "-O", str(gatk_output),
                    *annotator_args(work)])
        if gatk.returncode != 0:
            print(json.dumps({"status": "fail", "stage": "gatk",
                              "stderr": gatk.stderr[-2000:]}))
            return 1
        native = run([str(NATIVE), "-V", str(work / "input.vcf"),
                      "-O", str(native_output), *annotator_args(work)])
        if native.returncode != 0:
            print(json.dumps({"status": "fail", "stage": "native",
                              "stderr": native.stderr[-2000:]}))
            return 1

        expected = without_provenance(gatk_output)
        observed = without_provenance(native_output)
        rows = [line for line in expected if not line.startswith("#")]
        mismatches = [
            {"line": index, "gatk": line, "native": observed[index]
                if index < len(observed) else None}
            for index, line in enumerate(expected)
            if index >= len(observed) or observed[index] != line
        ]
        payload = {
            "status": "pass" if not mismatches and len(observed) == len(expected)
                      else "fail",
            "rows": len(rows),
            "header_lines": len(expected) - len(rows),
            "mismatches": mismatches,
            "oracle": "GATK-4.6.2.0",
            "covered": ["comp", "dbsnp", "resource-expression",
                        "expression-allele-mapping", "canonical-header-lines"],
        }
        print(json.dumps(payload, sort_keys=True))
        return 0 if payload["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
