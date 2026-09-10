#!/usr/bin/env python3
"""Deep oracle: --alleles crossed with the allele-count limits and gVCF.

Cases where GenotypeGivenAlleles competes with the downstream allele budget:

  * --max-alternate-alleles 1 / 2 with three forced SNP events
    (GenotypingEngine.calculateGenotypes -> calculateMostLikelyAlleles)
  * --max-genotype-count 2 with three forced events
    (HaplotypeCallerGenotypingEngine.removeAltAllelesIfTooManyGenotypes ->
     AlleleScoredByHaplotypeScores; injected haplotypes carry a NaN score)
  * a forced homopolymer insertion under -ERC GVCF
    (ReferenceConfidenceModel must keep the forced ALT beside <NON_REF>)

Shares the reference/BAM fixture of verify_hc_alleles_deep_boundary.py.
Non-zero exit on divergence; --case NAME runs a single case.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from verify_hc_alleles_deep_boundary import (  # noqa: E402
    HOMOPOLYMER_INTERVAL,
    HOMOPOLYMER_START,
    build_reference,
    compare,
    records,
    write_feature,
    write_inputs,
)

SNP_POSITIONS = (620, 640, 660)


def case_definitions(reference_text: str) -> list[dict]:
    rows = []
    for position in SNP_POSITIONS:
        ref_base = reference_text[position - 1]
        rows.append(("chr1", position, ref_base,
                     next(base for base in "ACGT" if base != ref_base)))
    return [
        {
            "name": "max-alt-alleles-1",
            "interval": HOMOPOLYMER_INTERVAL,
            "rows": rows,
            "extra": ["--max-alternate-alleles", "1"],
            "why": "three forced SNPs but only one alt allele may be emitted",
        },
        {
            "name": "max-alt-alleles-2",
            "interval": HOMOPOLYMER_INTERVAL,
            "rows": rows,
            "extra": ["--max-alternate-alleles", "2"],
            "why": "three forced SNPs with a two-alt budget",
        },
        {
            "name": "max-genotype-count-3",
            "interval": HOMOPOLYMER_INTERVAL,
            "rows": rows,
            "extra": ["--max-genotype-count", "3"],
            "why": "three forced SNPs with a 3-genotype budget (Java-accepted value)",
        },
        {
            "name": "max-genotype-count-2",
            "interval": HOMOPOLYMER_INTERVAL,
            "rows": rows,
            "extra": ["--max-genotype-count", "2"],
            "why": "forced haplotypes carry NaN scores; trimming must still keep them",
        },
        {
            "name": "gvcf-homopolymer-insertion",
            "interval": HOMOPOLYMER_INTERVAL,
            "rows": [("chr1", HOMOPOLYMER_START, "A", "AA")],
            "extra": ["-ERC", "GVCF"],
            "why": "forced homopolymer insertion emitted beside <NON_REF> in gVCF",
        },
        {
            "name": "gvcf-max-alt-alleles-1",
            "interval": HOMOPOLYMER_INTERVAL,
            "rows": rows,
            "extra": ["-ERC", "GVCF", "--max-alternate-alleles", "1"],
            "why": "gVCF reference-confidence blocks under a forced-alt budget",
        },
    ]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", default=None)
    parser.add_argument("--native", default=os.environ.get("FASTGATK_HC_BINARY"))
    parser.add_argument("--drop-alleles", action="store_true",
                        help="control run: omit --alleles entirely")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    native = Path(args.native) if args.native else root / "fastgatk-native/build/fastgatk-hc-call"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    java = os.environ.get("JAVA", str(root / "third_party/jdk17/bin/java"))
    if not native.exists() or not gatk.exists():
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE"):
            raise SystemExit("missing HC binary or pinned GATK 4.6.2.0 jar")
        print(json.dumps({"status": "skipped", "reason": "GATK jar or native binary absent"}))
        return 0

    results = []
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-alleles-limits-") as directory:
        work = Path(directory)
        reference, sam, reference_text = write_inputs(work)
        bam = work / "input.bam"
        subprocess.run([java, "-Xmx1g", "-jar", str(gatk), "CreateSequenceDictionary",
                        "-R", str(reference), "-O", str(work / "ref.dict"),
                        "--TRUNCATE_NAMES_AT_WHITESPACE", "true"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       timeout=300)
        subprocess.run([java, "-Xmx1g", "-jar", str(gatk), "SortSam", "-I", str(sam),
                        "-O", str(bam), "-SO", "coordinate", "--CREATE_INDEX", "true"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       timeout=300)

        for case in case_definitions(reference_text):
            if args.case and case["name"] not in args.case:
                continue
            feature = write_feature(work, case["name"], case["rows"], root)
            allele_args = [] if args.drop_alleles else ["--alleles", str(feature)]
            common = ["-R", str(reference), "-I", str(bam), "-L", case["interval"],
                      *allele_args, "--min-pruning", "1",
                      "--create-output-variant-index", "false",
                      "--add-output-vcf-command-line", "false", *case["extra"]]
            java_vcf = work / f"gatk.{case['name']}.vcf"
            native_vcf = work / f"native.{case['name']}.vcf"
            java_run = subprocess.run(
                [java, "-Xmx1g", "-jar", str(gatk), "HaplotypeCaller", *common,
                 "-O", str(java_vcf)],
                capture_output=True, text=True, timeout=300)
            native_run = subprocess.run(
                [str(native), *common, "--threads", "2", "-O", str(native_vcf)],
                capture_output=True, text=True, timeout=300)
            entry = {"case": case["name"], "why": case["why"],
                     "extra": case["extra"],
                     "java_exit": java_run.returncode, "native_exit": native_run.returncode}
            if java_run.returncode != 0 or native_run.returncode != 0:
                entry.update({"status": "error",
                              "java_stderr": java_run.stderr[-400:],
                              "native_stderr": native_run.stderr[-400:]})
                results.append(entry)
                continue
            java_rows, native_rows = records(java_vcf), records(native_vcf)
            entry.update(compare(java_rows, native_rows))
            entry["match"] = "first_diff_field" not in entry
            entry["status"] = "match" if entry["match"] else "diverge"
            if not entry["match"]:
                entry["java_rows_all"] = java_rows
                entry["native_rows_all"] = native_rows
            results.append(entry)

    diverged = [entry["case"] for entry in results if entry["status"] != "match"]
    print(json.dumps({
        "status": "pass" if not diverged else "diverge",
        "release": "GATK 4.6.2.0",
        "binary": str(native),
        "cases": results,
        "diverged_cases": diverged,
    }, indent=2, sort_keys=True))
    return 0 if not diverged else 1


if __name__ == "__main__":
    raise SystemExit(main())
