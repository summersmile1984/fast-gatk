#!/usr/bin/env python3
"""Pinned GATK oracle for VariantFiltration genotype no-call materialization."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path

from verify_variant_filtration_gatk_oracle import read_vcf


VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
##FORMAT=<ID=FT,Number=1,Type=String,Description=Genotype filter>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t10\trsLow\tA\tG\t25\tPASS\tDP=5\tGT:GQ:DP:AD:PL\t0/1:10:5:3,2:20,0,20\t1/1:50:20:0,20:50,20,0
chr1\t20\trsGood\tC\tT\t60\tPASS\tDP=20\tGT:GQ:DP:AD:PL\t0/0:50:20:20,0:0,30,60\t0/1:5:12:8,4:60,0,60
chr1\t30\trsPrior\tG\tA\t60\tPASS\tDP=30\tGT:GQ:DP:AD:PL:FT\t0/1:40:30:15,15:0,30,60:.\t1/1:50:30:0,30:50,20,0:Prior
"""


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get(
        "FASTGATK_VARIANT_FILTRATION_BINARY",
        root / "fastgatk-native/build/fastgatk-variant-filtration",
    ))
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    required = (binary, java, jar)
    if not all(path.is_file() for path in required):
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("VariantFiltration set-no-call oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-variant-filtration-set-nocall-oracle-") as directory:
        work = Path(directory)
        source = work / "input.vcf"
        source.write_text(VCF, encoding="utf-8")

        def gatk_run(output: Path, args: list[str]) -> subprocess.CompletedProcess[str]:
            return run([
                str(java), "-jar", str(jar), "VariantFiltration", "-V", str(source),
                "-O", str(output), "--create-output-variant-index", "false", *args,
            ])

        def native_run(output: Path, args: list[str]) -> subprocess.CompletedProcess[str]:
            return run([
                str(binary), "-V", str(source), "-O", str(output),
                "--create-output-variant-index=false", *args,
            ])

        def compare(label: str, args: list[str]) -> tuple[list[tuple[object, ...]], dict]:
            gatk_output = work / f"gatk-{label}.vcf"
            native_output = work / f"native-{label}.vcf"
            gatk_result = gatk_run(gatk_output, args)
            native_result = native_run(native_output, args)
            if gatk_result.returncode != 0 or native_result.returncode != 0:
                raise AssertionError(
                    f"{label} failed: GATK={gatk_result.returncode}, native={native_result.returncode}; "
                    f"native stderr={native_result.stderr[-2000:]}")
            gatk_rows = read_vcf(gatk_output)
            native_rows = read_vcf(native_output)
            if native_rows != gatk_rows:
                raise AssertionError(f"{label} native/GATK rows differ:\nnative={native_rows!r}\nGATK={gatk_rows!r}")
            summary_lines = native_result.stdout.splitlines()
            if not summary_lines:
                raise AssertionError(f"{label} native summary missing")
            return native_rows, json.loads(summary_lines[-1])

        filter_args = [
            "--genotype-filter-expression", "GQ < 20",
            "--genotype-filter-name", "LowGQ",
        ]
        false_rows, false_summary = compare("false", [*filter_args, "--set-filtered-genotype-to-no-call", "false"])
        true_rows, true_summary = compare("true", [*filter_args, "--set-filtered-genotype-to-no-call", "true"])
        bare_rows, bare_summary = compare("bare", [*filter_args, "--set-filtered-genotype-to-no-call"])
        if true_rows != bare_rows:
            raise AssertionError("bare and explicit true no-call rows differ")
        if false_rows == true_rows:
            raise AssertionError("set-no-call true did not change filtered genotypes")
        if false_summary.get("genotype_set_to_no_call") != 0:
            raise AssertionError(f"false no-call summary is nonzero: {false_summary}")
        if true_summary.get("genotype_set_to_no_call") != 3 or bare_summary.get("genotype_set_to_no_call") != 3:
            raise AssertionError(f"unexpected no-call count: true={true_summary}, bare={bare_summary}")

        # Existing FT state is also a filtered genotype in GATK's helper.  No
        # new genotype expression is required for this chained-stage case.
        prior_rows, prior_summary = compare(
            "existing-ft", ["--set-filtered-genotype-to-no-call", "true"])
        if prior_summary.get("genotype_set_to_no_call") != 1:
            raise AssertionError(f"existing FT was not converted: {prior_summary}")
        if "./." not in repr(prior_rows[2]):
            raise AssertionError(f"existing FT genotype was not no-called: {prior_rows[2]!r}")

        # Barclay rejects an embedded '=' spelling for this Boolean.  Native
        # must fail closed at the owning option instead of silently accepting
        # a command line Java would reject.
        invalid_gatk = gatk_run(work / "invalid-gatk.vcf", [
            *filter_args, "--set-filtered-genotype-to-no-call=true",
        ])
        invalid_native = native_run(work / "invalid-native.vcf", [
            *filter_args, "--set-filtered-genotype-to-no-call=true",
        ])
        if invalid_gatk.returncode == 0 or invalid_native.returncode == 0:
            raise AssertionError("embedded '=' no-call Boolean was accepted")

    print(json.dumps({
        "status": "pass",
        "gatk_version": "4.6.2.0",
        "set_filtered_genotype_to_no_call_exact": True,
        "recomputed_ac_an_af_exact": True,
        "existing_ft_conversion_exact": True,
        "bare_and_separated_true_exact": True,
        "embedded_equals_fail_closed": True,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
