#!/usr/bin/env python3
"""Pin FilterMutectCalls' ``-OVI`` alias to GATK 4.6.2.0.

Barclay exposes ``--create-output-variant-index,-OVI <Boolean>``.  The short
alias is used by existing GATK wrappers, so a native parser that only accepts
the long spelling is not a direct replacement even when its writer is
correct.  This oracle runs both implementations with the exact GATK syntax
and checks the observable index sidecar boundary for false and true.
"""

from __future__ import annotations

import gzip
import json
import os
import subprocess
import tempfile
from pathlib import Path
import oracle_guard


def run(command: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, text=True, capture_output=True, check=False)


def records(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        return [line for line in stream if line and not line.startswith("#")]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    binary = Path(os.environ.get(
        "FASTGATK_FILTER_MUTECT_BINARY",
        str(root / "fastgatk-native/build/fastgatk-filter-mutect-calls"),
    ))
    required = (java, gatk, reference, binary)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_filter_mutect_variant_index_alias_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK FilterMutectCalls -OVI oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK -OVI oracle unavailable"}))
        return 0

    header = """##fileformat=VCFv4.2
##contig=<ID=17,length=1000000>
##INFO=<ID=TLOD,Number=A,Type=Float,Description=Tumor log odds>
##INFO=<ID=NLOD,Number=A,Type=Float,Description=Normal log odds>
##INFO=<ID=POPAF,Number=A,Type=Float,Description=Population allele frequency>
##INFO=<ID=MBQ,Number=R,Type=Integer,Description=Median base quality>
##INFO=<ID=MMQ,Number=R,Type=Integer,Description=Median mapping quality>
##INFO=<ID=MPOS,Number=A,Type=Integer,Description=Median alternate read position>
##INFO=<ID=MFRL,Number=R,Type=Integer,Description=Median fragment length>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=AF,Number=A,Type=Float,Description=Allele fraction>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR
"""
    body = "17\t69000\t.\tA\tG\t.\tPASS\tTLOD=8;NLOD=0;POPAF=2;MBQ=30,30;MMQ=60,60;MPOS=5;MFRL=100,100\tGT:AD:AF\t0/1:90,10:0.1\n"

    with tempfile.TemporaryDirectory(prefix="fastgatk-filter-ovi-oracle-") as directory:
        work = Path(directory)
        input_path = work / "input.vcf"
        input_path.write_text(header + body, encoding="utf-8")
        index_run = run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(input_path)])
        assert index_run.returncode == 0, index_run.stderr
        stats = work / "calls.stats"
        stats.write_text(
            "#<METADATA>threshold=0.1\nstatistic\tvalue\ncallable\t1000000\n",
            encoding="utf-8",
        )

        help_run = run([str(java), "-jar", str(gatk), "FilterMutectCalls", "--help"])
        assert help_run.returncode == 0, help_run.stderr
        # GATK's launcher prints Barclay help on stderr even on a successful
        # `--help` exit, so retain both streams as the release evidence.
        help_text = help_run.stdout + help_run.stderr
        assert "--create-output-variant-index,-OVI <Boolean>" in help_text

        sidecars: dict[str, dict[str, bool]] = {}
        for value in ("false", "true"):
            gatk_output = work / f"gatk-{value}.vcf.gz"
            native_output = work / f"native-{value}.vcf.gz"
            native_common = [
                "--stats", str(stats),
                "--threshold-strategy", "CONSTANT", "--initial-threshold", "0.1",
                "-OVI", value,
            ]
            gatk_common = [*native_common, "--lenient", "true"]
            gatk_run = run([
                str(java), "-jar", str(gatk), "FilterMutectCalls", "-R", str(reference),
                "-V", str(input_path), "-O", str(gatk_output), *gatk_common,
            ])
            assert gatk_run.returncode == 0, gatk_run.stderr
            native_run = run([
                str(binary), "-R", str(reference), "-V", str(input_path),
                "-O", str(native_output), *native_common,
            ])
            assert native_run.returncode == 0, native_run.stderr
            assert records(gatk_output) and records(native_output)
            assert len(records(gatk_output)) == len(records(native_output)) == 1
            expected_index = value == "true"
            gatk_index = Path(f"{gatk_output}.tbi").is_file()
            native_index = Path(f"{native_output}.tbi").is_file()
            assert gatk_index == native_index == expected_index
            sidecars[value] = {"gatk_tbi": gatk_index, "native_tbi": native_index}

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "alias": "-OVI",
            "separated_boolean_values": ["false", "true"],
            "sidecars": sidecars,
            "index_sidecar_shape_exact": True,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
