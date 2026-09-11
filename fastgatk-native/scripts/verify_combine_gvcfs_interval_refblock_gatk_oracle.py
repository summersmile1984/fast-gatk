#!/usr/bin/env python3
"""Pinned GATK oracle for CombineGVCFs interval/reference-block boundaries.

CombineGVCFs does not clip a REF/<NON_REF> block to an interval.  If ``-L``
ends in the middle of the block, the Java walker warns and drops that block;
when the block END is inside the interval it emits the original full block.
This is a small direct-replacement boundary that ordinary overlap predicates
miss, so exercise both the aggregate and bounded stream native paths.
"""
from __future__ import annotations

import gzip
import json
import os
from pathlib import Path
import subprocess
import tempfile
import oracle_guard


HEADER = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description=All filters passed>
##contig=<ID=17,length=1000000>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=END,Number=1,Type=Integer,Description=End position of a reference block>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Genotype quality>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Phred-scaled genotype likelihoods>
"""


def write(path: Path, sample: str) -> None:
    path.write_text(
        HEADER
        + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
        + "17\t69000\t.\tA\t<NON_REF>\t.\tPASS\tEND=69005\t"
        + "GT:GQ:DP:AD:PL\t0/0:60:10:10,0:0,99,99\n",
        encoding="utf-8",
    )


def run(command: list[str], label: str) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise AssertionError(f"{label} failed ({result.returncode}):\n{result.stderr}")
    return result


def body(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        return [line.rstrip("\r\n") for line in stream
                if line and not line.startswith("#")]


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    java = root / "third_party/jdk17/bin/java"
    gatk = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    native = Path(os.environ.get(
        "FASTGATK_COMBINE_GVCFS_BINARY",
        str(root / "fastgatk-native/build/fastgatk-combine-gvcfs"),
    ))
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    required = (java, gatk, native, reference)
    if not all(path.is_file() for path in required):
        oracle_guard.oracle_not_verified('verify_combine_gvcfs_interval_refblock_gatk_oracle.py', java, gatk)
        if os.environ.get("FASTGATK_REQUIRE_GATK_ORACLE") == "1":
            raise SystemExit("bundled GATK interval/reference-block oracle inputs are required")
        print(json.dumps({"status": "skip", "reason": "bundled GATK oracle unavailable"}))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-combine-gvcfs-interval-oracle-") as directory:
        work = Path(directory)
        first = work / "first.g.vcf"
        second = work / "second.g.vcf"
        write(first, "SAMPLE1")
        write(second, "SAMPLE2")
        for input_path in (first, second):
            run([str(java), "-jar", str(gatk), "IndexFeatureFile", "-I", str(input_path)],
                "GATK IndexFeatureFile")

        results: dict[str, dict[str, int | bool]] = {}
        # The first interval ends inside END=69005: Java drops the block.  The
        # second contains END: Java emits the original untrimmed block.
        for name, interval, expected_count in (
            ("cut", "17:69000-69003", 0),
            ("contained_end", "17:69003-69010", 1),
        ):
            gatk_output = work / f"gatk-{name}.vcf.gz"
            native_output = work / f"native-{name}.vcf.gz"
            native_stream = work / f"native-{name}-stream.vcf.gz"
            run([
                str(java), "-jar", str(gatk), "CombineGVCFs", "-R", str(reference),
                "-V", str(first), "-V", str(second), "-L", interval,
                "-O", str(gatk_output), "--create-output-variant-index", "false",
            ], "GATK CombineGVCFs")
            for output, extra in (
                (native_output, []),
                (native_stream, ["--stream-merge"]),
            ):
                run([
                    str(native), "-R", str(reference), "-V", str(first),
                    "-V", str(second), "-L", interval, "-O", str(output),
                    "--create-output-variant-index=false", *extra,
                ], f"native CombineGVCFs {name} {'stream' if extra else 'aggregate'}")
            gatk_rows = body(gatk_output)
            native_rows = body(native_output)
            stream_rows = body(native_stream)
            assert len(gatk_rows) == expected_count, (name, gatk_rows)
            assert native_rows == gatk_rows, {
                "name": name, "gatk": gatk_rows, "native": native_rows,
            }
            assert stream_rows == gatk_rows, {
                "name": name, "gatk": gatk_rows, "native_stream": stream_rows,
            }
            if name == "cut":
                assert not gatk_rows
            else:
                assert gatk_rows and "END=69005" in gatk_rows[0]
            results[name] = {
                "gatk_rows": len(gatk_rows),
                "aggregate_exact": True,
                "stream_exact": True,
            }

        print(json.dumps({
            "status": "pass",
            "gatk_version": "4.6.2.0",
            "reference_block_end": 69005,
            "interval_cut_drops_block": True,
            "contained_end_preserves_full_block": True,
            "cases": results,
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
