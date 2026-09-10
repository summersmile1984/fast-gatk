#!/usr/bin/env python3
"""File-boundary benchmarks for the native FASTA reference walkers."""

from __future__ import annotations

import json
import statistics
import subprocess
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
FASTA_BIN = ROOT / "fastgatk-native/build/fastgatk-fasta-reference-maker"
ALTERNATE_BIN = ROOT / "fastgatk-native/build/fastgatk-fasta-alternate-reference-maker"


def run(command: list[str]) -> tuple[float, Path]:
    output = Path(command[command.index("-O") + 1])
    started = time.perf_counter()
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, check=False)
    elapsed = time.perf_counter() - started
    if result.returncode != 0:
        raise RuntimeError(f"command failed ({result.returncode}): {' '.join(command)}\n{result.stderr}")
    return elapsed, output


def main() -> int:
    assert REF.exists() and REF.with_suffix(REF.suffix + ".fai").exists()
    assert FASTA_BIN.exists() and ALTERNATE_BIN.exists()
    with tempfile.TemporaryDirectory(prefix="fastgatk-fasta-bench-") as directory:
        work = Path(directory)
        sequence = "".join(line.strip() for line in REF.read_text().splitlines()
                            if not line.startswith(">"))
        vcf = work / "variants.vcf"
        pos = 1000
        ref_base = sequence[pos - 1]
        alt_base = {"A": "C", "C": "G", "G": "T", "T": "A"}[ref_base.upper()]
        vcf.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=17,length=1000000>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            f"17\t{pos}\t.\t{ref_base}\t{alt_base}\t50\tPASS\t.\n"
        )
        cases = []
        for label, binary, extra in (
            ("FastaReferenceMaker", FASTA_BIN, []),
            ("FastaAlternateReferenceMaker", ALTERNATE_BIN, ["-V", str(vcf)]),
        ):
            timings = []
            output_bytes = 0
            for repeat in range(4):
                output = work / f"{label}.{repeat}.fasta"
                elapsed, path = run([str(binary), "-R", str(REF), *extra,
                                     "-L", "17:1-1000000", "-O", str(output),
                                     "--line-width", "80", "--threads", "2"])
                timings.append(elapsed)
                output_bytes = path.stat().st_size
            cases.append({
                "tool": label,
                "bases": 1_000_000,
                "repetitions": len(timings),
                "warmup": timings[0],
                "p50_seconds": statistics.median(timings[1:]),
                "p95_seconds": max(timings[1:]),
                "bases_per_second_p50": 1_000_000 / statistics.median(timings[1:]),
                "output_bytes": output_bytes,
                "status": "pass",
            })
    print(json.dumps({"benchmark": "fasta-reference-tools-file-boundary",
                      "cases": cases, "schema_version": 1, "status": "pass"},
                     sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
