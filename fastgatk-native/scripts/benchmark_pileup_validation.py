#!/usr/bin/env python3
"""File-boundary benchmark for Get/GatherPileupSummaries and ValidateVariants.

The benchmark intentionally exercises the same HTSlib -> Kokkos -> output
boundary used by the contract tests.  It reports warmup-free wall time and
output/manifest sizes so a run is comparable across machines without claiming
Java/GATK equivalence from a synthetic fixture.
"""

from __future__ import annotations

import json
import pathlib
import statistics
import subprocess
import tempfile
import time


def p50_p95(samples: list[float]) -> tuple[float, float]:
    ordered = sorted(samples)
    index = min(len(ordered) - 1, int(0.95 * len(ordered)))
    return statistics.median(ordered), ordered[index]


def run_repeated(command: list[str], output: pathlib.Path, manifest: pathlib.Path,
                 repetitions: int = 6) -> dict[str, object]:
    samples: list[float] = []
    last_summary: dict[str, object] = {}
    for iteration in range(repetitions + 1):
        started = time.perf_counter()
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        elapsed = time.perf_counter() - started
        if completed.returncode != 0:
            raise SystemExit(completed.stderr or completed.stdout)
        if iteration > 0:
            samples.append(elapsed)
        if manifest.is_file():
            last_summary = json.loads(manifest.read_text(encoding="utf-8"))
    p50, p95 = p50_p95(samples)
    return {
        "p50_seconds": p50,
        "p95_seconds": p95,
        "output_bytes": output.stat().st_size,
        "manifest_bytes": manifest.stat().st_size if manifest.is_file() else 0,
        "telemetry": last_summary.get("telemetry", {}),
    }


def main() -> int:
    root = pathlib.Path(__file__).resolve().parents[2]
    get_pileups = root / "fastgatk-native/build/fastgatk-get-pileup-summaries"
    validate = root / "fastgatk-native/build/fastgatk-validate-variants"
    gather = root / "fastgatk-native/build/fastgatk-gather-pileup-summaries"
    if not all(path.is_file() for path in (get_pileups, validate, gather)):
        raise SystemExit("native pileup/validation binaries are required")

    with tempfile.TemporaryDirectory(prefix="fastgatk-pileup-validation-benchmark-") as directory:
        work = pathlib.Path(directory)
        reads = work / "reads.sam"
        read_lines = [
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:100000\n",
            "@RG\tID:rg1\tSM:TUMOR\n",
        ]
        for index in range(512):
            start = 1 + (index % 512)
            sequence = ("ACGT" * 26)[:100]
            read_lines.append(
                f"r{index}\t0\tchr1\t{start}\t60\t100M\t*\t0\t0\t"
                f"{sequence}\t{'I' * 100}\tRG:Z:rg1\n"
            )
        reads.write_text("".join(read_lines), encoding="ascii")

        sites = work / "sites.vcf"
        site_lines = [
            "##fileformat=VCFv4.2\n",
            "##contig=<ID=chr1,length=100000>\n",
            "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        ]
        for index in range(128):
            position = 2 + index * 8
            site_lines.append(
                f"chr1\t{position}\trs{index}\tC\tG\t50\tPASS\tAF=0.10\n"
            )
        sites.write_text("".join(site_lines), encoding="ascii")

        get_output = work / "pileups.table"
        get_manifest = work / "pileups.manifest.json"
        get_result = run_repeated([
            str(get_pileups), "-I", str(reads), "-V", str(sites), "-L", str(sites),
            "-O", str(get_output), "--output-manifest", str(get_manifest), "--threads", "1",
        ], get_output, get_manifest)

        shard_a = work / "shard-a.table"
        shard_b = work / "shard-b.table"
        dictionary = work / "reference.dict"
        dictionary.write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n", encoding="ascii")
        header = "#<METADATA>SAMPLE=TUMOR\ncontig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n"
        shard_a.write_text(header + "chr1\t2\t1\t1\t0\t0.1\n", encoding="ascii")
        shard_b.write_text(header + "chr1\t10\t1\t1\t0\t0.1\n", encoding="ascii")
        gather_output = work / "gathered.table"
        gather_manifest = work / "gathered.manifest.json"
        gather_result = run_repeated([
            str(gather), "-I", str(shard_a), "-I", str(shard_b),
            "-SD", str(dictionary),
            "-O", str(gather_output), "--output-manifest", str(gather_manifest),
        ], gather_output, gather_manifest)

        reference = work / "reference.fa"
        reference.write_text(
            ">chr1\n" + ("A" * 100000) + "\n", encoding="ascii"
        )
        (work / "reference.fa.fai").write_text(
            "chr1\t100000\t6\t100000\t100001\n", encoding="ascii"
        )
        validation_vcf = work / "validation.vcf"
        validation_lines = [
            "##fileformat=VCFv4.2\n",
            "##contig=<ID=chr1,length=100000>\n",
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=Allele count>\n",
            "##INFO=<ID=AN,Number=1,Type=Integer,Description=Allele number>\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
        ]
        for index in range(128):
            position = 1 + index * 8
            validation_lines.append(
                f"chr1\t{position}\trs{index}\tA\tC\t50\tPASS\tAC=1;AN=2\tGT\t0/1\n"
            )
        validation_vcf.write_text("".join(validation_lines), encoding="ascii")
        validation_output = work / "validation.report.tsv"
        validation_manifest = work / "validation.manifest.json"
        validation_result = run_repeated([
            str(validate), "-V", str(validation_vcf), "-R", str(reference),
            "-O", str(validation_output), "--output-manifest", str(validation_manifest),
        ], validation_output, validation_manifest)

        print(json.dumps({
            "schema_version": 1,
            "suite": "fastgatk-pileup-validation-file-boundary",
            "backend": "Kokkos",
            "status": "pass",
            "warmup": 1,
            "repetitions": 6,
            "cases": {
                "GetPileupSummaries": get_result,
                "GatherPileupSummaries": gather_result,
                "ValidateVariants": validation_result,
            },
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
