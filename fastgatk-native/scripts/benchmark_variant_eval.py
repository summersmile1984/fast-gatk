#!/usr/bin/env python3
"""File-boundary benchmark for VariantEval standard strata/concordance."""

from __future__ import annotations

import json
import os
import pathlib
import statistics
import subprocess
import tempfile
import time
import argparse


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##contig=<ID=chr2,length=1000000>
##FILTER=<ID=LowQual,Description=low quality>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
"""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gatk-report", action="store_true",
                        help="benchmark the GATKReport v1.1-compatible output mode")
    args = parser.parse_args()
    root = pathlib.Path(__file__).resolve().parents[2]
    binary = pathlib.Path(os.environ.get(
        "FASTGATK_VARIANT_EVAL_BINARY",
        str(root / "fastgatk-native/build/fastgatk-variant-eval")))
    if not binary.is_file():
        raise SystemExit("native VariantEval binary is required")
    with tempfile.TemporaryDirectory(prefix="fastgatk-variant-eval-benchmark-") as directory:
        work = pathlib.Path(directory)
        evaluation = work / "evaluation.vcf"
        comparison = work / "comparison.vcf"
        pedigree = work / "pedigree.ped"
        pedigree.write_text("fam S1 0 0 1 1\nfam S2 0 0 2 1\nfam S3 S1 S2 1 1\n",
                            encoding="ascii")
        eval_lines = [HEADER]
        comp_lines = [HEADER]
        for index in range(256):
            chrom = "chr1" if index % 2 == 0 else "chr2"
            pos = index * 10 + 1
            if index % 5 == 0:
                ref, alt, genotype, filt = "A", "AT", "0/1", "PASS"
            elif index % 7 == 0:
                ref, alt, genotype, filt = "C", "G", "1/1", "LowQual"
            else:
                ref, alt, genotype, filt = "A", "G", "0/1", "PASS"
            eval_lines.append(
                f"{chrom}\t{pos}\tvar{index}\t{ref}\t{alt}\t50\t{filt}\t.\tGT\t{genotype}\t0/0\t0/1\n"
            )
            if index % 3 == 0:
                comp_lines.append(
                    f"{chrom}\t{pos}\tvar{index}\t{ref}\t{alt}\t50\tPASS\t.\tGT\t{genotype}\t0/0\t0/1\n"
                )
        evaluation.write_text("".join(eval_lines), encoding="ascii")
        comparison.write_text("".join(comp_lines), encoding="ascii")
        output = work / "evaluation.report"
        manifest = work / "evaluation.manifest.json"
        command = [
            str(binary), "-eval", str(evaluation), "-comp", str(comparison),
            "-O", str(output), "--output-manifest", str(manifest),
            "-EV", "GenotypeConcordance", "-S", "Contig", "-S", "VariantType",
            "-EV", "VariantAFEvaluator",
            "-EV", "ThetaVariantEvaluator",
            "-EV", "MendelianViolationEvaluator", "--pedigree", str(pedigree),
        ]
        if args.gatk_report:
            command.append("--gatk-report")
        samples: list[float] = []
        for iteration in range(6):
            started = time.perf_counter()
            completed = subprocess.run(command, text=True, capture_output=True, check=False)
            elapsed = time.perf_counter() - started
            if completed.returncode != 0:
                raise SystemExit(completed.stderr or completed.stdout)
            if iteration:
                samples.append(elapsed)
        ordered = sorted(samples)
        p95 = ordered[min(len(ordered) - 1, int(0.95 * len(ordered)))]
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        telemetry = metadata["telemetry"]
        print(json.dumps({
            "schema_version": 1,
            "suite": "fastgatk-variant-eval-file-boundary",
            "tool": "VariantEval",
            "backend": "HTSlib+Host",
            "status": "pass",
            "warmup": 1,
            "repetitions": len(samples),
            "evaluation_records": 256,
            "comparison_records": 86,
            "p50_seconds": statistics.median(samples),
            "p95_seconds": p95,
            "input_bytes": evaluation.stat().st_size + comparison.stat().st_size,
            "output_bytes": output.stat().st_size,
            "compared_genotypes": telemetry["compared_genotypes"],
            "concordant_genotypes": telemetry["concordant_genotypes"],
            "sample_names": telemetry.get("sample_names", 0),
            "stratification_tables": telemetry["stratification_tables"],
            "variant_af_called_sites": telemetry.get("variant_af_called_sites", 0),
            "theta_num_sites": telemetry.get("theta_num_sites", 0),
            "mendelian_n_variants": telemetry.get("mendelian_n_variants", 0),
            "report_format": "gatk-v1.1" if args.gatk_report else "native-tabular",
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
