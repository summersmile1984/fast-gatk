#!/usr/bin/env python3
"""File-boundary benchmark for joint GenotypeGVCFs PL projection."""

from __future__ import annotations

import gzip
import json
import os
import pathlib
import statistics
import subprocess
import tempfile
import time
import argparse
import sys


HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length={contig_length}>
##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>
##INFO=<ID=DP,Number=1,Type=Integer,Description=Depth>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stream-by-locus", action="store_true",
                        help="benchmark the bounded joint-locus path")
    parser.add_argument(
        "--records-per-shard", type=int, default=128,
        help=("synthetic loci written to each input shard (default: 128; "
              "use >=4096 for a batch-scale workload; maximum: 100000)"),
    )
    parser.add_argument(
        "--repetitions", type=int, default=5,
        help="measured repetitions after one warmup (default: 5)",
    )
    args = parser.parse_args()
    if not 1 <= args.records_per_shard <= 100_000:
        parser.error("--records-per-shard must be between 1 and 100000")
    if not 1 <= args.repetitions <= 20:
        parser.error("--repetitions must be between 1 and 20")
    # This script intentionally does not invoke Java GATK, so its wall time
    # is a native regression/telemetry measurement, never a speedup claim.
    # Keep the warning visible in CI logs when the default micro workload is
    # used, since tiny per-locus calls are dominated by OpenMP launch costs.
    if args.records_per_shard < 4096:
        print(
            "NOTICE: records_per_shard<4096 is a micro workload; this run "
            "must not be used as a Java/GATK speedup claim (no Java baseline).",
            file=sys.stderr,
        )
    root = pathlib.Path(__file__).resolve().parents[2]
    binary = pathlib.Path(os.environ.get(
        "FASTGATK_GENOTYPE_BINARY", root / "fastgatk-native/build/fastgatk-genotype-gvcf"))
    hc_binary = pathlib.Path(os.environ.get(
        "FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    if not binary.is_file():
        raise SystemExit("native GenotypeGVCFs binary is required")
    with tempfile.TemporaryDirectory(prefix="fastgatk-genotype-gvcf-benchmark-") as directory:
        work = pathlib.Path(directory)
        shard_g = work / "shard-g.g.vcf.gz"
        shard_reverse = work / "shard-reverse.g.vcf.gz"
        contig_length = max(1_000_000, args.records_per_shard * 10 + 100)
        with gzip.open(shard_g, "wt", encoding="ascii") as first, gzip.open(
            shard_reverse, "wt", encoding="ascii"
        ) as second:
            first.write(HEADER.format(sample="S1", contig_length=contig_length))
            second.write(HEADER.format(sample="S2", contig_length=contig_length))
            for index in range(args.records_per_shard):
                position = index * 10 + 1
                first.write(
                    f"chr1\t{position}\t.\tA\tG,<NON_REF>\t.\tPASS\tDP=20\t"
                    "GT:DP:AD:PL\t0/1:20:12,8,0:50,0,60,70,80,90\n"
                )
                # Deliberately reverse the source ALT order.  The joint union
                # is G,T, while this shard is T,G; the Kokkos remap sorts the
                # mapped tuple before computing the source genotype rank.
                second.write(
                    f"chr1\t{position}\t.\tA\tT,G,<NON_REF>\t.\tPASS\tDP=21\t"
                    "GT:DP:AD:PL\t1/2:21:2,9,10,0:90,80,70,60,50,0\n"
                )
        output = work / "output.vcf.gz"
        manifest = work / "output.manifest.json"
        command = [
            str(binary), "-V", str(shard_g), "-V", str(shard_reverse),
            "-O", str(output), "--output-manifest", str(manifest),
        ]
        if args.stream_by_locus:
            command.append("--stream-by-locus")
        samples: list[float] = []
        for iteration in range(args.repetitions + 1):
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
        if telemetry["pl_remap_kernel_calls"] < args.records_per_shard:
            raise SystemExit("joint PL remap kernel was not exercised")
        if telemetry["pipeline_lifecycle"] != (
                "Host decode->bounded queue->Kokkos compute->encode->sink"):
            raise SystemExit("GenotypeGVCFs pipeline was not exercised")
        if not (telemetry["pipeline_decoded_items"] ==
                telemetry["pipeline_computed_items"] ==
                telemetry["pipeline_encoded_items"] == telemetry["output_records"] > 0):
            raise SystemExit("GenotypeGVCFs pipeline item counts are inconsistent")
        if not all(telemetry[key] > 0 for key in (
                "pipeline_peak_decoded_bytes", "pipeline_peak_computed_bytes",
                "pipeline_peak_encoded_bytes")):
            raise SystemExit("GenotypeGVCFs pipeline peaks were not recorded")

        # The synthetic shards above intentionally exercise the ALT-union
        # kernels without an index.  Add a small real BGZF/TBI workload to
        # make the decode-path optimization measurable as well.  The HC
        # fixture produces the index through the same writer used by the
        # production pipeline, so this benchmark does not require a separate
        # tabix executable in the worker image.
        indexed_p50 = None
        indexed_queries = 0
        indexed_inputs = 0
        if hc_binary.is_file():
            bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
            reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
            indexed_gvcf = work / "indexed.g.vcf.gz"
            if bam.is_file() and reference.is_file():
                generated = subprocess.run([
                    str(hc_binary), "-I", str(bam), "-R", str(reference),
                    "-L", "17:69000-69100", "-O", str(indexed_gvcf),
                    "-ERC", "GVCF", "--min-depth", "1", "--min-alt-support", "1",
                ], text=True, capture_output=True, check=False)
                if generated.returncode != 0:
                    raise SystemExit(generated.stderr or generated.stdout)
                indexed_output = work / "indexed-output.vcf.gz"
                indexed_manifest = work / "indexed-output.manifest.json"
                indexed_command = [
                    str(binary), "-V", str(indexed_gvcf), "-L", "17:69000-69100",
                    "-O", str(indexed_output), "--output-manifest", str(indexed_manifest),
                ]
                if args.stream_by_locus:
                    indexed_command.append("--stream-by-locus")
                indexed_samples: list[float] = []
                for iteration in range(4):
                    started = time.perf_counter()
                    completed = subprocess.run(indexed_command, text=True,
                                               capture_output=True, check=False)
                    elapsed = time.perf_counter() - started
                    if completed.returncode != 0:
                        raise SystemExit(completed.stderr or completed.stdout)
                    if iteration:
                        indexed_samples.append(elapsed)
                indexed_metadata = json.loads(indexed_manifest.read_text(encoding="utf-8"))
                indexed_telemetry = indexed_metadata["telemetry"]
                indexed_inputs = indexed_telemetry["indexed_inputs"]
                indexed_queries = indexed_telemetry["indexed_interval_queries"]
                if indexed_inputs != 1 or indexed_queries < 1 or indexed_telemetry["interval_skipped"] != 0:
                    raise SystemExit("indexed GenotypeGVCFs traversal was not exercised")
                indexed_p50 = statistics.median(indexed_samples)
        print(json.dumps({
            "schema_version": 1,
            "suite": "fastgatk-genotype-gvcf-file-boundary",
            "tool": "GenotypeGVCFs",
            "backend": "Kokkos",
            "status": "pass",
            "warmup": 1,
            "repetitions": len(samples),
            "samples": 2,
            "records_per_shard": args.records_per_shard,
            "workload_class": "micro" if args.records_per_shard < 4096 else "batch",
            "speedup_claim_allowed": False,
            "speedup_claim_blocked_reason": (
                "No matched Java/GATK baseline is run; native wall time is "
                "regression/telemetry only."
            ),
            "p50_seconds": statistics.median(samples),
            "p95_seconds": p95,
            "input_bytes": shard_g.stat().st_size + shard_reverse.stat().st_size,
            "output_bytes": output.stat().st_size,
            "pl_remap_kernel_calls": telemetry["pl_remap_kernel_calls"],
            "pl_remap_kernel_seconds": telemetry["pl_remap_kernel_seconds"],
            "allele_field_remap_kernel_calls": telemetry["allele_field_remap_kernel_calls"],
            "allele_field_remap_kernel_seconds": telemetry["allele_field_remap_kernel_seconds"],
            "genotype_kernel_calls": telemetry["genotype_kernel_calls"],
            "genotype_kernel_execution_space": telemetry["genotype_kernel_execution_space"],
            "genotype_prior_kernel_calls": telemetry["genotype_prior_kernel_calls"],
            "genotype_prior_kernel_seconds": telemetry["genotype_prior_kernel_seconds"],
            "genotype_prior_kernel_execution_space": telemetry["genotype_prior_kernel_execution_space"],
            "cohort_af_kernel_calls": telemetry["cohort_af_kernel_calls"],
            "orphan_spanning_deletion_loci": telemetry.get("orphan_spanning_deletion_loci", 0),
            "pipeline_decoded_items": telemetry["pipeline_decoded_items"],
            "pipeline_computed_items": telemetry["pipeline_computed_items"],
            "pipeline_encoded_items": telemetry["pipeline_encoded_items"],
            "pipeline_decoded_bytes": telemetry["pipeline_decoded_bytes"],
            "pipeline_computed_bytes": telemetry["pipeline_computed_bytes"],
            "pipeline_encoded_bytes": telemetry["pipeline_encoded_bytes"],
            "pipeline_peak_decoded_bytes": telemetry["pipeline_peak_decoded_bytes"],
            "pipeline_peak_computed_bytes": telemetry["pipeline_peak_computed_bytes"],
            "pipeline_peak_encoded_bytes": telemetry["pipeline_peak_encoded_bytes"],
            "pipeline_stage_capacity_bytes": telemetry["pipeline_stage_capacity_bytes"],
            "indexed_interval_traversal": indexed_inputs == 1,
            "indexed_inputs": indexed_inputs,
            "indexed_interval_queries": indexed_queries,
            "indexed_interval_p50_seconds": indexed_p50,
            "stream_by_locus": args.stream_by_locus,
            "streamed_loci": telemetry.get("streamed_loci", 0),
            "streamed_peak_host_bytes": telemetry.get("streamed_peak_host_bytes", 0),
            "stream_max_inflight_records": telemetry.get("stream_max_inflight_records", 0),
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
