#!/usr/bin/env python3
"""File-boundary HC benchmark for the broad assembly oracle fixture."""
from __future__ import annotations

import json
import hashlib
import os
import subprocess
import tempfile
import time
from pathlib import Path


def run_native(binary: Path, bam: Path, reference: Path, work: Path, repeat: int,
               *, graph: bool = False, stream_region: int | None = None,
               double_precision: bool | None = None,
               rcm_realignment: bool = False) -> dict:
    if graph:
        prefix = "native-graph"
    elif stream_region is not None:
        prefix = "native-region"
    elif double_precision is False:
        prefix = "native-float32"
    elif double_precision is True:
        prefix = "native-double"
    elif rcm_realignment:
        prefix = "native-rcm-realignment"
    else:
        prefix = "native"
    output = work / f"{prefix}-{repeat}.vcf"
    manifest = work / f"{prefix}-{repeat}.json"
    command = [
        str(binary), "-I", str(bam), "-R", str(reference), "-L", "17:69000-70000",
        "-O", str(output), "--output-manifest", str(manifest), "--threads", "2",
    ]
    if graph:
        command += ["--linked-de-bruijn-graph", "--adaptive-pruning",
                    "--num-pruning-samples", "1", "--min-pruning", "1",
                    "--allow-non-unique-kmers-in-ref"]
    if stream_region is not None:
        command += ["--stream-by-region", str(stream_region)]
    if double_precision is not None:
        command += ["--native-pair-hmm-use-double-precision",
                    "true" if double_precision else "false"]
    if rcm_realignment:
        command += ["--use-haplotype-realignment-for-rcm"]
    start = time.perf_counter()
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    elapsed = time.perf_counter() - start
    metadata = json.loads(manifest.read_text(encoding="utf-8"))
    telemetry = metadata["telemetry"]
    calls = sum(1 for line in output.read_text(encoding="utf-8").splitlines()
                if line and not line.startswith("#"))
    return {
        "seconds": elapsed,
        "output_bytes": output.stat().st_size,
        "calls": calls,
        "candidate_sites": telemetry["candidate_sites"],
        "candidate_softclip_suppressed": telemetry["candidate_softclip_suppressed"],
        "candidate_fragment_suppressed": telemetry["candidate_fragment_suppressed"],
        "candidate_low_support_suppressed": telemetry.get(
            "candidate_low_support_suppressed", 0),
        "rcm_haplotype_realignment_used": telemetry.get(
            "rcm_haplotype_realignment_used", False),
        "rcm_realigned_observations": telemetry.get(
            "rcm_realigned_observations", 0),
        "rcm_realignment_fallback_observations": telemetry.get(
            "rcm_realignment_fallback_observations", 0),
        "pairhmm_pairs": telemetry["pairhmm_pairs"],
        "pairhmm_precision": telemetry.get("pairhmm_effective_precision", ""),
        "pairhmm_pcr_indel_model": telemetry.get("pairhmm_pcr_indel_model", ""),
        "pairhmm_pcr_error_rate_factor": telemetry.get("pairhmm_pcr_error_rate_factor", 0.0),
        "pairhmm_pcr_adjusted_positions": telemetry.get("pairhmm_pcr_adjusted_positions", 0),
        "assembly_region_union_count": telemetry.get("assembly_region_union_count", 0),
        "assembly_region_union_merges": telemetry.get("assembly_region_union_merges", 0),
        "assembly_cross_region_rescued": telemetry.get("assembly_cross_region_rescued", 0),
        "read_haplotype_cigar_pairs": telemetry.get("read_haplotype_cigar_pairs", 0),
        "read_haplotype_cigar_filtered_pairs": telemetry.get(
            "read_haplotype_cigar_filtered_pairs", 0),
        "read_haplotype_softclip_filtered_pairs": telemetry.get(
            "read_haplotype_softclip_filtered_pairs", 0),
        "pairhmm_flow_reads": telemetry.get("pairhmm_flow_reads", 0),
        "pairhmm_flow_reads_clipped": telemetry.get("pairhmm_flow_reads_clipped", 0),
        "pairhmm_flow_haplotypes_clipped": telemetry.get(
            "pairhmm_flow_haplotypes_clipped", 0),
        "pairhmm_flow_haplotypes_collapsed": telemetry.get(
            "pairhmm_flow_haplotypes_collapsed", 0),
        "pairhmm_flow_haplotypes_uncollapsed": telemetry.get(
            "pairhmm_flow_haplotypes_uncollapsed", 0),
        "pairhmm_flow_haplotype_remaps": telemetry.get(
            "pairhmm_flow_haplotype_remaps", 0),
        "pairhmm_flow_identical_haplotype_groups": telemetry.get(
            "pairhmm_flow_identical_haplotype_groups", 0),
        "pairhmm_graph_snp_posterior_pairs": telemetry.get(
            "pairhmm_graph_snp_posterior_pairs", 0),
        "pairhmm_normalization_execution_space": telemetry.get(
            "pairhmm_normalization_execution_space", ""),
        "pairhmm_normalization_prepare_seconds": telemetry.get(
            "pairhmm_normalization_prepare_seconds", 0.0),
        "pairhmm_normalization_seconds": telemetry.get(
            "pairhmm_normalization_seconds", 0.0),
        "pairhmm_marginalization_execution_space": telemetry.get(
            "pairhmm_marginalization_execution_space", ""),
        "pairhmm_marginalization_prepare_seconds": telemetry.get(
            "pairhmm_marginalization_prepare_seconds", 0.0),
        "pairhmm_marginalization_seconds": telemetry.get(
            "pairhmm_marginalization_seconds", 0.0),
        "pairhmm_uncertainty_execution_space": telemetry.get(
            "pairhmm_uncertainty_execution_space", ""),
        "pairhmm_uncertainty_prepare_seconds": telemetry.get(
            "pairhmm_uncertainty_prepare_seconds", 0.0),
        "pairhmm_uncertainty_seconds": telemetry.get(
            "pairhmm_uncertainty_seconds", 0.0),
        "read_haplotype_uncertain_reads": telemetry.get("read_haplotype_uncertain_reads", 0),
        "read_haplotype_informative_reads": telemetry.get("read_haplotype_informative_reads", 0),
        "graph_haplotype_provenance_count": telemetry.get("graph_haplotype_provenance_count", 0),
        "graph_haplotype_cigar_signature": telemetry.get("graph_haplotype_cigar_signature", 0),
        "graph_haplotype_sw_used": telemetry.get("graph_haplotype_sw_used", False),
        "graph_haplotype_sw_execution_space": telemetry.get("graph_haplotype_sw_execution_space", ""),
        "graph_haplotype_sw_simd_width": telemetry.get("graph_haplotype_sw_simd_width", 1),
        "graph_haplotype_sw_simd_groups": telemetry.get("graph_haplotype_sw_simd_groups", 0),
        "graph_haplotype_sw_seconds": telemetry.get("graph_haplotype_sw_seconds", 0.0),
        "execution_space": telemetry["execution_space"],
        "reference_block_lookup_indexed": telemetry.get(
            "reference_block_lookup_indexed", False),
        "pipeline_lifecycle": telemetry.get("pipeline_lifecycle", "not-used"),
        "pipeline_decoded_items": telemetry.get("pipeline_decoded_items", 0),
        "pipeline_computed_items": telemetry.get("pipeline_computed_items", 0),
        "pipeline_encoded_items": telemetry.get("pipeline_encoded_items", 0),
        "pipeline_decoded_bytes": telemetry.get("pipeline_decoded_bytes", 0),
        "pipeline_computed_bytes": telemetry.get("pipeline_computed_bytes", 0),
        "pipeline_encoded_bytes": telemetry.get("pipeline_encoded_bytes", 0),
        "pipeline_peak_decoded_bytes": telemetry.get("pipeline_peak_decoded_bytes", 0),
        "pipeline_peak_computed_bytes": telemetry.get("pipeline_peak_computed_bytes", 0),
        "pipeline_peak_encoded_bytes": telemetry.get("pipeline_peak_encoded_bytes", 0),
        "stream_region": stream_region or 0,
        "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
    }


def run_gatk(java: Path, jar: Path, bam: Path, reference: Path, work: Path,
             repeat: int) -> dict:
    output = work / f"gatk-{repeat}.vcf"
    start = time.perf_counter()
    subprocess.run([
        str(java), "-Xmx1g", "-jar", str(jar), "HaplotypeCaller",
        "-R", str(reference), "-I", str(bam), "-L", "17:69000-70000",
        "-O", str(output), "--native-pair-hmm-threads", "2",
        "--create-output-variant-index", "false", "--seconds-between-progress-updates", "1",
    ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    elapsed = time.perf_counter() - start
    calls = sum(1 for line in output.read_text(encoding="utf-8").splitlines()
                if line and not line.startswith("#"))
    return {"seconds": elapsed, "output_bytes": output.stat().st_size, "calls": calls}


def summarize(samples: list[dict]) -> dict:
    timings = sorted(item["seconds"] for item in samples)
    return {
        "repetitions": len(samples),
        "warmup_seconds": timings[0],
        "p50_seconds": timings[len(timings) // 2],
        "p95_seconds": timings[-1],
        "samples": samples,
    }


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = Path(os.environ.get("FASTGATK_HC_BINARY", root / "fastgatk-native/build/fastgatk-hc-call"))
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    if not binary.exists() or not bam.exists() or not reference.exists():
        raise SystemExit("missing native broad HC benchmark assets")
    java = Path(os.environ.get("JAVA", root / "third_party/jdk17/bin/java"))
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-broad-benchmark-") as directory:
        work = Path(directory)
        native = [run_native(binary, bam, reference, work, repeat) for repeat in range(3)]
        native_float32 = [run_native(binary, bam, reference, work, repeat,
                                      double_precision=False) for repeat in range(2)]
        native_double = [run_native(binary, bam, reference, work, repeat,
                                     double_precision=True) for repeat in range(2)]
        native_graph = [run_native(binary, bam, reference, work, repeat, graph=True)
                        for repeat in range(2)]
        native_region = [run_native(binary, bam, reference, work, repeat, stream_region=500)
                         for repeat in range(2)]
        native_rcm_realignment = []
        if os.environ.get("FASTGATK_BENCHMARK_RCM_REALIGNMENT", "0") == "1":
            native_rcm_realignment = [run_native(
                binary, bam, reference, work, repeat, rcm_realignment=True)
                for repeat in range(2)]
        if not all(item["reference_block_lookup_indexed"]
                   for item in (*native, *native_graph, *native_region)):
            raise SystemExit("reference-block lookup index was not exercised")
        if not all(item["pipeline_lifecycle"] ==
                   "Host decode->bounded queue->Kokkos compute->encode->sink"
                   and item["pipeline_decoded_items"] == item["pipeline_computed_items"]
                   == item["pipeline_encoded_items"] > 1
                   for item in native_region):
            raise SystemExit("region streaming pipeline telemetry was not exercised")
        # Region tiles intentionally recompute activity/graph context per
        # core+halo.  Their per-site evidence can differ from the unbounded
        # aggregate path, so benchmark determinism within the tiled mode while
        # reporting (rather than asserting) aggregate SHA equality.
        region_sha_exact = all(item["output_sha256"] == native_region[0]["output_sha256"]
                               for item in native_region)
        if not all(item["pairhmm_precision"] == "float32" for item in native_float32):
            raise SystemExit("explicit Float32 HC benchmark did not select the Kokkos float path")
        if not all(item["pairhmm_precision"] == "double" for item in native_double):
            raise SystemExit("explicit double HC benchmark did not select the strict path")
        result = {"schema_version": 1, "status": "pass", "tool": "HaplotypeCaller",
                  "benchmark": "broad-assembly-file-boundary", "region": "17:69000-70000",
                  "native": summarize(native),
                  "native_float32": summarize(native_float32),
                  "native_double": summarize(native_double),
                  "native_graph": summarize(native_graph),
                  "native_region_streaming": summarize(native_region),
                  "region_vs_aggregate_sha256_exact": all(
                      item["output_sha256"] == native[0]["output_sha256"]
                      for item in native_region),
                  "region_repeat_sha256_exact": region_sha_exact,
                  "note": "wall time includes BAM/FASTA decode, Host projection, Kokkos assembly/PairHMM, VCF writing; not a Java speedup claim"}
        if native_rcm_realignment:
            result["native_rcm_realignment"] = summarize(native_rcm_realignment)
        if jar.exists() and java.exists() and os.environ.get("FASTGATK_SKIP_GATK_BENCHMARK", "0") != "1":
            gatk = [run_gatk(java, jar, bam, reference, work, repeat) for repeat in range(2)]
            result["gatk_java"] = summarize(gatk)
        print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
