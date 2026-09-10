#!/usr/bin/env python3
"""Small reproducible file-boundary benchmark for native HC/Mutect2 adapters.

This is deliberately separate from kernel-only benchmarks: it reports input
decode, Kokkos prepare/execute, VCF compression/index and total wall time in
one artifact, so prepare/IO cannot be hidden behind a kernel speedup claim.
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
import time
from pathlib import Path


def run_case(name: str, command: list[str], output: Path, work: Path,
             summary_json: bool = True) -> dict:
    env = os.environ.copy()
    env.setdefault("OMP_PROC_BIND", "true")
    env.setdefault("OMP_PLACES", "threads")
    timing = work / (name.lower().replace(" ", "-") + ".time")
    start = time.perf_counter()
    timed_command = ["/usr/bin/time", "-f", "%e\\t%M\\t%I\\t%O",
                     "-o", str(timing), *command]
    completed = subprocess.run(timed_command, check=False, text=True,
                               capture_output=True, env=env)
    if completed.returncode != 0:
        raise RuntimeError(
            f"{name} exited with status {completed.returncode}\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}")
    wall = time.perf_counter() - start
    values = timing.read_text(encoding="utf-8").strip().split("\t")
    measured = {
        "max_rss_kb": int(values[1]) if len(values) > 1 else 0,
        "filesystem_inputs": int(values[2]) if len(values) > 2 else 0,
        "filesystem_outputs": int(values[3]) if len(values) > 3 else 0,
    }
    result = {
        "name": name,
        "wall_seconds": wall,
        "output_bytes": output.stat().st_size if output.exists() else 0,
        "resource": measured,
    }
    if summary_json:
        result["summary"] = json.loads(completed.stdout.strip().splitlines()[-1])
    manifest = Path(str(output) + ".manifest.json")
    if manifest.exists():
        metadata = json.loads(manifest.read_text(encoding="utf-8"))
        telemetry = metadata.get("telemetry", {})
        result["kernel_telemetry"] = {
            key: telemetry[key] for key in (
                "pairhmm_pairs", "pairhmm_request_links", "pairhmm_haplotypes",
                "pairhmm_haplotype_combination_blocks",
                "pairhmm_haplotype_combination_chunked", "pairhmm_cache_hits",
                "pairhmm_assembly_region_groups", "pairhmm_unassigned_candidates",
                "pairhmm_assembly_region_partitioned", "assembly_unassigned_candidates",
                "assembly_cross_region_candidates", "assembly_region_partitioned",
                "pairhmm_reads_clipped", "pairhmm_reads_dropped_after_clipping",
                "downsampled_reads", "read_filter_max_reads_per_locus",
                "sw_pairs", "sw_simd_width", "sw_simd_groups",
                "stream_by_contig", "streamed_contigs", "streamed_peak_host_bytes",
                "stream_by_region", "streamed_regions", "streamed_region_splits",
                "stream_indexed",
            ) if key in telemetry
        }
    return result


def run_gatk_case(name: str, java: Path, gatk_jar: Path, arguments: list[str],
                  output: Path, work: Path) -> dict:
    result = run_case(name, [str(java), "-Xmx1g", "-jar", str(gatk_jar), *arguments],
                      output, work, summary_json=False)
    result["implementation"] = "gatk-java"
    return result


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    native = root / "fastgatk-native/build"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    with tempfile.TemporaryDirectory(prefix="fastgatk-e2e-benchmark-") as directory:
        work = Path(directory)
        hc_output = work / "hc.g.vcf.gz"
        hc_stream_output = work / "hc-stream.g.vcf.gz"
        hc_region_output = work / "hc-region-stream.vcf.gz"
        mutect_output = work / "mutect.vcf.gz"
        mutect_filtered_output = work / "mutect.filtered.vcf.gz"
        hc = run_case("HaplotypeCaller", [
            str(native / "fastgatk-hc-call"), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-70000", "-O", str(hc_output), "-ERC", "GVCF",
            "--min-depth", "1", "--min-alt-support", "1", "--threads", "4",
        ], hc_output, work)
        hc_stream = run_case("HaplotypeCaller-stream-by-contig", [
            str(native / "fastgatk-hc-call"), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-70000", "-O", str(hc_stream_output), "-ERC", "GVCF",
            "--min-depth", "1", "--min-alt-support", "1", "--threads", "4",
            "--stream-by-contig",
        ], hc_stream_output, work)
        hc_region_stream = run_case("HaplotypeCaller-stream-by-region", [
            str(native / "fastgatk-hc-call"), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-70000", "-O", str(hc_region_output),
            "--min-depth", "1", "--min-alt-support", "1", "--threads", "4",
            "--stream-by-region", "500",
        ], hc_region_output, work)
        mutect = run_case("Mutect2", [
            str(native / "fastgatk-mutect2"), "-I", str(bam), "-R", str(reference),
            "-L", "17:69000-70000", "-O", str(mutect_output), "--normal-input", str(bam),
            "--min-depth", "1", "--min-alt-support", "1", "--threads", "4",
        ], mutect_output, work)
        mutect_filtered = run_case("FilterMutectCalls", [
            str(native / "fastgatk-filter-mutect-calls"), "-V", str(mutect_output),
            "-O", str(mutect_filtered_output), "--min-tlod", "0",
        ], mutect_filtered_output, work)
        gatk_jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
        gatk_java = root / "third_party/jdk17/bin/java"
        gatk_cases = []
        if gatk_jar.exists() and gatk_java.exists() and os.environ.get(
                "FASTGATK_SKIP_GATK_BENCHMARK", "0") != "1":
            gatk_hc_output = work / "gatk.hc.g.vcf.gz"
            gatk_cases.append(run_gatk_case("GATK-HaplotypeCaller", gatk_java, gatk_jar, [
                "HaplotypeCaller", "-R", str(reference), "-I", str(bam),
                "-L", "17:69000-70000", "-O", str(gatk_hc_output), "-ERC", "GVCF",
                "--native-pair-hmm-threads", "4", "--create-output-variant-index", "true",
                "--seconds-between-progress-updates", "1",
            ], gatk_hc_output, work))
        artifact = {
            "schema_version": 4,
            "status": "pass",
            "input": str(bam),
            "reference": str(reference),
            "interval": "17:69000-70000",
            "cases": [hc, hc_stream, hc_region_stream, mutect, mutect_filtered],
            "baselines": gatk_cases,
            "note": "file-boundary wall time includes HTSlib decode, Kokkos preparation/execution, compression and Tabix indexing; resource fields come from /usr/bin/time and Java GATK is an optional same-fixture baseline",
        }
        print(json.dumps(artifact, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
