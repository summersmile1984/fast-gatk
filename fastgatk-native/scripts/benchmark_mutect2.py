#!/usr/bin/env python3
"""File-boundary benchmark for Mutect2's aggregate or region-streamed path."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
import time
import argparse
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stream-by-region", type=int, default=0,
                        help="benchmark indexed Mutect2 core/halo streaming with this tile size")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    native = root / "fastgatk-native/build/fastgatk-mutect2"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    repeats = max(1, int(os.environ.get("FASTGATK_MUTECT2_BENCH_REPEATS", "2")))
    env = os.environ.copy()
    env.setdefault("OMP_PROC_BIND", "true")
    env.setdefault("OMP_PLACES", "threads")
    durations: list[float] = []
    summaries: list[dict] = []
    with tempfile.TemporaryDirectory(prefix="fastgatk-mutect2-bench-") as directory:
        work = Path(directory)
        for repeat in range(repeats):
            output = work / f"calls-{repeat}.vcf.gz"
            stats = work / f"calls-{repeat}.stats.json"
            manifest = work / f"calls-{repeat}.manifest.json"
            sidecar = work / f"calls-{repeat}.f1r2.tar.gz"
            command = [
                str(native), "-I", str(bam), "-R", str(reference),
                "-L", "17:69000-70000", "-O", str(output),
                "--normal-input", str(bam),
                "--min-depth", "1", "--min-alt-support", "1",
                "--kmer-size", "7", "--min-kmer-count", "1",
                "--max-num-haplotypes-in-population", "8", "--max-haplotype-depth", "64",
                "--dont-increase-kmer-sizes-for-cycles",
                "--max-haplotype-combination-alleles", "5",
                "--stats", str(stats), "--f1r2-tar-gz", str(sidecar),
                "--output-manifest", str(manifest),
            ]
            if args.stream_by_region > 0:
                command += ["--stream-by-region", str(args.stream_by_region)]
            started = time.perf_counter()
            result = subprocess.run(command, text=True, env=env, capture_output=True, check=True)
            durations.append(time.perf_counter() - started)
            summary = json.loads(result.stdout.splitlines()[-1])
            summaries.append(summary)
            metadata = json.loads(manifest.read_text(encoding="utf-8"))
            assert metadata["compatibility"]["f1r2_standard_tar"] is True
            if args.stream_by_region > 0:
                assert metadata["compatibility"]["stream_by_region"] is True
                assert metadata["telemetry"]["pipeline_decoded_items"] == \
                    metadata["telemetry"]["streamed_regions"]
            assert sidecar.stat().st_size > 0 and output.stat().st_size > 0
        durations.sort()
        summary = summaries[-1]
        p50 = durations[len(durations) // 2]
        stats = json.loads((work / f"calls-{repeats - 1}.stats.json").read_text(encoding="utf-8"))
        print(json.dumps({
            "status": "pass", "tool": "Mutect2", "backend": "Kokkos",
            "f1r2_format": stats["f1r2_format"], "repeats": repeats,
            "tumor_calls": summary["tumor_calls"], "normal_calls": summary["normal_calls"],
            "wall_seconds": sum(durations) / len(durations), "p50_wall_seconds": p50,
            "calls_per_second": summary["tumor_calls"] / p50,
            "sidecar_bytes": (work / f"calls-{repeats - 1}.f1r2.tar.gz").stat().st_size,
            "somatic_evidence_groups": stats["somatic_evidence_groups"],
            "somatic_evidence_grouping": stats["somatic_evidence_grouping"],
            "somatic_likelihood_model": stats["somatic_likelihood_model"],
            "somatic_likelihood_execution_space": stats["somatic_likelihood_execution_space"],
            "somatic_posterior_execution_space": stats["somatic_posterior_execution_space"],
            "pcr_indel_model": stats.get("pcr_indel_model", ""),
            "pcr_error_rate_factor": stats.get("pcr_error_rate_factor", 0.0),
            "pairhmm_pcr_adjusted_positions": stats.get("pairhmm_pcr_adjusted_positions", 0),
            "stream_by_region": args.stream_by_region > 0,
            "streamed_regions": stats.get("streamed_regions", 0),
            "streamed_peak_host_bytes": stats.get("streamed_peak_host_bytes", 0),
            "pipeline_decoded_items": stats.get("pipeline_decoded_items", 0),
            "pipeline_peak_computed_bytes": stats.get("pipeline_peak_computed_bytes", 0),
        }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
