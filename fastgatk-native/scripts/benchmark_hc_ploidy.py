#!/usr/bin/env python3
"""File-boundary benchmark for reference-backed HC ploidy modes."""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
import time
from pathlib import Path


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    binary = root / "fastgatk-native/build/fastgatk-hc-call"
    bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
    reference = root / "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    os.environ.setdefault("OMP_PROC_BIND", "true")
    os.environ.setdefault("OMP_PLACES", "threads")
    cases = []
    with tempfile.TemporaryDirectory(prefix="fastgatk-hc-ploidy-benchmark-") as directory:
        work = Path(directory)
        for ploidy in (2, 3, 8):
            timings = []
            output_bytes = 0
            manifest_data = {}
            for repeat in range(3):
                output = work / f"ploidy-{ploidy}-{repeat}.vcf"
                manifest = work / f"ploidy-{ploidy}-{repeat}.manifest.json"
                start = time.perf_counter()
                subprocess.run([
                    str(binary), "-I", str(bam), "-R", str(reference), "-L", "17:69000-69100",
                    "-O", str(output), "--threads", "2", "--min-depth", "1",
                    "--min-alt-support", "1", "--sample-ploidy", str(ploidy),
                    "--output-manifest", str(manifest),
                ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                timings.append(time.perf_counter() - start)
                output_bytes = output.stat().st_size
                manifest_data = json.loads(manifest.read_text(encoding="utf-8"))
            timings.sort()
            cases.append({"sample_ploidy": ploidy, "warmup_seconds": timings[0],
                          "mode": "VCF",
                          "p50_seconds": timings[1], "p95_seconds": timings[-1],
                          "output_bytes": output_bytes, "repetitions": len(timings),
                          "execution_space": manifest_data["telemetry"]["execution_space"]})
        for ploidy in (2, 3):
            timings = []
            output_bytes = 0
            manifest_data = {}
            for repeat in range(2):
                output = work / f"gvcf-ploidy-{ploidy}-{repeat}.g.vcf.gz"
                manifest = work / f"gvcf-ploidy-{ploidy}-{repeat}.manifest.json"
                start = time.perf_counter()
                subprocess.run([
                    str(binary), "-I", str(bam), "-R", str(reference), "-L", "17:69000-69100",
                    "-O", str(output), "--emit-ref-confidence", "GVCF", "--threads", "2",
                    "--min-depth", "1", "--min-alt-support", "1", "--sample-ploidy", str(ploidy),
                    "--output-manifest", str(manifest),
                ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                timings.append(time.perf_counter() - start)
                output_bytes = output.stat().st_size
                manifest_data = json.loads(manifest.read_text(encoding="utf-8"))
            timings.sort()
            cases.append({"sample_ploidy": ploidy, "warmup_seconds": timings[0],
                          "mode": "GVCF", "p50_seconds": timings[-1],
                          "p95_seconds": timings[-1], "output_bytes": output_bytes,
                          "repetitions": len(timings),
                          "execution_space": manifest_data["telemetry"]["execution_space"],
                          "reference_confidence_execution_space": manifest_data["telemetry"]["reference_confidence_execution_space"],
                          "genotype_priors": manifest_data["telemetry"]["genotype_priors"],
                          "genotype_assignment_method": manifest_data["telemetry"]["genotype_assignment_method"],
                          "joint_genotype_priors_used": manifest_data["telemetry"]["joint_genotype_priors_used"],
                          "gvcf_standard_fields": manifest_data["compatibility"]["gvcf_standard_fields"],
                          "gvcf_candidate_standard_fields": manifest_data["compatibility"]["gvcf_candidate_standard_fields"],
                          "gvcf_candidate_rebuild_after_genotyping": manifest_data["compatibility"]["gvcf_candidate_rebuild_after_genotyping"]})
        # The MNP distance is a graph-SW grouping policy, so benchmark it at
        # the same file boundary as the ploidy cases.  Recording the emitted
        # graph_mnp_candidates prevents a timing-only result from hiding a
        # no-op CLI path.
        for distance in (0, 1, 2):
            timings = []
            output_bytes = 0
            manifest_data = {}
            for repeat in range(2):
                output = work / f"mnp-{distance}-{repeat}.vcf"
                manifest = work / f"mnp-{distance}-{repeat}.manifest.json"
                start = time.perf_counter()
                subprocess.run([
                    str(binary), "-I", str(bam), "-R", str(reference), "-L", "17:69000-69100",
                    "-O", str(output), "--threads", "2", "--min-depth", "1",
                    "--min-alt-support", "1", "--max-mnp-distance", str(distance),
                    "--output-manifest", str(manifest),
                ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                timings.append(time.perf_counter() - start)
                output_bytes = output.stat().st_size
                manifest_data = json.loads(manifest.read_text(encoding="utf-8"))
            timings.sort()
            telemetry = manifest_data["telemetry"]
            cases.append({"max_mnp_distance": distance, "warmup_seconds": timings[0],
                          "mode": "MNP", "p50_seconds": timings[-1],
                          "p95_seconds": timings[-1], "output_bytes": output_bytes,
                          "graph_mnp_candidates": telemetry["graph_mnp_candidates"],
                          "execution_space": telemetry["execution_space"],
                          "repetitions": len(timings)})
    print(json.dumps({"schema_version": 1, "status": "pass", "tool": "HaplotypeCaller",
                      "benchmark": "ploidy-file-boundary", "cases": cases}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
