#!/usr/bin/env python3
"""Run the shared Kokkos kernel smoke across available execution spaces."""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path


def run_json(path: Path, *, benchmark: bool = False) -> dict:
    env = os.environ.copy()
    if benchmark:
        # Keep the matrix check short while exercising the exact same mixed
        # workload/checksum contract as the standalone benchmark.
        env.update({"FASTGATK_BENCH_WARMUP": "1", "FASTGATK_BENCH_REPEATS": "3"})
    result = subprocess.run([str(path)], text=True, capture_output=True, check=False, env=env)
    if result.returncode != 0:
        raise RuntimeError(f"{path}: {result.stderr}")
    return json.loads(result.stdout.strip().splitlines()[-1])


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    builds = [("openmp", root / "fastgatk-native/build"),
              ("serial", root / "fastgatk-native/build-serial")]
    reports = []
    for name, build in builds:
        api = build / "fastgatk-kernels" / "fastgatk-kernels-api-smoke"
        if not api.is_file():
            api = build / "fastgatk-kernels-api-smoke"
        if not api.is_file():
            continue
        report = run_json(api)
        if report.get("status") != "pass":
            raise RuntimeError(f"{name} backend smoke failed: {report}")
        entry = {"backend": name, "status": report["status"],
                        "execution_space": report.get("execution_space", name),
                        "sw_score": report.get("sw_kokkos_score", report.get("sw_score"))}
        benchmark = build / "fastgatk-kernels" / "fastgatk-kernels-benchmark"
        if not benchmark.is_file():
            benchmark = build / "fastgatk-kernels-benchmark"
        if benchmark.is_file():
            bench = run_json(benchmark, benchmark=True)
            if bench.get("status") != "pass":
                raise RuntimeError(f"{name} benchmark failed: {bench}")
            entry["benchmark"] = {
                "schema_version": bench.get("schema_version"),
                "simd_width": bench.get("simd_width"),
                "sw_simd_width": bench.get("sw_simd_width"),
                "pairhmm_checksum": bench.get("pairhmm_checksum"),
                "float_pairhmm_checksum": bench.get("float_pairhmm_checksum"),
                "float_pairhmm_max_abs_error": bench.get("float_pairhmm_max_abs_error"),
                "float_pairhmm_simd_width": bench.get("float_pairhmm_simd_width"),
                "sw_checksum": bench.get("sw_checksum"),
                "genotype_checksum": bench.get("genotype_checksum"),
            }
        reports.append(entry)
    if not reports:
        if os.environ.get("FASTGATK_REQUIRE_KOKKOS_BACKENDS") == "1":
            raise SystemExit("no Kokkos backend smoke artifact found")
        print(json.dumps({"status": "skip", "reason": "no backend smoke artifact"}))
        return 0
    if os.environ.get("FASTGATK_REQUIRE_KOKKOS_BACKENDS") == "1" and len(reports) < 2:
        raise SystemExit("required OpenMP+Serial backend matrix is incomplete")
    benchmark_reports = [entry["benchmark"] for entry in reports if "benchmark" in entry]
    if len(benchmark_reports) >= 2:
        reference = benchmark_reports[0]
        checksum_keys = ("pairhmm_checksum", "float_pairhmm_checksum",
                         "sw_checksum", "genotype_checksum")
        if any(any(item.get(key) != reference.get(key) for key in checksum_keys)
               for item in benchmark_reports[1:]):
            raise SystemExit("Kokkos backend benchmark checksums diverged")
    elif os.environ.get("FASTGATK_REQUIRE_KOKKOS_BACKENDS") == "1":
        raise SystemExit("required backend benchmark artifacts are incomplete")
    print(json.dumps({"status": "pass", "backends": reports,
                      "shared_api": True, "kokkos_only": True,
                      "benchmark_checksums_equal": len(benchmark_reports) >= 2}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
