#!/usr/bin/env python3
"""Benchmark one-source Kokkos PairHMM variants against GATK Java and GKL.

Independent pairs are compared with GATK's pure-Java LoglessPairHMM.  The
matrix8 workload is compared with GKL's real amount of work: GKL's Java API
reports only the diagonal record count even though each block computes n*n
read/haplotype combinations, so its reported throughput is corrected here.
"""
from __future__ import annotations

import argparse
import json
import os
import statistics
import subprocess
from pathlib import Path


def run(command: list[str], cpus: str, env: dict[str, str]) -> dict:
    output = subprocess.check_output(
        ["taskset", "-c", cpus, *command], text=True, env=env,
        stderr=subprocess.DEVNULL,
    )
    for line in reversed(output.splitlines()):
        if line.startswith("{") and line.endswith("}"):
            return json.loads(line)
    raise RuntimeError(f"no JSON result in output:\n{output}")


def median(results: list[dict], key: str) -> float:
    return statistics.median(float(item[key]) for item in results)


def matrix_pair_count(records: int, block: int = 8) -> int:
    return sum(min(block, records - base) ** 2 for base in range(0, records, block))


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", type=int, default=512,
                        help="number of read/haplotype records")
    parser.add_argument("--read-len", type=int, default=150)
    parser.add_argument("--hap-len", type=int, default=160)
    parser.add_argument("--iterations", type=int, default=10)
    parser.add_argument("--independent-iterations", type=int)
    parser.add_argument("--matrix-iterations", type=int)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--cores", default="1,2,4,8,16")
    parser.add_argument("--variants", default="scalar,avx2,avx512")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    variants = [item.strip() for item in args.variants.split(",") if item.strip()]
    cores_list = [int(item) for item in args.cores.split(",")]
    java = root / "third_party/jdk17/bin/java"
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    classes = root / "pairhmm-demo/java-classes"
    tables = root / "pairhmm-demo/results/gatk-tables.hex"
    binaries = {
        variant: root / f"pairhmm-demo/build-kokkos-{variant}/pairhmm-kokkos"
        for variant in variants
    }
    missing = [path for path in [java, jar, classes, tables, *binaries.values()] if not path.exists()]
    if missing:
        raise SystemExit("missing benchmark assets: " + ", ".join(str(path) for path in missing))

    env = os.environ.copy()
    env["FAST_GATK_PAIRHMM_TABLES"] = str(tables)
    env["OMP_PROC_BIND"] = "true"
    env["OMP_PLACES"] = "threads"
    independent_iterations = args.independent_iterations or args.iterations
    matrix_iterations = args.matrix_iterations or args.iterations
    base = [f"--pairs={args.pairs}", f"--read-len={args.read_len}",
            f"--hap-len={args.hap_len}", f"--seed={args.seed}"]
    independent_base = [*base, f"--iterations={independent_iterations}"]
    matrix_base = [*base, f"--iterations={matrix_iterations}"]
    java_base = [str(java), "-cp", f"{classes}:{jar}",
                 "org.broadinstitute.hellbender.utils.pairhmm.GatkPairHmmBenchmark"]
    matrix_pairs = matrix_pair_count(args.pairs)
    gkl_work_multiplier = matrix_pairs / args.pairs
    rows: list[dict] = []

    for cores in cores_list:
        cpus = f"0-{cores - 1}"
        for variant in variants:
            java_command = [*java_base, "--mode=java", f"--threads={cores}", *independent_base]
            independent_command = [str(binaries[variant]), "--workload=independent",
                                   f"--threads={cores}", *independent_base]
            independent_runs: list[dict] = []
            java_runs: list[dict] = []
            for repeat in range(args.repeats):
                if repeat % 2 == 0:
                    java_runs.append(run(java_command, cpus, env))
                    independent_runs.append(run(independent_command, cpus, env))
                else:
                    independent_runs.append(run(independent_command, cpus, env))
                    java_runs.append(run(java_command, cpus, env))
            java_rate = median(java_runs, "pairs_per_second")
            kokkos_rate = median(independent_runs, "pairs_per_second")
            prepare = median(independent_runs, "prepare_seconds")
            kernel_seconds = median(independent_runs, "seconds")
            rows.append({
                "workload": "independent",
                "cores": cores,
                "variant": variant,
                "simd_width": independent_runs[0]["simd_width"],
                "gatk_java_pairs_per_second": java_rate,
                "kokkos_kernel_pairs_per_second": kokkos_rate,
                "kokkos_amortized_pairs_per_second":
                    args.pairs * independent_iterations / (prepare + kernel_seconds),
                "kokkos_over_java": kokkos_rate / java_rate,
                "prepare_seconds": prepare,
                "checksum": independent_runs[0]["checksum"],
            })

            gkl_command = [*java_base, "--mode=gatk-avx-omp", f"--threads={cores}", *matrix_base]
            matrix_command = [str(binaries[variant]), "--workload=matrix8",
                              f"--threads={cores}", *matrix_base]
            matrix_runs: list[dict] = []
            gkl_runs: list[dict] = []
            for repeat in range(args.repeats):
                if repeat % 2 == 0:
                    gkl_runs.append(run(gkl_command, cpus, env))
                    matrix_runs.append(run(matrix_command, cpus, env))
                else:
                    matrix_runs.append(run(matrix_command, cpus, env))
                    gkl_runs.append(run(gkl_command, cpus, env))
            gkl_reported = median(gkl_runs, "pairs_per_second")
            gkl_actual = gkl_reported * gkl_work_multiplier
            kokkos_matrix_rate = median(matrix_runs, "pairs_per_second")
            matrix_prepare = median(matrix_runs, "prepare_seconds")
            matrix_kernel_seconds = median(matrix_runs, "seconds")
            rows.append({
                "workload": "matrix8",
                "records": args.pairs,
                "computed_pairs": matrix_pairs,
                "cores": cores,
                "variant": variant,
                "simd_width": matrix_runs[0]["simd_width"],
                "gkl_reported_diagonal_pairs_per_second": gkl_reported,
                "gkl_actual_matrix_pairs_per_second": gkl_actual,
                "gkl_work_multiplier": gkl_work_multiplier,
                "kokkos_kernel_pairs_per_second": kokkos_matrix_rate,
                "kokkos_amortized_pairs_per_second":
                    matrix_pairs * matrix_iterations / (matrix_prepare + matrix_kernel_seconds),
                "kokkos_over_gkl_actual": kokkos_matrix_rate / gkl_actual,
                "prepare_seconds": matrix_prepare,
                "checksum": matrix_runs[0]["checksum"],
            })

    result = {
        "input": {
            "records": args.pairs,
            "read_len": args.read_len,
            "hap_len": args.hap_len,
            "independent_iterations": independent_iterations,
            "matrix_iterations": matrix_iterations,
            "repeats": args.repeats,
            "seed": args.seed,
            "warmup_iterations": 1,
            "affinity": "physical logical CPUs 0..cores-1",
            "kokkos_version": "5.2.0",
        },
        "methodology": {
            "statistic": "median",
            "order": "alternating baseline/Kokkos",
            "kokkos_kernel_rate_excludes": "input packing and host/device deep copies",
            "kokkos_amortized_rate_includes": "one preparation per benchmark invocation",
            "gkl_rate_correction": "reported diagonal record rate multiplied by actual n*n/block work",
        },
        "rows": rows,
    }
    rendered = json.dumps(result, indent=2) + "\n"
    print(rendered, end="")
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
