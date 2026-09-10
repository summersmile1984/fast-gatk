#!/usr/bin/env python3
"""Verify all Kokkos architecture builds against a pinned strict GATK oracle."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
from pathlib import Path


def read_bits(path: Path) -> list[int]:
    return [int(line, 16) for line in path.read_text(encoding="utf-8").splitlines() if line]


def run(command: list[str], env: dict[str, str]) -> dict:
    output = subprocess.check_output(command, text=True, env=env, stderr=subprocess.DEVNULL)
    for line in reversed(output.splitlines()):
        if line.startswith("{") and line.endswith("}"):
            return json.loads(line)
    raise RuntimeError(f"no JSON result in output:\n{output}")


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", type=int, default=512)
    parser.add_argument("--read-len", type=int, default=150)
    parser.add_argument("--hap-len", type=int, default=160)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--workload", choices=("independent", "matrix8"), default="matrix8")
    parser.add_argument("--variants", default="scalar,avx2,avx512")
    args = parser.parse_args()

    java = root / "third_party/jdk17/bin/java"
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    classes = root / "pairhmm-demo/java-classes"
    tables = root / "pairhmm-demo/results/gatk-tables.hex"
    variants = [item.strip() for item in args.variants.split(",") if item.strip()]
    env = os.environ.copy()
    env["FAST_GATK_PAIRHMM_TABLES"] = str(tables)
    # Also exercise the production default.  The native build should find the
    # pinned Java table through FASTGATK_DEFAULT_PAIRHMM_TABLES_PATH without a
    # caller having to export a process-specific environment variable.
    default_env = env.copy()
    default_env.pop("FAST_GATK_PAIRHMM_TABLES", None)
    env["OMP_PROC_BIND"] = "true"
    env["OMP_PLACES"] = "threads"
    default_env["OMP_PROC_BIND"] = "true"
    default_env["OMP_PLACES"] = "threads"
    common = [f"--pairs={args.records}", f"--read-len={args.read_len}",
              f"--hap-len={args.hap_len}", "--threads=1", "--iterations=1",
              f"--seed={args.seed}", f"--workload={args.workload}"]
    java_base = [str(java), "-XX:+UnlockDiagnosticVMOptions", "-XX:DisableIntrinsic=_dlog10",
                 "-cp", f"{classes}:{jar}",
                 "org.broadinstitute.hellbender.utils.pairhmm.GatkPairHmmBenchmark",
                 "--mode=java", *common]

    report: list[dict] = []
    failed = False
    with tempfile.TemporaryDirectory(prefix="fast-gatk-kokkos-") as temp:
        temp_path = Path(temp)
        oracle_path = temp_path / "gatk-strict.hex"
        oracle_result = run([*java_base, f"--values-out={oracle_path}"], env)
        oracle = read_bits(oracle_path)
        for variant in variants:
            binary = root / f"pairhmm-demo/build-kokkos-{variant}/pairhmm-kokkos"
            output_path = temp_path / f"kokkos-{variant}.hex"
            actual_result = run(["taskset", "-c", "0", str(binary), *common,
                                 f"--values-out={output_path}"], env)
            actual = read_bits(output_path)
            differing = sum(left != right for left, right in zip(actual, oracle))
            differing += abs(len(actual) - len(oracle))
            default_output_path = temp_path / f"kokkos-{variant}-default.hex"
            default_result = run(["taskset", "-c", "0", str(binary), *common,
                                  f"--values-out={default_output_path}"], default_env)
            default_actual = read_bits(default_output_path)
            default_differing = sum(left != right for left, right in zip(default_actual, oracle))
            default_differing += abs(len(default_actual) - len(oracle))
            failed = failed or differing != 0 or default_differing != 0
            report.append({
                "variant": variant,
                "simd_width": actual_result["simd_width"],
                "values": len(actual),
                "bit_different": differing,
                "default_table_bit_different": default_differing,
                "default_table": "pinned-java-generated",
                "status": "pass" if differing == 0 and default_differing == 0 else "fail",
            })
    print(json.dumps({
        "oracle": "GATK 4.6.2.0 LoglessPairHMM + Java 17 StrictMath fdlibm",
        "java_flag": "-XX:DisableIntrinsic=_dlog10",
        "workload": args.workload,
        "oracle_values": oracle_result.get("computed_pairs", oracle_result.get("pairs")),
        "report": report,
    }, indent=2))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
