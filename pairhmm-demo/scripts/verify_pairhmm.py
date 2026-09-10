#!/usr/bin/env python3
"""Compare native PairHMM backends with the real GATK Java oracle bit-for-bit."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
from pathlib import Path


def read_bits(path: Path) -> list[int]:
    return [int(line.strip(), 16) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def run_json(command: list[str], env: dict[str, str]) -> dict:
    output = subprocess.check_output(command, text=True, env=env)
    for line in reversed(output.splitlines()):
        if line.startswith("{") and line.endswith("}"):
            return json.loads(line)
    raise RuntimeError(f"no JSON result in output:\n{output}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", type=int, default=512)
    parser.add_argument("--read-len", type=int, default=150)
    parser.add_argument("--hap-len", type=int, default=160)
    parser.add_argument("--iterations", type=int, default=3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--backends", default="scalar,avx2,avx512")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    native = root / "pairhmm-demo/pairhmm-demo"
    java = root / "third_party/jdk17/bin/java"
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    classes = root / "pairhmm-demo/java-classes"
    tables = root / "pairhmm-demo/results/gatk-tables.hex"
    if not native.exists() or not java.exists() or not jar.exists() or not tables.exists():
        raise SystemExit("missing native/GATK package/tables; build and export the pinned test assets first")

    env = os.environ.copy()
    env["FAST_GATK_PAIRHMM_TABLES"] = str(tables)
    common = [f"--pairs={args.pairs}", f"--read-len={args.read_len}", f"--hap-len={args.hap_len}",
              f"--iterations={args.iterations}", f"--seed={args.seed}"]
    with tempfile.TemporaryDirectory(prefix="fast-gatk-pairhmm-") as temp:
        temp_path = Path(temp)
        oracle_path = temp_path / "gatk.hex"
        run_json([str(java), "-cp", f"{classes}:{jar}",
                  "org.broadinstitute.hellbender.utils.pairhmm.GatkPairHmmBenchmark",
                  "--mode=java", *common, f"--values-out={oracle_path}"], env)
        oracle = read_bits(oracle_path)
        report = []
        failed = False
        for backend in [x.strip() for x in args.backends.split(",") if x.strip()]:
            native_path = temp_path / f"native-{backend}.hex"
            try:
                result = run_json([str(native), "--mode=simd" if backend != "scalar" else "--mode=scalar",
                                   f"--backend={backend}" if backend != "scalar" else "--backend=auto",
                                   *common, f"--values-out={native_path}"], env)
            except subprocess.CalledProcessError as exc:
                report.append({"backend": backend, "status": "unavailable", "returncode": exc.returncode})
                continue
            actual = read_bits(native_path)
            differing = sum(a != b for a, b in zip(actual, oracle)) + abs(len(actual) - len(oracle))
            report.append({"backend": backend, "status": "pass" if differing == 0 else "fail",
                           "bit_different": differing, "pairs_per_second": result["pairs_per_second"]})
            failed |= differing != 0
        print(json.dumps({"oracle": "GATK pure-Java LoglessPairHMM", "report": report}, indent=2))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
