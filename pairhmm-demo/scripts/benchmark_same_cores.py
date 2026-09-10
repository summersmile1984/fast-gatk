#!/usr/bin/env python3
"""Run the GATK Java and the portable Kokkos PairHMM with the same CPU set.

The command is intentionally process-level pinned with taskset. It does not use a
machine-wide CPU count or allow Java/native to choose different affinity masks.
The native command is always ``pairhmm-kokkos`` from the selected scalar,
ZEN3/AVX2, or ZEN4/AVX-512 build; the legacy raw-intrinsics demo is never used
by this benchmark.
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import statistics
import subprocess
from pathlib import Path


def run(command: list[str], cpus: str, env: dict[str, str]) -> dict:
    output = subprocess.check_output(["taskset", "-c", cpus, *command], text=True, env=env)
    for line in reversed(output.splitlines()):
        if line.startswith("{") and line.endswith("}"):
            return json.loads(line)
    raise RuntimeError(f"no JSON result in output:\n{output}")


def _cpu_flags() -> set[str]:
    """Return normalized host CPU flags without executing an ISA probe.

    Reading the kernel's advertised flags is deliberately conservative for a
    benchmark launcher: an AVX-512 binary must never be started merely because
    it exists on disk.  On non-x86 hosts the set is empty and the portable
    scalar Kokkos build is selected.
    """
    if platform.machine().lower() not in {"x86_64", "amd64", "i386", "i686"}:
        return set()
    candidates = [Path("/proc/cpuinfo"), Path("/sys/devices/system/cpu/cpu0/capability")]
    for path in candidates:
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for line in text.splitlines():
            if line.lower().startswith(("flags", "features")) and ":" in line:
                return {flag.lower() for flag in line.split(":", 1)[1].split()}
    return set()


def choose_native_variant(root: Path, requested: str) -> str:
    """Choose a runnable Kokkos ISA build and reject unsafe explicit choices."""
    available = {
        variant: (root / f"pairhmm-demo/build-kokkos-{variant}/pairhmm-kokkos").is_file()
        for variant in ("scalar", "avx2", "avx512")
    }
    flags = _cpu_flags()
    supported = {
        "scalar": True,
        "avx2": "avx2" in flags,
        "avx512": "avx512f" in flags,
    }
    if requested != "auto":
        if not available.get(requested, False):
            raise SystemExit(f"missing Kokkos PairHMM binary for --native-backend={requested}")
        if not supported[requested]:
            raise SystemExit(
                f"--native-backend={requested} is not supported by this host; "
                f"advertised flags include: {','.join(sorted(flags))}"
            )
        return requested
    for variant in ("avx512", "avx2", "scalar"):
        if available[variant] and supported[variant]:
            return variant
    raise SystemExit(
        "no runnable Kokkos PairHMM binary; build scalar, avx2, or avx512 "
        "variant appropriate for this host"
    )


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", type=int, default=512)
    parser.add_argument("--read-len", type=int, default=150)
    parser.add_argument("--hap-len", type=int, default=160)
    parser.add_argument("--iterations", type=int, default=10)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--cores", default="1,2,4,8,16")
    parser.add_argument("--native-backend", choices=("auto", "scalar", "avx2", "avx512"), default="auto")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    java = root / "third_party/jdk17/bin/java"
    jar = root / "third_party/gatk-package/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
    tables = root / "pairhmm-demo/results/gatk-tables.hex"
    java_classpath = root / "pairhmm-demo/java-classes"
    native_variant = choose_native_variant(root, args.native_backend)
    native = root / f"pairhmm-demo/build-kokkos-{native_variant}/pairhmm-kokkos"
    if not native.is_file():
        raise SystemExit(f"missing Kokkos PairHMM binary: {native}")
    base = [f"--pairs={args.pairs}", f"--read-len={args.read_len}", f"--hap-len={args.hap_len}",
            f"--iterations={args.iterations}", f"--seed={args.seed}"]
    env = os.environ.copy()
    env["FAST_GATK_PAIRHMM_TABLES"] = str(tables)
    # Keep Kokkos/OpenMP placement inside the exact taskset mask used by both
    # implementations.  This also avoids an implementation-specific default
    # that could silently oversubscribe a SLURM allocation.
    env.setdefault("OMP_PROC_BIND", "true")
    env.setdefault("OMP_PLACES", "threads")
    all_flags = _cpu_flags()
    host_flags = sorted(all_flags.intersection({
        "avx2", "avx512f", "avx512bw", "avx512vl", "avx512vnni", "fma", "sse4_2"
    }))
    rows = []
    for cores in [int(x) for x in args.cores.split(",")]:
        cpus = f"0-{cores - 1}"
        java_runs = []
        native_runs = []
        java_command = [str(java), "-cp", f"{java_classpath}:{jar}",
                        "org.broadinstitute.hellbender.utils.pairhmm.GatkPairHmmBenchmark",
                        "--mode=java", f"--threads={cores}", *base]
        native_command = [str(native), "--workload=independent", f"--threads={cores}", *base]
        for repeat in range(args.repeats):
            # Alternate order to reduce systematic boost/thermal bias from always
            # running one implementation immediately after the other.
            if repeat % 2 == 0:
                java_runs.append(run(java_command, cpus, env))
                native_runs.append(run(native_command, cpus, env))
            else:
                native_runs.append(run(native_command, cpus, env))
                java_runs.append(run(java_command, cpus, env))
        j = statistics.median(x["pairs_per_second"] for x in java_runs)
        n = statistics.median(x["pairs_per_second"] for x in native_runs)
        native_first = native_runs[0]
        rows.append({"cores": cores, "gatk_java_pairs_per_second": j,
                     "native_backend": native_variant,
                     "native_execution_space": native_first.get("execution_space", ""),
                     "native_simd_width": native_first.get("simd_width", 1),
                     "native_prepare_seconds": native_first.get("prepare_seconds", 0.0),
                     "native_pairs_per_second": n,
                     f"native_{native_variant}_pairs_per_second": n,
                     "speedup": n / j,
                     "checksum": native_runs[0]["checksum"]})
    result = {"input": {"pairs": args.pairs, "read_len": args.read_len, "hap_len": args.hap_len,
                         "iterations": args.iterations, "seed": args.seed,
                         "warmup_iterations": 1,
                         "host_arch": platform.machine(),
                         "host_cpu_flags": host_flags,
                         "requested_native_backend": args.native_backend,
                         "native_backend": native_variant,
                         "isa_selection_policy": "host-flags-and-binary-availability",
                         "native_api": "Kokkos::Experimental::simd<double>"}, "rows": rows}
    print(json.dumps(result, indent=2))
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
