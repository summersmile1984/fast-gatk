#!/usr/bin/env python3
"""Reproducible ReblockGVCF native-versus-Java file-boundary benchmark.

The benchmark deliberately builds one common ReblockGVCF argument vector and
injects only the implementation-specific executable prefix.  This prevents a
native run from silently using a different interval, output-index policy, or
input fixture than the Java GATK baseline.  Warmup, wall time, maximum RSS,
and filesystem I/O counters are collected with GNU ``time`` for both paths.

The default fixture is the pinned GATK ReblockGVCF input.  Its small size makes
the command contract cheap to check, while ``--repetitions`` can be increased
for a real file-boundary measurement.  Full VCF byte identity is intentionally
not a benchmark gate; use the dedicated GATK oracle for that.  A throughput
claim is enabled only when both implementations ran successfully and emitted
the same number of records from the same input scope.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shlex
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any


TIME_BIN = Path("/usr/bin/time")
SCHEMA_VERSION = 1


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def count_records(path: Path) -> int:
    count = 0
    with path.open("rt", encoding="utf-8") as stream:
        for line in stream:
            if line and not line.startswith("#"):
                count += 1
    return count


def remove_output(path: Path) -> None:
    """Remove one output and all conventional variant-index sidecars."""
    for candidate in (
        path,
        Path(str(path) + ".tbi"),
        Path(str(path) + ".csi"),
        Path(str(path) + ".idx"),
        Path(str(path) + ".md5"),
    ):
        try:
            candidate.unlink()
        except FileNotFoundError:
            pass


def parse_time_file(path: Path) -> dict[str, int | float]:
    fields = path.read_text(encoding="ascii").strip().split("\t")
    if len(fields) != 4:
        raise RuntimeError(f"unexpected GNU time output: {path.read_text()!r}")
    return {
        # %e is elapsed wall-clock seconds for the complete process.
        "wall_seconds": float(fields[0]),
        # %M is the maximum resident set size in KiB.
        "max_rss_kb": int(fields[1]),
        # %I/%O are filesystem input/output operations, not bytes.  Keeping
        # their units explicit makes Java/native measurements comparable.
        "filesystem_inputs": int(fields[2]),
        "filesystem_outputs": int(fields[3]),
    }


def run_once(command: list[str], output: Path, timing: Path,
             environment: dict[str, str]) -> dict[str, Any]:
    remove_output(output)
    timed = [
        str(TIME_BIN), "-f", "%e\\t%M\\t%I\\t%O", "-o", str(timing),
        *command,
    ]
    started = time.perf_counter()
    completed = subprocess.run(
        timed, text=True, capture_output=True, check=False, env=environment,
    )
    python_wall = time.perf_counter() - started
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError(
            f"command failed ({completed.returncode}): "
            f"{shlex.join(command)}\n{detail}"
        )
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"command did not create a non-empty output: {output}")
    measured = parse_time_file(timing)
    # GNU time's %e is only centisecond precision, which rounds very small
    # native runs to 0.00 and would make a speedup ratio undefined.  Use the
    # same high-resolution parent-clock measurement for both implementations
    # as the canonical wall field, while retaining %e for auditability.
    measured["gnu_time_wall_seconds"] = measured["wall_seconds"]
    measured["wall_seconds"] = python_wall
    measured["output_bytes"] = output.stat().st_size
    measured["output_records"] = count_records(output)
    return measured


def summarize_runs(runs: list[dict[str, Any]]) -> dict[str, Any]:
    wall = sorted(float(run["wall_seconds"]) for run in runs)
    middle = wall[len(wall) // 2]
    p95 = wall[min(len(wall) - 1, int(0.95 * len(wall)))]
    return {
        "warmup": 1,
        "repetitions": len(runs),
        "runs": runs,
        "p50_seconds": middle,
        "p95_seconds": p95,
        "max_rss_kb": max(int(run["max_rss_kb"]) for run in runs),
        "filesystem_inputs": sum(int(run["filesystem_inputs"]) for run in runs),
        "filesystem_outputs": sum(int(run["filesystem_outputs"]) for run in runs),
        "output_bytes": int(runs[-1]["output_bytes"]),
        "output_records": int(runs[-1]["output_records"]),
    }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Compare ReblockGVCF native and GATK Java with aligned resources"
    )
    result.add_argument("--source", type=Path,
                        help="single-sample GVCF; without it, generate one from the pinned HC fixture")
    result.add_argument("--reference", type=Path,
                        help="reference used by both ReblockGVCF commands")
    result.add_argument("--interval", default="17:69000-70000")
    result.add_argument("--repetitions", type=int, default=3)
    result.add_argument("--java-xmx", default="1g")
    result.add_argument("--minimum-records", type=int, default=1024,
                        help="minimum output records before a speedup claim is enabled (default: 1024)")
    result.add_argument("--native-binary", type=Path)
    result.add_argument("--java", type=Path)
    result.add_argument("--gatk-jar", type=Path)
    result.add_argument("--omp-threads", type=int, default=1)
    result.add_argument("--dry-run", action="store_true",
                        help="validate aligned commands without launching either implementation")
    return result


def main() -> int:
    args = parser().parse_args()
    if not 1 <= args.repetitions <= 20:
        raise SystemExit("--repetitions must be between 1 and 20")
    if args.omp_threads < 1:
        raise SystemExit("--omp-threads must be positive")
    if args.minimum_records < 1:
        raise SystemExit("--minimum-records must be positive")
    if not TIME_BIN.is_file() and not args.dry_run:
        raise SystemExit(f"GNU time is required: {TIME_BIN}")

    root = Path(__file__).resolve().parents[2]
    source = args.source
    native = args.native_binary or Path(os.environ.get(
        "FASTGATK_REBLOCK_BINARY", str(root / "fastgatk-native/build/fastgatk-reblock-gvcf")
    ))
    java = args.java or Path(os.environ.get(
        "JAVA", str(root / "third_party/jdk17/bin/java")
    ))
    gatk = args.gatk_jar or Path(os.environ.get(
        "FASTGATK_GATK_JAR",
        str(root / "third_party/gatk-package/gatk-4.6.2.0/"
            "gatk-package-4.6.2.0-local.jar"),
    ))
    reference = args.reference or root / (
        "gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
    )

    # One common tool argument vector is used by both implementations.  The
    # output placeholder is substituted only after the alignment assertions.
    common_args = [
        "-R", "{REFERENCE}", "-V", "{INPUT}", "-L", args.interval,
        "-O", "{OUTPUT}",
        "--create-output-variant-index", "false",
    ]
    native_template = [str(native), *common_args]
    java_template = [
        str(java), f"-Xmx{args.java_xmx}", "-jar", str(gatk), "ReblockGVCF",
        *common_args,
    ]
    native_algorithm_args = native_template[1:]
    java_algorithm_args = java_template[5:]
    if native_algorithm_args != java_algorithm_args:
        raise AssertionError("native and Java ReblockGVCF arguments diverged")
    if common_args[5] != args.interval:
        raise AssertionError("interval was not propagated to common arguments")
    if not args.dry_run and source is not None and not source.is_file():
        raise SystemExit(f"missing benchmark source: {source}")
    if not args.dry_run and not reference.is_file():
        raise SystemExit(f"missing benchmark reference: {reference}")
    if not args.dry_run and not native.is_file():
        raise SystemExit(f"missing native ReblockGVCF binary: {native}")
    if not args.dry_run and (not java.is_file() or not gatk.is_file()):
        raise SystemExit("pinned Java/GATK baseline is required; use --dry-run to inspect commands")

    source_meta: dict[str, Any] = {
        "path": str(source) if source else "<generated HaplotypeCaller fixture>",
        "bytes": source.stat().st_size if source and source.is_file() else None,
        "sha256": sha256(source) if source and source.is_file() else None,
    }
    command_scope = {
        "input": "{INPUT}",
        "interval": args.interval,
        "output_index": False,
        "reference": str(reference),
        "algorithm_args": common_args,
        "native_command_template": native_template,
        "java_command_template": java_template,
        "native_command": shlex.join(native_template),
        "java_command": shlex.join(java_template),
    }
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(args.omp_threads)
    environment.setdefault("OMP_PROC_BIND", "true")
    environment.setdefault("OMP_PLACES", "threads")

    if args.dry_run:
        print(json.dumps({
            "schema_version": SCHEMA_VERSION,
            "suite": "fastgatk-reblock-gvcf-java-baseline",
            "tool": "ReblockGVCF",
            "status": "pass",
            "run_mode": "dry-run",
            "parameter_alignment": {
                "same_input_placeholder": True,
                "same_interval": True,
                "same_output_index_policy": True,
                "same_algorithm_arguments": True,
                "same_omp_threads_environment": True,
            },
            "command_scope": command_scope,
            "resource_schema": [
                "wall_seconds", "max_rss_kb", "filesystem_inputs",
                "filesystem_outputs", "output_bytes", "output_records",
            ],
            "speedup_claim_allowed": False,
            "speedup_claim_blocked_reason": "dry-run does not measure a workload",
        }, sort_keys=True))
        return 0

    with tempfile.TemporaryDirectory(prefix="fastgatk-reblock-java-benchmark-") as directory:
        work = Path(directory)
        setup_command: list[str] | None = None
        if source is None:
            bam = root / "gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
            # A compressed, indexed source is required because GATK's
            # interval traversal rejects an unindexed feature input.  The
            # same source is then consumed by native and Java ReblockGVCF.
            source = work / "hc.g.vcf.gz"
            setup_command = [
                str(java), f"-Xmx{args.java_xmx}", "-jar", str(gatk),
                "HaplotypeCaller", "-R", str(reference), "-I", str(bam),
                "-L", args.interval, "-O", str(source), "-ERC", "GVCF",
                "--create-output-variant-index", "true",
            ]
            setup = subprocess.run(setup_command, text=True, capture_output=True,
                                   check=False, env=environment)
            if setup.returncode != 0:
                detail = setup.stderr.strip() or setup.stdout.strip()
                raise RuntimeError(f"HaplotypeCaller fixture generation failed:\n{detail}")
        source_meta = {
            "path": str(source),
            "bytes": source.stat().st_size,
            "sha256": sha256(source),
        }
        native_runs: list[dict[str, Any]] = []
        java_runs: list[dict[str, Any]] = []
        for implementation, template, runs in (
            ("native-kokkos", native_template, native_runs),
            ("gatk-java", java_template, java_runs),
        ):
            output = work / f"{implementation}.vcf"
            command = [
                item.replace("{INPUT}", str(source)).replace("{OUTPUT}", str(output))
                .replace("{REFERENCE}", str(reference))
                for item in template
            ]
            warmup_timing = work / f"{implementation}-warmup.time"
            run_once(command, output, warmup_timing, environment)
            for repeat in range(args.repetitions):
                timing = work / f"{implementation}-{repeat}.time"
                runs.append(run_once(command, output, timing, environment))

        native_summary = summarize_runs(native_runs)
        java_summary = summarize_runs(java_runs)
        counts_match = native_summary["output_records"] == java_summary["output_records"]
        record_count = min(int(native_summary["output_records"]),
                           int(java_summary["output_records"]))
        workload_sufficient = record_count >= args.minimum_records
        claim_allowed = bool(counts_match and workload_sufficient)
        reasons = []
        if not counts_match:
            reasons.append("native and Java emitted different record counts")
        if not workload_sufficient:
            reasons.append(
                f"workload has {record_count} records; minimum is {args.minimum_records}"
            )
        reason = "; ".join(reasons) or None
        native_p50 = float(native_summary["p50_seconds"])
        java_p50 = float(java_summary["p50_seconds"])
        artifact = {
            "schema_version": SCHEMA_VERSION,
            "suite": "fastgatk-reblock-gvcf-java-baseline",
            "tool": "ReblockGVCF",
            "status": "pass",
            "run_mode": "measured",
            "gatk_version": "4.6.2.0",
            "backend": "Kokkos",
            "input": source_meta,
            "fixture_setup": {
                "generated_by_gatk_haplotype_caller": setup_command is not None,
                "command": setup_command,
                "timed": False,
            },
            "command_scope": command_scope,
            "parameter_alignment": {
                "same_input": True,
                "same_input_sha256": True,
                "same_interval": True,
                "same_output_index_policy": True,
                "same_algorithm_arguments": True,
                "same_omp_threads_environment": True,
            },
            "environment": {
                "platform": platform.platform(),
                "machine": platform.machine(),
                "omp_threads": args.omp_threads,
                "omp_proc_bind": environment.get("OMP_PROC_BIND"),
                "omp_places": environment.get("OMP_PLACES"),
                "time_binary": str(TIME_BIN),
                "time_io_units": "filesystem_operations",
            },
            "implementations": {
                "native-kokkos": native_summary,
                "gatk-java": java_summary,
            },
            "output_record_counts_match": counts_match,
            "output_record_count_min": record_count,
            "minimum_records_for_speedup_claim": args.minimum_records,
            "workload_sufficient_for_speedup_claim": workload_sufficient,
            "speedup_claim_allowed": claim_allowed,
            "speedup_claim_blocked_reason": reason,
            "p50_speedup_java_over_native": java_p50 / native_p50 if native_p50 else None,
            "note": (
                "wall/RSS/filesystem-I/O fields are collected with the same GNU "
                "time format for both processes; full record/bit identity is "
                "covered by verify_reblock_gatk_oracle.py"
            ),
        }
        print(json.dumps(artifact, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
