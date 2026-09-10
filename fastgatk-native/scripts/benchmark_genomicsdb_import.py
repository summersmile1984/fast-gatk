#!/usr/bin/env python3
"""File-boundary benchmark for the resource-aware GenomicsDBImport adapter.

The sparse GenomicsDB/TileDB backend remains the configured GATK backend, so
this benchmark measures the native Host adapter contract: argument staging,
resource adaptation, workspace validation and deterministic input-index
publication.  It never reports the fake backend as GenomicsDB throughput.
"""

from __future__ import annotations

import json
import gzip
import os
import statistics
import subprocess
import tempfile
import time
from pathlib import Path


def run_once(binary: Path, backend: Path, input_vcf: Path, workspace: Path, manifest: Path,
             batch: int, threads: int, scratch: Path) -> tuple[float, dict, dict]:
    start = time.perf_counter()
    result = subprocess.run(
        [str(binary), "-V", str(input_vcf),
         "--genomicsdb-workspace-path", str(workspace),
         "--batch-size", str(batch), "--reader-threads", str(threads),
         "--tmp-dir", str(scratch), "--output-manifest", str(manifest)],
        text=True, capture_output=True, check=False,
        env={**os.environ, "FASTGATK_GATK_BINARY": str(backend),
             "SLURM_CPUS_PER_TASK": "4", "SLURM_MEM_PER_NODE": "64",
             "SLURM_TMPDIR": str(scratch), "SLURM_JOB_ID": "benchmark"},
    )
    elapsed = time.perf_counter() - start
    if result.returncode != 0:
        raise RuntimeError(result.stderr or result.stdout)
    summary = json.loads(result.stdout.strip().splitlines()[-1])
    metadata = json.loads(manifest.read_text(encoding="utf-8"))
    if summary["status"] != "adapter" or not metadata["outputs"][0]["complete"]:
        raise AssertionError("adapter did not publish a complete workspace")
    return elapsed, summary, metadata


def run_native_once(binary: Path, input_vcf: Path, workspace: Path, manifest: Path,
                    scratch: Path) -> tuple[float, dict, dict]:
    start = time.perf_counter()
    result = subprocess.run(
        [str(binary), "-V", str(input_vcf),
         "--genomicsdb-workspace-path", str(workspace),
         "--fastgatk-native-workspace", "--batch-size", "64",
         "--reader-threads", "2", "--tmp-dir", str(scratch),
         "--output-manifest", str(manifest)],
        text=True, capture_output=True, check=False,
        env={**os.environ, "SLURM_CPUS_PER_TASK": "2", "SLURM_MEM_PER_NODE": "64",
             "SLURM_TMPDIR": str(scratch), "SLURM_JOB_ID": "benchmark-native"},
    )
    elapsed = time.perf_counter() - start
    if result.returncode != 0:
        raise RuntimeError(result.stderr or result.stdout)
    summary = json.loads(result.stdout.strip().splitlines()[-1])
    metadata = json.loads(manifest.read_text(encoding="utf-8"))
    if summary["execution_mode"] != "native-sparse-index" or not metadata["outputs"][0]["complete"]:
        raise AssertionError("native adapter did not publish a complete sparse workspace")
    return elapsed, summary, metadata


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    build = Path(os.environ.get("FASTGATK_NATIVE_BUILD", root / "fastgatk-native/build"))
    binary = build / "fastgatk-genomicsdb-import"
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise SystemExit(f"missing executable: {binary}")

    with tempfile.TemporaryDirectory(prefix="fastgatk-genomicsdb-benchmark-") as temp:
        work = Path(temp)
        input_vcf = work / "sample.g.vcf.gz"
        # The external backend is intentionally tiny; the input is only used
        # to exercise local-path validation and deterministic indexing.
        input_vcf.write_bytes(b"placeholder-gvcf\n" * 4096)
        native_input = work / "native-sample.g.vcf.gz"
        with gzip.open(native_input, "wt", encoding="ascii") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            handle.write("##contig=<ID=chr1,length=100>\n")
            handle.write("##ALT=<ID=NON_REF,Description=Represents any possible alternate allele>\n")
            handle.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n")
            handle.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Read depth>\n")
            handle.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=Likelihoods>\n")
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
            handle.write("chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=10\tGT:DP:PL\t0/0:10:0,99,990\n")
        scratch = work / "scratch"
        scratch.mkdir()
        backend = work / "fake-gatk"
        backend.write_text(
            "#!/usr/bin/env python3\n"
            "import pathlib, sys\n"
            "args = sys.argv[1:]\n"
            "workspace = pathlib.Path(args[args.index('--genomicsdb-workspace-path') + 1])\n"
            "workspace.mkdir(parents=True, exist_ok=True)\n"
            "(workspace / 'callset.json').write_text('{\\\"records\\\":1}\\n')\n"
            "sys.exit(0)\n",
            encoding="utf-8",
        )
        backend.chmod(0o755)

        rows: list[dict] = []
        for batch, threads in ((1, 1), (999, 999)):
            times: list[float] = []
            summary: dict | None = None
            metadata: dict | None = None
            for repeat in range(4):
                workspace = work / f"workspace-{batch}-{threads}-{repeat}"
                manifest = work / f"manifest-{batch}-{threads}-{repeat}.json"
                elapsed, summary, metadata = run_once(
                    binary, backend, input_vcf, workspace, manifest,
                    batch, threads, scratch)
                if repeat > 0:
                    times.append(elapsed)
            assert summary is not None and metadata is not None
            p50 = statistics.median(times)
            p95 = max(times)
            rows.append({
                "requested_batch_size": batch,
                "requested_reader_threads": threads,
                "effective_batch_size": summary["effective_batch_size"],
                "effective_reader_threads": summary["effective_reader_threads"],
                "p50_seconds": p50,
                "p95_seconds": p95,
                "workspace_bytes": sum(
                    path.stat().st_size for path in workspace.rglob("*") if path.is_file()),
                "input_count": metadata["telemetry"]["workspace_input_count"],
                "execution_mode": metadata["execution_mode"],
            })

        native_times: list[float] = []
        native_summary: dict | None = None
        native_metadata: dict | None = None
        for repeat in range(4):
            workspace = work / f"native-workspace-{repeat}"
            manifest = work / f"native-manifest-{repeat}.json"
            elapsed, native_summary, native_metadata = run_native_once(
                binary, native_input, workspace, manifest, scratch)
            if repeat > 0:
                native_times.append(elapsed)
        assert native_summary is not None and native_metadata is not None
        rows.append({
            "requested_batch_size": 64,
            "requested_reader_threads": 2,
            "effective_batch_size": native_summary["effective_batch_size"],
            "effective_reader_threads": native_summary["effective_reader_threads"],
            "p50_seconds": statistics.median(native_times),
            "p95_seconds": max(native_times),
            "workspace_bytes": sum(
                path.stat().st_size for path in (work / "native-workspace-3").rglob("*")
                if path.is_file()),
            "input_count": native_metadata["telemetry"]["workspace_input_count"],
            "record_index": True,
            "execution_mode": native_metadata["execution_mode"],
        })

    print(json.dumps({
        "schema_version": 1,
        "status": "pass",
        "tool": "GenomicsDBImport",
        "benchmark": "host-adapter-file-boundary",
        "backend": "fake-external-contract-and-native-sparse-index",
        "warmup_runs": 1,
        "measured_repeats": 3,
        "workloads": rows,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
