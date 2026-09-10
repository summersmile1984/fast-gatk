#!/usr/bin/env python3
"""Small end-to-end benchmark for native IndexFeatureFile."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BINARY = Path(os.environ.get("FASTGATK_INDEX_FEATURE_FILE_BINARY", ROOT / "fastgatk-native/build/fastgatk-index-feature-file"))
INPUT = ROOT / "gatk-source/src/test/resources/org/broadinstitute/hellbender/tools/IndexFeatureFile/4featuresHG38Header.unindexed.vcf.gz"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--iterations", type=int, default=3)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--input", type=Path, default=INPUT,
                        help="feature file to index (defaults to the bundled BGZF VCF)")
    args = parser.parse_args()
    input_path = args.input if args.input.is_absolute() else ROOT / args.input
    if not BINARY.is_file() or not input_path.is_file():
        print(json.dumps({"status": "skip", "suite": "index-feature-file-benchmark"}))
        return 0
    timings: list[float] = []
    index_kinds: list[str] = []
    with tempfile.TemporaryDirectory(prefix="fastgatk-index-feature-file-bench-") as directory:
        plain_tribble_input = input_path.name.endswith((".g.vcf", ".vcf", ".bed"))
        output = Path(directory) / ("sample.idx" if plain_tribble_input else "sample.tbi")
        manifest = Path(directory) / "sample.manifest.json"
        for _ in range(max(1, args.iterations)):
            output.unlink(missing_ok=True)
            manifest.unlink(missing_ok=True)
            start = time.perf_counter()
            result = subprocess.run(
                [str(BINARY), "-I", str(input_path), "-O", str(output), "--threads", str(args.threads),
                 "--output-manifest", str(manifest)],
                text=True, capture_output=True,
            )
            if result.returncode != 0:
                raise RuntimeError(result.stderr)
            timings.append(time.perf_counter() - start)
            index_kinds.append(json.loads(manifest.read_text())["index_type"])
    timings.sort()
    p50 = timings[len(timings) // 2]
    print(json.dumps({
        "status": "pass", "suite": "index-feature-file-benchmark", "input": str(input_path),
        "input_bytes": input_path.stat().st_size,
        "threads": args.threads, "iterations": len(timings), "p50_seconds": p50,
        "p95_seconds": timings[min(len(timings) - 1, int(len(timings) * 0.95))],
        "bytes_per_second": input_path.stat().st_size / p50 if p50 else None,
        "compressed_bytes_per_second": input_path.stat().st_size / p50 if p50 else None,
        "index_kind": index_kinds[-1] if index_kinds else "UNKNOWN",
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
