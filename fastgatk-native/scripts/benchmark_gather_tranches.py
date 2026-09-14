#!/usr/bin/env python3
"""File-boundary benchmark for GatherTranches report merging."""

from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import tempfile
import time


ROOT = pathlib.Path(__file__).resolve().parents[2]
BINARY = ROOT / "fastgatk-native" / "build" / "fastgatk-gather-tranches"
HEADER = (
    "requestedVQSLOD,numKnown,numNovel,knownTiTv,novelTiTv,minVQSLod,"
    "filterName,model,accessibleTruthSites,callsAtTruthSites,truthSensitivity\n"
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shards", type=int, default=32)
    parser.add_argument("--tranches", type=int, default=100)
    args = parser.parse_args()
    if args.shards < 1 or args.tranches < 1:
        raise SystemExit("--shards and --tranches must be positive")
    with tempfile.TemporaryDirectory(prefix="fastgatk-gather-tranches-bench-") as temporary:
        work = pathlib.Path(temporary)
        inputs: list[pathlib.Path] = []
        for shard in range(args.shards):
            path = work / f"shard-{shard}.tranches"
            with path.open("w", encoding="utf-8") as handle:
                handle.write("# Variant quality score tranches file\n# Version number 6\n" + HEADER)
                for index in range(args.tranches):
                    lod = 10.0 - index * 0.1
                    calls = max(1, 1000 - index * 7)
                    handle.write(
                        f"100.0,10,20,2.0,1.5,{lod:.4f},filter,SNP,1000,{calls},{calls / 1000:.4f}\n"
                    )
            inputs.append(path)
        output = work / "gathered.tranches"
        manifest = work / "gathered.manifest.json"
        begin = time.perf_counter()
        result = subprocess.run(
            [str(BINARY), *sum((["-I", str(path)] for path in inputs), []),
             "--mode", "SNP",
             "--truth-sensitivity-tranche", "99", "-O", str(output),
             "--output-manifest", str(manifest)],
            text=True, capture_output=True, check=True,
        )
        elapsed = time.perf_counter() - begin
        records = sum(1 for line in output.read_text(encoding="utf-8").splitlines()
                      if line and not line.startswith("#"))
        manifest_payload = json.loads(manifest.read_text(encoding="utf-8"))
        if manifest_payload["primary_output"] != str(output):
            raise AssertionError(manifest_payload)
        if not manifest_payload["outputs"][0]["complete"]:
            raise AssertionError(manifest_payload)
        if manifest_payload["telemetry"]["merged_tranches"] != args.tranches:
            raise AssertionError(manifest_payload)
        print(json.dumps({
            "status": "pass",
            "tool": "GatherTranches",
            "shards": args.shards,
            "input_tranches_per_shard": args.tranches,
            "merged_tranches": args.tranches,
            "output_records": records,
            "wall_seconds": elapsed,
            "rows_per_second": args.shards * args.tranches / elapsed if elapsed else 0.0,
            "output_bytes": output.stat().st_size,
            "manifest_bytes": manifest.stat().st_size,
            "manifest_contract": True,
            "summary": json.loads(result.stdout),
        }, sort_keys=True))


if __name__ == "__main__":
    main()
