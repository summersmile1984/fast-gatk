#!/usr/bin/env python3
"""Build M0 cluster-evidence.json from a real scatter-gather run.

The point of this script is to write an evidence file that proves the
M0 infrastructure (slurm_smoke.sh wrapper, scatter_gather_smoke.sh driver,
manifest contract, gather contract) works end-to-end on the *same* code
path that production will run.  It is NOT a substitute for a real SLURM
allocation; the cluster_kind field disambiguates the two.

Invocation:
    python3 fastgatk-native/scripts/build_cluster_evidence.py \
        --shard-dir /tmp/m0-shards \
        --output fastgatk-native/evidence/cluster-evidence.json \
        --cluster-kind local-infrastructure-proof

Exit codes:
    0 = evidence written; schema_valid
    1 = shard manifest count != 8 or any shard manifest incomplete
    2 = bad CLI arguments
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shard-dir", type=pathlib.Path, required=True,
                        help="Output dir of scatter_gather_smoke.sh")
    parser.add_argument("--output", type=pathlib.Path, required=True,
                        help="Path to write cluster-evidence.json")
    parser.add_argument("--cluster-kind", default="local-infrastructure-proof",
                        choices=("slurm", "local-infrastructure-proof"),
                        help="How the shards were run")
    parser.add_argument("--queue", default="compute")
    parser.add_argument("--account", default="biology")
    parser.add_argument("--cpus-per-task", type=int, default=2)
    parser.add_argument("--mem-per-cpu", default="4G")
    return parser.parse_args()


def load_manifest(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    args = parse_args()
    shard_dir: pathlib.Path = args.shard_dir
    if not shard_dir.is_dir():
        print(f"shard dir does not exist: {shard_dir}", file=sys.stderr)
        return 2
    scatter_manifest = shard_dir / "scatter-gather.manifest.json"
    if not scatter_manifest.is_file():
        print(f"scatter manifest missing: {scatter_manifest}", file=sys.stderr)
        return 2
    scatter = load_manifest(scatter_manifest)
    gather_manifest = shard_dir / "gathered.vcf.gz.manifest.json"
    if not gather_manifest.is_file():
        print(f"gather manifest missing: {gather_manifest}", file=sys.stderr)
        return 2
    gather = load_manifest(gather_manifest)
    shard_manifests = sorted(shard_dir.glob("calls-*.vcf.manifest.json"))
    if len(shard_manifests) != 8:
        print(f"M0 requires 8 shard manifests; got {len(shard_manifests)}", file=sys.stderr)
        return 1
    shards: list[dict[str, Any]] = []
    for manifest_path in shard_manifests:
        data = load_manifest(manifest_path)
        outputs = data.get("outputs", [])
        if not outputs or not outputs[0].get("complete"):
            print(f"shard manifest incomplete: {manifest_path}", file=sys.stderr)
            return 1
        telemetry = data.get("telemetry", {})
        shards.append({
            "interval": data.get("interval"),
            "manifest": str(manifest_path),
            "output": data.get("primary_output"),
            "wall_seconds": telemetry.get("execute_seconds"),
            "reads": telemetry.get("reads"),
            "execution_space": telemetry.get("execution_space"),
            "sacct": {},  # no real SLURM; populated only by verify_cluster_smoke.sh on a real cluster
        })
    evidence = {
        "schema_version": 1,
        "status": "pass",
        "cluster_kind": args.cluster_kind,
        "queue": args.queue,
        "account": args.account,
        "cpus_per_task": args.cpus_per_task,
        "mem_per_cpu": args.mem_per_cpu,
        "scatter": scatter,
        "gather": gather,
        "shards": shards,
        "workflow_signature": scatter.get("telemetry", {}).get("workflow_signature"),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(evidence, sort_keys=True, indent=2) + "\n", encoding="utf-8",
    )
    print(json.dumps({"status": "pass",
                      "shards": len(shards),
                      "cluster_kind": args.cluster_kind,
                      "evidence": str(args.output)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())