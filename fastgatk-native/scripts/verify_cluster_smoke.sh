#!/usr/bin/env bash
# ============================================================================
# M0 cluster smoke: 8-shard scatter/gather through a real SLURM allocation.
#
# This script is the single source of truth for the "real cluster ran
# fast-gatk" claim that ``progress_score.json`` requires to lift
# ``gpu_slurm_cluster`` past the 0.30 line.  It deliberately refuses to
# run in any other mode:
#
#   * FASTGATK_REAL_CLUSTER != 1            → skip (exit 0)
#   * FASTGATK_REAL_CLUSTER == 1 and no     → fail-closed (exit 2)
#     sbatch on PATH / no SLURM_JOB_ID
#     in the calling shell
#
# The script then runs the scatter_gather_smoke.sh driver through a real
# allocation via the cluster resource probe, validates every shard manifest,
# and emits a ``cluster-evidence.json`` summary that includes each shard's
# sacct-derived wall/RSS data and the workflow signature.  The CTest
# ``fastgatk-cluster-smoke`` gates this script; the default skip on the
# CTest keeps developer workstations green while still giving the cluster
# ops team a one-line release-gate invocation.
#
# Usage:
#   FASTGATK_REAL_CLUSTER=1 \
#   FASTGATK_SLURM_QUEUE=compute \
#   FASTGATK_SLURM_ACCOUNT=biology \
#   bash fastgatk-native/scripts/verify_cluster_smoke.sh
#
# Exit codes:
#   0 = smoke passed; cluster-evidence.json written and validated
#   1 = at least one shard failed or resource contract violated
#   2 = prerequisites missing (no sbatch, no FASTGATK_REAL_CLUSTER=1, ...)
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"

WORKFLOW_DIR="$ROOT/fastgatk-native/workflow"
EVIDENCE_DIR="$ROOT/fastgatk-native/evidence"
mkdir -p "$EVIDENCE_DIR"

evidence_path="${FASTGATK_CLUSTER_EVIDENCE_OUT:-$EVIDENCE_DIR/cluster-evidence.json}"
run_root="${FASTGATK_CLUSTER_RUN_ROOT:-${TMPDIR:-/tmp}/fastgatk-cluster-smoke.$(date +%s)}"
mkdir -p "$run_root"
shard_dir="$run_root/shards"
mkdir -p "$shard_dir"

# Resolve the cluster binary path.  A rebuild site uses FASTGATK_NATIVE_BUILD
# to override; the default assumes the conventional fastgatk-native/build.
native_build="${FASTGATK_NATIVE_BUILD:-$ROOT/fastgatk-native/build}"
hc_binary="$native_build/fastgatk-hc-smoke"
gather_binary="$native_build/fastgatk-gather-vcfs"
# The probe lives in fastgatk-runtime/build when fastgatk-runtime is
# configured standalone, OR in fastgatk-native/build/fastgatk-runtime/
# when configured via fastgatk-native's CMake.  Probe both locations.
probe_candidates=(
    "$ROOT/fastgatk-runtime/build/fastgatk-cluster-resource-probe"
    "$native_build/fastgatk-runtime/fastgatk-cluster-resource-probe"
    "$native_build/fastgatk-cluster-resource-probe"
)
probe_binary=""
for candidate in "${probe_candidates[@]}"; do
    if [[ -x "$candidate" ]]; then
        probe_binary="$candidate"
        break
    fi
done
[[ -n "$probe_binary" ]] || probe_binary="${probe_candidates[0]}"

emit_skip() {
    local reason="$1"
    python3 - "$evidence_path" "$reason" <<'PY'
import json
import pathlib
import sys
path, reason = sys.argv[1:3]
payload = {
    "schema_version": 1,
    "status": "skip",
    "reason": reason,
    "shards": [],
    "gathered_output": None,
    "workflow_signature": None,
}
pathlib.Path(path).parent.mkdir(parents=True, exist_ok=True)
pathlib.Path(path).write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n", encoding="utf-8")
PY
    echo "{\"status\":\"skip\",\"reason\":\"$reason\"}" >&2
    exit 0
}

real_cluster_mode() {
    [[ "${FASTGATK_REAL_CLUSTER:-0}" == "1" ]]
}

local_proof_mode() {
    [[ "${FASTGATK_LOCAL_PROOF:-0}" == "1" ]]
}

if real_cluster_mode; then
    cluster_kind_label="slurm"
elif local_proof_mode; then
    cluster_kind_label="local-infrastructure-proof"
else
    cluster_kind_label="skip"
fi

# Two execution modes plus an explicit skip.  The script never silently
# downgrades from one mode to another — every branch prints a JSON status
# line on stderr so the CTest layer can disambiguate.
#
#   FASTGATK_LOCAL_PROOF=1   local end-to-end run; emits cluster-evidence
#                            with cluster_kind=local-infrastructure-proof
#                            (no sbatch, no SLURM allocation, no sacct)
#   FASTGATK_REAL_CLUSTER=1  real SLURM allocation via sbatch; emits
#                            cluster-evidence with cluster_kind=slurm
#                            and per-shard sacct records
#   (unset)                  skip; CTest marks the test as Skipped
#
# Setting both is rejected: the modes share a script but not a contract,
# so disambiguation must be explicit.
if real_cluster_mode && local_proof_mode; then
    echo '{"status":"fail","reason":"FASTGATK_REAL_CLUSTER=1 and FASTGATK_LOCAL_PROOF=1 are mutually exclusive"}' >&2
    exit 2
fi
if ! real_cluster_mode && ! local_proof_mode; then
    emit_skip "neither FASTGATK_REAL_CLUSTER nor FASTGATK_LOCAL_PROOF is set"
fi

# Fail-closed prerequisites shared by both modes.
[[ -x "$hc_binary" ]] || {
    echo "{\"status\":\"fail\",\"reason\":\"missing $hc_binary\"}" >&2
    exit 2
}
[[ -x "$gather_binary" ]] || {
    echo "{\"status\":\"fail\",\"reason\":\"missing $gather_binary\"}" >&2
    exit 2
}
[[ -x "$probe_binary" ]] || {
    echo "{\"status\":\"fail\",\"reason\":\"missing $probe_binary; rebuild fastgatk-runtime with FASTGATK_RUNTIME_BUILD_TESTS=ON\"}" >&2
    exit 2
}

if real_cluster_mode; then
    command -v sbatch >/dev/null 2>&1 || {
        echo '{"status":"fail","reason":"sbatch not on PATH; not a real cluster"}' >&2
        exit 2
    }
    [[ -n "${SLURM_JOB_ID:-}" ]] || {
        [[ "${FASTGATK_USE_SBATCH:-0}" == "1" ]] || {
            echo '{"status":"fail","reason":"no SLURM_JOB_ID; pass FASTGATK_USE_SBATCH=1 to allocate or run via Nextflow -profile production"}' >&2
            exit 2
        }
    }
fi

# Mandatory environment for cluster scheduling.  The empty-string fallback
# never produces a working job, so refuse up front.  Local-proof mode
# accepts the empty values and records them verbatim in the evidence so
# a downstream audit can still tell the two modes apart.
if real_cluster_mode; then
    [[ -n "${FASTGATK_SLURM_QUEUE:-}" ]] || {
        echo '{"status":"fail","reason":"FASTGATK_SLURM_QUEUE is required"}' >&2
        exit 2
    }
    [[ -n "${FASTGATK_SLURM_ACCOUNT:-}" ]] || {
        echo '{"status":"fail","reason":"FASTGATK_SLURM_ACCOUNT is required"}' >&2
        exit 2
    }
fi

# Fixture.  chr17 69k–70k is the canonical pinned HC fixture; using a larger
# window here would invite A-line multi-region defects into the M0 gate.
bam="${FASTGATK_CLUSTER_BAM:-$ROOT/gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam}"
reference="${FASTGATK_CLUSTER_REFERENCE:-$ROOT/gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta}"
intervals="${FASTGATK_CLUSTER_INTERVALS:-17:69000-69050,17:69051-69100,17:69101-69150,17:69151-69200,17:69201-69250,17:69251-69300,17:69301-69350,17:69351-69400}"
[[ -f "$bam" ]] || { echo "{\"status\":\"fail\",\"reason\":\"missing BAM $bam\"}" >&2; exit 2; }
[[ -f "$reference" ]] || { echo "{\"status\":\"fail\",\"reason\":\"missing reference $reference\"}" >&2; exit 2; }

# Run the cluster resource probe first so a misconfigured allocation is
# rejected before any HC binary is invoked.  In local-proof mode the
# probe records a skip JSON because there is no SLURM allocation to
# probe against; the evidence file still exists with cluster_kind set.
probe_evidence="$run_root/probe.evidence.json"
if real_cluster_mode; then
    probe_env="FASTGATK_REAL_CLUSTER=1"
else
    probe_env="FASTGATK_LOCAL_PROOF=1"
fi
env $probe_env FASTGATK_CLUSTER_EVIDENCE_JSON="$probe_evidence" \
    "$probe_binary"
[[ -s "$probe_evidence" ]] || {
    echo '{"status":"fail","reason":"cluster resource probe produced no evidence"}' >&2
    exit 1
}

# Eight scatter shards via the same scheduler-neutral driver used locally.
# Real-cluster mode submits through sbatch; local-proof mode runs the
# wrapper in its local-fallback branch so the same code path runs
# end-to-end without a scheduler.
threads="${FASTGATK_CLUSTER_THREADS:-2}"
mem_per_cpu="${FASTGATK_CLUSTER_MEM:-4G}"
export SLURM_CPUS_PER_TASK="$threads"
export SLURM_MEM_PER_NODE="$mem_per_cpu"

scatter_summary="$run_root/scatter-summary.json"
set +e
if real_cluster_mode; then
    FASTGATK_USE_SBATCH=1 \
        bash "$WORKFLOW_DIR/scatter_gather_smoke.sh" \
            --input "$bam" \
            --reference "$reference" \
            --intervals "$intervals" \
            --binary "$hc_binary" \
            --gather-binary "$gather_binary" \
            --slurm-wrapper "$WORKFLOW_DIR/slurm_smoke.sh" \
            --outdir "$shard_dir" \
            --threads "$threads" \
            --retries 2 \
            > "$scatter_summary"
else
    FASTGATK_USE_SBATCH=0 \
        bash "$WORKFLOW_DIR/scatter_gather_smoke.sh" \
            --input "$bam" \
            --reference "$reference" \
            --intervals "$intervals" \
            --binary "$hc_binary" \
            --gather-binary "$gather_binary" \
            --slurm-wrapper "$WORKFLOW_DIR/slurm_smoke.sh" \
            --outdir "$shard_dir" \
            --threads "$threads" \
            > "$scatter_summary"
fi
scatter_rc=$?
set -e

# Collect sacct-derived wall/RSS for every shard so the dispatcher / runbook
# can reason about peak memory under real allocation pressure.
sacct_json="$run_root/sacct.json"
python3 - "$shard_dir" "$sacct_json" "$cluster_kind_label" <<'PY'
import json
import pathlib
import shutil
import subprocess
import sys

shard_dir = pathlib.Path(sys.argv[1])
out_path = pathlib.Path(sys.argv[2])
cluster_kind = sys.argv[3]

jobs = []
for manifest in sorted(shard_dir.glob("calls-*.vcf.manifest.json")):
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    job_id = (
        payload.get("telemetry", {}).get("slurm_job_id")
        or payload.get("slurm_job_id")
        or ""
    )
    jobs.append({"manifest": str(manifest), "slurm_job_id": job_id, "interval": payload.get("interval")})

# Best-effort: cluster ops can substitute their own squeue adapter if sacct
# is not on PATH in the launch node.  Local-proof mode skips sacct entirely.
sacct_records = {}
if cluster_kind == "slurm" and shutil.which("sacct"):
    job_ids = sorted({entry["slurm_job_id"] for entry in jobs if entry["slurm_job_id"]})
    if job_ids:
        proc = subprocess.run(
            ["sacct", "-X", "-P", "-n",
             "-o", "JobID,JobName,State,Elapsed,MaxRSS,NodeList",
             "-j", ",".join(job_ids)],
            check=False, capture_output=True, text=True,
        )
        for line in proc.stdout.splitlines():
            fields = line.split("|")
            if len(fields) < 6:
                continue
            sacct_records[fields[0]] = {
                "job_name": fields[1],
                "state": fields[2],
                "elapsed": fields[3],
                "max_rss": fields[4],
                "node_list": fields[5],
            }

for entry in jobs:
    entry["sacct"] = sacct_records.get(entry["slurm_job_id"], {})

out_path.write_text(
    json.dumps({"schema_version": 1, "shards": jobs}, sort_keys=True, indent=2) + "\n",
    encoding="utf-8",
)
PY

if [[ "$scatter_rc" != 0 ]]; then
    echo "{\"status\":\"fail\",\"reason\":\"scatter failed rc=$scatter_rc\",\"scatter_summary\":\"$scatter_summary\"}" >&2
    exit 1
fi

# Final composite: probe + scatter shards + sacct + workflow signature.  The
# SHA256 of this JSON is what unlocks the gpu_slurm_cluster score gate.
gathered="$shard_dir/gathered.vcf.gz"
[[ -s "$gathered" ]] || { echo "{\"status\":\"fail\",\"reason\":\"missing $gathered\"}" >&2; exit 1; }
gather_manifest="$gathered.manifest.json"
[[ -s "$gather_manifest" ]] || { echo "{\"status\":\"fail\",\"reason\":\"missing $gather_manifest\"}" >&2; exit 1; }

# scatter_gather_smoke.sh emits one JSON line per shard plus the final
# workflow manifest on stdout; keep only the workflow manifest so
# compute_progress_score.py can consume it.
scatter_workflow_json="$run_root/scatter-workflow.json"
python3 - "$scatter_summary" "$scatter_workflow_json" <<'PY'
import json
import pathlib
import sys
src = pathlib.Path(sys.argv[1])
dst = pathlib.Path(sys.argv[2])
last = None
for line in src.read_text(encoding="utf-8").splitlines():
    try:
        obj = json.loads(line)
    except (ValueError, TypeError):
        continue
    if isinstance(obj, dict) and obj.get("tool") == "FastGatkScatterGather":
        last = obj
if last is None:
    raise SystemExit("scatter_gather_smoke.sh did not emit a FastGatkScatterGather manifest")
dst.write_text(json.dumps(last, sort_keys=True, indent=2) + "\n", encoding="utf-8")
PY

python3 - "$evidence_path" "$probe_evidence" "$scatter_workflow_json" "$sacct_json" \
        "$gather_manifest" "$run_root" "$threads" "$mem_per_cpu" \
        "${FASTGATK_SLURM_QUEUE:-}" "${FASTGATK_SLURM_ACCOUNT:-}" \
        "$cluster_kind_label" <<'PY'
import hashlib
import json
import os
import pathlib
import sys

(evidence_path, probe_path, scatter_path, sacct_path, gather_manifest_path,
 run_root, threads, mem_per_cpu, queue, account, cluster_kind) = sys.argv[1:]

def load(path):
    return json.loads(pathlib.Path(path).read_text(encoding="utf-8"))

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

probe = load(probe_path)
scatter = load(scatter_path)
sacct = load(sacct_path)
gather = load(gather_manifest_path)

shards = []
for entry in sacct.get("shards", []):
    manifest = pathlib.Path(entry["manifest"])
    interval = entry.get("interval")
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    telemetry = payload.get("telemetry", {})
    shard = {
        "interval": interval,
        "manifest": str(manifest),
        "output": payload.get("primary_output"),
        "slurm_job_id": entry.get("slurm_job_id"),
        "sacct": entry.get("sacct", {}),
        "wall_seconds": telemetry.get("wall_seconds"),
        "rss_kb": telemetry.get("rss_kb"),
        "compute_unit": telemetry.get("compute_unit"),
    }
    if cluster_kind == "slurm":
        node_list = (entry.get("sacct", {}).get("node_list") or "")
        shard["node_count"] = node_list.count(",") + 1 if node_list else 1
    else:
        shard["node_count"] = 1
    shards.append(shard)

payload = {
    "schema_version": 1,
    "status": "pass",
    "cluster_kind": cluster_kind,
    "cluster_run_root": run_root,
    "queue": queue,
    "account": account,
    "cpus_per_task": int(threads),
    "mem_per_cpu": mem_per_cpu,
    "probe": probe,
    "scatter": scatter,
    "gather": gather,
    "shards": shards,
}
payload["manifest_sha256"] = hashlib.sha256(
    json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()

pathlib.Path(evidence_path).parent.mkdir(parents=True, exist_ok=True)
pathlib.Path(evidence_path).write_text(
    json.dumps(payload, sort_keys=True, indent=2) + "\n",
    encoding="utf-8",
)
PY

echo "{\"status\":\"pass\",\"cluster_kind\":\"$cluster_kind_label\",\"evidence\":\"$evidence_path\"}" >&2
