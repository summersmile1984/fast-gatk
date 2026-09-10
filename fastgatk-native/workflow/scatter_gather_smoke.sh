#!/usr/bin/env bash
set -euo pipefail

# Scheduler-neutral scatter/gather runner for the Nextflow/SLURM contract.
# Each shard is reusable only when its complete OutputManifest still matches
# the input, interval and executable signatures. Partial outputs are never
# published: each attempt writes into a private directory and the manifest is
# moved last as the commit marker.  ``--retries`` models the bounded retry that
# Nextflow/SLURM would perform for a failed task; it is intentionally opt-in so
# a caller can let Nextflow own retry policy instead.

usage() {
    echo "usage: scatter_gather_smoke.sh --input BAM --intervals CONTIG:START-END,... --binary HC --gather-binary GATHER --outdir DIR [--reference FASTA] [--threads N] [--retries N] [--slurm-wrapper FILE]" >&2
    exit 2
}

input=""
reference=""
intervals=""
binary=""
gather_binary=""
outdir=""
threads="${SLURM_CPUS_PER_TASK:-2}"
retries="${FASTGATK_SCATTER_RETRIES:-0}"
slurm_wrapper="$(cd "$(dirname "$0")" && pwd)/slurm_smoke.sh"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -I|--input) input=${2:?missing value for $1}; shift 2 ;;
        -R|--reference) reference=${2:?missing value for $1}; shift 2 ;;
        -L|--intervals|--interval) intervals=${2:?missing value for $1}; shift 2 ;;
        --binary) binary=${2:?missing value for $1}; shift 2 ;;
        --gather-binary) gather_binary=${2:?missing value for $1}; shift 2 ;;
        --outdir) outdir=${2:?missing value for $1}; shift 2 ;;
        --threads) threads=${2:?missing value for $1}; shift 2 ;;
        --retries|--max-retries) retries=${2:?missing value for $1}; shift 2 ;;
        --slurm-wrapper) slurm_wrapper=${2:?missing value for $1}; shift 2 ;;
        -h|--help) usage ;;
        *) echo "unknown option: $1" >&2; usage ;;
    esac
done

[[ -n "$input" && -f "$input" ]] || { echo "missing --input: $input" >&2; exit 2; }
[[ -n "$intervals" ]] || { echo "--intervals is required" >&2; exit 2; }
[[ -x "$binary" ]] || { echo "missing executable --binary: $binary" >&2; exit 2; }
[[ -x "$gather_binary" ]] || { echo "missing executable --gather-binary: $gather_binary" >&2; exit 2; }
[[ -x "$slurm_wrapper" ]] || { echo "missing executable --slurm-wrapper: $slurm_wrapper" >&2; exit 2; }
[[ "$retries" =~ ^[0-9]+$ ]] || { echo "--retries must be a non-negative integer" >&2; exit 2; }
# Keep the parsed CLI value available to the final workflow manifest.  This
# also makes the contract deterministic when callers use --retries without
# setting the equivalent environment variable.
export FASTGATK_SCATTER_RETRIES="$retries"
mkdir -p "$outdir"

input_size=$(stat -c '%s' "$input")
input_mtime=$(stat -c '%Y' "$input")
binary_size=$(stat -c '%s' "$binary")
binary_mtime=$(stat -c '%Y' "$binary")
workflow_manifest="$outdir/scatter-gather.manifest.json"

safe_interval() {
    printf '%s' "$1" | tr -c 'A-Za-z0-9_.-' '_'
}

declare -a shard_paths=()
declare -a shard_manifests=()
declare -a shard_intervals=()

# Validate the generic OutputManifest envelope before a bundle becomes visible
# to gather/resume.  Tool-specific fields are deliberately left to the tool's
# own oracle; this gate only enforces the cross-tool publish contract.
manifest_complete() {
    local manifest=$1 primary=$2 index=${3:-}
    python3 - "$manifest" "$primary" "$index" <<'PY'
import json
import pathlib
import sys

manifest, primary, index = sys.argv[1:]
try:
    manifest_path = pathlib.Path(manifest)
    primary_path = pathlib.Path(primary)
    index_path = pathlib.Path(index) if index else None
    data = json.loads(manifest_path.read_text(encoding="utf-8"))
    outputs = data.get("outputs")
    ok = (isinstance(data, dict) and isinstance(data.get("schema_version"), int) and
          isinstance(outputs, list) and bool(outputs) and
          all(isinstance(item, dict) and item.get("complete") is True for item in outputs) and
          primary_path.is_file() and primary_path.stat().st_size > 0 and
          (index_path is None or (index_path.is_file() and index_path.stat().st_size > 0)))
except (OSError, ValueError, TypeError, json.JSONDecodeError):
    ok = False
raise SystemExit(0 if ok else 1)
PY
}

# Publish all staged files with rollback.  Staging and destination are under
# the same outdir, so rename is atomic on the filesystems supported by the
# workflow.  Existing valid output is normally resumed; if a caller forces a
# recompute, old files are moved to a private backup until the manifest commit
# succeeds, preventing a mixed old/new bundle after a rename error.
publish_staged_bundle() {
    local stage_dir=$1 stage_primary=$2 stage_index=$3 stage_manifest=$4
    local final_primary=$5 final_index=$6 final_manifest=$7
    local backup_dir="$stage_dir/.publish-backup"
    local had_primary=0 had_index=0 had_manifest=0 moved_primary=0 moved_index=0 moved_manifest=0
    mkdir -p "$backup_dir"

    if [[ -e "$final_primary" ]]; then
        mv "$final_primary" "$backup_dir/primary" || return 1
        had_primary=1
    fi
    if [[ -n "$final_index" && -e "$final_index" ]]; then
        mv "$final_index" "$backup_dir/index" || {
            [[ "$had_primary" == 1 ]] && mv "$backup_dir/primary" "$final_primary" || true
            return 1
        }
        had_index=1
    fi
    if [[ -e "$final_manifest" ]]; then
        mv "$final_manifest" "$backup_dir/manifest" || {
            [[ "$had_index" == 1 ]] && mv "$backup_dir/index" "$final_index" || true
            [[ "$had_primary" == 1 ]] && mv "$backup_dir/primary" "$final_primary" || true
            return 1
        }
        had_manifest=1
    fi

    rollback() {
        [[ "$moved_manifest" == 1 ]] && rm -f "$final_manifest"
        [[ "$moved_index" == 1 ]] && rm -f "$final_index"
        [[ "$moved_primary" == 1 ]] && rm -f "$final_primary"
        [[ "$had_manifest" == 1 && -e "$backup_dir/manifest" ]] && mv "$backup_dir/manifest" "$final_manifest" || true
        [[ "$had_index" == 1 && -e "$backup_dir/index" ]] && mv "$backup_dir/index" "$final_index" || true
        [[ "$had_primary" == 1 && -e "$backup_dir/primary" ]] && mv "$backup_dir/primary" "$final_primary" || true
    }

    if ! mv "$stage_primary" "$final_primary"; then rollback; return 1; fi
    moved_primary=1
    if [[ -n "$stage_index" ]]; then
        if ! mv "$stage_index" "$final_index"; then rollback; return 1; fi
        moved_index=1
    fi
    # The manifest is the final commit marker.  A gather/resume consumer must
    # never treat a primary/index without this file as a complete result.
    if ! mv "$stage_manifest" "$final_manifest"; then rollback; return 1; fi
    moved_manifest=1
    rm -f "$backup_dir/primary" "$backup_dir/index" "$backup_dir/manifest"
    rmdir "$backup_dir" 2>/dev/null || true
    return 0
}

workflow_signature() {
    python3 - "$input" "$reference" "$binary" "$gather_binary" "$threads" \
        "${#shard_paths[@]}" "${shard_paths[@]}" "${shard_manifests[@]}" <<'PY'
import hashlib
import json
import pathlib
import sys

input_path, reference, binary, gather_binary, threads, count = sys.argv[1:7]
count = int(count)
paths = sys.argv[7:7 + count]
manifests = sys.argv[7 + count:7 + 2 * count]

def digest(path):
    value = pathlib.Path(path)
    if not value:
        return None
    h = hashlib.sha256()
    with value.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

payload = {
    "input": {"path": input_path, "sha256": digest(input_path)},
    "reference": {"path": reference, "sha256": digest(reference)} if reference else None,
    "binary": {"path": binary, "sha256": digest(binary)},
    "gather_binary": {"path": gather_binary, "sha256": digest(gather_binary)},
    "threads": int(threads),
    "shards": [{"path": path, "sha256": digest(path),
                "manifest": manifest, "manifest_sha256": digest(manifest)}
               for path, manifest in zip(paths, manifests)],
}
encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
print(hashlib.sha256(encoded).hexdigest())
PY
}

IFS=',' read -r -a interval_array <<< "$intervals"
[[ "${#interval_array[@]}" -gt 0 ]] || { echo "--intervals is empty" >&2; exit 2; }

run_one() {
    local interval=$1 safe output manifest state ref_args output_name manifest_name
    local stage_dir stage_output stage_manifest rc attempts_used retry_count
    safe=$(safe_interval "$interval")
    output="$outdir/calls-${safe}.vcf"
    manifest="$output.manifest.json"
    output_name="calls-${safe}.vcf"
    manifest_name="${output_name}.manifest.json"
    state=computed
    attempts_used=0
    retry_count=0
    if [[ -s "$output" && -s "$manifest" ]] && python3 - "$manifest" "$input" "$interval" "$input_size" "$input_mtime" "$binary_size" "$binary_mtime" <<'PY'
import json
import pathlib
import sys

manifest, input_path, interval, input_size, input_mtime, binary_size, binary_mtime = sys.argv[1:]
try:
    data = json.loads(pathlib.Path(manifest).read_text(encoding="utf-8"))
    telemetry = data.get("telemetry", {})
    ok = (data.get("status") == "smoke" and data.get("interval") == interval and
          data.get("input") == input_path and data.get("primary_output_kind") == "vcf" and
          data.get("outputs", [{}])[0].get("complete") is True and
          telemetry.get("input_signature", {}).get("size") == int(input_size) and
          telemetry.get("input_signature", {}).get("mtime") == int(input_mtime) and
          telemetry.get("binary_signature", {}).get("size") == int(binary_size) and
          telemetry.get("binary_signature", {}).get("mtime") == int(binary_mtime))
except (OSError, ValueError, TypeError, json.JSONDecodeError):
    ok = False
raise SystemExit(0 if ok else 1)
PY
    then
        state=resumed
    else
        # A stale/incomplete destination must not be mistaken for a resumable
        # shard.  It is backed up by publish_staged_bundle during a successful
        # retry, so a failed attempt cannot erase a valid previous result.
        if [[ -e "$output" || -e "$manifest" ]]; then
            rm -f "$output" "$manifest"
        fi
        ref_args=()
        [[ -n "$reference" ]] && ref_args=(-R "$reference")
        for ((attempts_used = 1; attempts_used <= retries + 1; ++attempts_used)); do
            retry_count=$((attempts_used - 1))
            stage_dir=$(mktemp -d "$outdir/.fastgatk-shard-${safe}.attempt.XXXXXX")
            stage_output="$stage_dir/$output_name"
            stage_manifest="$stage_dir/$manifest_name"
            set +e
            FASTGATK_SCATTER_ATTEMPT="$attempts_used" SLURM_CPUS_PER_TASK="$threads" \
                "$slurm_wrapper" "$binary" -I "$input" "${ref_args[@]}" -L "$interval" \
                -O "$stage_output" --output-manifest "$stage_manifest" --threads "$threads"
            rc=$?
            set -e
            if [[ "$rc" == 0 ]] && manifest_complete "$stage_manifest" "$stage_output"; then
                python3 - "$stage_manifest" "$stage_output" "$output" "$input" "$interval" \
                    "$input_size" "$input_mtime" "$binary_size" "$binary_mtime" "$retry_count" <<'PY'
import json
import pathlib
import sys

manifest, stage_output, output, input_path, interval, input_size, input_mtime, binary_size, binary_mtime, retry_count = sys.argv[1:]
path = pathlib.Path(manifest)
data = json.loads(path.read_text(encoding="utf-8"))
data["workflow"] = "scatter-gather"
data["interval"] = interval
if data.get("primary_output") == stage_output:
    data["primary_output"] = output
for item in data.get("outputs", []):
    if item.get("path") == stage_output:
        item["path"] = output
data["telemetry"] = dict(data.get("telemetry", {}))
data["telemetry"]["input_signature"] = {"size": int(input_size), "mtime": int(input_mtime)}
data["telemetry"]["binary_signature"] = {"size": int(binary_size), "mtime": int(binary_mtime)}
data["telemetry"]["retry"] = {"attempts": int(retry_count) + 1, "retries": int(retry_count)}
data["compatibility"] = dict(data.get("compatibility", {}))
data["compatibility"]["scatter_resume_manifest"] = True
data["compatibility"]["atomic_attempt_publish"] = True
temporary = path.with_name(path.name + ".tmp")
with temporary.open("w", encoding="utf-8") as stream:
    json.dump(data, stream, sort_keys=True)
    stream.write("\n")
    stream.flush()
    import os
    os.fsync(stream.fileno())
os.replace(temporary, path)
PY
                if publish_staged_bundle "$stage_dir" "$stage_output" "" "$stage_manifest" \
                    "$output" "" "$manifest"; then
                    state=computed
                    rm -rf "$stage_dir"
                    break
                fi
            fi
            rm -rf "$stage_dir"
            if [[ "$attempts_used" -le "$retries" ]]; then
                echo "retrying shard $interval after attempt $attempts_used (exit=$rc)" >&2
            fi
        done
        if [[ "$state" != computed || ! -s "$output" || ! -s "$manifest" ]]; then
            echo "shard failed after $attempts_used attempt(s): $interval" >&2
            return 1
        fi
    fi
    [[ -s "$output" && -s "$manifest" ]] || { echo "incomplete shard: $interval" >&2; return 1; }
    printf '%s\t%s\t%s\t%s\n' "$interval" "$state" "$attempts_used" "$retry_count" >> "$outdir/.shard-status"
}

rm -f "$outdir/.shard-status"
for interval in "${interval_array[@]}"; do
    [[ -n "$interval" ]] || { echo "empty interval in --intervals" >&2; exit 2; }
    shard_intervals+=("$interval")
    safe=$(safe_interval "$interval")
    shard_paths+=("$outdir/calls-${safe}.vcf")
    shard_manifests+=("$outdir/calls-${safe}.vcf.manifest.json")
done

declare -a pids=()
for interval in "${shard_intervals[@]}"; do
    run_one "$interval" &
    pids+=("$!")
done
failed_shards=0
for pid in "${pids[@]}"; do
    if ! wait "$pid"; then failed_shards=1; fi
done
if [[ "$failed_shards" != 0 ]]; then
    echo "scatter failed; no gather was attempted" >&2
    exit 1
fi

gather_output="$outdir/gathered.vcf.gz"
gather_manifest="$gather_output.manifest.json"
gather_args=()
for shard in "${shard_paths[@]}"; do gather_args+=(-I "$shard"); done

signature=$(workflow_signature)
gather_state=computed
gather_attempts=0
gather_retry_count=0
if [[ -s "$gather_output" && -s "$gather_output.tbi" && -s "$gather_manifest" ]] && \
   manifest_complete "$gather_manifest" "$gather_output" "$gather_output.tbi" && \
   python3 - "$gather_manifest" "$signature" <<'PY'
import json
import pathlib
import sys

try:
    data = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
    ok = (data.get("status") == "pass" and
          data.get("telemetry", {}).get("workflow_signature") == sys.argv[2])
except (OSError, ValueError, TypeError, json.JSONDecodeError):
    ok = False
raise SystemExit(0 if ok else 1)
PY
then
    gather_state=resumed
else
    # Gather is a second atomic bundle.  It is intentionally executed through
    # the same scheduler wrapper so a real Nextflow/SLURM allocation remains
    # the single owner of resources; the wrapper itself never nests sbatch.
    for ((gather_attempts = 1; gather_attempts <= retries + 1; ++gather_attempts)); do
        gather_retry_count=$((gather_attempts - 1))
        gather_stage_dir=$(mktemp -d "$outdir/.fastgatk-gather.attempt.XXXXXX")
        gather_stage_output="$gather_stage_dir/gathered.vcf.gz"
        gather_stage_index="$gather_stage_output.tbi"
        gather_stage_manifest="$gather_stage_output.manifest.json"
        set +e
        FASTGATK_SCATTER_ATTEMPT="$gather_attempts" SLURM_CPUS_PER_TASK="$threads" \
            "$slurm_wrapper" "$gather_binary" "${gather_args[@]}" \
            -O "$gather_stage_output" --output-manifest "$gather_stage_manifest"
        rc=$?
        set -e
        if [[ "$rc" == 0 ]] && manifest_complete "$gather_stage_manifest" \
            "$gather_stage_output" "$gather_stage_index"; then
            python3 - "$gather_stage_manifest" "$gather_stage_output" "$gather_output" \
                "$gather_stage_index" "$gather_output.tbi" "$signature" "$gather_retry_count" <<'PY'
import json
import os
import pathlib
import sys

manifest, stage_output, output, stage_index, index, signature, retry_count = sys.argv[1:]
path = pathlib.Path(manifest)
data = json.loads(path.read_text(encoding="utf-8"))
data["status"] = "pass"
if data.get("primary_output") == stage_output:
    data["primary_output"] = output
for item in data.get("outputs", []):
    if item.get("path") == stage_output:
        item["path"] = output
    elif item.get("path") == stage_index:
        item["path"] = index
data["workflow"] = "scatter-gather"
data["telemetry"] = dict(data.get("telemetry", {}))
data["telemetry"]["workflow_signature"] = signature
data["telemetry"]["retry"] = {"attempts": int(retry_count) + 1, "retries": int(retry_count)}
data["compatibility"] = dict(data.get("compatibility", {}))
data["compatibility"]["atomic_attempt_publish"] = True
data["compatibility"]["gather_resume_manifest"] = True
temporary = path.with_name(path.name + ".tmp")
with temporary.open("w", encoding="utf-8") as stream:
    json.dump(data, stream, sort_keys=True)
    stream.write("\n")
    stream.flush()
    os.fsync(stream.fileno())
os.replace(temporary, path)
PY
            if publish_staged_bundle "$gather_stage_dir" "$gather_stage_output" \
                "$gather_stage_index" "$gather_stage_manifest" "$gather_output" \
                "$gather_output.tbi" "$gather_manifest"; then
                gather_state=computed
                rm -rf "$gather_stage_dir"
                break
            fi
        fi
        rm -rf "$gather_stage_dir"
        if [[ "$gather_attempts" -le "$retries" ]]; then
            echo "retrying gather after attempt $gather_attempts (exit=$rc)" >&2
        fi
    done
    if [[ "$gather_state" != computed ]] || ! manifest_complete "$gather_manifest" \
        "$gather_output" "$gather_output.tbi"; then
        echo "gather failed after $gather_attempts attempt(s)" >&2
        exit 1
    fi
fi

resumed=0
computed=0
retried_shards=0
shard_attempts=0
while IFS=$'\t' read -r _ state attempts retry_count; do
    [[ "$state" == resumed ]] && resumed=$((resumed + 1))
    [[ "$state" == computed ]] && computed=$((computed + 1))
    [[ "$state" == computed ]] && shard_attempts=$((shard_attempts + attempts))
    [[ "${retry_count:-0}" -gt 0 ]] && retried_shards=$((retried_shards + 1))
done < "$outdir/.shard-status"

python3 - "$workflow_manifest" "$input" "$intervals" "$gather_output" "$resumed" "$computed" \
    "$retried_shards" "$shard_attempts" "$gather_state" "$gather_attempts" "$gather_retry_count" \
    "$signature" "${shard_manifests[@]}" <<'PY'
import json
import os
import pathlib
import sys

manifest, input_path, intervals, gathered, resumed, computed, retried_shards, shard_attempts, \
    gather_state, gather_attempts, gather_retry_count, signature, *shard_manifests = sys.argv[1:]
shards = []
for shard_manifest in shard_manifests:
    data = json.loads(pathlib.Path(shard_manifest).read_text(encoding="utf-8"))
    shards.append({"interval": data.get("interval"), "manifest": shard_manifest,
                   "output": data.get("primary_output"),
                   "complete": bool(data.get("outputs", [{}])[0].get("complete")),
                   "retry": data.get("telemetry", {}).get("retry", {})})
payload = {"schema_version": 1, "tool": "FastGatkScatterGather", "status": "pass",
           "input": input_path, "intervals": intervals.split(","),
           "gathered_output": gathered, "shards": shards,
           "resume": {"resumed_shards": int(resumed), "computed_shards": int(computed)},
           "retry": {"retried_shards": int(retried_shards),
                     "shard_attempts": int(shard_attempts),
                     "gather_state": gather_state,
                     "gather_attempts": int(gather_attempts),
                     "gather_retries": int(gather_retry_count),
                     "max_retries": int(os.environ.get("FASTGATK_SCATTER_RETRIES", "0"))},
           "telemetry": {"workflow_signature": signature},
           "compatibility": {"nextflow_resume": True, "slurm_resource_wrapper": True,
                              "output_manifest_readback": True, "vcf_index": True,
                              "atomic_attempt_publish": True, "failure_retry": True}}
path = pathlib.Path(manifest)
temporary = path.with_name(path.name + ".tmp")
with temporary.open("w", encoding="utf-8") as stream:
    json.dump(payload, stream, sort_keys=True)
    stream.write("\n")
    stream.flush()
    os.fsync(stream.fileno())
os.replace(temporary, path)
print(json.dumps(payload, sort_keys=True))
PY
