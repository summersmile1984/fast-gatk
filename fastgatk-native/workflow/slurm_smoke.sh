#!/usr/bin/env bash
set -euo pipefail

# SLURM contract wrapper. On a cluster set FASTGATK_USE_SBATCH=1 to submit a
# real allocation. The default is a local simulation, preserving the resource
# environment used by the native dispatcher and needing no cluster.
if [[ $# -lt 1 ]]; then
    echo "usage: slurm_smoke.sh FASTGATK_BINARY [arguments...]" >&2
    exit 2
fi
binary=$1
shift

inside_allocation=0
if [[ -n "${SLURM_JOB_ID:-}" || -n "${SLURM_STEP_ID:-}" ]]; then
    inside_allocation=1
fi

# A Nextflow process submitted with the SLURM executor already runs inside a
# job allocation.  Never submit a nested sbatch in that case, even when a
# site-wide profile happens to export FASTGATK_USE_SBATCH=1.  The explicit
# override variables are useful for a launcher invoked outside Nextflow (and
# make local wrapper tests deterministic) while preserving SLURM's exported
# values inside an existing allocation.
if [[ "${FASTGATK_USE_SBATCH:-0}" == "1" && "$inside_allocation" == "0" && -n "$(command -v sbatch || true)" ]]; then
    cpus=${FASTGATK_SLURM_CPUS_PER_TASK:-${SLURM_CPUS_PER_TASK:-2}}
    # Keep resource translation at the workflow boundary.  Native code reads
    # the allocation that SLURM exports; it must not submit a nested job or
    # silently exceed the cgroup.  SLURM accepts either an integer MiB value
    # or a suffixed value such as 4G, so preserve the caller's spelling.
    job_name=${FASTGATK_SLURM_JOB_NAME:-fastgatk-smoke}
    sbatch_args=(--wait --job-name="$job_name" --cpus-per-task="$cpus")
    memory=${FASTGATK_SLURM_MEMORY:-${SLURM_MEM_PER_NODE:-}}
    memory_per_cpu=${SLURM_MEM_PER_CPU:-}
    if [[ -n "$memory" ]]; then
        sbatch_args+=(--mem="$memory")
    elif [[ -n "$memory_per_cpu" ]]; then
        sbatch_args+=(--mem-per-cpu="$memory_per_cpu")
    fi
    gres=${FASTGATK_SLURM_GRES:-}
    gpu_count=${FASTGATK_SLURM_GPU_COUNT:-${FASTGATK_GPU_COUNT:-}}
    gpu_type=${FASTGATK_SLURM_GPU_TYPE:-${FASTGATK_GPU_TYPE:-}}
    if [[ -n "$gres" ]]; then
        sbatch_args+=(--gres="$gres")
    elif [[ -n "$gpu_count" ]]; then
        [[ "$gpu_count" =~ ^[1-9][0-9]*$ ]] || {
            echo "FASTGATK_SLURM_GPU_COUNT must be a positive integer" >&2
            exit 2
        }
        if [[ -n "$gpu_type" ]]; then
            sbatch_args+=(--gpus="${gpu_type}:${gpu_count}")
        else
            sbatch_args+=(--gpus="$gpu_count")
        fi
    fi
    [[ -n "${FASTGATK_SLURM_PARTITION:-}" ]] && sbatch_args+=(--partition="$FASTGATK_SLURM_PARTITION")
    [[ -n "${FASTGATK_SLURM_ACCOUNT:-}" ]] && sbatch_args+=(--account="$FASTGATK_SLURM_ACCOUNT")
    [[ -n "${FASTGATK_SLURM_TIME:-}" ]] && sbatch_args+=(--time="$FASTGATK_SLURM_TIME")
    wrap="SLURM_CPUS_PER_TASK=$(printf '%q' "$cpus")"
    [[ -n "$memory" ]] && wrap+=" SLURM_MEM_PER_NODE=$(printf '%q' "$memory")"
    [[ -n "$memory_per_cpu" ]] && wrap+=" SLURM_MEM_PER_CPU=$(printf '%q' "$memory_per_cpu")"
    [[ -n "$gres" ]] && wrap+=" SLURM_JOB_GRES=$(printf '%q' "$gres")"
    # SLURM normally creates CUDA_VISIBLE_DEVICES itself after allocating a
    # GPU.  Keep an explicit launcher override available for deterministic
    # site wrappers and local fake-sbatch tests, but never invent a device id
    # merely because --gpus was requested.
    cuda_devices=${FASTGATK_CUDA_VISIBLE_DEVICES:-}
    [[ -n "$cuda_devices" ]] && wrap+=" CUDA_VISIBLE_DEVICES=$(printf '%q' "$cuda_devices")"
    if [[ -n "${FASTGATK_SLURM_TMPDIR:-}" ]]; then
        wrap+=" SLURM_TMPDIR=$(printf '%q' "$FASTGATK_SLURM_TMPDIR")"
    fi
    wrap+=" $(printf '%q ' "$binary" "$@")"
    sbatch "${sbatch_args[@]}" --wrap="$wrap"
else
    export SLURM_JOB_ID=${SLURM_JOB_ID:-local-fastgatk-smoke}
    export SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-${FASTGATK_SLURM_CPUS_PER_TASK:-2}}
    export SLURM_MEM_PER_NODE=${SLURM_MEM_PER_NODE:-${FASTGATK_SLURM_MEMORY:-4096}}
    if [[ -n "${SLURM_MEM_PER_CPU:-}" ]]; then
        export SLURM_MEM_PER_CPU
    fi
    export SLURM_TMPDIR=${SLURM_TMPDIR:-${FASTGATK_SLURM_TMPDIR:-${TMPDIR:-/tmp}}}
    if [[ -n "${FASTGATK_SLURM_GRES:-}" && -z "${SLURM_JOB_GRES:-}" ]]; then
        export SLURM_JOB_GRES="$FASTGATK_SLURM_GRES"
    fi
    if [[ -n "${FASTGATK_CUDA_VISIBLE_DEVICES:-}" && -z "${CUDA_VISIBLE_DEVICES:-}" ]]; then
        export CUDA_VISIBLE_DEVICES="$FASTGATK_CUDA_VISIBLE_DEVICES"
    fi
    exec "$binary" "$@"
fi
