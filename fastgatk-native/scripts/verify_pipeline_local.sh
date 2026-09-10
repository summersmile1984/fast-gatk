#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/../.." && pwd)
build_dir="${FASTGATK_WORKFLOW_BUILD_DIR:-$root/fastgatk-native/build}"
binary="$build_dir/fastgatk-hc-smoke"
bam="$root/gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
run_dir=$(mktemp -d "${TMPDIR:-/tmp}/fastgatk-pipeline-smoke.XXXXXX")
trap 'rm -rf "$run_dir"' EXIT

[[ -x "$binary" ]] || { echo "missing $binary; run build_native.sh first" >&2; exit 2; }
[[ -f "$bam" ]] || { echo "missing BAM fixture: $bam" >&2; exit 2; }

compat="$root/fastgatk-native/scripts/gatk_compat.sh"
slurm="$root/fastgatk-native/workflow/slurm_smoke.sh"
scatter="$root/fastgatk-native/workflow/scatter_gather_smoke.sh"
gather_binary="$build_dir/fastgatk-gather-vcfs"
[[ -x "$gather_binary" ]] || { echo "missing $gather_binary; run build_native.sh first" >&2; exit 2; }
split_binary="$build_dir/fastgatk-split-intervals"
split_reference="$root/gatk-source/src/test/resources/hg19micro.fasta"
[[ -x "$split_binary" ]] || { echo "missing $split_binary; run build_native.sh first" >&2; exit 2; }
[[ -f "$split_reference" && -f "$split_reference.fai" ]] || {
    echo "missing SplitIntervals reference fixture: $split_reference" >&2
    exit 2
}

# Nextflow process-equivalent local execution. This is also the fallback when
# the nextflow executable is not installed in a developer/CI image.
FASTGATK_HC_BINARY="$binary" "$compat" --java-options "-Xmx1g" HaplotypeCaller \
    -I "$bam" -L 17:69000-69100 -O "$run_dir/nextflow-calls.vcf" \
    --output-manifest "$run_dir/nextflow-calls.vcf.manifest.json" --threads 2

# The same command through the SLURM resource adapter, locally simulated.
SLURM_CPUS_PER_TASK=2 "$slurm" "$binary" \
    -I "$bam" -L 17:69000-69100 -O "$run_dir/slurm-calls.vcf" \
    --output-manifest "$run_dir/slurm-calls.vcf.manifest.json" \
    --native-pair-hmm-threads 2

# Multi-interval scatter/gather through the same SLURM resource wrapper.  The
# second invocation is the restart contract: complete shard manifests are
# reused and only the final deterministic gather is repeated.
scatter_dir="$run_dir/scatter-gather"
FASTGATK_HC_BINARY="$binary" "$scatter" \
    -I "$bam" -L "17:69000-69050,17:69051-69100" \
    --binary "$binary" --gather-binary "$gather_binary" \
    --outdir "$scatter_dir" --threads 2 > "$run_dir/scatter-first.json"
FASTGATK_HC_BINARY="$binary" "$scatter" \
    -I "$bam" -L "17:69000-69050,17:69051-69100" \
    --binary "$binary" --gather-binary "$gather_binary" \
    --outdir "$scatter_dir" --threads 2 > "$run_dir/scatter-resume.json"

# SplitIntervals is the interval producer for the same scatter boundary.  Run
# it through the SLURM wrapper so the generated shard names/manifests are
# directly consumable by a Nextflow or sbatch process.
split_dir="$run_dir/split-intervals"
SLURM_CPUS_PER_TASK=2 "$slurm" "$split_binary" \
    -R "$split_reference" -O "$split_dir" --scatter-count 3 \
    --output-manifest "$run_dir/split-intervals.manifest.json" > "$run_dir/split-intervals.json"

python3 - "$run_dir" <<'PY'
import json
import pathlib
import sys

run_dir = pathlib.Path(sys.argv[1])
for stem in ("nextflow-calls", "slurm-calls"):
    vcf = run_dir / f"{stem}.vcf"
    manifest = run_dir / f"{stem}.vcf.manifest.json"
    text = vcf.read_text()
    assert "##fileformat=VCFv4.2" in text
    assert "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" in text
    data = json.loads(manifest.read_text())
    assert data["status"] == "smoke"
    assert data["primary_output_kind"] == "vcf"
    assert data["compatibility"]["gatk_parameter_aliases"] is True
    assert data["telemetry"]["reads"] == 7
scatter_manifest = json.loads((run_dir / "scatter-gather" / "scatter-gather.manifest.json").read_text())
assert scatter_manifest["status"] == "pass"
assert scatter_manifest["resume"] == {"resumed_shards": 2, "computed_shards": 0}
assert scatter_manifest["compatibility"]["nextflow_resume"] is True
assert (run_dir / "scatter-gather" / "gathered.vcf.gz.tbi").stat().st_size > 0
assert all(item["complete"] for item in scatter_manifest["shards"])
split_manifest = json.loads((run_dir / "split-intervals.manifest.json").read_text())
assert split_manifest["tool"] == "SplitIntervals"
assert split_manifest["output_shards"] == 3
assert split_manifest["kernel"]["policy"] == "RangePolicy"
assert split_manifest["kernel"]["execution_space"] in {"OpenMP", "Serial"}
assert all(pathlib.Path(item["path"]).is_file() and item["complete"]
           for item in split_manifest["outputs"])
print(json.dumps({"status": "pass", "mode": "local-nextflow-equivalent+local-slurm",
                  "split_intervals": {"shards": split_manifest["output_shards"],
                                      "kernel_policy": split_manifest["kernel"]["policy"]}}))
PY

# If Nextflow is available (or explicitly supplied), execute the actual DSL2
# process as a second check.  A bundled JDK can be selected without changing
# the host PATH, which makes this contract reproducible in CI images.
nextflow_bin="${FASTGATK_NEXTFLOW_BIN:-$(command -v nextflow || true)}"
nextflow_java="${FASTGATK_NEXTFLOW_JAVA_HOME:-$root/third_party/jdk17}"
nextflow_cmd=()
if [[ -n "$nextflow_bin" && -x "$nextflow_bin" ]]; then
    nextflow_cmd=("$nextflow_bin")
else
    owner_root=$(cd "$root/../.." && pwd)
    nextflow_jar="${FASTGATK_NEXTFLOW_JAR:-$(find "$owner_root/.nextflow/framework" -type f -name 'nextflow-*-one.jar' 2>/dev/null | sort | tail -1)}"
    if [[ -n "$nextflow_jar" && -f "$nextflow_jar" ]]; then
        java_bin="$nextflow_java/bin/java"
        [[ -x "$java_bin" ]] || java_bin="$(command -v java || true)"
        [[ -n "$java_bin" ]] || { echo "Nextflow jar found but no Java runtime" >&2; exit 2; }
        nextflow_cmd=("$java_bin" -jar "$nextflow_jar")
    fi
fi
if [[ "${#nextflow_cmd[@]}" -gt 0 ]]; then
    JAVA_HOME="$nextflow_java" "${nextflow_cmd[@]}" run "$root/fastgatk-native/workflow/nextflow_smoke.nf" \
        --input "$bam" --region 17:69000-69100 --binary "$binary" \
        --outdir "$run_dir/nextflow-real" --threads 2 -resume
    [[ -f "$run_dir/nextflow-real/calls.vcf" ]]
    [[ -f "$run_dir/nextflow-real/calls.vcf.manifest.json" ]]
    JAVA_HOME="$nextflow_java" "${nextflow_cmd[@]}" run "$root/fastgatk-native/workflow/nextflow_scatter_gather.nf" \
        --input "$bam" --intervals "17:69000-69050,17:69051-69100" \
        --binary "$binary" --gather_binary "$gather_binary" \
        --outdir "$run_dir/nextflow-scatter-real" --threads 2 -resume
    [[ -f "$run_dir/nextflow-scatter-real/gathered.vcf.gz" ]]
    [[ -f "$run_dir/nextflow-scatter-real/gathered.vcf.gz.tbi" ]]
    # Generic direct-replacement process: the same DSL2 boundary invokes the
    # drop-in `gatk` launcher for a non-HC registry tool and publishes both
    # primary output and OutputManifest.  This catches regressions where a
    # workflow accidentally bypasses dispatcher validation/fallback.
    JAVA_HOME="$nextflow_java" "${nextflow_cmd[@]}" run "$root/fastgatk-native/workflow/nextflow_gatk_compat.nf" \
        --tool CountReads --input "$bam" --output custom-counts.txt \
        --launcher "$root/fastgatk-native/dispatcher/gatk" \
        --extra_args '--batch-records 37 --threads 2' \
        --outdir "$run_dir/nextflow-generic-count" -resume
    [[ "$(tr -d '[:space:]' < "$run_dir/nextflow-generic-count/custom-counts.txt")" == "493" ]]
    [[ -s "$run_dir/nextflow-generic-count/custom-counts.txt.manifest.json" ]]
    echo '{"nextflow": "executed"}'
else
    echo '{"nextflow": "not-installed; local process simulation passed"}'
fi
