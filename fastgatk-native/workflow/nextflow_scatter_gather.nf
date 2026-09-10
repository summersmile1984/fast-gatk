nextflow.enable.dsl=2

/*
 * GATK-compatible DSL2 scatter/gather contract.  Nextflow owns task caching
 * and `-resume`; each native task also writes an OutputManifest so the same
 * shard boundary can be inspected or retried by a SLURM-only runner.
 *
 * Example:
 *   nextflow run nextflow_scatter_gather.nf -resume \
 *     --input reads.bam --intervals '17:69000-69050,17:69051-69100' \
 *     --binary ../build/fastgatk-hc-smoke \
 *     --gather_binary ../build/fastgatk-gather-vcfs \
 *     --max_retries 1 --error_strategy retry
 */

params.input = null
params.reference = null
params.intervals = null
params.binary = "${projectDir}/../build/fastgatk-hc-smoke"
params.gather_binary = "${projectDir}/../build/fastgatk-gather-vcfs"
params.slurm_wrapper = "${projectDir}/slurm_smoke.sh"
params.outdir = "nextflow-scatter-gather-results"
params.threads = 2
// Nextflow owns process retry policy.  The default is fail-fast; production
// profiles can set --error_strategy retry and --max_retries N.  The native
// wrapper itself remains single-attempt inside an allocation.
params.max_retries = 0
params.error_strategy = 'terminate'

process FASTGATK_SCATTER {
    tag "${interval}"
    cpus params.threads as int
    errorStrategy { params.error_strategy ?: 'terminate' }
    maxRetries params.max_retries as int
    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    tuple val(interval), path(reads)

    output:
    tuple val(interval), path("calls-*.vcf"), path("calls-*.vcf.manifest.json"), emit: shards

    script:
    def safe = interval.replaceAll(/[^A-Za-z0-9_.-]/, "_")
    def output_name = "calls-${safe}.vcf"
    def manifest_name = "${output_name}.manifest.json"
    def ref_arg = params.reference ? "-R '${params.reference}'" : ""
    """
    '${params.slurm_wrapper}' '${params.binary}' -I '${reads}' ${ref_arg} \
      -L '${interval}' -O '${output_name}' --output-manifest '${manifest_name}' \
      --threads ${task.cpus}
    """
}

process FASTGATK_GATHER {
    tag "gathered-vcfs"
    cpus 1
    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    path shard_vcfs

    output:
    path "gathered.vcf.gz"
    path "gathered.vcf.gz.tbi"
    path "gathered.vcf.gz.manifest.json"

    script:
    def shard_args = shard_vcfs.collect { "-I '${it}'" }.join(' ')
    """
    '${params.gather_binary}' ${shard_args} -O gathered.vcf.gz \
      --output-manifest gathered.vcf.gz.manifest.json
    """
}

workflow {
    if (!params.input) error "--input is required"
    if (!params.intervals) error "--intervals is required (comma-separated)"
    def interval_list = params.intervals instanceof List ? params.intervals :
        params.intervals.toString().split(',').collect { it.trim() }.findAll { it }
    if (!interval_list) error "--intervals must contain at least one interval"

    def reads = file(params.input, checkIfExists: true)
    scatter_input = Channel.fromList(interval_list).map { interval -> tuple(interval, reads) }
    scattered = FASTGATK_SCATTER(scatter_input)
    shard_vcfs = scattered.shards.map { item -> item[1] }.collect()
    FASTGATK_GATHER(shard_vcfs)
}
