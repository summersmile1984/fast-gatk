nextflow.enable.dsl=2

/*
 * Generic DSL2 process for replacing `gatk TOOL ...` with the fastgatk
 * launcher.  The process deliberately keeps the tool argv at the workflow
 * boundary: Nextflow still owns staging, cache/resume and the executor
 * allocation, while the dispatcher owns registry validation, native-vs-Java
 * fallback and OutputManifest publication.
 *
 * Example (single input):
 *   nextflow run nextflow_gatk_compat.nf -resume \
 *     --tool CountReads --input reads.bam --output result.txt \
 *     --launcher ../dispatcher/gatk --extra_args '--batch-records 4096'
 */

params.tool = 'HaplotypeCaller'
params.input = null
params.reference = null
params.intervals = null
params.output = 'result.out'
params.launcher = "${projectDir}/../dispatcher/gatk"
params.extra_args = ''
params.outdir = 'nextflow-gatk-compat-results'
params.threads = 2

process FASTGATK_COMPAT_TOOL {
    tag "${params.tool}:${reads.simpleName}"
    cpus params.threads as int
    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    path reads

    output:
    // Keep the caller's output stem visible to Nextflow.  The previous
    // hard-coded `result*` declaration happened to work for the smoke
    // fixture, but made a normal GATK process fail after the native command
    // succeeded whenever `--output` used another filename.  GATK tools also
    // emit sidecars (manifest/index/metrics), so the stem wildcard preserves
    // those artifacts for publishDir and -resume without changing names.
    path("${params.output}*"), emit: outputs

    script:
    def ref_arg = params.reference ? "-R '${params.reference}'" : ''
    def interval_values = params.intervals instanceof List ? params.intervals :
        (params.intervals ? params.intervals.toString().split(',').collect { it.trim() }.findAll { it } : [])
    def interval_arg = interval_values.collect { "-L '${it}'" }.join(' ')
    def extra = params.extra_args ?: ''
    """
    '${params.launcher}' '${params.tool}' -I '${reads}' ${ref_arg} ${interval_arg} \
      -O '${params.output}' --output-manifest '${params.output}.manifest.json' \
      --threads ${task.cpus} ${extra}
    """
}

workflow {
    if (!params.input) error '--input is required'
    if (!params.tool) error '--tool is required'
    if (!params.output) error '--output is required'
    reads_ch = Channel.fromPath(params.input, checkIfExists: true)
    FASTGATK_COMPAT_TOOL(reads_ch)
}
