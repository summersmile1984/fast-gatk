nextflow.enable.dsl=2

/*
 * Local Nextflow contract test. It exercises process staging, -resume-compatible
 * output names, and GATK-style -I/-R/-L/-O arguments. The default executable is
 * the next-stage native assembly/likelihood/genotyping prototype; callers can
 * pass fastgatk-hc-smoke explicitly when only file-boundary validation is desired.
 */
params.input = null
params.reference = null
params.region = null
params.binary = "${projectDir}/../build/fastgatk-hc-call"
params.outdir = "nextflow-smoke-results"
params.threads = 2

process FASTGATK_HAPLOTYPECALLER {
    tag "${reads.baseName}"
    cpus params.threads as int
    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    path reads

    output:
    path "calls.vcf"
    path "calls.vcf.manifest.json"

    script:
    def ref_arg = params.reference ? "-R ${params.reference}" : ""
    def interval_arg = params.region ? "-L ${params.region}" : ""
    """
    ${params.binary} -I ${reads} ${ref_arg} ${interval_arg} -O calls.vcf \
      --output-manifest calls.vcf.manifest.json --threads ${task.cpus}
    """
}

workflow {
    if (!params.input) error "--input is required"
    reads_ch = Channel.fromPath(params.input, checkIfExists: true)
    FASTGATK_HAPLOTYPECALLER(reads_ch)
}
