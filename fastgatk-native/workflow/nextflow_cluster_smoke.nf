nextflow.enable.dsl=2

/*
 * M0 cluster smoke.  Eight scatter shards + one gather through the real
 * SLURM allocation via the `production` profile in nextflow.config.  This is
 * the only workflow that exercises the cluster-side resource contract
 * end-to-end:
 *
 *   1. Each FASTGATK_SCATTER task runs the cluster resource probe first so
 *      a misconfigured allocation (missing SLURM_JOB_ID, cgroup pressure,
 *      budget overflow) fails the workflow before any variant caller runs.
 *   2. The native HC binary is then invoked through slurm_smoke.sh, which
 *      detects the existing allocation and never submits a nested sbatch.
 *   3. The gather task pulls the eight VCFs and writes a final manifest
 *      whose `workflow_signature` is pinned in `cluster_evidence.json` so
 *      the dispatcher / runbook can prove provenance.
 *
 * Defaults to the chr17 69k–70k fixture already pinned in the tree.  Pass
 * `--input`/`--intervals`/`--reference`/`--output` to point at a different
 * shard boundary.
 *
 * Invocation:
 *   nextflow run nextflow_cluster_smoke.nf -profile production \
 *       --container fastgatk-native:latest --cluster_account biology
 */

params.input = "${projectDir}/../../gatk-source/src/test/resources/NA12878.chr17_69k_70k.dictFix.bam"
params.reference = "${projectDir}/../../gatk-source/src/test/resources/human_g1k_v37.chr17_1Mb.fasta"
params.intervals = "17:69000-69050,17:69051-69100,17:69101-69150,17:69151-69200,17:69201-69250,17:69251-69300,17:69301-69350,17:69351-69400"
params.binary = "${projectDir}/../build/fastgatk-hc-smoke"
params.gather_binary = "${projectDir}/../build/fastgatk-gather-vcfs"
params.probe_binary = "${projectDir}/../../fastgatk-runtime/build/fastgatk-cluster-resource-probe"
params.slurm_wrapper = "${projectDir}/slurm_smoke.sh"
params.outdir = "nextflow-cluster-smoke-results"
params.threads = 2
params.evidence_json = "cluster-evidence.json"

process FASTGATK_CLUSTER_PROBE {
    tag "allocation-${SLURM_JOB_ID ?: 'none'}"
    cpus 1
    memory '256 MB'
    publishDir params.outdir, mode: 'copy', overwrite: true

    output:
    path "probe.evidence.json"

    script:
    def evidence = "${params.outdir}/probe-${SLURM_JOB_ID ?: 'unknown'}.evidence.json"
    """
    mkdir -p ${params.outdir}
    FASTGATK_REAL_CLUSTER=1 FASTGATK_CLUSTER_EVIDENCE_JSON='${evidence}' \
        '${params.probe_binary}'
    test -s '${evidence}'
    cp '${evidence}' probe.evidence.json
    """
}

process FASTGATK_SCATTER {
    tag "${interval}"
    cpus params.threads as int
    errorStrategy 'retry'
    maxRetries 2
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
    FASTGATK_REAL_CLUSTER=1 '${params.slurm_wrapper}' '${params.binary}' \
        -I '${reads}' ${ref_arg} -L '${interval}' \
        -O '${output_name}' --output-manifest '${manifest_name}' \
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
    FASTGATK_REAL_CLUSTER=1 '${params.slurm_wrapper}' '${params.gather_binary}' ${shard_args} \
        -O gathered.vcf.gz --output-manifest gathered.vcf.gz.manifest.json
    """
}

workflow {
    FASTGATK_CLUSTER_PROBE()

    def reads = file(params.input, checkIfExists: true)
    def interval_list = params.intervals instanceof List ? params.intervals :
        params.intervals.toString().split(',').collect { it.trim() }.findAll { it }
    if (!interval_list) error "--intervals must contain at least one interval"
    if (interval_list.size() != 8) {
        log.warn "M0 smoke expects 8 intervals; got ${interval_list.size()}"
    }
    scatter_input = Channel.fromList(interval_list).map { interval -> tuple(interval, reads) }
    scattered = FASTGATK_SCATTER(scatter_input)
    shard_vcfs = scattered.shards.map { item -> item[1] }.collect()
    FASTGATK_GATHER(shard_vcfs)
}
