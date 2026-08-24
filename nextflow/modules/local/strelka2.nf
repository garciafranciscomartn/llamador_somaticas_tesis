// ============================================================
// modules/strelka2.nf
// Strelka2 somatic calling processes (pairwise).
// Default: STRELKA2_PAIR_DIRECT (no Manta).
// Set params.strelka2_run_manta = true to use MANTA_PAIR first.
// ============================================================

process PREP_STRELKA_BED {
    tag "prep_bed"
    label 'process_low'

    input:
    path bed

    output:
    tuple path("${bed}.gz"), path("${bed}.gz.tbi"), emit: bed_gz

    script:
    """
    bgzip -c ${bed} > ${bed}.gz
    tabix -p bed ${bed}.gz
    """
}


process MANTA_PAIR {
    tag "${pair_name}"
    label 'strelka2'

    publishDir "${params.outdir}/variant_calling/strelka2/${pair_name}/manta/results/variants",
        mode: 'copy', pattern: 'candidateSmallIndels.vcf.gz*'

    input:
    tuple val(pair_name), val(tumor_id), path(tumor_bam), path(tumor_bai),
          val(normal_id), path(normal_bam), path(normal_bai)
    path genome
    tuple path(bed_gz), path(bed_tbi)

    output:
    tuple val(pair_name),
          path("candidateSmallIndels.vcf.gz"),
          path("candidateSmallIndels.vcf.gz.tbi"), emit: indels
    path "logs/${pair_name}.manta.log",             emit: log

    script:
    def manta_bin = params.manta_bin
    """
    mkdir -p logs
    rm -rf manta_run

    ${manta_bin}/configManta.py \\
        --normalBam ${normal_bam} \\
        --tumorBam  ${tumor_bam} \\
        --referenceFasta ${genome} \\
        --callRegions ${bed_gz} \\
        --exome \\
        --runDir manta_run \\
        > logs/${pair_name}.manta.log 2>&1

    manta_run/runWorkflow.py \\
        -m local -j ${task.cpus} \\
        >> logs/${pair_name}.manta.log 2>&1

    cp manta_run/results/variants/candidateSmallIndels.vcf.gz .
    cp manta_run/results/variants/candidateSmallIndels.vcf.gz.tbi .
    """
}


process STRELKA2_PAIR {
    // Strelka2 WITH Manta indel candidates (params.strelka2_run_manta = true)
    tag "${pair_name}"
    label 'strelka2'

    input:
    tuple val(pair_name), val(tumor_id), path(tumor_bam), path(tumor_bai),
          val(normal_id), path(normal_bam), path(normal_bai),
          path(indels_vcf), path(indels_tbi)
    path genome
    tuple path(bed_gz), path(bed_tbi)

    output:
    tuple val(pair_name),
          path("strelka_run/results/variants/somatic.snvs.vcf.gz"),
          path("strelka_run/results/variants/somatic.snvs.vcf.gz.tbi"),
          path("strelka_run/results/variants/somatic.indels.vcf.gz"),
          path("strelka_run/results/variants/somatic.indels.vcf.gz.tbi"), emit: results
    path "logs/${pair_name}.strelka2.log",                                 emit: log

    script:
    def strelka_bin = params.strelka2_bin
    """
    mkdir -p logs
    rm -rf strelka_run

    ${strelka_bin}/configureStrelkaSomaticWorkflow.py \\
        --normalBam ${normal_bam} \\
        --tumorBam  ${tumor_bam} \\
        --referenceFasta ${genome} \\
        --indelCandidates ${indels_vcf} \\
        --callRegions ${bed_gz} \\
        --exome \\
        --runDir strelka_run \\
        > logs/${pair_name}.strelka2.log 2>&1

    strelka_run/runWorkflow.py \\
        -m local -j ${task.cpus} \\
        >> logs/${pair_name}.strelka2.log 2>&1
    """
}


process STRELKA2_PAIR_DIRECT {
    // Strelka2 WITHOUT Manta (default, as used in the original analysis)
    tag "${pair_name}"
    label 'strelka2'

    input:
    tuple val(pair_name), val(tumor_id), path(tumor_bam), path(tumor_bai),
          val(normal_id), path(normal_bam), path(normal_bai)
    path genome
    tuple path(bed_gz), path(bed_tbi)

    output:
    tuple val(pair_name),
          path("strelka_run/results/variants/somatic.snvs.vcf.gz"),
          path("strelka_run/results/variants/somatic.snvs.vcf.gz.tbi"),
          path("strelka_run/results/variants/somatic.indels.vcf.gz"),
          path("strelka_run/results/variants/somatic.indels.vcf.gz.tbi"), emit: results
    path "logs/${pair_name}.strelka2.log",                                 emit: log

    script:
    def strelka_bin = params.strelka2_bin
    """
    mkdir -p logs
    rm -rf strelka_run

    ${strelka_bin}/configureStrelkaSomaticWorkflow.py \\
        --normalBam ${normal_bam} \\
        --tumorBam  ${tumor_bam} \\
        --referenceFasta ${genome} \\
        --callRegions ${bed_gz} \\
        --exome \\
        --runDir strelka_run \\
        > logs/${pair_name}.strelka2.log 2>&1

    strelka_run/runWorkflow.py \\
        -m local -j ${task.cpus} \\
        >> logs/${pair_name}.strelka2.log 2>&1
    """
}


process STRELKA2_MERGE_PAIR {
    tag "${pair_name}"
    label 'process_low'

    publishDir "${params.outdir}/variant_calling/strelka2/${pair_name}",
        mode: 'copy', pattern: '*.{vcf.gz,vcf.gz.tbi}'

    input:
    tuple val(pair_name),
          path(snvs), path(snvs_tbi),
          path(indels), path(indels_tbi)

    output:
    tuple val(pair_name),
          path("${pair_name}.strelka2.pass.vcf.gz"),
          path("${pair_name}.strelka2.pass.vcf.gz.tbi"), emit: vcf
    path "logs/${pair_name}.strelka2_merge.log",          emit: log

    script:
    """
    mkdir -p logs
    (
        bcftools concat -a ${snvs} ${indels} 2>/dev/null \\
            | bcftools sort 2>/dev/null \\
            | bcftools view -f .,PASS -Oz -o ${pair_name}.strelka2.pass.vcf.gz
        tabix -p vcf ${pair_name}.strelka2.pass.vcf.gz
    ) > logs/${pair_name}.strelka2_merge.log 2>&1
    """
}
