// ============================================================
// modules/varscan2.nf
// VarScan2 somatic calling processes (pairwise).
// ============================================================

process VARSCAN_MPILEUP_PAIR {
    tag "${pair_name}"
    label 'varscan2'

    input:
    tuple val(pair_name), val(tumor_id), path(tumor_bam), path(tumor_bai),
          val(normal_id), path(normal_bam), path(normal_bai)
    path genome
    path bed

    output:
    tuple val(pair_name),
          path("${pair_name}.normal.pileup"),
          path("${pair_name}.tumor.pileup"), emit: pileups
    path "logs/${pair_name}.mpileup.log",    emit: log

    script:
    """
    mkdir -p logs
    samtools mpileup -q 1 -Q 20 -f ${genome} \\
        -l ${bed} \\
        ${normal_bam} > ${pair_name}.normal.pileup 2> logs/${pair_name}.mpileup.log

    samtools mpileup -q 1 -Q 20 -f ${genome} \\
        -l ${bed} \\
        ${tumor_bam} > ${pair_name}.tumor.pileup 2>> logs/${pair_name}.mpileup.log
    """
}


process VARSCAN_SOMATIC_PAIR {
    tag "${pair_name}"
    label 'varscan2'

    input:
    tuple val(pair_name), path(normal_pileup), path(tumor_pileup)

    output:
    tuple val(pair_name),
          path("${pair_name}.varscan2.snp.vcf"),
          path("${pair_name}.varscan2.indel.vcf"), emit: vcfs
    path "logs/${pair_name}.varscan_somatic.log",  emit: log

    script:
    """
    mkdir -p logs
    varscan somatic \\
        ${normal_pileup} \\
        ${tumor_pileup} \\
        ${pair_name}.varscan2 \\
        --min-coverage ${params.varscan_min_coverage} \\
        --min-var-freq ${params.varscan_min_var_freq} \\
        --p-value ${params.varscan_p_value} \\
        --somatic-p-value ${params.varscan_somatic_p_value} \\
        --strand-filter 1 \\
        --output-vcf 1 \\
        > logs/${pair_name}.varscan_somatic.log 2>&1
    """
}


process VARSCAN_MERGE_PAIR {
    // Merge SNP+indel, keep SS=2 (somatic only) to exclude germline and LOH.
    tag "${pair_name}"
    label 'process_low'

    publishDir "${params.outdir}/variant_calling/varscan2/${pair_name}",
        mode: 'copy', pattern: '*.{vcf.gz,vcf.gz.tbi}'

    input:
    tuple val(pair_name), path(snp_vcf), path(indel_vcf)

    output:
    tuple val(pair_name),
          path("${pair_name}.varscan2.pass.vcf.gz"),
          path("${pair_name}.varscan2.pass.vcf.gz.tbi"), emit: vcf
    path "logs/${pair_name}.varscan_merge.log",           emit: log

    script:
    """
    mkdir -p logs
    (
        bgzip -c ${snp_vcf}   > ${snp_vcf}.gz
        tabix -p vcf ${snp_vcf}.gz
        bgzip -c ${indel_vcf} > ${indel_vcf}.gz
        tabix -p vcf ${indel_vcf}.gz

        bcftools concat -a ${snp_vcf}.gz ${indel_vcf}.gz \\
            | bcftools view -i 'INFO/SS=2' \\
            | bcftools sort -Oz -o ${pair_name}.varscan2.pass.vcf.gz
        tabix -p vcf ${pair_name}.varscan2.pass.vcf.gz

        rm -f ${snp_vcf}.gz ${snp_vcf}.gz.tbi \\
              ${indel_vcf}.gz ${indel_vcf}.gz.tbi
    ) > logs/${pair_name}.varscan_merge.log 2>&1
    """
}
