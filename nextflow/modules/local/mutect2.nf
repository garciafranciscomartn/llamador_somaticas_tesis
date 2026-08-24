// ============================================================
// modules/mutect2.nf
// GATK Mutect2 processes for pairwise somatic calling.
// ============================================================

process MUTECT2_PAIR {
    tag "${pair_name}"
    label 'mutect2'

    publishDir "${params.outdir}/variant_calling/mutect2/${pair_name}",
        mode: 'copy', pattern: '*.{vcf.gz,vcf.gz.tbi,tar.gz,stats}'

    input:
    tuple val(pair_name), val(tumor_id), path(tumor_bam), path(tumor_bai),
          val(normal_id), path(normal_bam), path(normal_bai)
    path genome

    output:
    tuple val(pair_name), path("${pair_name}.mutect2.vcf.gz"),
          path("${pair_name}.mutect2.vcf.gz.tbi"),
          path("${pair_name}.mutect2.vcf.gz.stats"),
          path("${pair_name}.f1r2.tar.gz"),      emit: vcf
    path "logs/${pair_name}.mutect2.log",        emit: log

    script:
    def pon_arg      = params.panel_of_normals  ? "--panel-of-normals ${params.panel_of_normals}"  : ''
    def germline_arg = params.germline_resource ? "--germline-resource ${params.germline_resource}" : ''
    def interv_arg   = params.intervals         ? "--intervals ${params.intervals}"                 : ''
    def java_mem     = task.memory.toMega().intValue()
    """
    mkdir -p logs
    TNAME=\$(samtools view -H ${tumor_bam} \\
        | awk -F'\\t' '/^@RG/{ for(i=1;i<=NF;i++) if(\$i~/^SM:/){ print substr(\$i,4); exit } }')
    [ -z "\$TNAME" ] && TNAME="${tumor_id}"

    NNAME=\$(samtools view -H ${normal_bam} \\
        | awk -F'\\t' '/^@RG/{ for(i=1;i<=NF;i++) if(\$i~/^SM:/){ print substr(\$i,4); exit } }')
    [ -z "\$NNAME" ] && NNAME="${normal_id}"

    gatk --java-options '-Xmx${java_mem}m' Mutect2 \\
        --input ${tumor_bam} \\
        --input ${normal_bam} \\
        --normal-sample "\$NNAME" \\
        --tumor-sample  "\$TNAME" \\
        --reference ${genome} \\
        ${pon_arg} \\
        ${germline_arg} \\
        ${interv_arg} \\
        --f1r2-tar-gz ${pair_name}.f1r2.tar.gz \\
        --output ${pair_name}.mutect2.vcf.gz \\
        ${params.mutect2_extra} \\
        > logs/${pair_name}.mutect2.log 2>&1
    """
}


process LEARN_ORIENTATION_PAIR {
    tag "${pair_name}"
    label 'gatk_default'

    publishDir "${params.outdir}/variant_calling/mutect2/${pair_name}",
        mode: 'copy', pattern: '*.tar.gz'

    input:
    tuple val(pair_name), path(f1r2)

    output:
    tuple val(pair_name), path("${pair_name}.read_orientation_model.tar.gz"), emit: model
    path "logs/${pair_name}.orientation.log",                                  emit: log

    script:
    def java_mem = task.memory.toMega().intValue()
    """
    mkdir -p logs
    gatk --java-options '-Xmx${java_mem}m' LearnReadOrientationModel \\
        --input ${f1r2} \\
        --output ${pair_name}.read_orientation_model.tar.gz \\
        > logs/${pair_name}.orientation.log 2>&1
    """
}


process PILEUP_TUMOR_PAIR {
    tag "${pair_name}"
    label 'pileup'

    publishDir "${params.outdir}/variant_calling/mutect2/${pair_name}",
        mode: 'copy', pattern: '*.table'

    input:
    tuple val(pair_name), val(tumor_id), path(tumor_bam), path(tumor_bai),
          val(normal_id), path(normal_bam), path(normal_bai)
    path germline_resource
    path germline_tbi
    path intervals

    output:
    tuple val(pair_name), path("${pair_name}.tumor.pileup_summaries.table"), emit: table
    path "logs/${pair_name}.pileup.tumor.log",                               emit: log

    script:
    def java_mem = task.memory.toMega().intValue()
    """
    mkdir -p logs
    gatk --java-options '-Xmx${java_mem}m' GetPileupSummaries \\
        --input ${tumor_bam} \\
        --variant ${germline_resource} \\
        --intervals ${intervals} \\
        --output ${pair_name}.tumor.pileup_summaries.table \\
        > logs/${pair_name}.pileup.tumor.log 2>&1
    """
}


process PILEUP_NORMAL_PAIR {
    tag "${pair_name}"
    label 'pileup'

    publishDir "${params.outdir}/variant_calling/mutect2/${pair_name}",
        mode: 'copy', pattern: '*.table'

    input:
    tuple val(pair_name), val(tumor_id), path(tumor_bam), path(tumor_bai),
          val(normal_id), path(normal_bam), path(normal_bai)
    path germline_resource
    path germline_tbi
    path intervals

    output:
    tuple val(pair_name), path("${pair_name}.normal.pileup_summaries.table"), emit: table
    path "logs/${pair_name}.pileup.normal.log",                               emit: log

    script:
    def java_mem = task.memory.toMega().intValue()
    """
    mkdir -p logs
    gatk --java-options '-Xmx${java_mem}m' GetPileupSummaries \\
        --input ${normal_bam} \\
        --variant ${germline_resource} \\
        --intervals ${intervals} \\
        --output ${pair_name}.normal.pileup_summaries.table \\
        > logs/${pair_name}.pileup.normal.log 2>&1
    """
}


process CALCULATE_CONTAMINATION_PAIR {
    tag "${pair_name}"
    label 'gatk_default'

    publishDir "${params.outdir}/variant_calling/mutect2/${pair_name}",
        mode: 'copy', pattern: '*.table'

    input:
    tuple val(pair_name), path(tumor_table), path(normal_table)

    output:
    tuple val(pair_name),
          path("${pair_name}.contamination.table"),
          path("${pair_name}.segments.table"),       emit: tables
    path "logs/${pair_name}.contamination.log",      emit: log

    script:
    def java_mem = task.memory.toMega().intValue()
    """
    mkdir -p logs
    gatk --java-options '-Xmx${java_mem}m' CalculateContamination \\
        --input ${tumor_table} \\
        --matched-normal ${normal_table} \\
        --tumor-segmentation ${pair_name}.segments.table \\
        --output ${pair_name}.contamination.table \\
        > logs/${pair_name}.contamination.log 2>&1
    """
}


process FILTER_MUTECT_PAIR {
    tag "${pair_name}"
    label 'gatk_default'

    publishDir "${params.outdir}/variant_calling/mutect2/${pair_name}",
        mode: 'copy', pattern: '*.{vcf.gz,vcf.gz.tbi}'

    input:
    tuple val(pair_name),
          path(vcf), path(tbi), path(stats),
          path(contamination), path(segments),
          path(orientation_model)
    path genome

    output:
    tuple val(pair_name),
          path("${pair_name}.mutect2.filtered.vcf.gz"),
          path("${pair_name}.mutect2.filtered.vcf.gz.tbi"), emit: vcf
    path "logs/${pair_name}.filter.log",                    emit: log

    script:
    def java_mem = task.memory.toMega().intValue()
    """
    mkdir -p logs
    gatk --java-options '-Xmx${java_mem}m' FilterMutectCalls \\
        --variant ${vcf} \\
        --reference ${genome} \\
        --contamination-table ${contamination} \\
        --tumor-segmentation ${segments} \\
        --orientation-bias-artifact-priors ${orientation_model} \\
        --stats ${stats} \\
        --output ${pair_name}.mutect2.filtered.vcf.gz \\
        > logs/${pair_name}.filter.log 2>&1
    """
}


process SELECT_PASS_PAIR {
    tag "${pair_name}"
    label 'gatk_default'

    publishDir "${params.outdir}/variant_calling/mutect2/${pair_name}",
        mode: 'copy', pattern: '*.pass.*'

    input:
    tuple val(pair_name), path(vcf), path(tbi)
    path genome

    output:
    tuple val(pair_name),
          path("${pair_name}.mutect2.pass.vcf.gz"),
          path("${pair_name}.mutect2.pass.vcf.gz.tbi"), emit: vcf
    path "logs/${pair_name}.select_pass.log",            emit: log

    script:
    def java_mem = task.memory.toMega().intValue()
    """
    mkdir -p logs
    gatk --java-options '-Xmx${java_mem}m' SelectVariants \\
        --variant ${vcf} \\
        --reference ${genome} \\
        --exclude-filtered \\
        --output ${pair_name}.mutect2.pass.vcf.gz \\
        > logs/${pair_name}.select_pass.log 2>&1
    """
}
