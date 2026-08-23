// ============================================================
// modules/vardict.nf
// VarDict somatic calling process (pairwise).
// ============================================================

process VARDICT_PAIR {
    tag "${pair_name}"
    label 'vardict'

    publishDir "${params.outdir}/variant_calling/vardict/${pair_name}",
        mode: 'copy', pattern: '*.{vcf.gz,vcf.gz.tbi}'

    input:
    tuple val(pair_name), val(tumor_id), path(tumor_bam), path(tumor_bai),
          val(normal_id), path(normal_bam), path(normal_bai)
    path genome
    path bed

    output:
    tuple val(pair_name),
          path("${pair_name}.vardict.pass.vcf.gz"),
          path("${pair_name}.vardict.pass.vcf.gz.tbi"), emit: vcf
    path "logs/${pair_name}.vardict.log",               emit: log

    script:
    def min_af = params.vardict_min_af
    """
    mkdir -p logs
    (
        vardict-java \\
            -G ${genome} \\
            -f ${min_af} \\
            -N ${tumor_id} \\
            -b "${tumor_bam}|${normal_bam}" \\
            -c 1 -S 2 -E 3 \\
            --nosv \\
            -th ${task.cpus} \\
            ${bed} \\
        | testsomatic.R \\
        | var2vcf_paired.pl \\
            -N "${tumor_id}|${normal_id}" \\
            -f ${min_af} \\
        | bcftools view \\
            -i 'INFO/STATUS="StrongSomatic" || INFO/STATUS="LikelySomatic"' \\
        | bcftools view -f .,PASS -Oz -o ${pair_name}.vardict.pass.vcf.gz

        tabix -p vcf ${pair_name}.vardict.pass.vcf.gz
    ) > logs/${pair_name}.vardict.log 2>&1
    """
}
