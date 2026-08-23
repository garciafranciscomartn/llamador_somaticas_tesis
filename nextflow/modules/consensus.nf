// ============================================================
// modules/consensus.nf
// Normalize caller VCFs then find consensus with bcftools isec.
// ============================================================

process NORM_VCF {
    tag "${id}:${caller}"
    label 'process_low'

    publishDir "${params.outdir}/variant_calling/norm/${id}",
        mode: 'copy', pattern: '*.{vcf.gz,vcf.gz.tbi}'

    input:
    tuple val(id), val(caller), path(vcf), path(tbi)
    path genome

    output:
    tuple val(id), val(caller),
          path("${id}.${caller}.norm.vcf.gz"),
          path("${id}.${caller}.norm.vcf.gz.tbi"), emit: vcf
    path "logs/${id}.${caller}.norm.log",           emit: log

    script:
    """
    mkdir -p logs
    (
        bcftools norm -m -both -f ${genome} ${vcf} \\
            | bcftools sort -Oz -o ${id}.${caller}.norm.vcf.gz
        tabix -p vcf ${id}.${caller}.norm.vcf.gz
    ) > logs/${id}.${caller}.norm.log 2>&1
    """
}


process CONSENSUS_ISEC {
    /*
     * bcftools isec -n+min_callers to find variants agreed on by
     * at least consensus_min_callers of the enabled callers.
     */
    tag "${id}"
    label 'process_medium'

    publishDir "${params.outdir}/variant_calling/consensus/${id}",
        mode: 'copy', pattern: '*.{vcf.gz,vcf.gz.tbi,sites.txt}'

    input:
    tuple val(id), path(vcfs), path(tbis)

    output:
    tuple val(id),
          path("${id}.consensus.vcf.gz"),
          path("${id}.consensus.vcf.gz.tbi"), emit: vcf
    path "${id}.sites.txt",                   emit: sites
    path "logs/${id}.consensus.log",           emit: log

    script:
    def min_callers = params.consensus_min_callers
    def n_callers   = vcfs instanceof List ? vcfs.size() : 1
    """
    mkdir -p logs isec_tmp
    (
        bcftools isec -n+${min_callers} \\
            ${vcfs.join(' ')} \\
            -p isec_tmp

        cp isec_tmp/sites.txt ${id}.sites.txt

        vcf_files=""
        for i in \$(seq 0 \$(( ${n_callers} - 1 ))); do
            f=\$(printf "isec_tmp/%04d.vcf" "\$i")
            g=\$(printf "isec_tmp/%04d.sites.vcf.gz" "\$i")
            if [ -f "\$f" ]; then
                bcftools view -G "\$f" -Oz -o "\$g"
                tabix -p vcf "\$g"
                vcf_files="\$vcf_files \$g"
            fi
        done

        if [ -n "\$vcf_files" ]; then
            bcftools concat -a \$vcf_files \\
                | bcftools sort \\
                | bcftools norm -d exact -Oz -o ${id}.consensus.vcf.gz
        else
            echo '##fileformat=VCFv4.2' | bgzip -c > ${id}.consensus.vcf.gz
        fi

        tabix -p vcf ${id}.consensus.vcf.gz
    ) > logs/${id}.consensus.log 2>&1
    """
}
