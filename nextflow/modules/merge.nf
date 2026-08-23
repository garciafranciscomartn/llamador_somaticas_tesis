// ============================================================
// modules/merge.nf
// Per-tumor merge of cross-group pairwise consensus VCFs,
// then frequency filter.
// ============================================================

process MERGE_PAIR_VCFS {
    /*
     * Merge all (cross-group) pairwise consensus VCFs for one tumor.
     * Adds PAIR_COUNT, PAIR_FREQ, PAIR_NAMES INFO fields.
     * The denominator (n_pairs) equals the number of cross-group pairs
     * so intra-group comparisons do not lower PAIR_FREQ.
     */
    tag "${tumor_id}"
    label 'process_medium'

    publishDir "${params.outdir}/merge/${tumor_id}",
        mode: 'copy', pattern: '*.{vcf.gz,vcf.gz.tbi}'

    input:
    tuple val(tumor_id), path(vcfs), path(tbis)

    output:
    tuple val(tumor_id),
          path("${tumor_id}.all_pairs.vcf.gz"),
          path("${tumor_id}.all_pairs.vcf.gz.tbi"), emit: vcf
    path "logs/${tumor_id}.merge.log",               emit: log

    script:
    def n_pairs  = vcfs instanceof List ? vcfs.size() : 1
    def vcf_list = vcfs instanceof List ? vcfs.join(' ') : vcfs.toString()
    def scripts  = "${projectDir}/../snakemake/scripts"
    """
    mkdir -p logs
    python3 ${scripts}/merge_pair_vcfs.py \\
        --vcfs ${vcf_list} \\
        --output ${tumor_id}.all_pairs.vcf.gz.tmp.gz \\
        --tumor-id ${tumor_id} \\
        --n-pairs ${n_pairs} \\
        > logs/${tumor_id}.merge.log 2>&1

    zcat ${tumor_id}.all_pairs.vcf.gz.tmp.gz | bcftools sort -Oz -o ${tumor_id}.all_pairs.vcf.gz
    rm -f ${tumor_id}.all_pairs.vcf.gz.tmp.gz
    tabix -p vcf ${tumor_id}.all_pairs.vcf.gz >> logs/${tumor_id}.merge.log 2>&1
    """
}


process FREQUENCY_FILTER {
    /*
     * Keep only variants with PAIR_FREQ >= params.min_frequency (default 0.75).
     * True somatic variants appear consistently across cross-group comparisons.
     */
    tag "${tumor_id}"
    label 'process_low'

    publishDir "${params.outdir}/merge/${tumor_id}",
        mode: 'copy', pattern: '*.{vcf.gz,vcf.gz.tbi}'

    input:
    tuple val(tumor_id), path(vcf), path(tbi)

    output:
    tuple val(tumor_id),
          path("${tumor_id}.freq_filtered.vcf.gz"),
          path("${tumor_id}.freq_filtered.vcf.gz.tbi"), emit: vcf
    path "logs/${tumor_id}.freq_filter.log",             emit: log

    script:
    def scripts  = "${projectDir}/../snakemake/scripts"
    def min_freq = params.min_frequency
    """
    mkdir -p logs
    python3 ${scripts}/frequency_filter.py \\
        --input ${vcf} \\
        --output ${tumor_id}.freq_filtered.vcf.gz.tmp.gz \\
        --min-frequency ${min_freq} \\
        > logs/${tumor_id}.freq_filter.log 2>&1

    zcat ${tumor_id}.freq_filtered.vcf.gz.tmp.gz | bcftools sort -Oz -o ${tumor_id}.freq_filtered.vcf.gz
    rm -f ${tumor_id}.freq_filtered.vcf.gz.tmp.gz
    tabix -p vcf ${tumor_id}.freq_filtered.vcf.gz >> logs/${tumor_id}.freq_filter.log 2>&1
    """
}
