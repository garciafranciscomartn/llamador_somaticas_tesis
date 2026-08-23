# ============================================================
# Merge & Frequency Filter Rules
# For each tumor sample:
#   1. Merge all CROSS-GROUP pairwise consensus VCFs into a
#      combined VCF, annotating each variant with:
#        PAIR_COUNT  – number of cross-group pairs that called it
#        PAIR_FREQ   – PAIR_COUNT / total cross-group pairs
#        PAIR_NAMES  – list of pair names that called it
#   2. Filter: keep only variants with PAIR_FREQ >= min_frequency
#
# Cross-group logic: the denominator only counts pairs where the
# pseudo-normal belongs to a DIFFERENT group than the tumor.
# This prevents intra-group comparisons (where a somatic variant
# looks shared/germline and is not called by Mutect2) from lowering
# PAIR_FREQ and causing real group-specific somatic calls to be
# discarded by the frequency filter.
# ============================================================


rule merge_pair_vcfs:
    """Merge cross-group consensus VCFs for a tumor; annotate PAIR_FREQ."""
    input:
        vcfs=get_pair_consensus_vcfs,   # cross-group pairs only
    output:
        all_vcf="{outdir}/merge/{tumor}/{tumor}.all_pairs.vcf.gz",
        all_tbi="{outdir}/merge/{tumor}/{tumor}.all_pairs.vcf.gz.tbi",
    params:
        n_pairs=lambda wc: len(CROSS_GROUP_PAIRS[wc.tumor]),
        tumor_id=lambda wc: wc.tumor,
    log:
        "{outdir}/logs/merge/{tumor}.merge.log",
    shell:
        """
        python3 snakemake/scripts/merge_pair_vcfs.py \
            --vcfs {input.vcfs} \
            --output {output.all_vcf}.tmp.gz \
            --tumor-id {params.tumor_id} \
            --n-pairs {params.n_pairs} \
            > {log} 2>&1

        zcat {output.all_vcf}.tmp.gz | bcftools sort -Oz -o {output.all_vcf}
        rm -f {output.all_vcf}.tmp.gz
        tabix -p vcf {output.all_vcf} >> {log} 2>&1
        """


rule frequency_filter:
    """Keep variants with PAIR_FREQ >= min_frequency (default 0.75)."""
    input:
        vcf="{outdir}/merge/{tumor}/{tumor}.all_pairs.vcf.gz",
    output:
        vcf="{outdir}/merge/{tumor}/{tumor}.freq_filtered.vcf.gz",
        tbi="{outdir}/merge/{tumor}/{tumor}.freq_filtered.vcf.gz.tbi",
    params:
        min_freq=MIN_FREQ,
    log:
        "{outdir}/logs/merge/{tumor}.freq_filter.log",
    shell:
        """
        python3 snakemake/scripts/frequency_filter.py \
            --input {input.vcf} \
            --output {output.vcf}.tmp.gz \
            --min-frequency {params.min_freq} \
            > {log} 2>&1

        zcat {output.vcf}.tmp.gz | bcftools sort -Oz -o {output.vcf}
        rm -f {output.vcf}.tmp.gz
        tabix -p vcf {output.vcf} >> {log} 2>&1
        """
