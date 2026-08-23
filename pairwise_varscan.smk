# ============================================================
# Pairwise VarScan2 Somatic Variant Calling
# samtools mpileup → VarScan somatic → merge SNP+indel → PASS
#
# IMPORTANT: VarScan marks both germline (SS=1) and LOH (SS=3) as
# FILTER=PASS. The bcftools -i 'INFO/SS=2' filter keeps only true
# somatic calls (SS=2) so non-somatic variants do not inflate
# consensus counts in the downstream isec step.
# ============================================================


rule varscan_mpileup_pair:
    """Generate samtools pileups for tumor and normal BAMs."""
    input:
        tumor_bam=get_tumor_bam,
        normal_bam=get_normal_bam,
        ref=config["reference"]["genome"],
        bed=config["reference"]["bed"],
    output:
        normal_pileup=temp("{outdir}/variant_calling/varscan2/{pair}/{pair}.normal.pileup"),
        tumor_pileup=temp("{outdir}/variant_calling/varscan2/{pair}/{pair}.tumor.pileup"),
    log:
        "{outdir}/logs/varscan2/{pair}.mpileup.log",
    threads: 2
    shell:
        """
        samtools mpileup -q 1 -Q 20 -f {input.ref} \
            -l {input.bed} \
            {input.normal_bam} > {output.normal_pileup} 2> {log}
        samtools mpileup -q 1 -Q 20 -f {input.ref} \
            -l {input.bed} \
            {input.tumor_bam} > {output.tumor_pileup} 2>> {log}
        """


rule varscan_somatic_pair:
    """Call somatic variants with VarScan2."""
    input:
        normal_pileup="{outdir}/variant_calling/varscan2/{pair}/{pair}.normal.pileup",
        tumor_pileup="{outdir}/variant_calling/varscan2/{pair}/{pair}.tumor.pileup",
    output:
        snp="{outdir}/variant_calling/varscan2/{pair}/{pair}.varscan2.snp.vcf",
        indel="{outdir}/variant_calling/varscan2/{pair}/{pair}.varscan2.indel.vcf",
    params:
        output_prefix="{outdir}/variant_calling/varscan2/{pair}/{pair}.varscan2",
        min_coverage=config.get("varscan2", {}).get("min_coverage", 10),
        min_var_freq=config.get("varscan2", {}).get("min_var_freq", 0.01),
        p_value=config.get("varscan2", {}).get("p_value", 0.99),
        somatic_p_value=config.get("varscan2", {}).get("somatic_p_value", 0.05),
    log:
        "{outdir}/logs/varscan2/{pair}.somatic.log",
    shell:
        """
        varscan somatic \
            {input.normal_pileup} \
            {input.tumor_pileup} \
            {params.output_prefix} \
            --min-coverage {params.min_coverage} \
            --min-var-freq {params.min_var_freq} \
            --p-value {params.p_value} \
            --somatic-p-value {params.somatic_p_value} \
            --strand-filter 1 \
            --output-vcf 1 \
            > {log} 2>&1
        """


rule varscan_merge_pair:
    """Merge VarScan SNPs and indels, keep SS=2 (somatic only), index."""
    input:
        snp="{outdir}/variant_calling/varscan2/{pair}/{pair}.varscan2.snp.vcf",
        indel="{outdir}/variant_calling/varscan2/{pair}/{pair}.varscan2.indel.vcf",
    output:
        vcf="{outdir}/variant_calling/varscan2/{pair}/{pair}.varscan2.pass.vcf.gz",
        tbi="{outdir}/variant_calling/varscan2/{pair}/{pair}.varscan2.pass.vcf.gz.tbi",
    log:
        "{outdir}/logs/varscan2/{pair}.merge.log",
    shell:
        """
        (
            bgzip -c {input.snp}   > {input.snp}.gz
            tabix -p vcf {input.snp}.gz
            bgzip -c {input.indel} > {input.indel}.gz
            tabix -p vcf {input.indel}.gz

            bcftools concat -a {input.snp}.gz {input.indel}.gz \
                | bcftools view -i 'INFO/SS=2' \
                | bcftools sort -Oz -o {output.vcf}
            tabix -p vcf {output.vcf}

            rm -f {input.snp}.gz {input.snp}.gz.tbi \
                  {input.indel}.gz {input.indel}.gz.tbi
        ) > {log} 2>&1
        """
