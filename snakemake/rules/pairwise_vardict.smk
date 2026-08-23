# ============================================================
# Pairwise VarDict Somatic Variant Calling
# vardict-java | testsomatic.R | var2vcf_paired.pl → PASS filter
# ============================================================


rule vardict_pair:
    """Call somatic variants with VarDict (tumor-normal)."""
    input:
        tumor_bam=get_tumor_bam,
        normal_bam=get_normal_bam,
        ref=config["reference"]["genome"],
        bed=config["reference"]["bed"],
    output:
        vcf="{outdir}/variant_calling/vardict/{pair}/{pair}.vardict.pass.vcf.gz",
        tbi="{outdir}/variant_calling/vardict/{pair}/{pair}.vardict.pass.vcf.gz.tbi",
    params:
        tumor_name=get_tumor_id,
        normal_name=get_normal_id,
        min_af=config.get("vardict", {}).get("min_allele_freq", 0.01),
    log:
        "{outdir}/logs/vardict/{pair}.log",
    threads: config.get("resources", {}).get("vardict", {}).get("threads", 4)
    shell:
        """
        (
            vardict-java \
                -G {input.ref} \
                -f {params.min_af} \
                -N {params.tumor_name} \
                -b "{input.tumor_bam}|{input.normal_bam}" \
                -c 1 -S 2 -E 3 \
                --nosv \
                -th {threads} \
                {input.bed} \
            | testsomatic.R \
            | var2vcf_paired.pl \
                -N "{params.tumor_name}|{params.normal_name}" \
                -f {params.min_af} \
            | bcftools view \
                -i 'INFO/STATUS="StrongSomatic" || INFO/STATUS="LikelySomatic"' \
            | bcftools view -f .,PASS -Oz -o {output.vcf}

            tabix -p vcf {output.vcf}
        ) > {log} 2>&1
        """
