# ============================================================
# Pairwise Mutect2 Calling Rules
# Runs Mutect2 in tumor-normal mode for each (tumor, normal) pair.
# BAM inputs, PAIR_MAP, and sample-name helpers are defined in Snakefile.
# ============================================================


def mutect2_pon_arg():
    pon = config.get("reference", {}).get("panel_of_normals", "")
    return f"--panel-of-normals {pon}" if pon else ""


def mutect2_germline_arg():
    gr = config.get("reference", {}).get("germline_resource", "")
    return f"--germline-resource {gr}" if gr else ""


def mutect2_intervals_arg():
    iv = config.get("mutect2", {}).get("intervals", "")
    return f"--intervals {iv}" if iv else ""


rule mutect2_pair:
    """Call somatic variants with Mutect2 (tumor vs pseudo-normal)."""
    input:
        tumor_bam=get_tumor_bam,
        normal_bam=get_normal_bam,
        ref=config["reference"]["genome"],
    output:
        vcf="{outdir}/variant_calling/mutect2/{pair}/{pair}.mutect2.vcf.gz",
        tbi="{outdir}/variant_calling/mutect2/{pair}/{pair}.mutect2.vcf.gz.tbi",
        stats="{outdir}/variant_calling/mutect2/{pair}/{pair}.mutect2.vcf.gz.stats",
        f1r2="{outdir}/variant_calling/mutect2/{pair}/{pair}.f1r2.tar.gz",
    params:
        tumor_name=get_tumor_sample_name,
        normal_name=get_normal_sample_name,
        pon=mutect2_pon_arg(),
        germline=mutect2_germline_arg(),
        intervals=mutect2_intervals_arg(),
        extra=config.get("mutect2", {}).get("extra", ""),
        java_opts="-Xmx{}m".format(config["resources"]["mutect2"]["mem_mb"]),
    log:
        "{outdir}/logs/mutect2/{pair}.log",
    threads: config["resources"]["mutect2"]["threads"]
    shell:
        """
        TNAME=$(samtools view -H {input.tumor_bam} \
            | awk -F'\t' '/^@RG/ {{ for(i=1;i<=NF;i++) if($i ~ /^SM:/) {{print substr($i,4); exit}} }}')
        [ -z "$TNAME" ] && TNAME="{params.tumor_name}"

        NNAME=$(samtools view -H {input.normal_bam} \
            | awk -F'\t' '/^@RG/ {{ for(i=1;i<=NF;i++) if($i ~ /^SM:/) {{print substr($i,4); exit}} }}')
        [ -z "$NNAME" ] && NNAME="{params.normal_name}"

        gatk --java-options '{params.java_opts}' Mutect2 \
            --input {input.tumor_bam} \
            --input {input.normal_bam} \
            --normal-sample "$NNAME" \
            --reference {input.ref} \
            --tumor-sample "$TNAME" \
            {params.pon} \
            {params.germline} \
            {params.intervals} \
            --f1r2-tar-gz {output.f1r2} \
            --output {output.vcf} \
            {params.extra} \
            > {log} 2>&1
        """


rule learn_read_orientation_pair:
    """Learn the read orientation model for OxoG artifact correction."""
    input:
        f1r2="{outdir}/variant_calling/mutect2/{pair}/{pair}.f1r2.tar.gz",
    output:
        model="{outdir}/variant_calling/mutect2/{pair}/{pair}.read_orientation_model.tar.gz",
    params:
        java_opts="-Xmx{}m".format(config["resources"]["default"]["mem_mb"]),
    log:
        "{outdir}/logs/learn_read_orientation/{pair}.log",
    shell:
        """
        gatk --java-options '{params.java_opts}' LearnReadOrientationModel \
            --input {input.f1r2} \
            --output {output.model} \
            > {log} 2>&1
        """


rule pileup_summaries_tumor_pair:
    """Pileup summaries for the tumor BAM."""
    input:
        bam=get_tumor_bam,
        germline_resource=config["reference"]["germline_resource"],
        intervals=config["mutect2"]["intervals"],
    output:
        table="{outdir}/variant_calling/mutect2/{pair}/{pair}.tumor.pileup_summaries.table",
    params:
        java_opts="-Xmx{}m".format(config["resources"]["pileup"]["mem_mb"]),
    threads: config["resources"]["pileup"]["threads"]
    log:
        "{outdir}/logs/pileup_summaries/{pair}.tumor.log",
    shell:
        """
        gatk --java-options '{params.java_opts}' GetPileupSummaries \
            --input {input.bam} \
            --variant {input.germline_resource} \
            --intervals {input.intervals} \
            --output {output.table} \
            > {log} 2>&1
        """


rule pileup_summaries_normal_pair:
    """Pileup summaries for the normal BAM."""
    input:
        bam=get_normal_bam,
        germline_resource=config["reference"]["germline_resource"],
        intervals=config["mutect2"]["intervals"],
    output:
        table="{outdir}/variant_calling/mutect2/{pair}/{pair}.normal.pileup_summaries.table",
    params:
        java_opts="-Xmx{}m".format(config["resources"]["pileup"]["mem_mb"]),
    threads: config["resources"]["pileup"]["threads"]
    log:
        "{outdir}/logs/pileup_summaries/{pair}.normal.log",
    shell:
        """
        gatk --java-options '{params.java_opts}' GetPileupSummaries \
            --input {input.bam} \
            --variant {input.germline_resource} \
            --intervals {input.intervals} \
            --output {output.table} \
            > {log} 2>&1
        """


rule calculate_contamination_pair:
    """Estimate cross-sample contamination."""
    input:
        tumor_table="{outdir}/variant_calling/mutect2/{pair}/{pair}.tumor.pileup_summaries.table",
        normal_table="{outdir}/variant_calling/mutect2/{pair}/{pair}.normal.pileup_summaries.table",
    output:
        contamination="{outdir}/variant_calling/mutect2/{pair}/{pair}.contamination.table",
        segments="{outdir}/variant_calling/mutect2/{pair}/{pair}.segments.table",
    params:
        java_opts="-Xmx{}m".format(config["resources"]["default"]["mem_mb"]),
    log:
        "{outdir}/logs/calculate_contamination/{pair}.log",
    shell:
        """
        gatk --java-options '{params.java_opts}' CalculateContamination \
            --input {input.tumor_table} \
            --matched-normal {input.normal_table} \
            --tumor-segmentation {output.segments} \
            --output {output.contamination} \
            > {log} 2>&1
        """


rule filter_mutect_calls_pair:
    """Apply FilterMutectCalls to the raw Mutect2 output."""
    input:
        vcf="{outdir}/variant_calling/mutect2/{pair}/{pair}.mutect2.vcf.gz",
        ref=config["reference"]["genome"],
        contamination="{outdir}/variant_calling/mutect2/{pair}/{pair}.contamination.table",
        segments="{outdir}/variant_calling/mutect2/{pair}/{pair}.segments.table",
        orientation_model="{outdir}/variant_calling/mutect2/{pair}/{pair}.read_orientation_model.tar.gz",
        stats="{outdir}/variant_calling/mutect2/{pair}/{pair}.mutect2.vcf.gz.stats",
    output:
        vcf="{outdir}/variant_calling/mutect2/{pair}/{pair}.mutect2.filtered.vcf.gz",
        tbi="{outdir}/variant_calling/mutect2/{pair}/{pair}.mutect2.filtered.vcf.gz.tbi",
    params:
        java_opts="-Xmx{}m".format(config["resources"]["filter"]["mem_mb"]),
    log:
        "{outdir}/logs/filter_mutect_calls/{pair}.log",
    shell:
        """
        gatk --java-options '{params.java_opts}' FilterMutectCalls \
            --variant {input.vcf} \
            --reference {input.ref} \
            --contamination-table {input.contamination} \
            --tumor-segmentation {input.segments} \
            --orientation-bias-artifact-priors {input.orientation_model} \
            --stats {input.stats} \
            --output {output.vcf} \
            > {log} 2>&1
        """


rule select_pass_pair:
    """Keep only PASS variants from filtered Mutect2 VCF."""
    input:
        vcf="{outdir}/variant_calling/mutect2/{pair}/{pair}.mutect2.filtered.vcf.gz",
        ref=config["reference"]["genome"],
    output:
        vcf="{outdir}/variant_calling/mutect2/{pair}/{pair}.mutect2.pass.vcf.gz",
        tbi="{outdir}/variant_calling/mutect2/{pair}/{pair}.mutect2.pass.vcf.gz.tbi",
    params:
        java_opts="-Xmx{}m".format(config["resources"]["default"]["mem_mb"]),
    log:
        "{outdir}/logs/select_pass/{pair}.log",
    shell:
        """
        gatk --java-options '{params.java_opts}' SelectVariants \
            --variant {input.vcf} \
            --reference {input.ref} \
            --exclude-filtered \
            --output {output.vcf} \
            > {log} 2>&1
        """
