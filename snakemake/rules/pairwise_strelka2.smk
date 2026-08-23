# ============================================================
# Pairwise Strelka2 Somatic Variant Calling
# Default: run Strelka2 directly without Manta indel candidates
#          (set strelka2.run_manta: true in config to enable Manta)
# ============================================================

STRELKA2_BIN = config.get("strelka2", {}).get("bin",
    "/home/andres/herramientas/strelka-2.9.10.centos6_x86_64/bin")
MANTA_BIN = config.get("strelka2", {}).get("manta_bin",
    "/home/andres/herramientas/manta-1.6.0.centos6_x86_64/bin")


rule prep_strelka_bed:
    """bgzip + tabix the BED file required by Strelka2/Manta --callRegions."""
    input:
        bed=config["reference"]["bed"],
    output:
        bed_gz=config["reference"]["bed"] + ".gz",
        bed_tbi=config["reference"]["bed"] + ".gz.tbi",
    log:
        OUTDIR + "/logs/strelka2/prep_bed.log",
    shell:
        """
        (
            bgzip -c {input.bed} > {output.bed_gz}
            tabix -p bed {output.bed_gz}
        ) > {log} 2>&1
        """


rule manta_pair:
    """Run Manta to generate indel candidates for Strelka2 (optional)."""
    input:
        tumor_bam=get_tumor_bam,
        normal_bam=get_normal_bam,
        ref=config["reference"]["genome"],
        bed_gz=config["reference"]["bed"] + ".gz",
        bed_tbi=config["reference"]["bed"] + ".gz.tbi",
    output:
        candidate_indels="{outdir}/variant_calling/strelka2/{pair}/manta/results/variants/candidateSmallIndels.vcf.gz",
    params:
        manta_dir="{outdir}/variant_calling/strelka2/{pair}/manta",
        manta_bin=MANTA_BIN,
    log:
        "{outdir}/logs/strelka2/{pair}.manta.log",
    threads: config.get("resources", {}).get("strelka2", {}).get("threads", 4)
    shell:
        """
        (
            rm -rf {params.manta_dir}
            {params.manta_bin}/configManta.py \
                --normalBam {input.normal_bam} \
                --tumorBam  {input.tumor_bam} \
                --referenceFasta {input.ref} \
                --callRegions {input.bed_gz} \
                --exome \
                --runDir {params.manta_dir}
            {params.manta_dir}/runWorkflow.py -m local -j {threads}
        ) > {log} 2>&1
        """


def strelka2_input(wildcards):
    """Build Strelka2 inputs depending on whether Manta is enabled."""
    base = {
        "tumor_bam":  get_tumor_bam(wildcards),
        "normal_bam": get_normal_bam(wildcards),
        "ref":        config["reference"]["genome"],
        "bed_gz":     config["reference"]["bed"] + ".gz",
        "bed_tbi":    config["reference"]["bed"] + ".gz.tbi",
    }
    run_manta = config.get("strelka2", {}).get("run_manta", False)
    if run_manta:
        base["candidate_indels"] = (
            f"{OUTDIR}/variant_calling/strelka2/{wildcards.pair}/manta"
            "/results/variants/candidateSmallIndels.vcf.gz"
        )
    return base


rule strelka2_pair:
    """Call somatic SNVs and indels with Strelka2 (tumor-normal)."""
    input:
        unpack(strelka2_input),
    output:
        snvs="{outdir}/variant_calling/strelka2/{pair}/strelka_run/results/variants/somatic.snvs.vcf.gz",
        indels="{outdir}/variant_calling/strelka2/{pair}/strelka_run/results/variants/somatic.indels.vcf.gz",
    params:
        strelka_dir="{outdir}/variant_calling/strelka2/{pair}/strelka_run",
        strelka_bin=STRELKA2_BIN,
        run_manta=config.get("strelka2", {}).get("run_manta", False),
    log:
        "{outdir}/logs/strelka2/{pair}.strelka2.log",
    threads: config.get("resources", {}).get("strelka2", {}).get("threads", 4)
    shell:
        """
        (
            rm -rf {params.strelka_dir}

            IND_ARG=""
            if [ "{params.run_manta}" = "True" ] || [ "{params.run_manta}" = "true" ]; then
                CAND="{params.strelka_dir}/../manta/results/variants/candidateSmallIndels.vcf.gz"
                [ -f "$CAND" ] && IND_ARG="--indelCandidates $CAND"
            fi

            {params.strelka_bin}/configureStrelkaSomaticWorkflow.py \
                --normalBam {input.normal_bam} \
                --tumorBam  {input.tumor_bam} \
                --referenceFasta {input.ref} \
                $IND_ARG \
                --callRegions {input.bed_gz} \
                --exome \
                --runDir {params.strelka_dir}

            {params.strelka_dir}/runWorkflow.py -m local -j {threads}
        ) > {log} 2>&1
        """


rule strelka2_merge_pair:
    """Merge Strelka2 SNVs + indels, keep PASS variants, index."""
    input:
        snvs="{outdir}/variant_calling/strelka2/{pair}/strelka_run/results/variants/somatic.snvs.vcf.gz",
        indels="{outdir}/variant_calling/strelka2/{pair}/strelka_run/results/variants/somatic.indels.vcf.gz",
    output:
        vcf="{outdir}/variant_calling/strelka2/{pair}/{pair}.strelka2.pass.vcf.gz",
        tbi="{outdir}/variant_calling/strelka2/{pair}/{pair}.strelka2.pass.vcf.gz.tbi",
    log:
        "{outdir}/logs/strelka2/{pair}.merge.log",
    shell:
        """
        (
            bcftools concat -a {input.snvs} {input.indels} 2>/dev/null \
                | bcftools sort 2>/dev/null \
                | bcftools view -f .,PASS -Oz -o {output.vcf}
            tabix -p vcf {output.vcf}
        ) > {log} 2>&1
        """
