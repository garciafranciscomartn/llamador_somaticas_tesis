# ============================================================
# Pairwise Consensus Rules
# 1) Normalize each caller's PASS VCF (split multi-allelic,
#    left-align indels) so all callers use the same representation.
# 2) bcftools isec -n+MIN_CALLERS to find variants agreed on by
#    at least MIN_CALLERS of the enabled callers.
# 3) Concatenate isec outputs, deduplicate → final consensus VCF.
# ============================================================

# CALLERS and MIN_CALLERS are defined in Snakefile


def caller_pass_vcf(pair, caller):
    """Return the PASS VCF path for a given pair and caller."""
    if caller == "mutect2":
        return f"{OUTDIR}/variant_calling/mutect2/{pair}/{pair}.mutect2.pass.vcf.gz"
    return f"{OUTDIR}/variant_calling/{caller}/{pair}/{pair}.{caller}.pass.vcf.gz"


rule norm_pair:
    """Normalize a caller's PASS VCF for consensus comparison."""
    wildcard_constraints:
        caller="mutect2|varscan2|vardict|strelka2",
    input:
        vcf=lambda wc: caller_pass_vcf(wc.pair, wc.caller),
        ref=config["reference"]["genome"],
    output:
        vcf="{outdir}/variant_calling/norm/{pair}/{pair}.{caller}.norm.vcf.gz",
        tbi="{outdir}/variant_calling/norm/{pair}/{pair}.{caller}.norm.vcf.gz.tbi",
    log:
        "{outdir}/logs/norm/{pair}.{caller}.log",
    shell:
        """
        (
            bcftools norm -m -both -f {input.ref} {input.vcf} \
                | bcftools sort -Oz -o {output.vcf}
            tabix -p vcf {output.vcf}
        ) > {log} 2>&1
        """


rule consensus_pair:
    """Find variants called by >= MIN_CALLERS, produce consensus VCF."""
    input:
        vcfs=lambda wc: expand(
            "{outdir}/variant_calling/norm/{pair}/{pair}.{caller}.norm.vcf.gz",
            outdir=wc.outdir, pair=wc.pair, caller=CALLERS,
        ),
        tbis=lambda wc: expand(
            "{outdir}/variant_calling/norm/{pair}/{pair}.{caller}.norm.vcf.gz.tbi",
            outdir=wc.outdir, pair=wc.pair, caller=CALLERS,
        ),
    output:
        vcf="{outdir}/variant_calling/consensus/{pair}/{pair}.consensus.vcf.gz",
        tbi="{outdir}/variant_calling/consensus/{pair}/{pair}.consensus.vcf.gz.tbi",
        sites="{outdir}/variant_calling/consensus/{pair}/{pair}.sites.txt",
    params:
        isec_dir="{outdir}/variant_calling/consensus/{pair}/isec",
        min_callers=MIN_CALLERS,
        n_callers=len(CALLERS),
    log:
        "{outdir}/logs/consensus/{pair}.log",
    shell:
        """
        (
            mkdir -p {params.isec_dir}

            bcftools isec -n+{params.min_callers} \
                {input.vcfs} \
                -p {params.isec_dir}

            cp {params.isec_dir}/sites.txt {output.sites}

            vcf_files=""
            for i in $(seq 0 $(({params.n_callers} - 1))); do
                f=$(printf "%s/%04d.vcf" "{params.isec_dir}" "$i")
                g=$(printf "%s/%04d.sites.vcf.gz" "{params.isec_dir}" "$i")
                if [ -f "$f" ]; then
                    bcftools view -G "$f" -Oz -o "$g"
                    tabix -p vcf "$g"
                    vcf_files="$vcf_files $g"
                fi
            done

            if [ -n "$vcf_files" ]; then
                bcftools concat -a $vcf_files \
                    | bcftools sort \
                    | bcftools norm -d exact -Oz -o {output.vcf}
            else
                echo '##fileformat=VCFv4.2' | bgzip -c > {output.vcf}
            fi

            tabix -p vcf {output.vcf}
        ) > {log} 2>&1
        """
