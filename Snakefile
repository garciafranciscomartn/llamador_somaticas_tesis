"""
Pairwise Multi-Caller Somatic Variant Calling Pipeline
=======================================================
Starts directly from BAM files (pre-processed, duplicate-marked, BQSR).

Per pair (all-vs-all permutations):
  1. Each of 4 callers produces a PASS VCF:
       Mutect2, VarScan2, VarDict, Strelka2
  2. Normalize + bcftools isec → consensus VCF (≥ MIN_CALLERS agree)

Per tumor:
  3. Merge consensus VCFs from cross-group pairs → count PAIR_FREQ
  4. Frequency filter: keep variants found in >= min_frequency of cross-group pairs
     (cross-group = pseudo-normal from a DIFFERENT group than the tumor)
     This prevents true group-specific somatic variants from being discarded.

Usage:
  snakemake --snakefile snakemake/Snakefile \
      --configfile snakemake/config/config.yaml \
      --cores <N>
"""

import pandas as pd
from itertools import permutations

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
configfile: "snakemake/config/config.yaml"

OUTDIR = config["outdir"]

CALLERS = config.get("consensus", {}).get("callers", ["mutect2", "varscan2", "vardict", "strelka2"])
MIN_CALLERS = config.get("consensus", {}).get("min_callers", 2)

# Load sample sheet: required columns = sample_id, bam_path
# Optional column: group  (used for cross-group frequency denominator)
samples_df = pd.read_table(config["samples"], dtype=str).set_index("sample_id", drop=False)
SAMPLE_IDS = list(samples_df["sample_id"])
BAM_MAP    = dict(zip(samples_df["sample_id"], samples_df["bam_path"]))

# All pairwise (tumor, normal) permutations
ALL_PAIRS  = [(t, n) for t, n in permutations(SAMPLE_IDS, 2)]
PAIR_NAMES = [f"{t}_vs_{n}" for t, n in ALL_PAIRS]
PAIR_MAP   = {f"{t}_vs_{n}": (t, n) for t, n in ALL_PAIRS}

# Per-tumor pair list (all pairs where this sample is the tumor)
TUMOR_PAIRS = {t: [f"{t}_vs_{n}" for n in SAMPLE_IDS if n != t] for t in SAMPLE_IDS}

# Build group map: sample_id -> group label
# If the sample sheet lacks a 'group' column every sample is its own group
# (falls back to using ALL pairs as the denominator – same as original behaviour).
if "group" in samples_df.columns:
    GROUP_MAP = dict(zip(samples_df["sample_id"], samples_df["group"].fillna("")))
else:
    GROUP_MAP = {s: s for s in SAMPLE_IDS}

# Cross-group pairs per tumor: only pairs where pseudo-normal ∈ a different group.
# Using ONLY these pairs in the frequency denominator prevents intra-group
# comparisons (where somatic variants look shared/germline and are not called)
# from artificially lowering PAIR_FREQ and causing real somatic calls to be
# discarded by the frequency filter.
CROSS_GROUP_PAIRS: dict = {}
for t in SAMPLE_IDS:
    t_group = GROUP_MAP.get(t, "")
    CROSS_GROUP_PAIRS[t] = [
        f"{t}_vs_{n}" for n in SAMPLE_IDS
        if n != t and GROUP_MAP.get(n, "") != t_group
    ]
    if not CROSS_GROUP_PAIRS[t]:          # edge case: all samples in same group
        CROSS_GROUP_PAIRS[t] = TUMOR_PAIRS[t]

# Frequency filter threshold
MIN_FREQ = config.get("frequency_filter", {}).get("min_frequency", 0.75)


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
def get_tumor_bam(wildcards):
    return BAM_MAP[PAIR_MAP[wildcards.pair][0]]

def get_normal_bam(wildcards):
    return BAM_MAP[PAIR_MAP[wildcards.pair][1]]

def get_tumor_id(wildcards):
    return PAIR_MAP[wildcards.pair][0]

def get_normal_id(wildcards):
    return PAIR_MAP[wildcards.pair][1]


import subprocess

def _bam_sample_name(bam_path):
    """Return SM tag from BAM header; fallback to RG ID; None if unreadable."""
    samtools_cmds = [
        "/home/andres/herramientas/samtools/bin/samtools",
        "/usr/bin/samtools",
        "samtools",
    ]
    out = None
    for cmd in samtools_cmds:
        try:
            out = subprocess.check_output([cmd, "view", "-H", bam_path], text=True)
            break
        except Exception:
            continue
    if not out:
        return None
    for line in out.splitlines():
        if line.startswith("@RG"):
            for part in line.split("\t"):
                if part.startswith("SM:"):
                    return part[3:]
    for line in out.splitlines():
        if line.startswith("@RG"):
            for part in line.split("\t"):
                if part.startswith("ID:"):
                    return part[3:]
    return None

def get_tumor_sample_name(wildcards):
    tid = PAIR_MAP[wildcards.pair][0]
    name = _bam_sample_name(BAM_MAP.get(tid))
    return name if name else tid

def get_normal_sample_name(wildcards):
    nid = PAIR_MAP[wildcards.pair][1]
    name = _bam_sample_name(BAM_MAP.get(nid))
    return name if name else nid

def get_pair_consensus_vcfs(wildcards):
    """VCFs for cross-group pairs only (used as frequency denominator)."""
    pairs = CROSS_GROUP_PAIRS[wildcards.tumor]
    return expand(
        "{outdir}/variant_calling/consensus/{pair}/{pair}.consensus.vcf.gz",
        outdir=OUTDIR, pair=pairs,
    )


# -----------------------------------------------------------------------------
# Include rule modules
# -----------------------------------------------------------------------------
include: "rules/pairwise_calling.smk"

if "varscan2" in CALLERS:
    include: "rules/pairwise_varscan.smk"

if "vardict" in CALLERS:
    include: "rules/pairwise_vardict.smk"

if "strelka2" in CALLERS:
    include: "rules/pairwise_strelka2.smk"

include: "rules/consensus.smk"
include: "rules/merge_frequency.smk"


# -----------------------------------------------------------------------------
# Target rule
# -----------------------------------------------------------------------------
rule all:
    input:
        # Per-pair consensus VCFs
        expand(
            "{outdir}/variant_calling/consensus/{pair}/{pair}.consensus.vcf.gz",
            outdir=OUTDIR, pair=PAIR_NAMES,
        ),
        # Per-tumor all-pairs merged VCF
        expand(
            "{outdir}/merge/{tumor}/{tumor}.all_pairs.vcf.gz",
            outdir=OUTDIR, tumor=SAMPLE_IDS,
        ),
        # Per-tumor frequency-filtered VCF (final output)
        expand(
            "{outdir}/merge/{tumor}/{tumor}.freq_filtered.vcf.gz",
            outdir=OUTDIR, tumor=SAMPLE_IDS,
        ),
