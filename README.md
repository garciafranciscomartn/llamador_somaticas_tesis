# Pairwise Somatic Variant Calling Pipeline

Pairwise multi-caller somatic variant calling pipeline that starts from **BAM files** and produces per-tumor frequency-filtered consensus VCFs.

## Pipeline Overview

For every (tumor, normal) permutation of the input samples:

1. **4 callers** run in tumor-normal mode:
   - **Mutect2** (GATK) – with OxoG correction, contamination estimation, and `FilterMutectCalls`
   - **VarScan2** – pileup-based calling; only SS=2 (somatic) variants retained
   - **VarDict** – StrongSomatic + LikelySomatic variants retained
   - **Strelka2** – SNVs and indels merged; PASS variants retained

2. Each caller's PASS VCF is **normalized** (multi-allelic split, left-aligned indels).

3. **Consensus** via `bcftools isec`: keep variants supported by ≥ 2 callers (configurable).

Per tumor:

4. **Merge** consensus VCFs from **cross-group pairs only** (pseudo-normal from a different group than the tumor). Each variant is annotated with `PAIR_COUNT`, `PAIR_FREQ`, and `PAIR_NAMES`.

5. **Frequency filter**: keep variants with `PAIR_FREQ >= 0.75` (found in ≥ 75% of cross-group comparisons). This threshold ensures true somatic variants — which appear consistently regardless of the pseudo-normal — are retained, while noise and artifacts are discarded.

   > **Why cross-group?** Intra-group pairs share biology. Variants that are somatic within a group may look germline/shared to Mutect2 when comparing two samples from the same group, causing them to be missed. Using only cross-group pairs as the denominator prevents these group-specific somatic calls from being penalized.

## Outputs

```
results/
├── variant_calling/
│   ├── mutect2/{pair}/        – raw + filtered + PASS Mutect2 VCFs
│   ├── varscan2/{pair}/       – PASS VarScan2 VCFs
│   ├── vardict/{pair}/        – PASS VarDict VCFs
│   ├── strelka2/{pair}/       – PASS Strelka2 VCFs
│   ├── norm/{pair}/           – normalized per-caller VCFs
│   └── consensus/{pair}/      – per-pair consensus VCFs
└── merge/{tumor}/
    ├── {tumor}.all_pairs.vcf.gz       – merged cross-group VCF (PAIR_FREQ annotated)
    └── {tumor}.freq_filtered.vcf.gz   – final output (PAIR_FREQ >= 0.75)
```

---

## Sample Sheet

Both versions use a TSV sample sheet with the following columns:

| Column      | Required | Description                                            |
|-------------|----------|--------------------------------------------------------|
| `sample_id` | Yes      | Unique sample identifier                               |
| `bam_path`  | Yes      | Absolute path to the preprocessed BAM (must have .bai)|
| `group`     | No       | Group label (e.g. `EBV+`, `EBV-`). Used to restrict the frequency filter denominator to cross-group pairs. If omitted every sample is treated as its own group. |

Example (`config/samples.tsv`):

```
sample_id    bam_path                          group
SampleA      /data/bams/SampleA.bqsr.bam       EBV+
SampleB      /data/bams/SampleB.bqsr.bam       EBV-
SampleC      /data/bams/SampleC.bqsr.bam       EBV-
```

---

## Running the Snakemake Version

### Prerequisites

- Snakemake ≥ 7
- Python ≥ 3.8 with `pandas`
- `gatk` (GATK 4), `samtools`, `bcftools`, `bgzip`, `tabix` on `PATH`
- `varscan` on `PATH`
- `vardict-java`, `testsomatic.R`, `var2vcf_paired.pl` on `PATH`
- Strelka2 and Manta binaries at the paths configured in `config.yaml`

### Configuration

Edit `snakemake/config/config.yaml`:

- Set `samples:` to the path of your sample sheet.
- Update all paths under `reference:` to point to your reference files.
- Adjust `outdir:` (default: `results`).
- Tune `frequency_filter.min_frequency` if needed (default: `0.75`).

### Run

```bash
# Dry-run (preview jobs)
snakemake --snakefile snakemake/Snakefile \
          --configfile snakemake/config/config.yaml \
          --cores 64 \
          --dry-run

# Full run
snakemake --snakefile snakemake/Snakefile \
          --configfile snakemake/config/config.yaml \
          --cores 64

# Resume after failure
snakemake --snakefile snakemake/Snakefile \
          --configfile snakemake/config/config.yaml \
          --cores 64 \
          --rerun-incomplete
```

### Override parameters on the command line

```bash
# Use a different sample sheet
snakemake --snakefile snakemake/Snakefile \
          --configfile snakemake/config/config.yaml \
          --config samples=snakemake/config/samples_subset.tsv \
          --cores 64

# Change the frequency filter threshold
snakemake ... --config frequency_filter='{"min_frequency": 0.5}'
```

### Running in a screen session (recommended for long runs)

```bash
screen -S somatic
snakemake --snakefile snakemake/Snakefile \
          --configfile snakemake/config/config.yaml \
          --cores 64 \
          --keep-going \
          2>&1 | tee snakemake_run.log
# Detach: Ctrl+A D
```

---

## Running the Nextflow Version

### Prerequisites

- Nextflow ≥ 23
- Java ≥ 11
- Same tool dependencies as Snakemake version (all must be on `PATH`)

### Configuration

Edit `nextflow/nextflow.config`:

- Update all `params.*` paths (genome, germline_resource, bed, intervals, strelka2_bin, etc.).
- Adjust `params.min_frequency` if needed (default: `0.75`).
- Update `params.samples` or pass it at runtime with `--samples`.

### Run

```bash
# From the repo root
nextflow run nextflow/main.nf \
    -profile local \
    --samples nextflow/config/samples.tsv \
    --outdir results

# Resume after failure
nextflow run nextflow/main.nf \
    -profile local \
    --samples nextflow/config/samples.tsv \
    --outdir results \
    -resume
```

### Override parameters at runtime

```bash
nextflow run nextflow/main.nf \
    -profile local \
    --samples /path/to/samples.tsv \
    --outdir /path/to/output \
    --min_frequency 0.5 \
    --consensus_min_callers 3
```

### Running in a screen session (recommended)

```bash
screen -S somatic_nf
nextflow run nextflow/main.nf \
    -profile local \
    --samples nextflow/config/samples.tsv \
    --outdir results \
    -resume \
    2>&1 | tee nextflow_run.log
# Detach: Ctrl+A D
```

### Execution reports

Nextflow generates HTML reports automatically:

- `results/pipeline_report.html` – resource usage per process
- `results/pipeline_timeline.html` – execution timeline

---

## Key Parameters Reference

| Parameter                   | Snakemake (`config.yaml`)              | Nextflow (`nextflow.config`)          | Default          |
|-----------------------------|----------------------------------------|---------------------------------------|------------------|
| Sample sheet                | `samples:`                             | `params.samples`                      | `config/samples.tsv` |
| Reference genome            | `reference.genome:`                    | `params.genome`                       | hg38             |
| Germline resource           | `reference.germline_resource:`         | `params.germline_resource`            | gnomAD af-only   |
| Panel of normals            | `reference.panel_of_normals:`          | `params.panel_of_normals`             | 1000g PON        |
| Target BED                  | `reference.bed:`                       | `params.bed`                          | bedfile.bed      |
| Mutect2 intervals           | `mutect2.intervals:`                   | `params.intervals`                    | bedfile.interval_list |
| Mutect2 extra args          | `mutect2.extra:`                       | `params.mutect2_extra`                | `--af-of-alleles-not-in-resource 0.0000025` |
| VarScan2 min coverage       | `varscan2.min_coverage:`               | `params.varscan_min_coverage`         | `10`             |
| VarScan2 min var freq       | `varscan2.min_var_freq:`               | `params.varscan_min_var_freq`         | `0.01`           |
| VarDict min AF              | `vardict.min_allele_freq:`             | `params.vardict_min_af`               | `0.01`           |
| Callers                     | `consensus.callers:`                   | `params.consensus_callers`            | all 4            |
| Min callers for consensus   | `consensus.min_callers:`               | `params.consensus_min_callers`        | `2`              |
| Frequency filter threshold  | `frequency_filter.min_frequency:`      | `params.min_frequency`                | `0.75`           |
| Run Manta before Strelka2   | `strelka2.run_manta:`                  | `params.strelka2_run_manta`           | `false`          |

---

## Folder Structure

```
somatic_calling_pipeline/
├── README.md
├── snakemake/
│   ├── Snakefile                    # Main workflow entry point
│   ├── config/
│   │   ├── config.yaml              # All configuration parameters
│   │   └── samples.tsv              # Sample sheet template
│   ├── rules/
│   │   ├── pairwise_calling.smk     # Mutect2 pairwise rules
│   │   ├── pairwise_varscan.smk     # VarScan2 pairwise rules
│   │   ├── pairwise_vardict.smk     # VarDict pairwise rules
│   │   ├── pairwise_strelka2.smk    # Strelka2 pairwise rules
│   │   ├── consensus.smk            # Normalize + bcftools isec
│   │   └── merge_frequency.smk      # Merge + frequency filter
│   └── scripts/
│       ├── merge_pair_vcfs.py       # Annotates PAIR_COUNT / PAIR_FREQ
│       └── frequency_filter.py      # Filters by PAIR_FREQ threshold
└── nextflow/
    ├── main.nf                      # Nextflow entry point
    ├── nextflow.config              # Global configuration and profiles
    ├── config/
    │   └── samples.tsv              # Sample sheet template
    ├── workflows/
    │   └── pairwise.nf              # Main pairwise workflow
    └── modules/
        ├── mutect2.nf               # GATK Mutect2 processes
        ├── strelka2.nf              # Strelka2 / Manta processes
        ├── varscan2.nf              # VarScan2 processes
        ├── vardict.nf               # VarDict process
        ├── consensus.nf             # Normalize + isec processes
        └── merge.nf                 # Merge + frequency filter processes
```
