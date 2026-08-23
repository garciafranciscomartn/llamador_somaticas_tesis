#!/usr/bin/env nextflow
// ============================================================
// Pairwise Multi-Caller Somatic Variant Calling Pipeline
// Nextflow DSL2
// ============================================================
// Starts from BAM files and runs all-vs-all pairwise calling
// with 4 callers (Mutect2, VarScan2, VarDict, Strelka2).
// Per pair: consensus via bcftools isec (>= min_callers agree).
// Per tumor: merge cross-group consensus VCFs → frequency filter.
//
// Run:
//   nextflow run nextflow/main.nf \
//       -params-file nextflow/config/params.yml \
//       -profile local
// ============================================================

nextflow.enable.dsl = 2

include { PAIRWISE_WORKFLOW } from './workflows/pairwise.nf'

log.info """
╔══════════════════════════════════════════════════════════╗
║   Pairwise Somatic Variant Calling  (Nextflow DSL2)      ║
║   Samples : ${params.samples.padRight(46)}║
║   Outdir  : ${params.outdir.padRight(46)}║
╚══════════════════════════════════════════════════════════╝
""".stripIndent()

workflow {
    PAIRWISE_WORKFLOW()
}

workflow.onComplete {
    log.info (workflow.success
        ? "\n✓ Pipeline finished successfully [${workflow.duration}]"
        : "\n✗ Pipeline failed – check .nextflow.log for details")
}
