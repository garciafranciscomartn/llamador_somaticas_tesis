// ============================================================
// workflows/pairwise.nf
// ============================================================
// All-vs-all pairwise calling with 4 callers, consensus, and
// cross-group frequency filter.
// Mirrors Snakefile_pairwise logic exactly.
// ============================================================

nextflow.enable.dsl = 2

include { MUTECT2_PAIR;
          LEARN_ORIENTATION_PAIR;
          PILEUP_TUMOR_PAIR;
          PILEUP_NORMAL_PAIR;
          CALCULATE_CONTAMINATION_PAIR;
          FILTER_MUTECT_PAIR;
          SELECT_PASS_PAIR      } from '../modules/mutect2.nf'

include { PREP_STRELKA_BED;
          MANTA_PAIR;
          STRELKA2_PAIR;
          STRELKA2_PAIR_DIRECT;
          STRELKA2_MERGE_PAIR   } from '../modules/strelka2.nf'

include { VARSCAN_MPILEUP_PAIR;
          VARSCAN_SOMATIC_PAIR;
          VARSCAN_MERGE_PAIR    } from '../modules/varscan2.nf'

include { VARDICT_PAIR          } from '../modules/vardict.nf'

include { NORM_VCF;
          CONSENSUS_ISEC        } from '../modules/consensus.nf'

include { MERGE_PAIR_VCFS;
          FREQUENCY_FILTER      } from '../modules/merge.nf'


// ── Helper ─────────────────────────────────────────────────
def enabled_callers() {
    return params.consensus_callers.tokenize(',').collect { it.trim() }
}


// ── Main workflow ───────────────────────────────────────────
workflow PAIRWISE_WORKFLOW {

    // ── 1. Load sample sheet ─────────────────────────────
    // Columns: sample_id, bam_path, [group]
    ch_samples = Channel
        .fromPath(params.samples)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def bam  = file(row.bam_path)
            def bai  = file("${row.bam_path}.bai")
            def grp  = row.containsKey('group') ? row.group : row.sample_id
            tuple(row.sample_id, bam, bai, grp)
        }
        .collect()

    // ── 2. Generate all (tumor, normal) pairs ────────────
    ch_pairs = ch_samples.flatMap { samples ->
        def result = []
        samples.eachWithIndex { t, i ->
            samples.eachWithIndex { n, j ->
                if (i != j) {
                    result << tuple(
                        "${t[0]}_vs_${n[0]}",   // pair_name
                        t[0], t[1], t[2],        // tumor_id, bam, bai
                        n[0], n[1], n[2]         // normal_id, bam, bai
                    )
                }
            }
        }
        result
    }

    // ── 3. Reference channels ────────────────────────────
    ch_genome            = file(params.genome)
    ch_germline_resource = file(params.germline_resource)
    ch_germline_tbi      = file("${params.germline_resource}.tbi")
    ch_intervals         = file(params.intervals)
    ch_bed               = file(params.bed)

    // ── 4. Mutect2 ───────────────────────────────────────
    MUTECT2_PAIR(ch_pairs, ch_genome)

    LEARN_ORIENTATION_PAIR(
        MUTECT2_PAIR.out.vcf.map { pair, vcf, tbi, stats, f1r2 ->
            tuple(pair, f1r2) }
    )

    PILEUP_TUMOR_PAIR(ch_pairs, ch_germline_resource, ch_germline_tbi, ch_intervals)
    PILEUP_NORMAL_PAIR(ch_pairs, ch_germline_resource, ch_germline_tbi, ch_intervals)

    ch_pileup_combined = PILEUP_TUMOR_PAIR.out.table
        .join(PILEUP_NORMAL_PAIR.out.table, by: 0)
        .map { pair, tumor_tbl, normal_tbl -> tuple(pair, tumor_tbl, normal_tbl) }

    CALCULATE_CONTAMINATION_PAIR(ch_pileup_combined)

    ch_for_filter = MUTECT2_PAIR.out.vcf
        .join(CALCULATE_CONTAMINATION_PAIR.out.tables, by: 0)
        .join(LEARN_ORIENTATION_PAIR.out.model, by: 0)
        .map { pair, vcf, tbi, stats, f1r2, contam, seg, orient ->
               tuple(pair, vcf, tbi, stats, contam, seg, orient) }

    FILTER_MUTECT_PAIR(ch_for_filter, ch_genome)
    SELECT_PASS_PAIR(FILTER_MUTECT_PAIR.out.vcf, ch_genome)
    ch_mutect2_pass = SELECT_PASS_PAIR.out.vcf   // (pair, vcf, tbi)

    // ── 5. VarScan2 ──────────────────────────────────────
    ch_varscan_pass = Channel.empty()
    if ('varscan2' in enabled_callers()) {
        VARSCAN_MPILEUP_PAIR(ch_pairs, ch_genome, ch_bed)
        VARSCAN_SOMATIC_PAIR(VARSCAN_MPILEUP_PAIR.out.pileups)
        VARSCAN_MERGE_PAIR(VARSCAN_SOMATIC_PAIR.out.vcfs)
        ch_varscan_pass = VARSCAN_MERGE_PAIR.out.vcf
    }

    // ── 6. VarDict ───────────────────────────────────────
    ch_vardict_pass = Channel.empty()
    if ('vardict' in enabled_callers()) {
        VARDICT_PAIR(ch_pairs, ch_genome, ch_bed)
        ch_vardict_pass = VARDICT_PAIR.out.vcf
    }

    // ── 7. Strelka2 ──────────────────────────────────────
    ch_strelka2_pass = Channel.empty()
    if ('strelka2' in enabled_callers()) {
        PREP_STRELKA_BED(ch_bed)
        ch_bed_gz = PREP_STRELKA_BED.out.bed_gz

        if (params.strelka2_run_manta) {
            MANTA_PAIR(ch_pairs, ch_genome, ch_bed_gz)
            ch_with_manta = ch_pairs
                .join(MANTA_PAIR.out.indels, by: 0)
                .map { pair, t, tb, tbi, n, nb, nbi, ivcf, itbi ->
                       tuple(pair, t, tb, tbi, n, nb, nbi, ivcf, itbi) }
            STRELKA2_PAIR(ch_with_manta, ch_genome, ch_bed_gz)
            STRELKA2_MERGE_PAIR(STRELKA2_PAIR.out.results)
        } else {
            // Default: no Manta (as used in the original analysis)
            STRELKA2_PAIR_DIRECT(ch_pairs, ch_genome, ch_bed_gz)
            STRELKA2_MERGE_PAIR(STRELKA2_PAIR_DIRECT.out.results)
        }
        ch_strelka2_pass = STRELKA2_MERGE_PAIR.out.vcf
    }

    // ── 8. Normalize each caller's PASS VCF ─────────────
    ch_to_norm = ch_mutect2_pass
        .map { pair, vcf, tbi -> tuple(pair, 'mutect2', vcf, tbi) }

    if ('varscan2' in enabled_callers()) {
        ch_to_norm = ch_to_norm.mix(
            ch_varscan_pass.map { pair, vcf, tbi -> tuple(pair, 'varscan2', vcf, tbi) }
        )
    }
    if ('vardict' in enabled_callers()) {
        ch_to_norm = ch_to_norm.mix(
            ch_vardict_pass.map { pair, vcf, tbi -> tuple(pair, 'vardict', vcf, tbi) }
        )
    }
    if ('strelka2' in enabled_callers()) {
        ch_to_norm = ch_to_norm.mix(
            ch_strelka2_pass.map { pair, vcf, tbi -> tuple(pair, 'strelka2', vcf, tbi) }
        )
    }

    NORM_VCF(ch_to_norm, ch_genome)

    // ── 9. Consensus per pair ────────────────────────────
    ch_norm_grouped = NORM_VCF.out.vcf
        .map { pair, caller, vcf, tbi -> tuple(pair, vcf, tbi) }
        .groupTuple(by: 0)
        .map { pair, vcfs, tbis -> tuple(pair, vcfs, tbis) }

    CONSENSUS_ISEC(ch_norm_grouped)

    // ── 10. Per-tumor merge using cross-group pairs ──────
    // Build group map from sample sheet so we can restrict the
    // frequency denominator to cross-group comparisons.
    ch_group_map = Channel
        .fromPath(params.samples)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def grp = row.containsKey('group') ? row.group : row.sample_id
            tuple(row.sample_id, grp)
        }
        .collect()
        .map { rows -> rows.collectEntries { sid, grp -> [sid, grp] } }

    // For each consensus VCF decide if it is a cross-group pair,
    // then groupTuple by tumor_id for those pairs only.
    ch_consensus_by_tumor = CONSENSUS_ISEC.out.vcf
        .combine(ch_group_map)
        .map { pair_name, vcf, tbi, grp_map ->
            def (tumor_id, normal_id) = pair_name.split('_vs_')
            def t_grp = grp_map.getOrDefault(tumor_id,  tumor_id)
            def n_grp = grp_map.getOrDefault(normal_id, normal_id)
            // Include only cross-group pairs in the denominator
            if (t_grp != n_grp) {
                tuple(tumor_id, vcf, tbi)
            } else {
                null
            }
        }
        .filter { it != null }
        .groupTuple(by: 0)
        .map { tumor_id, vcfs, tbis -> tuple(tumor_id, vcfs, tbis) }

    MERGE_PAIR_VCFS(ch_consensus_by_tumor)

    // ── 11. Frequency filter ─────────────────────────────
    FREQUENCY_FILTER(MERGE_PAIR_VCFS.out.vcf)
}
