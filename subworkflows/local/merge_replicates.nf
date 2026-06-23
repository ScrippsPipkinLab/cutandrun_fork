/*
 * Pool biological replicates by averaging their spike-in-normalised tracks,
 * then call peaks and build a coverage track on the pooled "merged library".
 *
 * The bedGraphs produced by PREPARE_PEAKCALLING are already on a common
 * per-cell scale (bedtools genomecov -scale s_i with s_i = C / E_i). Averaging
 * these already-scaled tracks position-by-position preserves the spike-in
 * normalisation while reducing per-replicate noise (variance ~1/n), giving
 * atacseq-style merged-signal peak calling without ever merging raw BAMs (which
 * would discard the normalisation).
 *
 * Merged-signal peak calling is performed with SEACR in "non" mode because the
 * pooled inputs are already normalised - "norm" would re-scale and undo the
 * spike-in correction.
 */

include { BEDTOOLS_UNIONBEDG_MEAN as MERGE_TARGET_BEDGRAPH  } from "../../modules/local/bedtools/unionbedg/main"
include { BEDTOOLS_UNIONBEDG_MEAN as MERGE_CONTROL_BEDGRAPH } from "../../modules/local/bedtools/unionbedg/main"
include { UCSC_BEDGRAPHTOBIGWIG as MERGE_BEDGRAPHTOBIGWIG    } from "../../modules/nf-core/ucsc/bedgraphtobigwig/main"
include { SEACR_CALLPEAK as MERGE_SEACR_IGG                  } from "../../modules/nf-core/seacr/callpeak/main"
include { SEACR_CALLPEAK as MERGE_SEACR_NOIGG                } from "../../modules/nf-core/seacr/callpeak/main"

workflow MERGE_REPLICATES {
    take:
    ch_bedgraph    // channel: [ val(meta), [ bedgraph ] ] (spike-in scaled, all samples)
    ch_chrom_sizes // channel: [ sizes ]
    use_control    // value:   boolean

    main:
    ch_versions = Channel.empty()

    /*
     * CHANNEL: Split into target / control tracks
     */
    ch_bedgraph.branch { it ->
        target:  it[0].is_control == false
        control: it[0].is_control == true
    }
    .set { ch_bedgraph_split }

    /*
     * CHANNEL: Collect target replicate bedgraphs per group
     */
    ch_bedgraph_split.target
        .map { meta, bedgraph -> [ meta.group, meta, bedgraph ] }
        .groupTuple(by: 0)
        .map { group, metas, bedgraphs ->
            def meta_new = [:]
            meta_new.id             = group
            meta_new.group          = group
            meta_new.target         = metas[0].target
            meta_new.control_group  = metas[0].control_group
            meta_new.is_control     = false
            meta_new.single_end     = metas[0].single_end
            meta_new.num_replicates = bedgraphs.size()
            [ meta_new, bedgraphs ]
        }
        .set { ch_target_grouped }
    // EXAMPLE CHANNEL STRUCT: [[id:h3k27me3, group:h3k27me3, target:h3k27me3, ...], [BEDGRAPH1, BEDGRAPH2]]

    /*
     * MODULE: Average target replicate tracks
     */
    MERGE_TARGET_BEDGRAPH ( ch_target_grouped )
    ch_versions = ch_versions.mix(MERGE_TARGET_BEDGRAPH.out.versions)

    /*
     * MODULE: Build a pooled bigWig for the merged target track
     */
    MERGE_BEDGRAPHTOBIGWIG ( MERGE_TARGET_BEDGRAPH.out.bedgraph, ch_chrom_sizes )
    ch_versions = ch_versions.mix(MERGE_BEDGRAPHTOBIGWIG.out.versions)

    ch_merged_peaks = Channel.empty()
    if (use_control) {
        /*
         * CHANNEL: Collect control replicate bedgraphs per control group
         */
        ch_bedgraph_split.control
            .map { meta, bedgraph -> [ meta.group, meta, bedgraph ] }
            .groupTuple(by: 0)
            .map { group, metas, bedgraphs ->
                def meta_new = [:]
                meta_new.id         = group
                meta_new.group      = group
                meta_new.is_control = true
                [ meta_new, bedgraphs ]
            }
            .set { ch_control_grouped }

        /*
         * MODULE: Average control replicate tracks
         */
        MERGE_CONTROL_BEDGRAPH ( ch_control_grouped )
        ch_versions = ch_versions.mix(MERGE_CONTROL_BEDGRAPH.out.versions)

        /*
         * CHANNEL: Pair each merged target with its merged control.
         * The per-replicate control_group is "<control>_<rep>" so a target group
         * matches a merged control whose group is the "<control>" prefix.
         */
        MERGE_TARGET_BEDGRAPH.out.bedgraph
            .combine ( MERGE_CONTROL_BEDGRAPH.out.bedgraph )
            .filter { target_meta, target_bg, control_meta, control_bg ->
                target_meta.control_group == control_meta.group ||
                    ( target_meta.control_group != null && target_meta.control_group.startsWith(control_meta.group + "_") )
            }
            .map { target_meta, target_bg, control_meta, control_bg ->
                [ target_meta, target_bg, control_bg ]
            }
            .set { ch_merged_paired }
        // EXAMPLE CHANNEL STRUCT: [[META], TARGET_BEDGRAPH, CONTROL_BEDGRAPH]

        MERGE_SEACR_IGG (
            ch_merged_paired,
            params.seacr_peak_threshold
        )
        ch_merged_peaks = MERGE_SEACR_IGG.out.bed
        ch_versions     = ch_versions.mix(MERGE_SEACR_IGG.out.versions)
    }
    else {
        /*
         * CHANNEL: Add fake control channel
         */
        MERGE_TARGET_BEDGRAPH.out.bedgraph
            .map { meta, bedgraph -> [ meta, bedgraph, [] ] }
            .set { ch_merged_target_fctrl }

        MERGE_SEACR_NOIGG (
            ch_merged_target_fctrl,
            params.seacr_peak_threshold
        )
        ch_merged_peaks = MERGE_SEACR_NOIGG.out.bed
        ch_versions     = ch_versions.mix(MERGE_SEACR_NOIGG.out.versions)
    }

    emit:
    bedgraph = MERGE_TARGET_BEDGRAPH.out.bedgraph // channel: [ val(meta), [ bedgraph ] ]
    bigwig   = MERGE_BEDGRAPHTOBIGWIG.out.bigwig  // channel: [ val(meta), [ bigwig ]   ]
    peaks    = ch_merged_peaks                    // channel: [ val(meta), [ bed ]      ]
    versions = ch_versions                        // channel: [ versions.yml ]
}
