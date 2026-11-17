/*
========================================================================================
    gauravkp/brbseq: main workflow
========================================================================================
*/

include { FASTQC as FASTQC_R1 }        from '../modules/local/fastqc/main'
include { FASTQC as FASTQC_R2 }        from '../modules/local/fastqc/main'
include { PREPARE_BARCODES }           from '../modules/local/prepare_barcodes/main'
include { STAR_INDEX }                 from '../modules/local/star/main'
include { STARSOLO_ALIGN }             from '../modules/local/starsolo/main'
include { RELABEL_MATRIX }             from '../modules/local/relabel_matrix/main'
include { MERGE_MATRICES }             from '../modules/local/merge_matrices/main'
include { EXTRACT_READ_COUNT_MATRIX }  from '../modules/local/extract_read_count_matrix/main'
include { DESEQ2_DIFFERENTIAL }        from '../modules/local/deseq2_differential/main'
include { GENE_ENRICHMENT }            from '../modules/local/gene_enrichment/main'
include { MULTIQC }                    from '../modules/local/multiqc/main'

workflow BRBSEQ {

    take:
    ch_runs   // channel: [ val(meta), path(fastq_1), path(fastq_2), path(barcode_map) ]

    main:

    ch_versions = Channel.empty()

    // ---- raw QC on both reads ----
    FASTQC_R1(ch_runs.map { meta, f1, f2, bc -> [meta, f1] })
    ch_versions = ch_versions.mix(FASTQC_R1.out.versions.first())

    FASTQC_R2(ch_runs.map { meta, f1, f2, bc -> [meta, f2] })
    ch_versions = ch_versions.mix(FASTQC_R2.out.versions.first())

    // ---- normalize the user's barcode/sample map into a STARsolo whitelist ----
    PREPARE_BARCODES(ch_runs.map { meta, f1, f2, bc -> [meta, bc] })
    ch_versions = ch_versions.mix(PREPARE_BARCODES.out.versions.first())

    // ---- STAR index: build it, or reuse a prebuilt one ----
    if (params.star_index) {
        star_index_ch = Channel.fromPath(params.star_index, checkIfExists: true).first()
    } else {
        STAR_INDEX(
            Channel.fromPath(params.genome_fasta, checkIfExists: true),
            Channel.fromPath(params.gtf, checkIfExists: true)
        )
        star_index_ch = STAR_INDEX.out.index.first()
        ch_versions = ch_versions.mix(STAR_INDEX.out.versions)
    }

    // ---- joint alignment + demultiplexing (STARsolo) ----
    ch_for_align = ch_runs
        .map { meta, f1, f2, bc -> [meta, f1, f2] }
        .join(PREPARE_BARCODES.out.whitelist)
        .map { meta, f1, f2, wl -> tuple(meta, f1, f2, wl) }

    STARSOLO_ALIGN(ch_for_align, star_index_ch)
    ch_versions = ch_versions.mix(STARSOLO_ALIGN.out.versions.first())

    // ---- relabel each run's raw barcode matrix to sample names ----
    ch_for_relabel = STARSOLO_ALIGN.out.raw_matrix
        .join(PREPARE_BARCODES.out.map)
        .map { meta, raw, bcmap -> tuple(meta, raw, bcmap) }

    RELABEL_MATRIX(ch_for_relabel)
    ch_versions = ch_versions.mix(RELABEL_MATRIX.out.versions.first())

    // ---- merge all runs into one gene x sample matrix ----
    MERGE_MATRICES(RELABEL_MATRIX.out.counts.map { meta, counts -> counts }.collect())
    ch_versions = ch_versions.mix(MERGE_MATRICES.out.versions)

    // ---- optional: raw (non-UMI) read-count matrix via BRB-seqTools ----
    ch_readcounts = Channel.empty()
    if (params.extract_read_counts) {
        ch_for_extract = STARSOLO_ALIGN.out.bam
            .join(PREPARE_BARCODES.out.map)
            .map { meta, bam, bcmap -> tuple(meta, bam, bcmap) }

        // Staged by Nextflow itself (https:// path input) so the plain JRE
        // container never needs curl/wget - see modules/local/extract_read_count_matrix/main.nf
        ch_brbseqtools_jar = Channel.fromPath(params.brbseqtools_jar_url, checkIfExists: true).first()

        EXTRACT_READ_COUNT_MATRIX(
            ch_for_extract,
            Channel.fromPath(params.gtf, checkIfExists: true).first(),
            ch_brbseqtools_jar
        )
        ch_versions = ch_versions.mix(EXTRACT_READ_COUNT_MATRIX.out.versions.first())
        ch_readcounts = EXTRACT_READ_COUNT_MATRIX.out.readcounts
    }

    // ---- optional: differential expression (DESeq2) + enrichment, once a
    // sample_meta + contrasts design is available. Both are per-pipeline-run
    // steps (not per-sequencing-run), so they take the final merged matrix. ----
    ch_deseq_results = Channel.empty()
    ch_ora_results   = Channel.empty()
    ch_gsea_results  = Channel.empty()

    if (params.sample_meta && params.contrasts) {
        DESEQ2_DIFFERENTIAL(
            MERGE_MATRICES.out.merged_counts,
            Channel.fromPath(params.sample_meta, checkIfExists: true),
            Channel.fromPath(params.contrasts, checkIfExists: true),
            params.deseq2_alpha,
            params.deseq2_min_count
        )
        ch_versions = ch_versions.mix(DESEQ2_DIFFERENTIAL.out.versions)
        ch_deseq_results = DESEQ2_DIFFERENTIAL.out.results

        // One result TSV per contrast is emitted as a single bundled list from
        // this one-shot process - flatten it back into a [contrast_id, file]
        // channel (id = filename with the .deseq2_results.tsv suffix stripped,
        // matching the `id` column of --contrasts) so enrichment can run per
        // contrast, same as any other per-sample/per-unit nf-core channel.
        ch_per_contrast_results = ch_deseq_results
            .flatten()
            .map { f -> tuple(f.name.replace('.deseq2_results.tsv', ''), f) }

        if (params.gmt_file) {
            GENE_ENRICHMENT(
                ch_per_contrast_results,
                Channel.fromPath(params.gmt_file, checkIfExists: true).first(),
                params.enrichment_alpha,
                params.enrichment_rank_metric,
                params.enrichment_min_set_size,
                params.enrichment_max_set_size,
                params.enrichment_nperm,
                params.enrichment_top_n
            )
            ch_versions = ch_versions.mix(GENE_ENRICHMENT.out.versions.first())
            ch_ora_results  = GENE_ENRICHMENT.out.ora_results
            ch_gsea_results = GENE_ENRICHMENT.out.gsea_results
        }
    }

    // ---- aggregate QC ----
    ch_multiqc_files = Channel.empty()
        .mix(FASTQC_R1.out.reports.map { meta, reports -> reports })
        .mix(FASTQC_R2.out.reports.map { meta, reports -> reports })
        .mix(STARSOLO_ALIGN.out.log_final.map { meta, log -> log })
        .mix(STARSOLO_ALIGN.out.summary.map { meta, summary -> summary })
        .mix(ch_versions.collectFile(name: 'collated_versions.yml'))
        .collect()

    MULTIQC(
        ch_multiqc_files,
        Channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true).first(),
        params.multiqc_title ?: ''
    )

    emit:
    merged_counts  = MERGE_MATRICES.out.merged_counts
    readcounts     = ch_readcounts
    deseq_results  = ch_deseq_results
    ora_results    = ch_ora_results
    gsea_results   = ch_gsea_results
    multiqc_report = MULTIQC.out.report
    versions       = ch_versions
}
