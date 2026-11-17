#!/usr/bin/env nextflow
/*
========================================================================================
    gauravkp/brbseq
========================================================================================
    Github : https://github.com/gauravkp/brbseq
    BRB-seq (Bulk RNA Barcoding and sequencing) analysis using STARsolo
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { BRBSEQ } from './workflows/brbseq'

workflow {

    main:

    // validate inputs
    if (!params.input) {
        error "Please provide --input samplesheet.csv (columns: run,fastq_1,fastq_2,barcode_map)"
    }
    if (!params.gtf) {
        error "Please provide --gtf /path/to/annotation.gtf"
    }
    if (!params.genome_fasta && !params.star_index) {
        error "Please provide --genome_fasta (to build a STAR index) or --star_index (a prebuilt one)"
    }
    if (!params.outdir) {
        error "Please provide --outdir /path/to/results"
    }

    // Downstream analysis is opt-in and layered: DESeq2 needs BOTH sample_meta and
    // contrasts (one without the other can't define a design), and enrichment needs
    // a completed DESeq2 run to enrich on top of.
    if ((params.sample_meta && !params.contrasts) || (!params.sample_meta && params.contrasts)) {
        error "Please provide BOTH --sample_meta and --contrasts to run differential expression, or neither to skip it"
    }
    if (params.gmt_file && !(params.sample_meta && params.contrasts)) {
        error "Please also provide --sample_meta and --contrasts - enrichment runs on each DESeq2 contrast's results, so DESeq2 must run first"
    }


    /*
    ========================================================================================
        RUN MAIN WORKFLOW
    ========================================================================================
    */

    def REQUIRED_COLUMNS = ['run', 'fastq_1', 'fastq_2', 'barcode_map']


    ch_sheet = Channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def missing = REQUIRED_COLUMNS.findAll { !row.containsKey(it) || !row[it]?.trim() }
            if (missing) {
                error "Samplesheet row is missing required value(s) ${missing}: ${row}"
            }

            def fastq_1     = file(row.fastq_1, checkIfExists: true)
            def fastq_2     = file(row.fastq_2, checkIfExists: true)
            def barcode_map = file(row.barcode_map, checkIfExists: true)

            tuple([id: row.run], fastq_1, fastq_2, barcode_map)
        }

    // fail fast on duplicate run IDs, which would silently clobber each other's outputs
    ch_sheet
        .map { meta, f1, f2, bc -> meta.id }
        .toList()
        .subscribe { ids ->
            def dupes = ids.countBy { it }.findAll { k, v -> v > 1 }.keySet()
            if (dupes) {
                error "Duplicate 'run' value(s) in samplesheet: ${dupes} - each run must be unique"
            }
        }

    BRBSEQ(ch_sheet)

    onComplete:
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Results are in: ${params.outdir}"
}
