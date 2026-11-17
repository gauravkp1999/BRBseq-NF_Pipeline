process GENE_ENRICHMENT {
    tag "$contrast_id"
    label 'low'

    // Only needs base R + ggplot2 (no DESeq2 dependency - it reads an already-computed
    // results TSV). Reuses the same verified Bioconductor-DESeq2 container as
    // DESEQ2_DIFFERENTIAL (which is a strict superset: R + ggplot2 are both present)
    // rather than introducing a second, separately-pinned container image for one
    // more package - one less container to keep in sync for a script this small.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/67/6713fae254ef549c08f00905da36abadefa357892de60180b7e2021a87d80c5d/data' :
        'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-limma_bioconductor-apeglm:7fccb58b6c45d084' }"

    input:
    tuple val(contrast_id), path(deseq_results)
    path gmt_file
    val  alpha
    val  rank_metric
    val  min_set_size
    val  max_set_size
    val  nperm
    val  top_n

    output:
    tuple val(contrast_id), path("*.ora_results.tsv"),  emit: ora_results
    tuple val(contrast_id), path("*.ora_dotplot.png"),  emit: ora_plot, optional: true
    tuple val(contrast_id), path("*.gsea_results.tsv"), emit: gsea_results
    tuple val(contrast_id), path("*.gsea_dotplot.png"), emit: gsea_plot, optional: true
    path "versions.yml",                                emit: versions

    script:
    """
    gene_enrichment.R \\
        --deseq_results ${deseq_results} \\
        --gmt_file ${gmt_file} \\
        --id ${contrast_id} \\
        --outdir . \\
        --alpha ${alpha} \\
        --rank_metric ${rank_metric} \\
        --min_set_size ${min_set_size} \\
        --max_set_size ${max_set_size} \\
        --nperm ${nperm} \\
        --top_n ${top_n}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$( R --version | head -n1 | sed 's/R version //; s/ .*//' )
    END_VERSIONS
    """

    stub:
    """
    touch ${contrast_id}.ora_results.tsv
    touch ${contrast_id}.ora_dotplot.png
    touch ${contrast_id}.gsea_results.tsv
    touch ${contrast_id}.gsea_dotplot.png
    touch versions.yml
    """
}
