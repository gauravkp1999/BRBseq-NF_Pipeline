process DESEQ2_DIFFERENTIAL {
    label 'medium'

    // This custom script deliberately depends on nothing but DESeq2 itself
    // (DESeq2's own DESCRIPTION lists ggplot2 as a real Import, so it's
    // guaranteed present here too) - no optparse/ggrepel/pheatmap, which are
    // NOT guaranteed in this container and would silently break under
    // docker/singularity while still working under -profile conda. Container
    // reused as-is from nf-core/modules' own verified deseq2/differential
    // module (https://github.com/nf-core/modules/tree/master/modules/nf-core/deseq2/differential).
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/67/6713fae254ef549c08f00905da36abadefa357892de60180b7e2021a87d80c5d/data' :
        'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-limma_bioconductor-apeglm:7fccb58b6c45d084' }"

    input:
    path counts
    path sample_meta
    path contrasts
    val  alpha
    val  min_count

    output:
    path "*.deseq2_results.tsv",       emit: results
    path "*.volcano.png",              emit: volcano_plots
    path "*.ma_plot.png",              emit: ma_plots
    path "normalized_counts.tsv",      emit: normalized_counts
    path "pca.png",                    emit: pca_plot
    path "sample_distance_heatmap.png", emit: sample_distance_plot
    path "dds.rds",                    emit: dds_rds
    path "versions.yml",               emit: versions

    script:
    """
    deseq2_differential.R \\
        --counts ${counts} \\
        --sample_meta ${sample_meta} \\
        --contrasts ${contrasts} \\
        --alpha ${alpha} \\
        --min_count ${min_count} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$( R --version | head -n1 | sed 's/R version //; s/ .*//' )
        bioconductor-deseq2: \$( Rscript -e 'cat(as.character(packageVersion("DESeq2")))' )
    END_VERSIONS
    """

    stub:
    """
    touch example_contrast.deseq2_results.tsv
    touch example_contrast.volcano.png
    touch example_contrast.ma_plot.png
    touch normalized_counts.tsv
    touch pca.png
    touch sample_distance_heatmap.png
    touch dds.rds
    touch versions.yml
    """
}
