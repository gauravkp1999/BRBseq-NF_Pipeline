process RELABEL_MATRIX {
    tag "$meta.id"
    label 'low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(raw_matrix_dir), path(barcode_map)

    output:
    tuple val(meta), path("${meta.id}.counts.tsv"), emit: counts
    path "versions.yml",                            emit: versions

    script:
    """
    relabel_matrix.py \\
        --matrix ${raw_matrix_dir}/matrix.mtx \\
        --barcodes ${raw_matrix_dir}/barcodes.tsv \\
        --features ${raw_matrix_dir}/features.tsv \\
        --barcode-map ${barcode_map} \\
        --run ${meta.id} \\
        --min-reads ${params.min_barcode_reads} \\
        -o ${meta.id}.counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python3 --version | sed 's/Python //g' )
        pandas: \$( python3 -c "import pandas; print(pandas.__version__)" )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.counts.tsv
    touch versions.yml
    """
}
