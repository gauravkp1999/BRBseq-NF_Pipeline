process PREPARE_BARCODES {
    tag "$meta.id"
    label 'low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'biocontainers/python:3.10' }"

    input:
    tuple val(meta), path(barcode_map)

    output:
    tuple val(meta), path("whitelist.txt"),            emit: whitelist
    tuple val(meta), path("barcode_sample_map.tsv"),   emit: map
    path "versions.yml",                               emit: versions

    script:
    """
    parse_barcode_map.py ${barcode_map} -o .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python3 --version | sed 's/Python //g' )
    END_VERSIONS
    """

    stub:
    """
    touch whitelist.txt
    touch barcode_sample_map.tsv
    touch versions.yml
    """
}
