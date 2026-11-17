process MERGE_MATRICES {
    label 'low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    path matrices

    output:
    path "brbseq.merged_counts.tsv", emit: merged_counts
    path "versions.yml",             emit: versions

    script:
    """
    merge_matrices.py ${matrices} -o brbseq.merged_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python3 --version | sed 's/Python //g' )
        pandas: \$( python3 -c "import pandas; print(pandas.__version__)" )
    END_VERSIONS
    """

    stub:
    """
    touch brbseq.merged_counts.tsv
    touch versions.yml
    """
}
