process MULTIQC {
    label 'low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.25.1--pyhdfd78af_0' :
        'biocontainers/multiqc:1.25.1--pyhdfd78af_0' }"

    input:
    path multiqc_files, stageAs: "?/*"
    path multiqc_config
    val  multiqc_title

    output:
    path "*multiqc_report.html", emit: report
    path "*_data",                emit: data
    path "versions.yml",          emit: versions

    script:
    def title = multiqc_title ? "--title \"${multiqc_title}\"" : ''
    def config = multiqc_config ? "--config ${multiqc_config}" : ''
    """
    multiqc \\
        --force \\
        ${title} \\
        ${config} \\
        -n multiqc_report.html \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html
    mkdir multiqc_report_data
    touch versions.yml
    """
}
