process STAR_INDEX {
    label 'star'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.11b--h5ca1c30_0' :
        'biocontainers/star:2.7.11b--h5ca1c30_0' }"

    input:
    path fasta
    path gtf

    output:
    path "star_index",   emit: index
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p star_index
    STAR \\
        --runMode genomeGenerate \\
        ${args} \\
        --genomeDir star_index \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
        --runThreadN ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$( STAR --version | sed -e "s/STAR_//g" )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p star_index
    touch star_index/dummy
    touch versions.yml
    """
}
