process EXTRACT_READ_COUNT_MATRIX {
    tag "$meta.id"
    label 'medium'

    // BRB-seqTools has no bioconda/biocontainers package (and is marked deprecated
    // upstream in favour of STARsolo - see docs/usage.md). Its jar is pulled in as a
    // plain `path` input (brbseqtools_jar, staged by Nextflow itself from the pinned
    // GitHub release URL in workflows/brbseq.nf - no curl/wget needed inside this
    // plain JRE container). Off by default (params.extract_read_counts); only
    // produces an EXTRA raw read-count matrix alongside STARsolo's own UMI matrix,
    // which is the pipeline's primary output.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/openjdk:17--haa801f8_0' :
        'biocontainers/openjdk:17--haa801f8_0' }"

    input:
    tuple val(meta), path(bam), path(barcode_sample_map)
    path gtf
    path brbseqtools_jar

    output:
    tuple val(meta), path("${meta.id}_readcounts"), emit: readcounts
    path "versions.yml",                            emit: versions

    script:
    """
    # BRB-seqTools expects its own mapping file as Name<TAB>Barcode (reverse of our
    # internal barcode<TAB>sample format), with a header line matching their example.
    echo -e "Name\\tB1" > brbseqtools_barcodes.txt
    awk -F'\\t' '{print \$2"\\t"\$1}' ${barcode_sample_map} >> brbseqtools_barcodes.txt

    mkdir -p ${meta.id}_readcounts

    java -jar ${brbseqtools_jar} ExtractReadCountMatrix \\
        -b ${bam} \\
        -c brbseqtools_barcodes.txt \\
        -gtf ${gtf} \\
        -o ${meta.id}_readcounts

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$( java -version 2>&1 | head -n1 | sed -e 's/.*version "//' -e 's/".*//' )
        brbseqtools: 1.6.1
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}_readcounts
    touch ${meta.id}_readcounts/placeholder.txt
    touch versions.yml
    """
}
