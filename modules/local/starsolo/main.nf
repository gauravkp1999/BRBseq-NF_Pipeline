process STARSOLO_ALIGN {
    tag "$meta.id"
    label 'star'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.11b--h5ca1c30_0' :
        'biocontainers/star:2.7.11b--h5ca1c30_0' }"

    input:
    tuple val(meta), path(fastq_1), path(fastq_2), path(whitelist)
    path index

    output:
    tuple val(meta), path("${meta.id}_Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(meta), path("${meta.id}_Log.final.out"),                 emit: log_final
    tuple val(meta), path("${meta.id}_Solo.out"),                      emit: solo_out
    tuple val(meta), path("${meta.id}_Solo.out/Gene/raw"),             emit: raw_matrix
    tuple val(meta), path("${meta.id}_Solo.out/Gene/Summary.csv"),     emit: summary
    path "versions.yml",                                               emit: versions

    script:
    // BRB-seq R1 = barcode(+UMI), R2 = cDNA fragment. STARsolo's --readFilesIn
    // convention is cDNA-read first, barcode-read second (i.e. R2 then R1) -
    // this is NOT a typo, it's how STARsolo expects the arguments.
    def read_cmd = fastq_2.toString().endsWith('.gz') ? '--readFilesCommand zcat' : ''
    """
    STAR \\
        --runMode alignReads \\
        --genomeDir ${index} \\
        --soloType CB_UMI_Simple \\
        --soloCBstart ${params.solo_cb_start} \\
        --soloCBlen ${params.solo_cb_len} \\
        --soloUMIstart ${params.solo_umi_start} \\
        --soloUMIlen ${params.solo_umi_len} \\
        --soloCBwhitelist ${whitelist} \\
        --soloCBmatchWLtype ${params.solo_cb_match_wl_type} \\
        --soloCellFilter None \\
        --soloFeatures Gene \\
        --soloStrand Forward \\
        --soloUMIdedup 1MM_Directional \\
        --clipAdapterType CellRanger4 \\
        --outFilterMultimapNmax 1 \\
        --outSAMmapqUnique 60 \\
        --outSAMunmapped Within \\
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \\
        --outSAMtype BAM SortedByCoordinate \\
        --quantMode GeneCounts \\
        ${read_cmd} \\
        --runThreadN ${task.cpus} \\
        --outFileNamePrefix ${meta.id}_ \\
        --readFilesIn ${fastq_2} ${fastq_1}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$( STAR --version | sed -e "s/STAR_//g" )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}_Solo.out/Gene/raw
    touch ${meta.id}_Aligned.sortedByCoord.out.bam
    touch ${meta.id}_Log.final.out
    touch ${meta.id}_Solo.out/Gene/Summary.csv
    touch ${meta.id}_Solo.out/Gene/raw/matrix.mtx
    touch ${meta.id}_Solo.out/Gene/raw/barcodes.tsv
    touch ${meta.id}_Solo.out/Gene/raw/features.tsv
    touch versions.yml
    """
}
