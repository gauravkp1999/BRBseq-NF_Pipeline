#!/bin/bash
#SBATCH --job-name=brbseq_nf
#SBATCH --output=nf_master_%j.log
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=standard

# This job just runs the Nextflow *head* process, which in turn submits
# one SLURM job per pipeline task (fastqc, STAR index, STARsolo, etc).

set -euo pipefail

nextflow run main.nf \
    -profile slurm,singularity \
    --input        samplesheet.csv \
    --genome_fasta reference/genome.fasta \
    --gtf          reference/genomic.gtf \
    --outdir       results \
    --solo_cb_start 1 --solo_cb_len 12 \
    --solo_umi_start 13 --solo_umi_len 16 \
    -resume
