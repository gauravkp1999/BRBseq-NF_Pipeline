#!/usr/bin/env bash
# Stand-alone check of the synthetic IPA dataset (no Nextflow needed): STAR index -> STARsolo (same flags as the
# brbseq STARSOLO_ALIGN module) -> per-sample BAMs for IGV -> per-site UMI counts -> baseline calls -> validation.
# Needs: STAR, samtools, python3 with numpy/scipy/pandas.   Usage: scripts/align_and_check.sh [DATA_DIR] [OUT_DIR] [THREADS]
set -euo pipefail
D=${1:-assets/test_data_ipa}; O=${2:-results_check}; T=${3:-2}; S=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$O/star_index" "$O/bam"; cd "$O"; D=$(cd - >/dev/null && cd "$D" && pwd)
STAR --runMode genomeGenerate --genomeSAindexNbases 8 --genomeDir star_index --genomeFastaFiles "$D/reference/genome.fasta" \
     --sjdbGTFfile "$D/reference/genes.gtf" --runThreadN "$T" --outFileNamePrefix idx_ >/dev/null
for p in plate1 plate2; do
  tail -n +2 "$D/barcode_map_$p.csv" | cut -d, -f2 > wl_$p.txt
  STAR --runMode alignReads --genomeDir star_index --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 12 \
       --soloUMIstart 13 --soloUMIlen 16 --soloCBwhitelist wl_$p.txt --soloCBmatchWLtype 1MM --soloCellFilter None \
       --soloFeatures Gene --soloStrand Forward --soloUMIdedup 1MM_Directional --clipAdapterType CellRanger4 \
       --outFilterMultimapNmax 1 --outSAMmapqUnique 60 --outSAMunmapped Within \
       --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outSAMtype BAM SortedByCoordinate \
       --quantMode GeneCounts --readFilesCommand zcat --runThreadN "$T" --outFileNamePrefix ${p}_ \
       --readFilesIn "$D/${p}_R2.fastq.gz" "$D/${p}_R1.fastq.gz" >/dev/null
  tail -n +2 "$D/barcode_map_$p.csv" | while IFS=, read -r sample bc; do          # split pooled run -> sample BAMs
    samtools view -b -F 4 -d "CB:$bc" ${p}_Aligned.sortedByCoord.out.bam > "bam/$sample.bam"
    samtools index "bam/$sample.bam"
  done
done
python3 "$S/count_site_umis.py" --pas "$D/reference/pas_reference.tsv" -o site_counts.tsv \
        --sample-bams $(for b in bam/*.bam; do n=$(basename "$b" .bam); echo "$n=$b"; done)
python3 "$S/baseline_ipa_caller.py" --counts site_counts.tsv --pas "$D/reference/pas_reference.tsv" \
        --meta "$D/sample_meta.csv" -o baseline_calls.tsv
python3 "$S/validate_ipa.py" baseline_calls.tsv "$D/expected/truth_sites.tsv"
