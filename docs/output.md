# gauravkp/brbseq: Output

## Introduction

This document describes the output produced by the pipeline, under the
directory given by `--outdir`.

## Pipeline overview

1. [FastQC](#fastqc) - QC on raw R1 (barcode/UMI) and R2 (cDNA) reads
2. [STARsolo](#starsolo) - joint alignment, demultiplexing and UMI counting, per run
3. [Relabelled matrix](#relabelled-matrix) - per-run barcode -> sample name matrix
4. [Merged matrix](#merged-matrix) - combined gene x sample matrix across all runs
5. [(Optional) Read-count matrix](#optional-read-count-matrix)
6. [(Optional) Differential expression](#optional-differential-expression) - DESeq2, one or more contrasts
7. [(Optional) Gene-set enrichment](#optional-gene-set-enrichment) - ORA + GSEA per contrast
8. [MultiQC](#multiqc) - aggregate QC report
9. [Pipeline information](#pipeline-information)

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/R1/` - FastQC report on each run's raw R1 (barcode/UMI read)
- `fastqc/R2/` - FastQC report on each run's raw R2 (cDNA read)

</details>

Note that standard FastQC metrics (GC content, sequence complexity, etc.)
are not very meaningful for R1, since it's a structured barcode+UMI
construct rather than biological sequence - R1 FastQC is mainly useful for
checking read length and overall quality scores.

### STARsolo

<details markdown="1">
<summary>Output files</summary>

- `starsolo/<run>/`
  - `<run>_Aligned.sortedByCoord.out.bam` - coordinate-sorted, CB/UB-tagged BAM
  - `<run>_Log.final.out` - alignment summary statistics
  - `<run>_Solo.out/Gene/Summary.csv` - STARsolo's own per-run summary (reads with valid barcode, saturation, etc.)
  - `<run>_Solo.out/Gene/raw/` - raw gene x **barcode** matrix (`matrix.mtx`, `barcodes.tsv`, `features.tsv`) - barcodes here are raw sequences, not yet sample names

</details>

[STARsolo](https://github.com/alexdobin/STAR) aligns the cDNA read (R2),
assigns each read to a sample barcode and gene using the barcode/UMI
construct in R1, and deduplicates by UMI.

### Relabelled matrix

<details markdown="1">
<summary>Output files</summary>

- `counts/per_run/<run>.counts.tsv` - gene x **sample name** matrix for that one run, with unmapped/unused barcodes dropped

</details>

### Merged matrix

<details markdown="1">
<summary>Output files</summary>

- `counts/brbseq.merged_counts.tsv` - **the main output**: combined gene x sample matrix across every run in the samplesheet (outer join on gene ID, missing values filled with 0)

</details>

### (Optional) Read-count matrix

Only produced if `--extract_read_counts true`.

<details markdown="1">
<summary>Output files</summary>

- `counts/per_run_readcounts/<run>_readcounts/` - BRB-seqTools `ExtractReadCountMatrix` output (raw, non-UMI-deduplicated counts). Contents are tool-generated; inspect this folder directly as the exact filename isn't pinned in this pipeline.

</details>

### (Optional) Differential expression

Only produced if both `--sample_meta` and `--contrasts` are provided. Runs
once on the final merged matrix (not per run).

<details markdown="1">
<summary>Output files</summary>

- `differential_expression/`
  - `<contrast_id>.deseq2_results.tsv` - one per row of `--contrasts`: `gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, significant`
  - `<contrast_id>.volcano.png` - log2FoldChange vs. -log10(padj), colored by direction, top 10 genes labelled
  - `<contrast_id>.ma_plot.png` - DESeq2's standard MA plot
  - `normalized_counts.tsv` - size-factor-normalized counts (shared across all contrasts, from one combined-design DESeq run)
  - `size_factors.tsv` - per-sample DESeq2 size factors
  - `pca.png` - PCA of variance-stabilized counts, colored by the first contrast variable
  - `sample_distance_heatmap.png` - sample-to-sample distance heatmap (VST counts)
  - `dds.rds` - the saved `DESeqDataSet` R object, for reuse/debugging in your own R session

</details>

`significant` = `padj < --deseq2_alpha` (default 0.05). Genes with fewer than
`--deseq2_min_count` total counts across all samples (default 10) are
excluded before fitting. Each contrast fits its own single-variable DESeq2
design (`~ variable`) - see [`docs/usage.md`](usage.md) if you need a
multi-covariate design.

### (Optional) Gene-set enrichment

Only produced if `--gmt_file` is also provided. Runs once per DESeq2
contrast.

<details markdown="1">
<summary>Output files</summary>

- `enrichment/<contrast_id>/`
  - `<contrast_id>.ora_results.tsv` - over-representation analysis: `gene_set, description, set_size, overlap_size, n_significant, universe_size, overlap_genes, pvalue, padj`
  - `<contrast_id>.ora_dotplot.png` - top gene sets by ORA p-value (point size = overlap size, color = padj)
  - `<contrast_id>.gsea_results.tsv` - preranked GSEA: `gene_set, description, size, ES, NES, pvalue, padj`
  - `<contrast_id>.gsea_dotplot.png` - top gene sets by normalized enrichment score (NES)

</details>

ORA tests genes significant at `--enrichment_alpha` (default 0.05) against
all tested genes as background (one-sided hypergeometric test). GSEA ranks
all tested genes by `--enrichment_rank_metric` (default `stat`) and uses a
gene-set-permutation null (`--enrichment_nperm`, default 1000 permutations).
Both are independent from-scratch implementations (no clusterProfiler/fgsea
dependency) - see the header comment in `bin/gene_enrichment.R` for the exact
statistics used.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/multiqc_report.html` - aggregated HTML report - **start here**
- `multiqc/multiqc_report_data/` - parsed data behind the report

</details>

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/` - `execution_report_*.html`, `execution_timeline_*.html`, `execution_trace_*.txt`, `pipeline_dag_*.html`

</details>
