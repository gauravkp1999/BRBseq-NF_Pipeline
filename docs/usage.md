# gauravkp/brbseq: Usage

## Introduction

`brbseq` analyzes BRB-seq (Bulk RNA Barcoding and sequencing) data: many
samples pooled into one or more sequencing runs, demultiplexed by an
in-read barcode, with a UMI for deduplication. It uses **STARsolo** to
jointly align and demultiplex each run in a single step, which is what the
[BRB-seqTools maintainers themselves now recommend](https://github.com/DeplanckeLab/BRB-seqTools)
over their own (deprecated) toolkit.

```
FastQC (R1, R2) → STARsolo (align + demultiplex + UMI-count) → relabel matrix (barcode -> sample) → merge across runs → [DESeq2 → ORA/GSEA enrichment] → MultiQC
```

## Before you run this: know your construct

BRB-seq oligo designs differ between protocol versions and labs. This
pipeline does NOT guess your barcode/UMI layout - you must set it explicitly
via these params (defaults match the published BRB-seq example construct:
12bp barcode at the start of R1, followed immediately by a 16bp UMI):

| Parameter | Meaning |
|---|---|
| `--solo_cb_start` | 1-based start position of the sample barcode in R1 |
| `--solo_cb_len` | Barcode length (bp) |
| `--solo_umi_start` | 1-based start position of the UMI in R1 |
| `--solo_umi_len` | UMI length (bp); set to 0 if your library has no UMI |

**Open a few R1 reads and check** (`zcat run_R1.fastq.gz \| head`) that these
positions actually match what you see, before trusting any results.

## Input: the samplesheet

`--input samplesheet.csv`, comma-separated, one row per **sequencing
run/pool** (not per biological sample - BRB-seq multiplexes many samples
into each run):

```csv
run,fastq_1,fastq_2,barcode_map
plate1,/data/plate1_R1.fastq.gz,/data/plate1_R2.fastq.gz,/data/plate1_barcodes.csv
plate2,/data/plate2_R1.fastq.gz,/data/plate2_R2.fastq.gz,/data/plate2_barcodes.csv
```

- `run` - unique ID for this pooled library/sequencing run
- `fastq_1` - R1 FASTQ: barcode + UMI reads
- `fastq_2` - R2 FASTQ: cDNA fragment reads (same read count/order as R1)
- `barcode_map` - a CSV or TSV mapping each barcode used in this run to a
  sample name (see below). **Sample names must be unique across the whole
  samplesheet**, even across different runs.

### Barcode map format

Comma- or tab-delimited, with a header row. The pipeline looks for a sample
column (`sample`, `samplename`, `sample_name`, or `name`) and a barcode
column (`barcode`, `bc`, `b1`, or `cb`), in any order:

```csv
sample,barcode
sampleA,AACCGGTTAACC
sampleB,TTGGCCAATTGG
sampleC,GGAATTCCGGAA
```

The pipeline derives both the STARsolo whitelist and the final matrix
column names from this one file - you don't need to prepare a separate
whitelist yourself.

## Running the pipeline

### 1. Local test run

```bash
nextflow run main.nf -profile test,docker --outdir results
```

### 2. Full run

```bash
nextflow run main.nf \
    -profile docker \
    --input        samplesheet.csv \
    --genome_fasta ref/genome.fa \
    --gtf          ref/annotation.gtf \
    --outdir       results \
    --solo_cb_start 1 --solo_cb_len 12 \
    --solo_umi_start 13 --solo_umi_len 16 \
    -resume
```

### 3. On a SLURM cluster

```bash
nextflow run main.nf \
    -profile slurm,singularity \
    --input        samplesheet.csv \
    --genome_fasta ref/genome.fa \
    --gtf          ref/annotation.gtf \
    --outdir       results \
    -resume
```

or submit `run_pipeline.sh` with `sbatch` (edit the paths/params inside first).

### 4. Reuse a prebuilt STAR index

```bash
--star_index /path/to/existing/star_index
```

### 5. Also get a raw (non-UMI) read-count matrix

```bash
--extract_read_counts true
```

Off by default. This runs BRB-seqTools' `ExtractReadCountMatrix` on the
STARsolo BAM to recover a non-UMI-deduplicated read count matrix (STARsolo
itself only outputs the UMI-deduplicated matrix). BRB-seqTools is marked
deprecated upstream, and this step stages its jar from the GitHub release URL
in `--brbseqtools_jar_url` (fetched by Nextflow itself, not by a tool inside
the container - override this param to point at a local jar if your execution
environment has no outbound network access). Treat its output as experimental
and inspect it before relying on it; see the comments in
`modules/local/extract_read_count_matrix/main.nf`.

### 6. Differential expression (DESeq2)

Off unless you provide **both** `--sample_meta` and `--contrasts`. Runs once,
on the final merged gene x sample matrix (`counts/brbseq.merged_counts.tsv`),
not per sequencing run.

```bash
nextflow run main.nf \
    -profile docker \
    --input samplesheet.csv --genome_fasta ref/genome.fa --gtf ref/annotation.gtf \
    --outdir results \
    --sample_meta sample_meta.csv \
    --contrasts   contrasts.csv \
    -resume
```

**`--sample_meta`** - one row per sample, a `sample` column plus one column
per covariate used in `--contrasts` (see
[`assets/sample_meta_example.csv`](../assets/sample_meta_example.csv)):

```csv
sample,condition
sampleA,control
sampleB,control
sampleC,treated
sampleD,treated
```

**`--contrasts`** - one row per comparison to run (see
[`assets/contrasts_example.csv`](../assets/contrasts_example.csv)):

```csv
id,variable,reference,target
treated_vs_control,condition,control,treated
```

- `variable` must be a column in `--sample_meta`; `reference`/`target` must be
  two of its values. `reference` is the baseline (denominator); `target` is
  compared against it (log2FoldChange > 0 means higher in `target`).
- Each contrast fits its own single-variable design (`~ variable`). If you
  need to control for a batch effect, edit `design_formula` in
  `bin/deseq2_differential.R` directly (multi-covariate designs aren't
  auto-detected from the samplesheet).
- Genes with fewer than `--deseq2_min_count` (default 10) total counts across
  all samples are dropped before fitting. Significance uses `--deseq2_alpha`
  (default 0.05) on `padj` (Benjamini-Hochberg).

Outputs land in `<outdir>/differential_expression/` - see
[`docs/output.md`](output.md).

### 7. Gene-set enrichment (ORA + GSEA)

Off unless you also provide `--gmt_file` (on top of `--sample_meta` /
`--contrasts` - enrichment runs per DESeq2 contrast, so DESeq2 must run
first):

```bash
--sample_meta sample_meta.csv --contrasts contrasts.csv --gmt_file genesets.gmt
```

`--gmt_file` is a standard GMT file - one gene set per line, tab-separated
(see [`assets/gmt_example.gmt`](../assets/gmt_example.gmt)):

```
SET_NAME<TAB>description<TAB>gene1<TAB>gene2<TAB>gene3...
```

**The gene IDs in your GMT must use the same namespace as your GTF's
`gene_id` attribute** (Ensembl IDs, gene symbols, etc. - whatever
featureCounts/STARsolo used to label rows in the count matrix). This pipeline
is organism-agnostic by design: it doesn't bundle or depend on any
species-specific annotation package (no `org.*.eg.db`), so you supply
whatever gene sets fit your organism and question - e.g. exported from
[MSigDB](https://www.gsea-msigdb.org/gsea/msigdb) (human/mouse),
[g:Profiler](https://biit.cs.ut.ee/gprofiler/gost) (many organisms), or your
own GO/KEGG term-to-gene mapping.

Two independently-implemented methods run on every contrast (not wrappers
around clusterProfiler/fgsea, to avoid a heavy Bioconductor annotation
dependency chain - see the comments in `bin/gene_enrichment.R`):

- **ORA** (over-representation analysis) - one-sided hypergeometric test on
  genes significant at `--enrichment_alpha` (default 0.05) vs. all tested
  genes as background.
- **GSEA** (preranked, weighted running-sum) - ranks all tested genes by
  `--enrichment_rank_metric` (`stat` by default, or `log2FoldChange`), with a
  gene-set-permutation null distribution (`--enrichment_nperm`, default
  1000).

Gene sets are restricted to `--enrichment_min_set_size`-`--enrichment_max_set_size`
genes (default 5-500) after overlapping with the tested gene universe.
Outputs land in `<outdir>/enrichment/<contrast_id>/`.

## Core Nextflow arguments

### `-profile`

Comma-separated, e.g. `-profile test,docker`.

- `docker` / `singularity` - containerized execution (recommended)
- `conda` - use an existing conda environment (`--conda_env`)
- `slurm` / `pbs` - submit each process as its own cluster job
- `test` - run the built-in minimal test dataset

### `-resume`

Resume from cached results. Always safe to add.

## Parameters

See [`nextflow_schema.json`](../nextflow_schema.json) for the full list with
defaults and descriptions.

## A note on what "sample" means here

STARsolo produces one matrix **per run**, with raw barcode sequences as
column headers - it has no concept of your sample names. This pipeline's
`RELABEL_MATRIX` step is what turns those barcode columns into your actual
sample names (and drops any barcode in STARsolo's output that isn't in your
map - i.e. unused/noise barcodes), and `MERGE_MATRICES` then combines every
run's relabeled matrix into one gene x sample table at
`<outdir>/counts/brbseq.merged_counts.tsv`.
