# gauravkp/brbseq

[![CI](https://github.com/gauravkp/brbseq/actions/workflows/ci.yml/badge.svg)](https://github.com/gauravkp/brbseq/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3eb049.svg)](https://docs.conda.io/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Introduction

**gauravkp/brbseq** is a Nextflow (DSL2) pipeline for **BRB-seq** (Bulk RNA
Barcoding and sequencing) analysis, built to nf-core pipeline conventions
(practice exercise toward an official nf-core submission).

BRB-seq pools many samples into one sequencing run using an in-read sample
barcode + UMI (on R1), with only the cDNA fragment (R2) needing alignment.
This pipeline uses **STARsolo** to align and demultiplex in one step, which
is what the [BRB-seqTools](https://github.com/DeplanckeLab/BRB-seqTools)
maintainers themselves now recommend over their own (deprecated) toolkit.

```
FastQC (R1, R2) → STARsolo (align + demultiplex + UMI-count) → relabel matrix (barcode → sample) → merge across runs → [DESeq2 → ORA/GSEA enrichment] → MultiQC
```

## Quick start

1. Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) (`>=23.04.0`)
2. Install [Docker](https://docs.docker.com/engine/installation/) or [Singularity](https://sylabs.io/guides/3.0/user-guide/installation.html) (or use `-profile conda`)
3. **Know your BRB-seq construct** - the barcode/UMI positions in R1 vary by
   oligo design/protocol version. Check `docs/usage.md` before trusting any
   output.
4. Test it runs on the built-in minimal dataset:

   ```bash
   nextflow run main.nf -profile test,docker --outdir results
   ```

5. Run with your own data:

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

6. Optionally add differential expression + enrichment on top:

   ```bash
   --sample_meta sample_meta.csv --contrasts contrasts.csv --gmt_file genesets.gmt
   ```

See [`docs/usage.md`](docs/usage.md) for full usage (including the
samplesheet, barcode-map, sample-metadata/contrasts and GMT formats) and
[`docs/output.md`](docs/output.md) for every output file.

## Samplesheet

One row per pooled **sequencing run**, not per biological sample:

```csv
run,fastq_1,fastq_2,barcode_map
plate1,/data/plate1_R1.fastq.gz,/data/plate1_R2.fastq.gz,/data/plate1_barcodes.csv
plate2,/data/plate2_R1.fastq.gz,/data/plate2_R2.fastq.gz,/data/plate2_barcodes.csv
```

`barcode_map` is a simple sample/barcode CSV or TSV (see
[`assets/barcode_map_example.csv`](assets/barcode_map_example.csv)) - the
pipeline derives the STARsolo whitelist and the final matrix's sample names
from this one file.

## Repository layout

```
brbseq/
├── main.nf                       # entry point + samplesheet validation
├── nextflow.config               # params, profiles, manifest
├── nextflow_schema.json          # parameter validation / auto-generated docs
├── conf/
│   ├── base.config
│   ├── modules.config
│   └── test.config
├── workflows/
│   └── brbseq.nf                 # main workflow logic
├── modules/local/
│   ├── fastqc/
│   ├── prepare_barcodes/         # normalizes user's barcode map -> whitelist + map
│   ├── star/                     # genome indexing
│   ├── starsolo/                 # joint align + demultiplex + UMI count
│   ├── relabel_matrix/           # barcode columns -> sample name columns
│   ├── merge_matrices/           # combine all runs into one matrix
│   ├── extract_read_count_matrix/ # optional: raw read-count matrix (BRB-seqTools)
│   ├── deseq2_differential/      # optional: DESeq2 differential expression
│   ├── gene_enrichment/          # optional: ORA + GSEA enrichment per contrast
│   └── multiqc/
├── bin/
│   ├── parse_barcode_map.py
│   ├── relabel_matrix.py
│   ├── merge_matrices.py
│   ├── deseq2_differential.R
│   └── gene_enrichment.R
├── assets/
│   ├── multiqc_config.yml
│   ├── samplesheet_example.csv
│   ├── barcode_map_example.csv
│   ├── sample_meta_example.csv
│   ├── contrasts_example.csv
│   ├── gmt_example.gmt
│   └── test_data/                # placeholder — add tiny test files here
├── docs/
│   ├── usage.md
│   └── output.md
├── .github/
│   ├── workflows/ci.yml
│   ├── workflows/linting.yml
│   ├── ISSUE_TEMPLATE/
│   └── PULL_REQUEST_TEMPLATE.md
├── .nf-core.yml
├── .editorconfig
├── .gitignore
├── LICENSE                       # MIT
├── CHANGELOG.md
├── CITATIONS.md
├── CODE_OF_CONDUCT.md
├── CONTRIBUTING.md
└── run_pipeline.sh               # example SLURM submission script
```

## Requirements

- Nextflow ≥ 23.04
- Docker, Singularity/Apptainer, **or** an existing conda environment (`--conda_env`) with `fastqc`, `star`, `pandas`, `numpy`, `multiqc` installed, plus `bioconductor-deseq2` (which brings `ggplot2`) if you'll run `--sample_meta`/`--contrasts` and/or `--gmt_file`

## Why STARsolo and not BRB-seqTools?

The BRB-seqTools authors themselves now mark their toolkit deprecated and
recommend STARsolo directly - it produces similar results in one command,
without BRB-seqTools' trimming/demultiplexing quirks (their own README warns
against feeding their Trim output into STARsolo). This pipeline follows that
recommendation for its primary (UMI) quantification path, and only calls out
to BRB-seqTools optionally, for the one thing STARsolo doesn't produce: a
non-UMI-deduplicated read-count matrix (`--extract_read_counts true`).

## Credits

Written by Gaurav Patil. See [`CITATIONS.md`](CITATIONS.md) for the tools
and the original BRB-seq method paper.

## What's still needed for a *real* nf-core submission

- Open a proposal on [nf-core/proposals](https://github.com/nf-core/proposals) **before** building further
- Add real CI test data under `assets/test_data/` (or host it on `nf-core/test-datasets`)
- Validate the samplesheet/barcode-map formats against nf-core's `nf-schema` plugin instead of the current inline Groovy checks
- Get `nf-core pipelines lint` fully passing
- Decide whether `extract_read_counts`/BRB-seqTools stays, given its deprecated upstream status and runtime jar download
- Community/core-team review

## License

[MIT](LICENSE)
