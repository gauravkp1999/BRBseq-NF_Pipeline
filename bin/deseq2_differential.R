#!/usr/bin/env Rscript
# Differential expression on a BRB-seq gene x sample count matrix using DESeq2.
#
# Deliberately depends on nothing but DESeq2 itself (which pulls in ggplot2 as
# a real Import - see its DESCRIPTION) plus base R: no optparse, ggrepel or
# pheatmap, so this runs in the same minimal DESeq2 container the module uses
# for both conda and docker/singularity profiles without any mismatch risk.
#
# Inputs:
#   --counts       merged gene x sample count matrix (TSV, first column gene_id) -
#                   i.e. the brbseq pipeline's counts/brbseq.merged_counts.tsv
#   --sample_meta  CSV/TSV: sample,<covariate columns...> - one row per sample
#                   name exactly matching a column in --counts
#   --contrasts    CSV/TSV: id,variable,reference,target - one row per comparison
#                   to run (variable/reference/target must match a column and two
#                   of its values in --sample_meta)
#   --outdir       output directory (created if missing)
#   --alpha        significance threshold on padj for the "significant" flag and
#                   volcano plot highlighting (default 0.05)
#   --min_count    drop genes with fewer than this many total counts across all
#                   samples before fitting (default 10)
#
# Outputs (per contrast <id>):
#   <id>.deseq2_results.tsv  - gene_id, baseMean, log2FoldChange, lfcSE, stat,
#                               pvalue, padj, significant
#   <id>.volcano.png
#   <id>.ma_plot.png
# Outputs (once, across all samples used by any contrast):
#   normalized_counts.tsv
#   size_factors.tsv
#   pca.png
#   sample_distance_heatmap.png
#   dds.rds                  - saved DESeqDataSet, for reuse/debugging

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

## ---- minimal base-R CLI parser (no optparse dependency) -------------------

arg_spec <- list(
  counts      = list(default = NULL,  type = "character"),
  sample_meta = list(default = NULL,  type = "character"),
  contrasts   = list(default = NULL,  type = "character"),
  outdir      = list(default = ".",   type = "character"),
  alpha       = list(default = 0.05,  type = "double"),
  min_count   = list(default = 10,    type = "integer")
)

parse_cli_args <- function(spec) {
  raw <- commandArgs(trailingOnly = TRUE)
  out <- lapply(spec, function(s) s$default)
  i <- 1
  while (i <= length(raw)) {
    token <- raw[i]
    if (!startsWith(token, "--")) stop(sprintf("Unexpected argument: %s", token), call. = FALSE)
    key <- sub("^--", "", token)
    if (grepl("=", key, fixed = TRUE)) {
      parts <- strsplit(key, "=", fixed = TRUE)[[1]]
      key <- parts[1]
      val <- paste(parts[-1], collapse = "=")
      i <- i + 1
    } else {
      if (i + 1 > length(raw)) stop(sprintf("Missing value for --%s", key), call. = FALSE)
      val <- raw[i + 1]
      i <- i + 2
    }
    if (!key %in% names(spec)) stop(sprintf("Unknown argument: --%s", key), call. = FALSE)
    type <- spec[[key]]$type
    out[[key]] <- switch(type,
      double = as.double(val), integer = as.integer(val), logical = as.logical(val), val
    )
  }
  out
}

opt <- parse_cli_args(arg_spec)
for (required in c("counts", "sample_meta", "contrasts")) {
  if (is.null(opt[[required]])) stop(sprintf("ERROR: --%s is required", required), call. = FALSE)
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

read_delim_auto <- function(path) {
  first_line <- readLines(path, n = 1)
  sep <- if (grepl("\t", first_line)) "\t" else ","
  read.delim(path, sep = sep, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

## ---- load inputs ----------------------------------------------------------

counts_raw <- read_delim_auto(opt$counts)
gene_col <- colnames(counts_raw)[1]
counts_mat <- as.matrix(counts_raw[, -1, drop = FALSE])
rownames(counts_mat) <- counts_raw[[gene_col]]
storage.mode(counts_mat) <- "integer"

sample_meta <- read_delim_auto(opt$sample_meta)
if (!"sample" %in% colnames(sample_meta)) {
  stop("ERROR: --sample_meta must have a 'sample' column", call. = FALSE)
}
rownames(sample_meta) <- sample_meta$sample

contrasts_df <- read_delim_auto(opt$contrasts)
required_contrast_cols <- c("id", "variable", "reference", "target")
missing_cols <- setdiff(required_contrast_cols, colnames(contrasts_df))
if (length(missing_cols) > 0) {
  stop(sprintf("ERROR: --contrasts is missing column(s): %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
}

## ---- align samples between counts and metadata ----------------------------

missing_samples <- setdiff(sample_meta$sample, colnames(counts_mat))
if (length(missing_samples) > 0) {
  stop(sprintf(
    "ERROR: sample(s) in --sample_meta not found as columns in --counts: %s",
    paste(missing_samples, collapse = ", ")
  ), call. = FALSE)
}

extra_samples <- setdiff(colnames(counts_mat), sample_meta$sample)
if (length(extra_samples) > 0) {
  message(sprintf(
    "NOTE: %d sample column(s) in --counts have no row in --sample_meta and will be dropped: %s",
    length(extra_samples), paste(extra_samples, collapse = ", ")
  ))
}

counts_mat <- counts_mat[, sample_meta$sample, drop = FALSE]

## ---- filter low-count genes once, shared by every contrast ---------------
## Each contrast's variable becomes a factor column; multi-variable designs
## (e.g. "~ batch + condition") are not auto-detected - this script fits a
## fresh single-variable design (~ variable) per contrast, which is the
## correct, safe default when you don't yet know your covariate structure.
## If you need to control for a batch effect, edit `design_formula` below.

keep_genes <- rowSums(counts_mat) >= opt$min_count
message(sprintf(
  "Keeping %d / %d genes with total count >= %d across all samples",
  sum(keep_genes), nrow(counts_mat), opt$min_count
))
counts_filt <- counts_mat[keep_genes, , drop = FALSE]

## ---- a single DESeq run (on the union of all contrast variables) for the
## shared QC outputs (PCA, sample-distance heatmap, normalized counts) ------

all_vars <- unique(contrasts_df$variable)
missing_vars <- setdiff(all_vars, colnames(sample_meta))
if (length(missing_vars) > 0) {
  stop(sprintf(
    "ERROR: --contrasts references variable(s) not present in --sample_meta: %s",
    paste(missing_vars, collapse = ", ")
  ), call. = FALSE)
}
for (v in all_vars) sample_meta[[v]] <- factor(sample_meta[[v]])

qc_design <- as.formula(paste("~", paste(all_vars, collapse = " + ")))
dds_qc <- DESeqDataSetFromMatrix(countData = counts_filt, colData = sample_meta, design = qc_design)
dds_qc <- DESeq(dds_qc)

norm_counts <- counts(dds_qc, normalized = TRUE)
norm_counts_out <- data.frame(gene_id = rownames(norm_counts), norm_counts, check.names = FALSE)
write.table(norm_counts_out, file.path(opt$outdir, "normalized_counts.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

size_factors_out <- data.frame(sample = colnames(dds_qc), size_factor = sizeFactors(dds_qc))
write.table(size_factors_out, file.path(opt$outdir, "size_factors.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

saveRDS(dds_qc, file.path(opt$outdir, "dds.rds"))

## ---- PCA plot (ggplot2, a real DESeq2 Import - always available) ----
vsd <- tryCatch(
  vst(dds_qc, blind = TRUE),
  error = function(e) varianceStabilizingTransformation(dds_qc, blind = TRUE)
)
pca_color_var <- all_vars[1]
pca_data <- plotPCA(vsd, intgroup = pca_color_var, returnData = TRUE)
pct_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = .data[[pca_color_var]], label = name)) +
  geom_point(size = 3) +
  geom_text(size = 3, vjust = -0.8, check_overlap = TRUE, show.legend = FALSE) +
  xlab(paste0("PC1 (", pct_var[1], "%)")) +
  ylab(paste0("PC2 (", pct_var[2], "%)")) +
  theme_bw() +
  ggtitle("PCA of variance-stabilized counts")
ggsave(file.path(opt$outdir, "pca.png"), p_pca, width = 6, height = 5, dpi = 150)

## ---- sample-distance heatmap (base R heatmap(), no pheatmap dependency) ----
sample_dist <- dist(t(assay(vsd)))
dist_mat <- as.matrix(sample_dist)
png(file.path(opt$outdir, "sample_distance_heatmap.png"), width = 900, height = 800, res = 150)
heatmap(dist_mat, symm = TRUE, main = "Sample-to-sample distance (VST counts)",
        margins = c(8, 8), col = heat.colors(32, rev = TRUE))
dev.off()

## ---- per-contrast differential expression ----------------------------------

run_one_contrast <- function(row) {
  id        <- row[["id"]]
  variable  <- row[["variable"]]
  reference <- row[["reference"]]
  target    <- row[["target"]]

  message(sprintf("Running contrast '%s': %s (%s vs %s)", id, variable, target, reference))

  levels_present <- levels(sample_meta[[variable]])
  if (!(reference %in% levels_present) || !(target %in% levels_present)) {
    stop(sprintf(
      "ERROR: contrast '%s' - reference/target ('%s'/'%s') not both found among values of '%s' (%s)",
      id, reference, target, variable, paste(levels_present, collapse = ", ")
    ), call. = FALSE)
  }

  design_formula <- as.formula(paste("~", variable))
  dds <- DESeqDataSetFromMatrix(countData = counts_filt, colData = sample_meta, design = design_formula)
  dds[[variable]] <- relevel(dds[[variable]], ref = reference)
  dds <- DESeq(dds)

  res <- results(dds, contrast = c(variable, target, reference), alpha = opt$alpha)
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  res_df$significant <- !is.na(res_df$padj) & res_df$padj < opt$alpha
  res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

  out_tsv <- file.path(opt$outdir, paste0(id, ".deseq2_results.tsv"))
  write.table(res_df, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  message(sprintf(
    "  -> %s: %d genes tested, %d significant at padj < %s",
    out_tsv, sum(!is.na(res_df$padj)), sum(res_df$significant, na.rm = TRUE), opt$alpha
  ))

  ## Volcano plot
  volcano_df <- res_df
  volcano_df$neg_log10_padj <- -log10(volcano_df$padj)
  volcano_df$direction <- "NS"
  volcano_df$direction[volcano_df$significant & volcano_df$log2FoldChange > 0] <- "Up"
  volcano_df$direction[volcano_df$significant & volcano_df$log2FoldChange < 0] <- "Down"
  volcano_df$direction <- factor(volcano_df$direction, levels = c("Down", "NS", "Up"))

  top_labels <- head(volcano_df[order(volcano_df$padj, na.last = TRUE), ], 10)

  p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_padj, color = direction)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c(Down = "#3366CC", NS = "grey70", Up = "#CC3333")) +
    geom_hline(yintercept = -log10(opt$alpha), linetype = "dashed", color = "grey40") +
    geom_text(data = top_labels, aes(label = gene_id), size = 2.8, vjust = -0.6,
              check_overlap = TRUE, color = "black") +
    theme_bw() +
    labs(title = sprintf("%s: %s vs %s", id, target, reference),
         x = "log2 fold change", y = "-log10(padj)", color = NULL)
  ggsave(file.path(opt$outdir, paste0(id, ".volcano.png")), p_volcano, width = 6.5, height = 5.5, dpi = 150)

  ## MA plot
  png(file.path(opt$outdir, paste0(id, ".ma_plot.png")), width = 900, height = 700, res = 150)
  plotMA(res, alpha = opt$alpha, main = sprintf("%s: %s vs %s", id, target, reference))
  dev.off()

  invisible(res_df)
}

for (i in seq_len(nrow(contrasts_df))) {
  run_one_contrast(contrasts_df[i, ])
}

message("Done.")
