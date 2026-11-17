#!/usr/bin/env Rscript
# Gene-set enrichment on a DESeq2 results table, against a user-supplied GMT
# gene-set file. Organism-agnostic by design (no org.*.eg.db / bioc annotation
# dependency) - works for any species as long as --gmt_file's gene IDs use the
# same namespace as --deseq_results' gene_id column (match your GTF's gene_id
# attribute, e.g. Ensembl IDs or gene symbols).
#
# Implements two standard, independently-coded methods rather than wrapping
# clusterProfiler/fgsea (kept out to avoid a heavy Bioconductor annotation
# dependency chain):
#   1. ORA  - one-sided hypergeometric (Fisher's exact) test per gene set,
#              on the significant genes (padj < --alpha), vs. all tested
#              genes as the background.
#   2. GSEA - classic weighted running-sum enrichment score (Subramanian et
#              al. 2005 "preranked" method), gene-set-label permutation for
#              the null distribution (same null model fgsea itself samples
#              from), ranking genes by --rank_metric.
#
# Inputs:
#   --deseq_results   <id>.deseq2_results.tsv from deseq2_differential.R
#   --gmt_file        GMT file: one gene set per line, tab-separated:
#                      <set_name>\t<description>\t<gene1>\t<gene2>\t...
#   --id              label used in output filenames/plot titles
#   --outdir          output directory
#   --alpha           padj threshold defining "significant" for ORA (default 0.05)
#   --rank_metric     column of --deseq_results to rank genes by for GSEA:
#                      "stat" (default) or "log2FoldChange"
#   --min_set_size    drop gene sets smaller than this after overlap with the
#                      tested gene universe (default 5)
#   --max_set_size    drop gene sets larger than this (default 500)
#   --nperm           permutations for the GSEA null distribution (default 1000)
#   --top_n           how many top gene sets to show in each plot (default 15)
#
# Outputs:
#   <id>.ora_results.tsv   gene_set, description, set_size, overlap_size,
#                           overlap_genes, pvalue, padj
#   <id>.ora_dotplot.png
#   <id>.gsea_results.tsv  gene_set, description, size, ES, NES, pvalue, padj
#   <id>.gsea_dotplot.png

suppressPackageStartupMessages({
  library(ggplot2)  # a real DESeq2 Import - always present in this module's container
})

## ---- minimal base-R CLI parser (no optparse dependency) -------------------

arg_spec <- list(
  deseq_results = list(default = NULL,       type = "character"),
  gmt_file      = list(default = NULL,       type = "character"),
  id            = list(default = "contrast", type = "character"),
  outdir        = list(default = ".",        type = "character"),
  alpha         = list(default = 0.05,       type = "double"),
  rank_metric   = list(default = "stat",     type = "character"),
  min_set_size  = list(default = 5,          type = "integer"),
  max_set_size  = list(default = 500,        type = "integer"),
  nperm         = list(default = 1000,       type = "integer"),
  top_n         = list(default = 15,         type = "integer"),
  seed          = list(default = 42,         type = "integer")
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

for (required in c("deseq_results", "gmt_file")) {
  if (is.null(opt[[required]])) stop(sprintf("ERROR: --%s is required", required), call. = FALSE)
}
if (!opt$rank_metric %in% c("stat", "log2FoldChange")) {
  stop("ERROR: --rank_metric must be 'stat' or 'log2FoldChange'", call. = FALSE)
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(opt$seed)

## ---- load DESeq2 results ---------------------------------------------------

res <- read.delim(opt$deseq_results, header = TRUE, stringsAsFactors = FALSE)
required_cols <- c("gene_id", "log2FoldChange", "stat", "padj")
missing_cols <- setdiff(required_cols, colnames(res))
if (length(missing_cols) > 0) {
  stop(sprintf("ERROR: --deseq_results is missing column(s): %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
}
res <- res[!is.na(res[[opt$rank_metric]]), ]
if (nrow(res) == 0) stop("ERROR: no genes with a non-NA rank_metric value in --deseq_results", call. = FALSE)

universe <- unique(res$gene_id)
sig_genes <- unique(res$gene_id[!is.na(res$padj) & res$padj < opt$alpha])
message(sprintf("Universe: %d genes tested, %d significant at padj < %s", length(universe), length(sig_genes), opt$alpha))

## ---- load GMT gene sets, restricted to the tested universe ----------------

read_gmt <- function(path) {
  lines <- readLines(path)
  lines <- lines[nzchar(trimws(lines))]
  sets <- lapply(lines, function(l) {
    parts <- strsplit(l, "\t")[[1]]
    list(name = parts[1], description = ifelse(length(parts) >= 2, parts[2], ""), genes = unique(parts[-c(1, 2)]))
  })
  names(sets) <- vapply(sets, function(s) s$name, character(1))
  sets
}

gmt <- read_gmt(opt$gmt_file)
if (length(gmt) == 0) stop(sprintf("ERROR: no gene sets parsed from %s", opt$gmt_file), call. = FALSE)

gmt_filtered <- lapply(gmt, function(s) {
  s$genes_in_universe <- intersect(s$genes, universe)
  s
})
set_sizes <- vapply(gmt_filtered, function(s) length(s$genes_in_universe), integer(1))
keep <- set_sizes >= opt$min_set_size & set_sizes <= opt$max_set_size
message(sprintf(
  "Gene sets: %d loaded, %d kept after restricting to the %d-gene universe and size filter [%d, %d]",
  length(gmt), sum(keep), length(universe), opt$min_set_size, opt$max_set_size
))
gmt_filtered <- gmt_filtered[keep]
if (length(gmt_filtered) == 0) {
  stop("ERROR: no gene sets survived the universe-overlap/size filter - check gene ID namespace matches between --gmt_file and --deseq_results", call. = FALSE)
}

## ---- 1. Over-representation analysis (hypergeometric test) ----------------
## For gene set S (size n_S, overlap k with the significant genes):
##   P(X >= k) under the hypergeometric model with population N = |universe|,
##   "successes in population" = n_S, "sample size" = |sig_genes|
## This is the one-sided test clusterProfiler's enricher()/DAVID/g:Profiler
## all compute under the hood.

run_ora <- function() {
  N <- length(universe)
  n_sig <- length(sig_genes)

  rows <- lapply(names(gmt_filtered), function(set_name) {
    s <- gmt_filtered[[set_name]]
    overlap_genes <- intersect(s$genes_in_universe, sig_genes)
    k <- length(overlap_genes)
    n_S <- length(s$genes_in_universe)
    ## P(X >= k) = 1 - P(X <= k-1)
    pval <- phyper(k - 1, n_S, N - n_S, n_sig, lower.tail = FALSE)
    data.frame(
      gene_set = set_name, description = s$description,
      set_size = n_S, overlap_size = k, n_significant = n_sig, universe_size = N,
      overlap_genes = paste(sort(overlap_genes), collapse = ","),
      pvalue = pval, stringsAsFactors = FALSE
    )
  })
  ora_df <- do.call(rbind, rows)
  ora_df$padj <- p.adjust(ora_df$pvalue, method = "BH")
  ora_df <- ora_df[order(ora_df$pvalue), ]
  ora_df
}

ora_df <- run_ora()
ora_out <- file.path(opt$outdir, paste0(opt$id, ".ora_results.tsv"))
write.table(ora_df, ora_out, sep = "\t", row.names = FALSE, quote = FALSE)
message(sprintf(
  "ORA: wrote %s (%d gene sets, %d significant at padj < %s)",
  ora_out, nrow(ora_df), sum(ora_df$padj < opt$alpha), opt$alpha
))

if (nrow(ora_df) > 0) {
  top_ora <- head(ora_df, opt$top_n)
  top_ora$gene_set <- factor(top_ora$gene_set, levels = rev(top_ora$gene_set))
  p_ora <- ggplot(top_ora, aes(x = -log10(pvalue), y = gene_set, size = overlap_size, color = padj)) +
    geom_point() +
    scale_color_gradient(low = "#CC3333", high = "#3366CC") +
    theme_bw() +
    labs(title = sprintf("%s: top gene sets (ORA)", opt$id),
         x = "-log10(p-value)", y = NULL, size = "Overlap size", color = "padj")
  ggsave(file.path(opt$outdir, paste0(opt$id, ".ora_dotplot.png")), p_ora,
         width = 7.5, height = max(3, 0.35 * nrow(top_ora) + 1.5), dpi = 150, limitsize = FALSE)
}

## ---- 2. GSEA (weighted running-sum, preranked, gene-set permutation) ------

compute_es <- function(hit_logical, rank_weights, p = 1) {
  ## hit_logical / rank_weights already ordered by decreasing rank_weights
  N <- length(hit_logical)
  n_hit <- sum(hit_logical)
  if (n_hit == 0 || n_hit == N) return(0)
  hit_score <- (abs(rank_weights))^p
  hit_score[!hit_logical] <- 0
  sum_hit <- sum(hit_score)
  if (sum_hit == 0) return(0)
  step_hit <- hit_score / sum_hit
  step_miss <- (!hit_logical) / (N - n_hit)
  running <- cumsum(step_hit - step_miss)
  running[which.max(abs(running))]
}

run_gsea <- function() {
  ranked <- res[order(-res[[opt$rank_metric]]), c("gene_id", opt$rank_metric)]
  ranked <- ranked[!duplicated(ranked$gene_id), ]
  rank_weights <- ranked[[opt$rank_metric]]
  gene_order <- ranked$gene_id
  N <- length(gene_order)

  rows <- lapply(names(gmt_filtered), function(set_name) {
    s <- gmt_filtered[[set_name]]
    hit_logical <- gene_order %in% s$genes_in_universe
    n_hit <- sum(hit_logical)
    es <- compute_es(hit_logical, rank_weights)

    ## Null distribution: permute which n_hit genes count as "hits", holding
    ## the ranked list fixed (the standard "gene-set permutation" preranked
    ## null - the same one fgsea's simple/non-analytic mode samples from).
    null_es <- replicate(opt$nperm, {
      perm_hit <- logical(N)
      perm_hit[sample.int(N, n_hit)] <- TRUE
      compute_es(perm_hit, rank_weights)
    })

    if (es >= 0) {
      pos_null <- null_es[null_es >= 0]
      denom <- if (length(pos_null) > 0) mean(pos_null) else NA_real_
      pval <- (sum(pos_null >= es) + 1) / (length(pos_null) + 1)
    } else {
      neg_null <- null_es[null_es < 0]
      denom <- if (length(neg_null) > 0) mean(abs(neg_null)) else NA_real_
      pval <- (sum(abs(neg_null) >= abs(es)) + 1) / (length(neg_null) + 1)
    }
    nes <- if (!is.na(denom) && denom > 0) es / denom else NA_real_

    data.frame(
      gene_set = set_name, description = s$description, size = n_hit,
      ES = es, NES = nes, pvalue = pval, stringsAsFactors = FALSE
    )
  })

  gsea_df <- do.call(rbind, rows)
  gsea_df$padj <- p.adjust(gsea_df$pvalue, method = "BH")
  gsea_df[order(gsea_df$pvalue), ]
}

gsea_df <- run_gsea()
gsea_out <- file.path(opt$outdir, paste0(opt$id, ".gsea_results.tsv"))
write.table(gsea_df, gsea_out, sep = "\t", row.names = FALSE, quote = FALSE)
message(sprintf(
  "GSEA: wrote %s (%d gene sets, %d significant at padj < %s, ranked by %s)",
  gsea_out, nrow(gsea_df), sum(gsea_df$padj < opt$alpha, na.rm = TRUE), opt$alpha, opt$rank_metric
))

if (nrow(gsea_df) > 0) {
  top_gsea <- head(gsea_df, opt$top_n)
  top_gsea$gene_set <- factor(top_gsea$gene_set, levels = rev(top_gsea$gene_set))
  p_gsea <- ggplot(top_gsea, aes(x = NES, y = gene_set, fill = padj)) +
    geom_col() +
    scale_fill_gradient(low = "#CC3333", high = "#3366CC") +
    geom_vline(xintercept = 0, color = "grey40") +
    theme_bw() +
    labs(title = sprintf("%s: top gene sets (GSEA, ranked by %s)", opt$id, opt$rank_metric),
         x = "Normalized enrichment score (NES)", y = NULL, fill = "padj")
  ggsave(file.path(opt$outdir, paste0(opt$id, ".gsea_dotplot.png")), p_gsea,
         width = 7.5, height = max(3, 0.35 * nrow(top_gsea) + 1.5), dpi = 150, limitsize = FALSE)
}

message("Done.")
