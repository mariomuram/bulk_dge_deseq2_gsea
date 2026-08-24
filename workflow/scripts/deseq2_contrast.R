#!/usr/bin/env Rscript
# DESeq2 contrast extraction: pull a single pairwise
# contrast's results out of a previously-fit DESeqDataSet.

library(DESeq2)
library(readr)

# Capture log
logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== DESeq2 Contrast Extraction ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Load fitted model -----------------------------------------------

dds_file <- snakemake@input[["dds_rds"]]
cat("Loading DESeqDataSet:", dds_file, "\n")
dds <- readRDS(dds_file)

contrast <- snakemake@params[["contrast"]]
str(contrast)
numerator <- contrast[["numerator"]]
denominator <- contrast[["denominator"]]
padj_threshold <- snakemake@params[["padj_threshold"]]
log2fc_threshold <- snakemake@params[["log2fc_threshold"]]

cat("Numerator:  ", numerator, "\n")
cat("Denominator:", denominator, "\n")

available_levels <- levels(factor(colData(dds)$condition))
missing_levels <- setdiff(c(numerator, denominator), available_levels)
if (length(missing_levels) > 0) {
  stop(
    "ERROR: condition level(s) not found in the fitted model: ",
    paste(missing_levels, collapse = ", "),
    ". Available levels: ", paste(available_levels, collapse = ", ")
  )
}

n_per_level <- table(colData(dds)$condition)[c(numerator, denominator)]
cat("\nReplicates in this contrast:\n")
print(n_per_level)
if (any(n_per_level < 2)) {
  cat(
    "WARNING: fewer than 2 replicates for this contrast. ",
    "Results may be unreliable or extraction may fail.\n"
  )
}

# --- Extract results ---------------------------------------------------

cat("\nExtracting results...\n")
res <- tryCatch(
  results(
    dds,
    contrast = c("condition", numerator, denominator),
    independentFiltering = TRUE,
    alpha = padj_threshold
  ),
  error = function(e) {
    stop("ERROR: results() failed for this contrast: ", e$message)
  }
)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res)
res_df <- res_df[, c(
  "gene_id", "baseMean", "log2FoldChange", "lfcSE",
  "stat", "pvalue", "padj"
)]
res_df <- res_df[order(res_df$padj), ]

# --- Report NAs so issues are visible, not silent -----------------------

n_outlier_na <- sum(is.na(res_df$pvalue))
n_filtered_na <- sum(is.na(res_df$padj) & !is.na(res_df$pvalue))

if (n_outlier_na > 0) {
  cat(
    "NOTE:", n_outlier_na,
    "gene(s) have NA pvalue/padj: flagged as Cook's distance ",
    "outliers (an extreme count in a single sample driving the ",
    "result). See the DESeq2 vignette section on count outliers.\n"
  )
}
if (n_filtered_na > 0) {
  cat(
    "NOTE:", n_filtered_na,
    "gene(s) have NA padj (but a pvalue): removed by ",
    "independent filtering of low mean-count genes to improve ",
    "power. This is expected behavior, not an error.\n"
  )
}

n_sig <- sum(
  res_df$padj < padj_threshold &
    abs(res_df$log2FoldChange) >= log2fc_threshold,
  na.rm = TRUE
)
cat(
  "\nSignificant genes (padj <", padj_threshold,
  "& |log2FC| >=", log2fc_threshold, "):", n_sig, "\n"
)

# --- Save (gzip --best compressed TSV) -----------------------------------

output_file <- snakemake@output[["results"]]
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
con <- gzfile(output_file, open = "wb", compression = 9)
write_tsv(res_df, con)
close(con)
cat("\nSaved:", output_file, "\n")

cat("\n=== DESeq2 Contrast Extraction Complete ===\n")
sink(type = "message")
sink()
close(logfile)
