#!/usr/bin/env Rscript
# DESeq2 contrast extraction: build a contrast as a difference
# of averaged model-matrix rows for arbitrary condition groups
# (per the DESeq2 vignette's approach for comparisons that
# aren't a single named coefficient). One mechanism covers
# pairwise contrasts, the interaction term, and simple-effect/
# combination contrasts alike -- the type of contrast is fully
# determined by the "groups" weights, nothing else changes.

library(DESeq2)
library(readr)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== DESeq2 Contrast Extraction ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

dds_file <- snakemake@input[["dds_rds"]]
cat("Loading DESeqDataSet:", dds_file, "\n")
dds <- readRDS(dds_file)

groups <- snakemake@params[["groups"]]  # named list: condition -> weight
padj_threshold <- snakemake@params[["padj_threshold"]]
log2fc_threshold <- snakemake@params[["log2fc_threshold"]]

cat("\nContrast groups (condition: weight):\n")
for (g in names(groups)) cat(" -", g, ":", groups[[g]], "\n")

available_conditions <- levels(factor(colData(dds)$condition))
missing_conditions <- setdiff(names(groups), available_conditions)
if (length(missing_conditions) > 0) {
  stop(
    "ERROR: condition(s) not found in the fitted model: ",
    paste(missing_conditions, collapse = ", "),
    ". Available: ", paste(available_conditions, collapse = ", ")
  )
}

n_per_condition <- table(colData(dds)$condition)[names(groups)]
cat("\nReplicates per group:\n")
print(n_per_condition)
if (any(n_per_condition < 2)) {
  cat(
    "WARNING: fewer than 2 replicates for a condition used in ",
    "this contrast. Results may be unreliable.\n"
  )
}

# --- Build the contrast as a difference of averaged model-matrix
# --- rows, weighted per condition group. ----------------------------

mod_mat <- model.matrix(design(dds), data = as.data.frame(colData(dds)))

contrast_vec <- rep(0, ncol(mod_mat))
names(contrast_vec) <- colnames(mod_mat)

for (cond in names(groups)) {
  rows <- colData(dds)$condition == cond
  group_mean_row <- colMeans(mod_mat[rows, , drop = FALSE])
  contrast_vec <- contrast_vec + groups[[cond]] * group_mean_row
}

cat("\nResolved coefficient weights:\n")
print(contrast_vec)

# --- Extract results -----------------------------------------------

cat("\nExtracting results...\n")
res <- tryCatch(
  results(
    dds,
    contrast = contrast_vec,
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

n_outlier_na <- sum(is.na(res_df$pvalue))
n_filtered_na <- sum(is.na(res_df$padj) & !is.na(res_df$pvalue))

if (n_outlier_na > 0) {
  cat(
    "NOTE:", n_outlier_na,
    "gene(s) have NA pvalue/padj: flagged as Cook's distance ",
    "outliers (an extreme count in a single sample driving the ",
    "result).\n"
  )
}
if (n_filtered_na > 0) {
  cat(
    "NOTE:", n_filtered_na,
    "gene(s) have NA padj (but a pvalue): removed by ",
    "independent filtering of low mean-count genes. This is ",
    "expected behavior, not an error.\n"
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

# --- Save (gzip --best compressed TSV) -----------------------------

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
