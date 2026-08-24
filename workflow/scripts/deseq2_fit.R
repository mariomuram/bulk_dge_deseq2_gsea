#!/usr/bin/env Rscript
# DESeq2 fit: build the count matrix, run normalization,
# dispersion estimation and model fitting for a single
# design formula. Produces one DESeqDataSet (RDS) that all
# contrasts belonging to this design extract results from.

library(DESeq2)
library(readr)

# Capture log
logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== DESeq2 Fit ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Read count matrix ---------------------------------------------

cat("Reading count matrix...\n")
counts_file <- snakemake@input[["counts"]]
counts_df <- read_tsv(counts_file, show_col_types = FALSE)

cat("Input file:", counts_file, "\n")
cat("Dimensions:", nrow(counts_df), "genes x",
    ncol(counts_df) - 2, "samples (raw file)\n")

gene_ids <- counts_df[[1]]
counts_matrix <- as.matrix(counts_df[, -c(1, 2)])
rownames(counts_matrix) <- gene_ids

# --- Build sample metadata ------------------------------------------

samples_info <- snakemake@params[["samples_info"]]

missing_samples <- setdiff(names(samples_info), colnames(counts_matrix))
if (length(missing_samples) > 0) {
  stop(
    "ERROR: sample(s) listed in config.yaml 'samples' are not ",
    "present in the count matrix: ",
    paste(missing_samples, collapse = ", ")
  )
}

cat("\nUsing", length(samples_info), "of",
    ncol(counts_matrix), "samples found in count matrix ",
    "(others are ignored, e.g. parental_clean).\n")

counts_matrix <- counts_matrix[, names(samples_info), drop = FALSE]

condition_labels <- unlist(samples_info, use.names = TRUE)
label_parts <- strsplit(condition_labels, "_", fixed = TRUE)

if (any(lengths(label_parts) != 2)) {
  bad <- names(condition_labels)[lengths(label_parts) != 2]
  stop(
    "ERROR: expected sample labels of the form ",
    "'<background>_<dox status>' (e.g. 'GATA2R396Q_dox'). ",
    "Could not parse: ", paste(bad, collapse = ", ")
  )
}

coldata <- data.frame(
  sample = colnames(counts_matrix),
  condition = condition_labels[colnames(counts_matrix)],
  background = vapply(label_parts, `[[`, character(1), 1)[
    colnames(counts_matrix)
  ],
  dox = vapply(label_parts, `[[`, character(1), 2)[
    colnames(counts_matrix)
  ],
  stringsAsFactors = FALSE
)
rownames(coldata) <- coldata$sample

cat("\nSample sheet:\n")
print(coldata)

if (anyNA(coldata$condition)) {
  stop(
    "ERROR: could not map all samples to a condition. ",
    "Check config.yaml 'samples' section."
  )
}

cat("\nReplicates per condition:\n")
print(table(coldata$condition))

n_per_condition <- table(coldata$condition)
if (any(n_per_condition < 2)) {
  cat(
    "WARNING: at least one condition has fewer than 2 ",
    "replicates. DESeq2 requires n >= 2 per group; results ",
    "for affected contrasts will likely fail or be unreliable.\n"
  )
} else if (any(n_per_condition == 2)) {
  cat(
    "NOTE: minimum replicate number is 2. Dispersion ",
    "estimates from n = 2 are noisy; treat results as ",
    "exploratory and interpret with caution.\n"
  )
}

# --- Build DESeqDataSet ----------------------------------------------

design_formula <- as.formula(snakemake@params[["design_formula"]])
cat("\nDesign formula:", deparse(design_formula), "\n")

cat("\nCreating DESeqDataSet...\n")
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_matrix),
  colData = coldata,
  design = design_formula
)

cat("DESeqDataSet created with", nrow(dds), "genes and",
    ncol(dds), "samples\n")

# Pre-filtering: remove genes with very low counts
cat("\nPre-filtering low-count genes (rowSums >= 10)...\n")
keep <- rowSums(counts(dds)) >= 10
cat("Genes removed:", sum(!keep), "\n")
dds <- dds[keep, ]
cat("Genes retained after filtering:", nrow(dds), "\n")

# --- Run DESeq2 -------------------------------------------------------

cat("\nRunning DESeq2 (normalize, estimate dispersions, fit)...\n")
dds <- tryCatch(
  DESeq(dds),
  error = function(e) {
    stop("ERROR: DESeq() failed to fit the model: ", e$message)
  }
)

cat("Normalization, dispersion estimation and model fitting ",
    "complete.\n")

if ("betaConv" %in% colnames(mcols(dds))) {
  n_not_converged <- sum(!mcols(dds)$betaConv, na.rm = TRUE)
  if (n_not_converged > 0) {
    cat(
      "WARNING:", n_not_converged,
      "gene(s) did not converge during beta (coefficient) ",
      "fitting. Their results should be treated with caution; ",
      "see mcols(dds)$betaConv for the affected genes.\n"
    )
  }
}

if ("dispOutlier" %in% colnames(mcols(dds))) {
  n_disp_outliers <- sum(mcols(dds)$dispOutlier, na.rm = TRUE)
  if (n_disp_outliers > 0) {
    cat(
      "NOTE:", n_disp_outliers,
      "gene(s) had dispersion estimates flagged as outliers ",
      "relative to the fitted mean-dispersion trend (not ",
      "shrunk). Expected for some genes; not necessarily a ",
      "problem, but worth a look for genes of interest.\n"
    )
  }
}

# --- Save --------------------------------------------------------------

output_dds <- snakemake@output[["dds_rds"]]
dir.create(dirname(output_dds), showWarnings = FALSE, recursive = TRUE)
saveRDS(dds, output_dds)
cat("\nSaved DESeqDataSet:", output_dds, "\n")

cat("\n=== DESeq2 Fit Complete ===\n")
sink(type = "message")
sink()
close(logfile)
