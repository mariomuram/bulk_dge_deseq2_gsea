#!/usr/bin/env Rscript
# DESeq2 fit: build the count matrix, run normalization,
# dispersion estimation and model fitting. One shared fit
# covers every contrast (pairwise and combination alike) via
# group-averaged model-matrix contrasts in deseq2_contrast.R.
# Sample fields (background, dox, ...) are parsed from labels
# in the Snakefile using config's sample_label_pattern, not
# hardcoded here.

library(DESeq2)
library(readr)

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

gene_ids <- counts_df[[1]]
counts_matrix <- as.matrix(counts_df[, -c(1, 2)])
rownames(counts_matrix) <- gene_ids

cat("Dimensions:", nrow(counts_df), "genes x",
    ncol(counts_matrix), "samples (raw file)\n")

# --- Build sample metadata (fields parsed in the Snakefile) --------

samples_info <- snakemake@params[["samples_info"]]    # sample -> label
sample_fields <- snakemake@params[["sample_fields"]]   # sample -> list

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

field_names <- names(sample_fields[[1]])
cat("Parsed sample fields:", paste(field_names, collapse = ", "), "\n")

coldata <- data.frame(
  sample = colnames(counts_matrix),
  condition = unlist(samples_info)[colnames(counts_matrix)],
  stringsAsFactors = FALSE
)
for (f in field_names) {
  coldata[[f]] <- vapply(
    sample_fields[colnames(counts_matrix)],
    function(x) x[[f]], character(1)
  )
}
rownames(coldata) <- coldata$sample

cat("\nSample sheet:\n")
print(coldata)

field_has_na <- vapply(coldata[field_names], anyNA, logical(1))
if (anyNA(coldata$condition) || any(field_has_na)) {
  stop(
    "ERROR: could not map all samples to their condition/",
    "fields. Check config.yaml 'samples' and ",
    "'sample_label_pattern'."
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

# --- Apply reference levels ------------------------------------------
# Cosmetic only: does not affect any contrast computed via the
# group-averaged model-matrix technique in deseq2_contrast.R.
# Kept so resultsNames(dds) and coefficient signs "make sense"
# for anyone loading this RDS into a fresh R session directly.

reference_levels <- snakemake@params[["reference_levels"]]
cat("\nReference levels (cosmetic; does not affect contrast\n",
    "values computed downstream):\n")
for (f in names(reference_levels)) {
  ref <- reference_levels[[f]]
  cat(" -", f, ":", ref, "\n")
  coldata[[f]] <- relevel(factor(coldata[[f]]), ref = ref)
}

# --- Build DESeqDataSet ------------------------------------------

design_formula <- as.formula(snakemake@params[["design_formula"]])
cat("\nDesign formula:", deparse(design_formula), "\n")

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
cat("\nModel coefficients (resultsNames):\n")
print(resultsNames(dds))

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

# --- Save --------------------------------------------------------

output_dds <- snakemake@output[["dds_rds"]]
dir.create(dirname(output_dds), showWarnings = FALSE, recursive = TRUE)
saveRDS(dds, output_dds)
cat("\nSaved DESeqDataSet:", output_dds, "\n")

cat("\n=== DESeq2 Fit Complete ===\n")
sink(type = "message")
sink()
close(logfile)
