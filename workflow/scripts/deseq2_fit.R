#!/usr/bin/env Rscript
# DESeq2 fit: build the count matrix, run normalization,
# dispersion estimation and model fitting.
#
# The design is the full crossing of every parsed sample field
# (all main effects and interactions), built and fit as an
# explicit model matrix rather than a formula string. This is
# necessary, not just tidier: if the fields don't form a
# complete factorial (e.g. a background x dox-status cell with
# zero samples, as will happen once 'clean' samples with no
# GATA2R396Q counterpart are added), the naive formula produces
# an all-zero model-matrix column and DESeq2 refuses to fit
# (checkFullRank error). Building the matrix directly lets us
# detect and drop exactly those structurally-empty columns --
# something a formula string cannot express -- following
# DESeq2's own documented workaround for this situation.
#
# The cleaned matrix is stashed in metadata(dds) so
# deseq2_contrast.R can reuse the exact matrix that was fit,
# rather than recomputing (and risking drift from) it.

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
    ncol(counts_matrix), "samples found in count matrix.\n")

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
    "WARNING: at least one condition has fewer than 2",
    "replicates. DESeq2 requires n >= 2 per group; results",
    "for affected contrasts will likely fail or be unreliable.\n"
  )
} else if (any(n_per_condition == 2)) {
  cat(
    "NOTE: minimum replicate number is 2. Dispersion",
    "estimates from n = 2 are noisy; treat results as",
    "exploratory and interpret with caution.\n"
  )
}

# --- Apply reference levels ------------------------------------------
# Affects resultsNames()/coefficient signs and which
# combinations end up as the "fully reference" corner of the
# full-crossing matrix below; does not affect fitted group means
# or dispersions.

reference_levels <- snakemake@params[["reference_levels"]]
cat("\nReference levels:\n")
for (f in names(reference_levels)) {
  ref <- reference_levels[[f]]
  cat(" -", f, ":", ref, "\n")
  coldata[[f]] <- relevel(factor(coldata[[f]]), ref = ref)
}

# --- Build the full-crossing design matrix, dropping any -----------
# --- structurally-empty (all-zero) columns. -------------------------

full_formula <- as.formula(
  paste("~", paste(field_names, collapse = "*"))
)
cat("\nFull-crossing formula (before rank check):",
    deparse(full_formula), "\n")

full_mat <- model.matrix(full_formula, data = coldata)
empty_cols <- colSums(full_mat != 0) == 0

if (any(empty_cols)) {
  cat(
    "\nNOTE: the parsed fields do not form a complete factorial",
    "-- the following coefficient(s) have no observations and",
    "are being dropped before fitting (expected, not an error,",
    "e.g. a construct-free condition existing for only one",
    "background):\n"
  )
  print(colnames(full_mat)[empty_cols])
  full_mat <- full_mat[, !empty_cols, drop = FALSE]
}

cat("\nFinal design matrix columns (", ncol(full_mat), "):\n", sep = "")
print(colnames(full_mat))

n_distinct_conditions <- length(unique(coldata$condition))
if (ncol(full_mat) != n_distinct_conditions) {
  stop(
    "ERROR: final design matrix has ", ncol(full_mat),
    " columns but there are ", n_distinct_conditions,
    " distinct observed conditions. The model is not saturated ",
    "relative to the observed groups -- this would break the ",
    "relevel-based contrast extraction used elsewhere in this ",
    "pipeline. This should not happen for a matrix built purely ",
    "from crossing the parsed fields; investigate before ",
    "trusting downstream results."
  )
}

# --- Build DESeqDataSet ------------------------------------------
# design is a placeholder (~condition, always full rank by
# construction) at object-creation time: DESeqDataSetFromMatrix
# would itself reject the full, uncleaned crossing formula via
# its own rank check. The actual (cleaned) design matrix is
# supplied to DESeq() below via `full =`.

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_matrix),
  colData = coldata,
  design = ~condition
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

cat("\nRunning DESeq2 (normalize, estimate dispersions, fit)",
    "with the cleaned full-crossing design matrix...\n")
dds <- tryCatch(
  DESeq(dds, full = full_mat, betaPrior = FALSE),
  error = function(e) {
    stop("ERROR: DESeq() failed to fit the model: ", e$message)
  }
)

cat("Normalization, dispersion estimation and model fitting",
    "complete.\n")
cat("\nModel coefficients (resultsNames):\n")
print(resultsNames(dds))

if ("betaConv" %in% colnames(mcols(dds))) {
  n_not_converged <- sum(!mcols(dds)$betaConv, na.rm = TRUE)
  if (n_not_converged > 0) {
    cat(
      "WARNING:", n_not_converged,
      "gene(s) did not converge during beta (coefficient)",
      "fitting. Their results should be treated with caution;",
      "see mcols(dds)$betaConv for the affected genes.\n"
    )
  }
}

if ("dispOutlier" %in% colnames(mcols(dds))) {
  n_disp_outliers <- sum(mcols(dds)$dispOutlier, na.rm = TRUE)
  if (n_disp_outliers > 0) {
    cat(
      "NOTE:", n_disp_outliers,
      "gene(s) had dispersion estimates flagged as outliers",
      "relative to the fitted mean-dispersion trend (not",
      "shrunk). Expected for some genes; not necessarily a",
      "problem, but worth a look for genes of interest.\n"
    )
  }
}

# --- Sync column names to resultsNames(dds) before stashing --------
# DESeq2 relabels raw model.matrix() coefficient names for its
# own resultsNames()/results(name=) interface -- notably ":" in
# interaction terms becomes "." (e.g. "backgroundGATA2R396Q:doxdox"
# -> "backgroundGATA2R396Q.doxdox"). This is a deliberate DESeq2
# display choice, not an R syntax requirement, and it means the
# raw model.matrix() colnames on full_mat do NOT match
# resultsNames(dds) for the identical coefficients. The pipeline
# itself resolves coefficients by position, not by name, so this
# mismatch was never a functional bug -- but syncing here means
# anyone loading this RDS for a manual, ad hoc check sees the
# same names in metadata(dds)$full_model_matrix as resultsNames()
# reports, and isn't caught out by it the way I was when I first
# encountered this.

if (length(resultsNames(dds)) != ncol(full_mat)) {
  stop(
    "ERROR: resultsNames(dds) has ", length(resultsNames(dds)),
    " entries but full_mat has ", ncol(full_mat), " columns; ",
    "cannot safely synchronize column names. Investigate before ",
    "trusting the stashed model matrix."
  )
}
colnames(full_mat) <- resultsNames(dds)

cat("\nStashed model matrix columns (synced to resultsNames):\n")
print(colnames(full_mat))

# --- Stash the exact fitted matrix for deseq2_contrast.R ------------

metadata(dds)$full_model_matrix <- full_mat

# --- Save --------------------------------------------------------

output_dds <- snakemake@output[["dds_rds"]]
dir.create(dirname(output_dds), showWarnings = FALSE, recursive = TRUE)
saveRDS(dds, output_dds)
cat("\nSaved DESeqDataSet:", output_dds, "\n")

cat("\n=== DESeq2 Fit Complete ===\n")
sink(type = "message")
sink()
close(logfile)
