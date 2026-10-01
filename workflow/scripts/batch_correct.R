#!/usr/bin/env Rscript
# Batch/confounder correction of the blind VST matrix.
#
# For VISUALIZATION AND DIAGNOSTICS ONLY (PCA, PC-field association). Never
# feed the result into DE testing: confounders are modelled in the DESeq2
# design instead. Uses the condition/confounder model matrices stashed by
# deseq2_fit.R, i.e. exactly the matrices the model was fit with.

library(DESeq2)
library(limma)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== Batch Correction ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

correction <- snakemake@wildcards[["correction"]]
cat("Correction method:", correction, "\n")

dds <- readRDS(snakemake@input[["dds_rds"]])
vst_mat <- readRDS(snakemake@input[["vst"]])
stopifnot(identical(colnames(vst_mat), colnames(dds)))

condition_mat <- metadata(dds)$condition_model_matrix
confounder_mat <- metadata(dds)$confounder_model_matrix
if (is.null(condition_mat) || is.null(confounder_mat)) {
  stop(
    "ERROR: metadata(dds) has no condition/confounder model matrix. Either ",
    "no confounders are configured (design.confounders is empty; there is ",
    "nothing to correct for and this rule should not have been requested), ",
    "or dds_object.rds predates the stash and deseq2_fit must be rerun."
  )
}
stopifnot(identical(rownames(condition_mat), colnames(vst_mat)),
          identical(rownames(confounder_mat), colnames(vst_mat)))
cat("Condition matrix columns (protected):",
    paste(colnames(condition_mat), collapse = ", "), "\n")
cat("Confounder matrix columns (removed):",
    paste(colnames(confounder_mat), collapse = ", "), "\n")

corrected <- switch(
  correction,
  # design = condition matrix: the biological signal to preserve.
  # covariates = confounder columns (dummy-coded if nominal): regressed out.
  limma = removeBatchEffect(vst_mat, design = condition_mat,
                            covariates = confounder_mat),
  stop("ERROR: unsupported correction method: ", correction)
)

output_file <- snakemake@output[["vst"]]
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(corrected, output_file)
cat("Saved:", output_file, "\n")

cat("\n=== Batch Correction Complete ===\n")
sink(type = "message")
sink()
close(logfile)
