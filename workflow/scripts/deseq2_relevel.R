#!/usr/bin/env Rscript
# Relevel + Wald retest: take the primary (interaction-design)
# DESeqDataSet, relevel `condition` to a new reference, and
# rerun only nbinomWaldTest(). Dispersions and size factors are
# reused as-is: the primary design is saturated (as many free
# parameters as distinct condition groups), so reparametrizing
# to ~condition changes which coefficient represents which
# comparison, but not the fitted group means or dispersions
# themselves. This makes every pairwise contrast against this
# reference extractable -- and apeglm-shrinkable -- as an
# ordinary named coefficient.

library(DESeq2)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== DESeq2 Relevel + Wald Retest ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

dds_file <- snakemake@input[["dds_rds"]]
reference_factor <- snakemake@params[["reference_factor"]]

cat("Loading primary DESeqDataSet:", dds_file, "\n")
dds <- readRDS(dds_file)

available_conditions <- levels(factor(colData(dds)$condition))
if (!(reference_factor %in% available_conditions)) {
  stop(
    "ERROR: reference_factor '", reference_factor, "' is not an ",
    "available condition. Available: ",
    paste(available_conditions, collapse = ", ")
  )
}

cat("Relevelling 'condition' to reference:", reference_factor, "\n")
colData(dds)$condition <- relevel(
  factor(colData(dds)$condition), ref = reference_factor
)

design(dds) <- ~condition

cat("\nRerunning nbinomWaldTest() only (dispersions and size",
    "factors reused from the primary fit)...\n")
dds <- tryCatch(
  nbinomWaldTest(dds),
  error = function(e) {
    stop("ERROR: nbinomWaldTest() failed: ", e$message)
  }
)

cat("\nModel coefficients (resultsNames) after relevel:\n")
print(resultsNames(dds))

output_dds <- snakemake@output[["dds_rds"]]
dir.create(dirname(output_dds), showWarnings = FALSE, recursive = TRUE)
saveRDS(dds, output_dds)
cat("\nSaved:", output_dds, "\n")

cat("\n=== DESeq2 Relevel + Wald Retest Complete ===\n")
sink(type = "message")
sink()
close(logfile)
