#!/usr/bin/env Rscript
# Blind variance-stabilizing transform of the fitted dds, stored once as a
# plain assay matrix. Every consumer (plot_pca, batch_correct) reads this file
# instead of recomputing the transform, so they all start from identical values.

library(DESeq2)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== VST (blind) ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

dds <- readRDS(snakemake@input[["dds_rds"]])
cat("Loaded dds:", ncol(dds), "samples,", nrow(dds), "genes\n")

# blind = TRUE: QC should not use the fitted design, so sample
# structure isn't partially "explained away" before we look at
# it. (blind = FALSE is appropriate later, for design-aware
# downstream use, not for QC.)
cat("\nRunning variance-stabilizing transform (blind = TRUE)...\n")
vst_data <- vst(dds, blind = TRUE)

output_file <- snakemake@output[["vst"]]
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(assay(vst_data), output_file)
cat("Saved:", output_file, "\n")

cat("\n=== VST (blind) Complete ===\n")
sink(type = "message")
sink()
close(logfile)
