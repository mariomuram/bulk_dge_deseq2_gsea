#!/usr/bin/env Rscript
# QC: PCA plot across all samples (not contrast-specific).
# Color = dox status, shape = background.

library(DESeq2)
library(ggplot2)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== QC: PCA Plot ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

dds <- readRDS(snakemake@input[["dds_rds"]])
cat("Loaded dds:", ncol(dds), "samples,", nrow(dds), "genes\n")

# blind = TRUE: QC should not use the fitted design, so sample
# structure isn't partially "explained away" before we look at
# it. (blind = FALSE is appropriate later, for design-aware
# downstream use, not for this exploratory plot.)
cat("\nRunning variance-stabilizing transform (blind = TRUE)...\n")
vst_data <- vst(dds, blind = TRUE)

pca_data <- plotPCA(
  vst_data,
  intgroup = c("dox", "background"),
  returnData = TRUE
)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = dox, shape = background)
) +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(
    title = "PCA (all samples)",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Dox status",
    shape = "Background"
  )

output_file <- snakemake@output[["pca"]]
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
ggsave(output_file, pca_plot, width = 8, height = 6, dpi = 300)
cat("Saved:", output_file, "\n")

cat("\n=== QC: PCA Plot Complete ===\n")
sink(type = "message")
sink()
close(logfile)
