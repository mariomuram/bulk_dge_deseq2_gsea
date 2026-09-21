#!/usr/bin/env Rscript
# QC: PCA plot across all samples (not contrast-specific).

library(DESeq2)
library(ggplot2)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== QC: PCA Plot ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

source(snakemake@input[["plot_theme_lib"]])

# Extract parameterized variables
color_var <- snakemake@params[["color"]]
shape_var <- snakemake@params[["shape"]]
color_label <- unname(color_var[[1]])
color_var <- names(color_var)[1]
shape_label <- unname(shape_var[[1]])
shape_var <- names(shape_var)[1]

cat("Plotting with: color =", color_var, ", shape =", shape_var, "\n\n")

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
  intgroup = c(color_var, shape_var),
  returnData = TRUE
)
percent_var <- round(100 * attr(pca_data, "percentVar"))

# TODO: Color scale! Play with scale.
pca_plot <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = .data[[color_var]],
      shape = .data[[shape_var]])
) +
  geom_point(size = default_point_size, alpha = 0.8) +
  labs(
    title = "PCA (all samples)",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = color_label, shape = shape_label
  )

output_file <- snakemake@output[["pca"]]
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
save_plot(output_file, pca_plot)
cat("Saved:", output_file, "\n")

cat("\n=== QC: PCA Plot Complete ===\n")
sink(type = "message")
sink()
close(logfile)
