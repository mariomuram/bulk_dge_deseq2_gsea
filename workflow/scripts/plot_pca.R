#!/usr/bin/env Rscript
# QC: PCA plot across all samples (not contrast-specific), for one correction
# stage ("uncorrected" = blind VST; otherwise the batch-corrected VST).
# The prcomp() object is stashed so pc_field_association doesn't recompute it.

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
correction <- snakemake@wildcards[["correction"]]
color_var <- snakemake@params[["color"]]
shape_var <- snakemake@params[["shape"]]
label_var <- snakemake@params[["label"]]
palette <- snakemake@params[["palette"]]
title_suffix <- snakemake@params[["title_suffix"]]
color_label <- unname(color_var[[1]])
color_var <- names(color_var)[1]
shape_label <- unname(shape_var[[1]])
shape_var <- names(shape_var)[1]

cat("Correction stage:", correction, "\n")
cat("Plotting with: color =", color_var, ", shape =", shape_var, "\n\n")

dds <- readRDS(snakemake@input[["dds_rds"]])
cat("Loaded dds:", ncol(dds), "samples,", nrow(dds), "genes\n")
vst_mat <- readRDS(snakemake@input[["assay"]])
stopifnot(identical(colnames(vst_mat), colnames(dds)))

# Same computation as DESeq2::plotPCA() with its defaults (which this script
# used before): top 500 genes by variance, centered, unscaled. Done by hand
# because plotPCA() returns only PC1/PC2 and not the prcomp object.
ntop <- 500
gene_var <- apply(vst_mat, 1, var)
top_genes <- order(gene_var, decreasing = TRUE)[
  seq_len(min(ntop, length(gene_var)))
]
pca <- prcomp(t(vst_mat[top_genes, ]))
percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2))

fields <- unique(c(color_var, shape_var, label_var))
pca_data <- cbind(data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2]),
                  as.data.frame(colData(dds))[rownames(pca$x), fields,
                                              drop = FALSE])

pca_plot <- ggplot(pca_data,
                   aes(x = PC1, y = PC2, color = .data[[color_var]],
                       shape = .data[[shape_var]])) +
  geom_point(size = default_point_size, alpha = 0.8) +
  scale_color_brewer(palette = palette) +
  labs(
    title = paste0("PCA (all samples)", title_suffix),
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = color_label, shape = shape_label
  )
if (!is.null(label_var)) {
  library(ggrepel)
  pca_plot <- pca_plot +
    geom_label_repel(
      aes(label = .data[[label_var]]), fill = "white", alpha = .5,
      size = default_text_size, fontface = "bold", box.padding = .5,
      point.padding = 0, point.size = default_point_size, force = 3,
      force_pull = .8, max.iter = 1e6, max.time = 60, max.overlaps = Inf,
      min.segment.length = 0, segment.color = "grey40", segment.size = 0.4,
      seed = 42L, show.legend = FALSE
    )
}

output_file <- snakemake@output[["pca"]]
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
save_plot(output_file, pca_plot)
cat("Saved:", output_file, "\n")

scores_file <- snakemake@output[["scores"]]
saveRDS(pca, scores_file)
cat("Saved:", scores_file, "\n")

cat("\n=== QC: PCA Plot Complete ===\n")
sink(type = "message")
sink()
close(logfile)
