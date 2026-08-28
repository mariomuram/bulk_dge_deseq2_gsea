#!/usr/bin/env Rscript
# QC: Volcano plot for a single contrast, unshrunk or shrunk LFC.

library(ggplot2)
library(dplyr)
library(readr)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== QC: Volcano Plot ===\n")
contrast_name <- snakemake@params[["contrast_name"]]
lfc_type <- snakemake@params[["lfc_type"]]
padj_threshold <- snakemake@params[["padj_threshold"]]
log2fc_threshold <- snakemake@params[["log2fc_threshold"]]
cat("Contrast:", contrast_name, "\n")
cat("LFC type:", lfc_type, "\n")
cat("padj threshold:", padj_threshold, "\n")
cat("log2FC threshold:", log2fc_threshold, "\n\n")

results_df <- read_tsv(snakemake@input[["results"]], show_col_types = FALSE)

lfc_col <- if (lfc_type == "shrunk") "log2FoldChange_shrunk" else "log2FoldChange"
if (!(lfc_col %in% colnames(results_df))) {
  stop("ERROR: column '", lfc_col, "' not found in results table.")
}
if (lfc_type == "shrunk" && all(is.na(results_df[[lfc_col]]))) {
  stop(
    "ERROR: log2FoldChange_shrunk is entirely NA for this ",
    "contrast -- it should not have been requested as a ",
    "shrinkable target. Check the Snakefile's ",
    "SHRINKABLE_CONTRASTS derivation."
  )
}

n_before <- nrow(results_df)
volcano_df <- results_df %>%
  filter(!is.na(padj)) %>%
  mutate(
    lfc_value = .data[[lfc_col]],
    neg_log10_padj = -log10(padj),
    color = case_when(
      padj < padj_threshold & lfc_value >= log2fc_threshold ~ "Up",
      padj < padj_threshold & lfc_value <= -log2fc_threshold ~ "Down",
      TRUE ~ "Not Significant"
    )
  )
cat(
  "Excluded", n_before - nrow(volcano_df),
  "gene(s) with NA padj (cannot plot -log10(NA)).\n"
)

lfc_label <- if (lfc_type == "shrunk") {
  "Log2 Fold Change (apeglm-shrunk)"
} else {
  "Log2 Fold Change (unshrunk MLE)"
}

volcano_plot <- ggplot(
  volcano_df,
  aes(x = lfc_value, y = neg_log10_padj, color = color)
) +
  geom_point(alpha = 0.6, size = 2) +
  geom_vline(
    xintercept = c(-log2fc_threshold, log2fc_threshold),
    linetype = "dashed", color = "gray50"
  ) +
  geom_hline(
    yintercept = -log10(padj_threshold),
    linetype = "dashed", color = "gray50"
  ) +
  scale_color_manual(
    values = c(
      "Up" = "#d73027", "Down" = "#4575b4",
      "Not Significant" = "#cccccc"
    )
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(
    title = paste0("Volcano Plot: ", contrast_name, " (", lfc_type, ")"),
    x = lfc_label,
    y = "-Log10(Adjusted p-value)",
    color = ""
  )

output_file <- snakemake@output[["volcano"]]
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
ggsave(output_file, volcano_plot, width = 8, height = 6, dpi = 300)
cat("Saved:", output_file, "\n")

cat("\n=== QC: Volcano Plot Complete ===\n")
sink(type = "message")
sink()
close(logfile)
