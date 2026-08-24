#!/usr/bin/env Rscript
# QC: Volcano plot for a single contrast.

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
padj_threshold <- snakemake@params[["padj_threshold"]]
log2fc_threshold <- snakemake@params[["log2fc_threshold"]]
cat("Contrast:", contrast_name, "\n")
cat("padj threshold:", padj_threshold, "\n")
cat("log2FC threshold:", log2fc_threshold, "\n\n")

results_df <- read_tsv(snakemake@input[["results"]], show_col_types = FALSE)

n_before <- nrow(results_df)
volcano_df <- results_df %>%
  filter(!is.na(padj)) %>%
  mutate(
    neg_log10_padj = -log10(padj),
    color = case_when(
      padj < padj_threshold &
        log2FoldChange >= log2fc_threshold ~ "Up",
      padj < padj_threshold &
        log2FoldChange <= -log2fc_threshold ~ "Down",
      TRUE ~ "Not Significant"
    )
  )
cat(
  "Excluded", n_before - nrow(volcano_df),
  "gene(s) with NA padj (cannot plot -log10(NA)).\n"
)

volcano_plot <- ggplot(
  volcano_df,
  aes(x = log2FoldChange, y = neg_log10_padj, color = color)
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
    title = paste("Volcano Plot:", contrast_name),
    x = "Log2 Fold Change",
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
