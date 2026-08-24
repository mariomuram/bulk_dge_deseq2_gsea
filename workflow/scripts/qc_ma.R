#!/usr/bin/env Rscript
# QC: MA plot for a single contrast.

library(ggplot2)
library(dplyr)
library(readr)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== QC: MA Plot ===\n")
contrast_name <- snakemake@params[["contrast_name"]]
padj_threshold <- snakemake@params[["padj_threshold"]]
cat("Contrast:", contrast_name, "\n")
cat("padj threshold:", padj_threshold, "\n\n")

results_df <- read_tsv(snakemake@input[["results"]], show_col_types = FALSE)

n_na <- sum(is.na(results_df$padj))
if (n_na > 0) {
  cat(
    "NOTE:", n_na, "gene(s) have NA padj (independent ",
    "filtering and/or Cook's outliers); shown as ",
    "'Not Significant' rather than dropped, so the low-count ",
    "filtering funnel stays visible on the plot.\n"
  )
}

ma_df <- results_df %>%
  mutate(
    color = ifelse(
      !is.na(padj) & padj < padj_threshold,
      "Significant", "Not Significant"
    )
  )

ma_plot <- ggplot(
  ma_df,
  aes(x = baseMean, y = log2FoldChange, color = color)
) +
  geom_point(alpha = 0.5, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_log10() +
  scale_color_manual(
    values = c("Significant" = "red", "Not Significant" = "gray")
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(
    title = paste("MA Plot:", contrast_name),
    x = "Mean Expression (baseMean, log10 scale)",
    y = "Log2 Fold Change",
    color = paste0("padj < ", padj_threshold)
  )

output_file <- snakemake@output[["ma"]]
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
ggsave(output_file, ma_plot, width = 8, height = 6, dpi = 300)
cat("Saved:", output_file, "\n")

cat("\n=== QC: MA Plot Complete ===\n")
sink(type = "message")
sink()
close(logfile)
