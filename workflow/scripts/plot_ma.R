#!/usr/bin/env Rscript
# Plot: MA plot for a single contrast, unshrunk or shrunk LFC.

library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== Plot: MA Plot ===\n")
contrast_name <- snakemake@params[["contrast_name"]]
lfc_type <- snakemake@params[["lfc_type"]]  # "unshrunk" or "shrunk"
padj_threshold <- snakemake@params[["padj_threshold"]]
cat("Contrast:", contrast_name, "\n")
cat("LFC type:", lfc_type, "\n")
cat("padj threshold:", padj_threshold, "\n\n")

source(snakemake@input[["plot_theme_lib"]])
source(snakemake@input[["plot_input_lib"]])

n_na <- sum(is.na(results_df$padj))
if (n_na > 0) {
  cat(
    "NOTE:", n_na, "gene(s) have NA padj (independent",
    "filtering and/or Cook's outliers); shown as",
    "'Not Significant' rather than dropped, so the low-count",
    "filtering funnel stays visible on the plot.\n"
  )
}

results_df <- results_df %>%
  mutate(
    lfc_value = .data[[lfc_col]],
    color = ifelse(
      !is.na(padj) & padj < padj_threshold,
      "Significant", "Not Significant"
    )
  )

n_na <- sum(is.na(results_df$padj))
if (n_na > 0) {
  cat(
    "NOTE:", n_na, "gene(s) have NA padj (independent",
    "filtering and/or Cook's outliers); shown as",
    "'Not Significant' rather than dropped, so the low-count",
    "filtering funnel stays visible on the plot.\n"
  )
}

results_df <- results_df %>%
  mutate(
    lfc_value = .data[[lfc_col]],
    color = ifelse(
      !is.na(padj) & padj < padj_threshold,
      "Significant", "Not Significant"
    )
  )

lfc_max <- max(abs(results_df$lfc_value), na.rm = TRUE)
ma_plot <- ggplot(
  results_df,
  aes(x = baseMean, y = lfc_value, color = color, label = label)
) +
  geom_point(alpha = 0.5, size = default_point_size) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_label_repel(fill="white", alpha=.5,
    size = default_text_size, fontface = "bold", box.padding = .5,
    point.padding = 0, point.size = default_point_size, force = 3,
    force_pull = .8, max.iter = 1e6, max.time = 60,
    max.overlaps = Inf, min.segment.length = 0, segment.color = "grey40",
    segment.size = 0.4, seed = 42L, show.legend = FALSE
  ) +
  scale_x_log10() +
  # +10% breathing room past the extremes
  coord_cartesian(ylim = c(-lfc_max, lfc_max) * 1.05) +
  scale_color_manual(
    values = c("Significant" = "red", "Not Significant" = "gray")
  ) +
  labs(
    title = paste0(contrast_name, " (", lfc_type, ")"),
    x = "Mean Expression (baseMean, log10 scale)",
    y = lfc_label,
    color = paste0("padj < ", padj_threshold)
  )

output_file <- snakemake@output[["ma"]]
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
save_plot(output_file, ma_plot)
cat("Saved:", output_file, "\n")

cat("\n=== Plot: MA Plot Complete ===\n")
sink(type = "message")
sink()
close(logfile)
