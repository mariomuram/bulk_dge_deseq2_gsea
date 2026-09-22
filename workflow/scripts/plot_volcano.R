#!/usr/bin/env Rscript
# Plot: Volcano plot for a single contrast, unshrunk or shrunk LFC.

library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== Plot: Volcano Plot ===\n")

source(snakemake@input[["plot_theme_lib"]])
source(snakemake@input[["plot_input_lib"]])
log2fc_threshold <- snakemake@params[["log2fc_threshold"]]
cat("log2FC threshold:", log2fc_threshold, "\n\n")


n_before <- nrow(results_df)
results_df <- results_df %>%
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
  "Excluded", n_before - nrow(results_df),
  "gene(s) with NA padj (cannot plot -log10(NA)).\n"
)

# TODO: Annotate up and down counts (within quadrants) + total n_genes (in title?)!
lfc_max <- max(abs(results_df$lfc_value), na.rm = TRUE)
volcano_plot <- ggplot(
  results_df,
  aes(x = lfc_value, y = neg_log10_padj, color = color, label = label)
) +
  geom_point(alpha = 0.6, size = default_point_size) +
  geom_vline(
    xintercept = c(-log2fc_threshold, log2fc_threshold),
    linetype = "dashed", color = "gray50", linewidth = default_line_width
  ) +
  geom_hline(
    yintercept = -log10(padj_threshold),
    linetype = "dashed", color = "gray50", linewidth = default_line_width
  ) +
  geom_label_repel(fill="white", alpha=.5,
    size = default_text_size, fontface = "bold", box.padding = .5,
    point.padding = 0, point.size = default_point_size, force = 3,
    force_pull = .8, max.iter = 1e6, max.time = 60,
    max.overlaps = Inf, min.segment.length = 0, segment.color = "grey40",
    segment.size = 0.4, seed = 42L, show.legend = FALSE
  ) +
  # +10% breathing room past the extremes
  coord_cartesian(xlim = c(-lfc_max, lfc_max) * 1.05) +
  scale_color_manual(
    values = c(
      # TODO: Use up and down colors!
      "Up" = "#d73027", "Down" = "#4575b4",
      "Not Significant" = "#cccccc"
    )
  ) +
  labs(
    title = paste0(contrast_name, " (", lfc_type, ")"),
    x = lfc_label,
    y = "-Log10(Adjusted p-value)",
    color = ""
  )

output_file <- snakemake@output[["volcano"]]
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
save_plot(output_file, volcano_plot)
cat("Saved:", output_file, "\n")

cat("\n=== Plot: Volcano Plot Complete ===\n")
sink(type = "message")
sink()
close(logfile)
