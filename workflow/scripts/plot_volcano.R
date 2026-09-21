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
  # TODO: Review repel params (see P_dox_vs_nodox)!
  geom_text_repel(
    max.overlaps = Inf,
    min.segment.length = 0,
    # TODO: Darken! (+check MA) [try `color_none`?]
    segment.color = "grey30",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      # TODO: Use up and down colors!
      "Up" = "#d73027", "Down" = "#4575b4",
      "Not Significant" = "#cccccc"
    )
  ) +
  labs(
    title = paste0("Volcano Plot: ", contrast_name, " (", lfc_type, ")"),
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
