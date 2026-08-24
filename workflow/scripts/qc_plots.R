#!/usr/bin/env Rscript
# QC Plots: PCA, MA plot, Volcano plot

library(DESeq2)
library(ggplot2)
library(dplyr)
library(readr)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== QC Plots ===\n")
cat("Contrast:", snakemake@params[["contrast_name"]], "\n\n")

# Load DDS object
dds <- readRDS(snakemake@input[["dds_rds"]])

# Load results
results_df <- read_tsv(snakemake@input[["results"]])

contrast_name <- snakemake@params[["contrast_name"]]

# === PCA Plot ===
cat("Generating PCA plot...\n")
vst_data <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vst_data,
                    intgroup = c("condition",
                                 "background"),
                    returnData = TRUE)

pca_plot <- ggplot(pca_data,
                   aes(x = PC1, y = PC2,
                       color = condition,
                       shape = background)) +
  geom_point(size = 4, alpha = 0.7) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(title = paste("PCA:", contrast_name),
       color = "Condition",
       shape = "Background")

ggsave(snakemake@output[["pca"]], pca_plot,
       width = 8, height = 6, dpi = 300)
cat("✓ PCA plot saved\n")

# === MA Plot ===
cat("Generating MA plot...\n")
ma_df <- results_df %>%
  mutate(
    mean_expr = baseMean,
    color = ifelse(padj < 0.05, "Significant",
                   "Not Significant")
  ) %>%
  filter(!is.na(padj))

ma_plot <- ggplot(ma_df,
                  aes(x = mean_expr,
                      y = log2FoldChange,
                      color = color)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50") +
  scale_color_manual(
    values = c("Significant" = "red",
               "Not Significant" = "gray")
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(
    title = paste("MA Plot:", contrast_name),
    x = "Mean Log2 Expression",
    y = "Log2 Fold Change",
    color = ""
  )

ggsave(snakemake@output[["ma"]], ma_plot,
       width = 8, height = 6, dpi = 300)
cat("✓ MA plot saved\n")

# === Volcano Plot ===
cat("Generating volcano plot...\n")
volcano_df <- results_df %>%
  mutate(
    neg_log10_padj = -log10(padj),
    color = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Up",
      padj < 0.05 & log2FoldChange < 0 ~ "Down",
      TRUE ~ "Not Significant"
    )
  ) %>%
  filter(!is.na(padj))

volcano_plot <- ggplot(volcano_df,
                       aes(x = log2FoldChange,
                           y = neg_log10_padj,
                           color = color)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_vline(xintercept = c(-1, 1), 
             linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "gray50") +
  scale_color_manual(
    values = c("Up" = "#d73027",
               "Down" = "#4575b4",
               "Not Significant" = "#cccccc")
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

ggsave(snakemake@output[["volcano"]],
       volcano_plot,
       width = 8, height = 6, dpi = 300)
cat("✓ Volcano plot saved\n")

cat("\n=== QC Plots Complete ===\n")
sink(type="message")
sink()
close(logfile)
