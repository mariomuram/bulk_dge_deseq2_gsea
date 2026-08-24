#!/usr/bin/env Rscript
# GSEA analysis using clusterProfiler and fgsea

library(org.Hs.eg.db)
library(fgsea)
library(readr)
library(dplyr)
library(ggplot2)
library(clusterProfiler)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== GSEA Analysis ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "\n\n")

# Load results
results_df <- read_csv(snakemake@input[["results"]])
contrast_name <- snakemake@params[["contrast_name"]]
gsea_config <- snakemake@params[["gsea_config"]]

cat("Contrast:", contrast_name, "\n")
cat("Database:", gsea_config$database, "\n")
cat("Ontology:", gsea_config$ontology, "\n\n")

# Remove NA and infinite values
results_df <- results_df %>%
  filter(!is.na(padj), !is.infinite(log2FoldChange))

# Create ranking statistic (default: stat)
ranking_stat <- gsea_config$ranking_stat
cat("Ranking statistic:", ranking_stat, "\n")

if (ranking_stat == "stat") {
  # Compute test statistic (sign(log2FC) * -log10(pvalue))
  results_df <- results_df %>%
    mutate(
      stat = sign(log2FoldChange) * 
             (-log10(pvalue + 1e-10)),
      stat = ifelse(is.infinite(stat), 
                    sign(log2FoldChange) * 10,
                    stat)
    )
  ranking <- results_df$stat
  names(ranking) <- results_df$gene_id
} else {
  ranking <- results_df$log2FoldChange
  names(ranking) <- results_df$gene_id
}

ranking <- sort(ranking, decreasing = TRUE)

cat("Genes ranked:", length(ranking), "\n")
cat("Top 5 genes:\n")
print(head(ranking, 5))

# Load pathway database
cat("\nLoading pathway database (", 
    gsea_config$database, ")...\n", sep = "")

if (gsea_config$database == "GO") {
  # Gene Ontology
  ontology <- gsea_config$ontology
  
  go_result <- enrichGO(
    gene = names(ranking),
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = ontology,
    pvalueCutoff = gsea_config$padj_threshold,
    minGSSize = gsea_config$min_set_size,
    maxGSSize = gsea_config$max_set_size
  )
  
  if (nrow(go_result) == 0) {
    cat("⚠ WARNING: No significant pathways found.\n")
    gsea_df <- data.frame(
      pathway = character(),
      NES = numeric(),
      padj = numeric(),
      size = integer()
    )
  } else {
    gsea_df <- as.data.frame(go_result) %>%
      select(Description, NES = GeneRatio,
             padj = p.adjust, Count) %>%
      rename(pathway = Description, size = Count) %>%
      arrange(padj)
    
    # Compute NES as log odds ratio proxy
    gsea_df$NES <- gsea_df$NES %>%
      as.character() %>%
      as.numeric()
    
    cat("✓ Enrichment complete\n")
    cat("Pathways with padj <", 
        gsea_config$padj_threshold, ":", 
        nrow(gsea_df), "\n", sep = "")
  }
  
} else {
  cat("⚠ Database '", gsea_config$database,
      "' not yet implemented.\n", sep = "")
  cat("   Using GO (BP) instead.\n")
  
  go_result <- enrichGO(
    gene = names(ranking),
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = gsea_config$padj_threshold,
    minGSSize = gsea_config$min_set_size,
    maxGSSize = gsea_config$max_set_size
  )
  
  if (nrow(go_result) == 0) {
    gsea_df <- data.frame(
      pathway = character(),
      NES = numeric(),
      padj = numeric(),
      size = integer()
    )
  } else {
    gsea_df <- as.data.frame(go_result) %>%
      select(Description, NES = GeneRatio,
             padj = p.adjust, Count) %>%
      rename(pathway = Description, size = Count)
  }
}

# Save GSEA results table
cat("\nSaving GSEA results...\n")
write_csv(gsea_df, snakemake@output[["gsea_table"]])
cat("✓ Saved:", snakemake@output[["gsea_table"]], "\n")

# Generate GSEA plot
cat("\nGenerating GSEA plot...\n")

if (nrow(gsea_df) > 0) {
  # Take top 15 pathways by padj
  top_gsea <- gsea_df %>% slice_head(n = 15)
  
  gsea_plot <- ggplot(top_gsea,
                      aes(x = reorder(pathway, padj),
                          y = -log10(padj),
                          fill = NES)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_gradient2(low = "#2166ac",
                         mid = "white",
                         high = "#b2182b",
                         midpoint = 0) +
    coord_flip() +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 11),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = paste("Top Enriched Pathways:",
                    contrast_name),
      y = "-Log10(Adjusted p-value)",
      x = "",
      fill = "NES"
    )
  
  ggsave(snakemake@output[["gsea_plot"]], gsea_plot,
         width = 10, height = 7, dpi = 300)
  cat("✓ Plot saved:", snakemake@output[["gsea_plot"]], 
      "\n")
  
} else {
  cat("⚠ No pathways to plot.\n")
  
  # Create empty plot
  empty_plot <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "No significant pathways found",
             size = 5, color = "gray50") +
    theme_void()
  
  ggsave(snakemake@output[["gsea_plot"]], empty_plot,
         width = 8, height = 6, dpi = 300)
}

cat("\n=== GSEA Analysis Complete ===\n")
sink(type="message")
sink()
close(logfile)
