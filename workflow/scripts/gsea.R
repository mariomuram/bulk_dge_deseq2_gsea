#!/usr/bin/env Rscript
# GSEA (preranked) using clusterProfiler::gseGO(),
# backed by fgsea. Ranking statistic comes directly from the
# DESeq2 results table (no reconstruction).

library(readr)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")

sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== GSEA ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Load DESeq2 results ---------------------------------------------

results_df <- read_tsv(snakemake@input[["results"]], show_col_types = FALSE)
contrast_name <- snakemake@params[["contrast_name"]]
gsea_config <- snakemake@params[["gsea_config"]]

cat("Contrast:", contrast_name, "\n")
cat("Analysis name:", snakemake@params[["gsea_name"]], "\n")
cat("Database:", gsea_config$database, "\n")
cat("Ontology:", gsea_config$ontology, "\n\n")

stat_col <- gsea_config$ranking_stat

cat("Ranking statistic:", stat_col, "\n")

if (!(stat_col %in% colnames(results_df))) {
  stop(
    "ERROR: gsea.ranking_stat = '", stat_col, "' is not a ",
    "column in the DESeq2 results table. Available columns: ",
    paste(colnames(results_df), collapse = ", ")
  )
}

# --- Build ranking -------------------------------------------------

# Filter on the ranking column itself, not padj: independent
# filtering leaves padj NA but stat/log2FoldChange intact, and
# preranked GSEA should use the full ranked transcriptome, not
# only genes that cleared the DE significance test.
n_before <- nrow(results_df)
results_df <- results_df %>%
  filter(
    !is.na(.data[[stat_col]]),
    !is.infinite(.data[[stat_col]])
  )

cat(
  "Dropped", n_before - nrow(results_df),
  "gene(s) with NA/infinite", stat_col, "\n"
)

if (any(duplicated(results_df$gene_id))) {

  n_dup <- sum(duplicated(results_df$gene_id))

  cat(
    "WARNING:", n_dup, "duplicate gene_id(s) found; keeping",
    "the entry with the largest |", stat_col, "| for each.",
    "fgsea requires unique gene names.\n"
  )

  results_df <- results_df %>%
    mutate(.abs_stat = abs(.data[[stat_col]])) %>%
    arrange(desc(.abs_stat)) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    select(-.abs_stat)
}

ranking <- results_df[[stat_col]]
names(ranking) <- results_df$gene_id
ranking <- sort(ranking, decreasing = TRUE)

cat("\nGenes ranked:", length(ranking), "\n")
cat("Top 5 genes:\n")
print(head(ranking, 5))
cat("Bottom 5 genes:\n")
print(tail(ranking, 5))

# --- Run GSEA -----------------------------------------------------

if (gsea_config$database == "GO") {

  library(org.Hs.eg.db)
  library(GOSemSim)

  cat(
    "\nRunning gseGO (database:", gsea_config$database,
    ", ontology:", gsea_config$ontology, ")...\n"
  )

  if (gsea_config$database != "GO") {

    cat(
      "WARNING: database '", gsea_config$database, "' is not",
      "yet implemented; running GO instead.\n", sep = ""
    )

  }

  gsea_result <- gseGO(
    geneList = ranking,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = gsea_config$ontology,
    minGSSize = gsea_config$min_set_size,
    maxGSSize = gsea_config$max_set_size,
    pvalueCutoff = gsea_config$padj_threshold,
    by = "fgsea",
    seed = TRUE
  )

} else if (gsea_config$database == "HALLMARK") {

  library(msigdbr)

  hallmark_sets <- msigdbr(
    species = "Homo sapiens",
    category = "H"
  )

  hallmark_list <- split(
    hallmark_sets$gene_symbol,
    hallmark_sets$gs_name
  )

  gsea_result <- GSEA(
    geneList = ranking,
    TERM2GENE = hallmark_sets[, c("gs_name", "gene_symbol")],
    minGSSize = gsea_config$min_set_size,
    maxGSSize = gsea_config$max_set_size,
    pvalueCutoff = gsea_config$padj_threshold,
    by = "fgsea",
    seed = TRUE
  )

} else {

  stop(
    "Unsupported database: ",
    gsea_config$database
  )

}

# --- Reduce GO term redundancy (GO only) ----------------------------

# GO:BP terms nest heavily (parent/child terms with highly
# overlapping gene membership); simplify() collapses semantically
# similar terms via GOSemSim, keeping the most significant
# representative from each redundant cluster. This runs before the
# results table and both plots are built, so all downstream outputs
# reflect the de-duplicated term set. Not applicable to HALLMARK,
# which isn't graph-structured the way GO is.
if (gsea_config$database == "GO") {

  n_terms_before <- nrow(as.data.frame(gsea_result))

  cat(
    "\nSimplifying redundant GO terms (cutoff = 0.7, by =",
    "'p.adjust', select_fun = min)...\n"
  )

  gsea_result <- simplify(
    gsea_result,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  n_terms_after <- nrow(as.data.frame(gsea_result))

  cat(
    "GO terms:", n_terms_before, "before simplification,",
    n_terms_after, "after (", n_terms_before - n_terms_after,
    "collapsed into a more significant parent/child term).\n"
  )

}

gsea_df_raw <- as.data.frame(gsea_result)

if (nrow(gsea_df_raw) == 0) {

  cat(
    "WARNING: no significantly enriched pathways found",
    "(padj <", gsea_config$padj_threshold, "). Check gene set",
    "size bounds, the ranking statistic used, and whether the",
    "underlying DE contrast had many significant genes to",
    "begin with.\n"
  )

  gsea_df <- data.frame(
    pathway = character(), NES = numeric(),
    padj = numeric(), size = integer()
  )

} else {

  gsea_df <- gsea_df_raw %>%
    transmute(
      pathway = Description,
      NES = NES,
      padj = p.adjust,
      size = setSize
    ) %>%
    arrange(padj)

  cat(
    "Pathways with padj <", gsea_config$padj_threshold, ":",
    nrow(gsea_df), "\n"
  )

}

# --- Save GSEA results (gzip --best compressed TSV) --------------------

cat("\nSaving GSEA results...\n")

output_table <- snakemake@output[["gsea_table"]]
dir.create(dirname(output_table), showWarnings = FALSE, recursive = TRUE)
con <- gzfile(output_table, open = "wb", compression = 9)
write_tsv(gsea_df, con)
close(con)

cat("Saved:", output_table, "\n")

# --- Plot top enriched pathways ------------------------------------

cat("\nGenerating GSEA plot...\n")

if (nrow(gsea_df) > 0) {

  top_gsea <- gsea_df %>% slice_head(n = 15)

  gsea_plot <- ggplot(
    top_gsea,
    aes(x = reorder(pathway, padj), y = -log10(padj), fill = NES)
  ) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0
    ) +
    coord_flip() +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 11),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = paste("Top Enriched Pathways:", contrast_name),
      y = "-Log10(Adjusted p-value)",
      x = "",
      fill = "NES"
    )

  ggsave(
    snakemake@output[["gsea_plot"]], gsea_plot,
    width = 10, height = 7, dpi = 300
  )

  cat("Plot saved:", snakemake@output[["gsea_plot"]], "\n")

} else {

  cat("No pathways to plot.\n")

  empty_plot <- ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5,
      label = "No significant pathways found",
      size = 5, color = "gray50"
    ) +
    theme_void()

  ggsave(
    snakemake@output[["gsea_plot"]], empty_plot,
    width = 8, height = 6, dpi = 300
  )

}

# --- Ridge plot of top pathways -------------------------------------

cat("\nGenerating GSEA ridge plot...\n")

if (nrow(gsea_df) > 0) {

  ridge_n <- min(15, nrow(gsea_df))

  gsea_ridge <- ridgeplot(gsea_result, showCategory = ridge_n) +
    labs(title = paste("GSEA Ridge Plot:", contrast_name)) +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(size = 12, face = "bold")
    )

  ggsave(
    snakemake@output[["gsea_ridge_plot"]], gsea_ridge,
    width = 10, height = 8, dpi = 300
  )

  cat("Ridge plot saved:", snakemake@output[["gsea_ridge_plot"]], "\n")
} else {

  cat("No pathways to plot (ridge).\n")
  ggsave(

    snakemake@output[["gsea_ridge_plot"]], empty,
    width = 8, height = 6, dpi = 300
  )

}

cat("\n=== GSEA Complete ===\n")

sink(type = "message")
sink()
close(logfile)
