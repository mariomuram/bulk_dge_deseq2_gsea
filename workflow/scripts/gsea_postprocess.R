#!/usr/bin/env Rscript
# Post-processes a live gseaResult object (produced by
# gsea_run_go.R or gsea_run_hallmark.R): GO-only redundancy
# reduction via simplify(), the results table(s), the top-pathway
# bar plot, and the ridge plot. Shared across both databases --
# collateral reruns here are cheap (seconds), so this stays one
# script with a database branch rather than being split further;
# the expensive fgsea run is what got isolated into gsea_run_*.R.
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

cat("=== GSEA postprocess ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

contrast_name <- snakemake@params[["contrast_name"]]
gsea_config <- snakemake@params[["gsea_config"]]
cat("Contrast:", contrast_name, "\n")
cat("Analysis name:", snakemake@params[["gsea_name"]], "\n")
cat("Database:", gsea_config$database, "\n\n")

# --- Load the enrichment result --------------------------------------
gsea_result <- readRDS(snakemake@input[["gsea_result_rds"]])
cat(
  "Loaded gseaResult:", nrow(as.data.frame(gsea_result)),
  "term(s) at padj <", gsea_config$padj_threshold, "\n"
)

# --- Helper: gseaResult -> tidy data frame ----------------------------
tidy_gsea <- function(x, padj_threshold) {
  raw <- as.data.frame(x)
  if (nrow(raw) == 0) {
    cat(
      "WARNING: no significantly enriched pathways found",
      "(padj <", padj_threshold, "). Check gene set size",
      "bounds, the ranking statistic used, and whether the",
      "underlying DE contrast had many significant genes to",
      "begin with.\n"
    )
    return(data.frame(
      pathway = character(), NES = numeric(),
      padj = numeric(), size = integer()
    ))
  }
  raw %>%
    transmute(
      pathway = Description,
      NES = NES,
      padj = p.adjust,
      size = setSize
    ) %>%
    arrange(padj)
}

# --- Save the raw (unfiltered-for-redundancy) table -------------------

# This is the canonical table: identical in shape/meaning for GO
# and HALLMARK, and the one scripts/gsea_summary.R reads across
# contrasts. It must NOT depend on GO-only redundancy reduction --
# which of two near-duplicate GO terms survives simplify() is
# decided by that contrast's own p-values, so using the simplified
# object here would let the SAME semantic signal show up under
# different term IDs in different contrasts, fragmenting pathway
# identity right where the comparison script's union needs it most.
gsea_df <- tidy_gsea(gsea_result, gsea_config$padj_threshold)
cat("\nPathways in raw table:", nrow(gsea_df), "\n")

cat("\nSaving raw GSEA table...\n")
output_table <- snakemake@output[["gsea_table"]]
dir.create(dirname(output_table), showWarnings = FALSE, recursive = TRUE)
con <- gzfile(output_table, open = "wb", compression = 9)
write_tsv(gsea_df, con)
close(con)
cat("Saved:", output_table, "\n")

# --- GO-only: redundancy reduction for human-facing table & plots -----

# GO:BP terms nest heavily (parent/child terms with highly
# overlapping gene membership); simplify() collapses semantically
# similar terms via GOSemSim, keeping the most significant
# representative from each redundant cluster. This is deliberately
# NOT what feeds gsea_table above -- but it IS what the plots and
# the simplified table are built from, since decluttering a single
# contrast's own display is exactly what simplify() is for.
plot_source <- gsea_result
plot_df <- gsea_df

if (gsea_config$database == "GO") {
  library(org.Hs.eg.db)
  library(GOSemSim)

  n_terms_before <- nrow(as.data.frame(gsea_result))
  cat(
    "\nSimplifying redundant GO terms (cutoff = 0.7, by =",
    "'p.adjust', select_fun = min)...\n"
  )
  gsea_result_simplified <- simplify(
    gsea_result,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  n_terms_after <- nrow(as.data.frame(gsea_result_simplified))
  cat(
    "GO terms:", n_terms_before, "before simplification,",
    n_terms_after, "after (", n_terms_before - n_terms_after,
    "collapsed into a more significant parent/child term).\n"
  )

  simplified_df <- tidy_gsea(
    gsea_result_simplified, gsea_config$padj_threshold
  )

  cat("\nSaving simplified GO table...\n")
  output_table_simplified <- snakemake@output[["gsea_table_simplified"]]
  dir.create(
    dirname(output_table_simplified),
    showWarnings = FALSE, recursive = TRUE
  )
  con <- gzfile(output_table_simplified, open = "wb", compression = 9)
  write_tsv(simplified_df, con)
  close(con)
  cat("Saved:", output_table_simplified, "\n")

  plot_source <- gsea_result_simplified
  plot_df <- simplified_df
} else {
  cat("\nNot GO -- no redundancy reduction, no simplified table.\n")
}

# --- Plot top enriched pathways ------------------------------------
cat("\nGenerating GSEA plot...\n")
if (nrow(plot_df) > 0) {
  top_gsea <- plot_df %>% slice_head(n = 15)
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
if (nrow(plot_df) > 0) {
  ridge_n <- min(15, nrow(plot_df))
  gsea_ridge <- ridgeplot(plot_source, showCategory = ridge_n) +
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
  empty_ridge <- ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5,
      label = "No significant pathways found",
      size = 5, color = "gray50"
    ) +
    theme_void()
  ggsave(
    snakemake@output[["gsea_ridge_plot"]], empty_ridge,
    width = 8, height = 6, dpi = 300
  )
}

cat("\n=== GSEA postprocess complete ===\n")
sink(type = "message")
sink()
close(logfile)
