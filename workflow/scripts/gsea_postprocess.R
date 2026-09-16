#!/usr/bin/env Rscript
# Post-processes a live gseaResult object (produced by
# gsea_run_go.R or gsea_run_hallmark.R). The run scripts now use
# pvalueCutoff = 1, so this script -- not the enrichment call --
# is where the padj_threshold significance decision is made. It
# writes three tables: every tested gene set (gsea_table_all, the
# input to scripts/gsea_summary.R), this contrast's significant
# subset (gsea_table_signif), and for GO only, that subset with
# redundant terms collapsed (gsea_table_simplified). It also builds
# the bar and ridge plots from the significant (GO: simplified) subset.
# Shared across both databases -- collateral reruns here are
# cheap, so this stays one script with a database branch.

GENE_SETS_OF_INTEREST_FALLBACK_N <- 15L

library(readr)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(forcats)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== GSEA postprocess ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

contrast_name <- snakemake@params[["contrast_name"]]
gsea_name <- snakemake@params[["gsea_name"]]
database <- snakemake@params[["database"]]
padj_threshold <- snakemake@params[["padj_threshold"]]
gene_sets_of_interest <- snakemake@params[["gene_sets_of_interest"]]

cat("Contrast:", contrast_name, "\n")
cat("Analysis name:", gsea_name, "\n")
cat("Database:", database, "\n\n")

# --- Load the enrichment result --------------------------------------
gsea_result <- readRDS(snakemake@input[["gsea_result_rds"]])
cat(
  "Loaded gseaResult:", nrow(as.data.frame(gsea_result)),
  "gene set(s) tested (all, unfiltered by significance)\n"
)

# --- Helper: gseaResult -> tidy data frame ----------------------------
tidy_gsea <- function(x) {
  raw <- as.data.frame(x)
  if (nrow(raw) == 0) {
    return(data.frame(
      id = character(), description = character(),
      NES = numeric(), padj = numeric(), size = integer()
    ))
  }
  raw %>%
    transmute(
      id = ID,
      description = Description,
      NES = NES,
      padj = p.adjust,
      size = setSize
    ) %>%
    arrange(padj)
}

# --- Save the "all tested gene sets" table -----------------------------

# This is the canonical table for cross-contrast comparison
# (scripts/gsea_summary.R). It is NOT restricted to this contrast's
# own significance call, so a gene set that's significant in one
# contrast and just short of threshold in another -- with a
# near-identical NES either way -- doesn't look like it's specific
# to the former. That's exactly the artifact a per-contrast
# pvalueCutoff would otherwise bake in.
gsea_df_all <- tidy_gsea(gsea_result)
if (nrow(gsea_df_all) == 0) {
  cat(
    "WARNING: no gene sets were tested at all -- check",
    "min_set_size/max_set_size against how many genes survived",
    "ranking. This is independent of padj_threshold: it means",
    "nothing was even in scope to test, not that nothing was",
    "significant.\n"
  )
}
cat("\nGene sets tested (all):", nrow(gsea_df_all), "\n")

cat("\nSaving all-gene-sets GSEA table...\n")
output_table_all <- snakemake@output[["gsea_table_all"]]
dir.create(dirname(output_table_all), showWarnings = FALSE, recursive = TRUE)
con <- gzfile(output_table_all, open = "wb", compression = 9)
write_tsv(gsea_df_all, con)
close(con)
cat("Saved:", output_table_all, "\n")

# --- Derive the significant-only table ---------------------------

# gsea_run_*.R runs with pvalueCutoff = 1, so this contrast's own
# significance call has to be applied explicitly here rather than
# inherited from the enrichment call. filter() has a registered
# method for gseaResult objects (returns a gseaResult, not a plain
# data frame), so simplify()/ridgeplot() below still work on it.
gsea_result_sig <- gsea_result %>%
  filter(p.adjust < padj_threshold)
gsea_df <- tidy_gsea(gsea_result_sig)
if (nrow(gsea_df) == 0) {
  cat(
    "WARNING: no significantly enriched pathways found",
    "(padj <", padj_threshold, ") among the", nrow(gsea_df_all),
    "gene set(s) tested. Check the ranking statistic used and",
    "whether the underlying DE contrast had many significant genes",
    "to begin with.\n"
  )
}
cat("\nPathways significant at padj <", padj_threshold, ":",
    nrow(gsea_df), "\n")

cat("\nSaving significant-only GSEA table...\n")
output_table <- snakemake@output[["gsea_table_signif"]]
dir.create(dirname(output_table), showWarnings = FALSE, recursive = TRUE)
con <- gzfile(output_table, open = "wb", compression = 9)
write_tsv(gsea_df, con)
close(con)
cat("Saved:", output_table, "\n")

# --- GO-only: redundancy reduction for human-facing table & plots -----

# GO:BP terms nest heavily (parent/child terms with highly
# overlapping gene membership); simplify() collapses semantically
# similar terms via GOSemSim, keeping the most significant
# representative from each redundant cluster. Runs on the
# significant subset, not the full tested universe -- clustering
# thousands of untested-as-relevant terms would be both wasteful
# and meaningless. This feeds the plots and the simplified table
# only; gsea_table_all and gsea_table_signif above are unaffected.
gsea_df_top_source <- gsea_df

if (database == "GO") {
  library(org.Hs.eg.db)
  library(GOSemSim)

  if (nrow(gsea_df) == 0) {
    cat(
      "\nNo significant GO terms to simplify; writing an empty",
      "simplified table.\n"
    )
    simplified_df <- gsea_df
  } else {
    n_terms_before <- nrow(gsea_df)
    cat(
      "\nSimplifying redundant GO terms (cutoff = 0.7, by =",
      "'p.adjust', select_fun = min)...\n"
    )
    gsea_result_simplified <- simplify(
      gsea_result_sig,
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
    simplified_df <- tidy_gsea(gsea_result_simplified)
  }

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

  gsea_df_top_source <- simplified_df
} else {
  cat("\nNot GO -- no redundancy reduction, no simplified table.\n")
}

# --- Resolve and validate gene sets of interest ---------------------
# Any configured ID not found is warned and omitted, and if that leaves
# nothing, this contrast's job stops rather than silently rendering
# an empty, misleadingly-labeled "no significant pathways" plot.
cat("\nResolving gene sets of interest...\n")
if (!is.null(gene_sets_of_interest) && length(gene_sets_of_interest) > 0) {
  tested_ids <- gsea_df_all$id
  missing_ids <- setdiff(gene_sets_of_interest, tested_ids)
  if (length(missing_ids) > 0) {
    cat(
      "WARNING:", length(missing_ids), "configured gene set ID(s)",
      "not found in this contrast's tested universe; omitting:",
      paste(missing_ids, collapse = ", "), "\n"
    )
  }
  gene_sets_of_interest <- intersect(gene_sets_of_interest, tested_ids)
  if (length(gene_sets_of_interest) == 0) {
    stop(
      "ERROR: none of the configured `gene_sets_of_interest` for '",
      gsea_name, "' were found among the ", nrow(gsea_df_all),
      " gene set(s) tested for contrast '", contrast_name, "'. ",
      "Check the ID(s) in `config[\"gsea\"][\"", gsea_name, "\"]",
      "[\"gene_sets_of_interest\"]` against the id column of ",
      "`gsea_table_all` for this contrast."
    )
  }
  cat(
    "Using", length(gene_sets_of_interest), "configured gene set(s)",
    "of interest.\n"
  )
} else {
  cat(
    "No gene sets of interest configured; falling back to the top",
    GENE_SETS_OF_INTEREST_FALLBACK_N, "by padj...\n"
  )
  gene_sets_of_interest <- gsea_df_top_source %>%
    slice_head(n = GENE_SETS_OF_INTEREST_FALLBACK_N) %>%
    pull(id)
}

# Unfiltered basis (all tested, not just this contrast's own
# significant subset): a gene set the person explicitly asked for
# should render with its real NES/padj even when it falls short of
# padj_threshold here -- the same principle gsea_summary.R already
# applies cross-contrast.
gsea_df_gene_sets_of_interest <- gsea_df_all %>%
  filter(id %in% gene_sets_of_interest)

# --- Leading-edge driver genes for pathways of interest ---------------
# core_enrichment lives only on the live gseaResult (tidy_gsea() never
# extracts it, since none of the on-disk tables want it) -- join it
# straight onto gsea_df_gene_sets_of_interest rather than keeping a
# second, parallel ID-filtered copy around. The join's left side is
# already restricted to gene_sets_of_interest, so it also does the
# filtering -- no separate filter() on gsea_result needed here.
cat("\nBuilding driver-gene table...\n")
gsea_df_gene_sets_of_interest <- gsea_df_gene_sets_of_interest %>%
  left_join(
    as.data.frame(gsea_result) %>% transmute(id = ID, core_enrichment),
    by = "id"
  )

core_enrichment_clean <- ifelse(
  is.na(gsea_df_gene_sets_of_interest$core_enrichment), "",
  gsea_df_gene_sets_of_interest$core_enrichment
)
driver_gene_lists <- strsplit(core_enrichment_clean, "/", fixed = TRUE)
gsea_df_gene_sets_of_interest$count <- lengths(driver_gene_lists)
n_empty <- sum(gsea_df_gene_sets_of_interest$count == 0)
if (n_empty > 0 && nrow(gsea_df_gene_sets_of_interest) > 0) {
  cat(
    "WARNING:", n_empty, "gene set(s) of interest have no",
    "leading-edge genes recorded (empty core_enrichment).\n"
  )
}

if (nrow(gsea_df_gene_sets_of_interest) > 0) {
  driver_genes_df <- data.frame(
    id = rep(
      gsea_df_gene_sets_of_interest$id,
      times = gsea_df_gene_sets_of_interest$count
    ),
    description = rep(
      gsea_df_gene_sets_of_interest$description,
      times = gsea_df_gene_sets_of_interest$count
    ),
    gene = unlist(driver_gene_lists, use.names = FALSE),
    stringsAsFactors = FALSE
  )
} else {
  driver_genes_df <- data.frame(
    id = character(), description = character(), gene = character()
  )
}
cat("Driver-gene rows:", nrow(driver_genes_df), "\n")

output_driver_genes <- snakemake@output[["gsea_driver_genes"]]
dir.create(dirname(output_driver_genes), showWarnings = FALSE,
           recursive = TRUE)
con <- gzfile(output_driver_genes, open = "wb", compression = 9)
write_tsv(driver_genes_df, con)
close(con)
cat("Saved:", output_driver_genes, "\n")

# --- Plot pathways of interest ------------------------------------
cat("\nGenerating GSEA plot...\n")
if (nrow(gsea_df_gene_sets_of_interest) > 0) {
  gsea_plot <- ggplot(
    gsea_df_gene_sets_of_interest,
    aes(x = reorder(description, padj), y = -log10(padj), fill = NES)
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

# --- Ridge plot of pathways of interest -----------------------------
cat("\nGenerating GSEA ridge plot...\n")
if (nrow(gsea_df_gene_sets_of_interest) > 0) {
  gsea_ridge <- ridgeplot(gsea_result, showCategory = gene_sets_of_interest) +
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

# --- Lollipop plot of pathways of interest ---------------------------
cat("\nGenerating GSEA lollipop plot...\n")
if (nrow(gsea_df_gene_sets_of_interest) > 0) {
  gsea_lollipop <- gsea_df_gene_sets_of_interest %>%
    ggplot(aes(NES, fct_reorder(description, NES))) +
    geom_segment(aes(xend = 0, yend = description)) +
    geom_point(aes(color = padj, size = count)) +
    scale_color_viridis_c(guide = guide_colorbar(reverse = TRUE)) +
    scale_size_continuous(range = c(1, 7)) +
    theme_minimal() +
    xlab("Normalized Enrichment Score") +
    ylab(NULL)
  ggsave(
    snakemake@output[["gsea_lollipop_plot"]], gsea_lollipop,
    width = 10, height = 8, dpi = 300
  )
  cat("Lollipop plot saved:", snakemake@output[["gsea_lollipop_plot"]], "\n")
} else {
  cat("No pathways to plot (lollipop).\n")
  empty_lollipop <- ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5,
      label = "No significant pathways found",
      size = 5, color = "gray50"
    ) +
    theme_void()
  ggsave(
    snakemake@output[["gsea_lollipop_plot"]], empty_lollipop,
    width = 8, height = 6, dpi = 300
  )
}

cat("\n=== GSEA postprocess complete ===\n")
sink(type = "message")
sink()
close(logfile)
