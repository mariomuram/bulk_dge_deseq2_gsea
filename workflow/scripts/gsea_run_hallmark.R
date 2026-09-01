#!/usr/bin/env Rscript
# GSEA (preranked) using clusterProfiler::GSEA() against MSigDB
# HALLMARK gene sets, backed by fgsea. Ranking statistic comes
# directly from the DESeq2 results table (no reconstruction). This
# script only runs the expensive enrichment step -- it saves the
# live gseaResult object to RDS for scripts/gsea_postprocess.R to
# pick up (tables and plots happen there instead, so tweaking a
# plot never reruns fgsea).
library(readr)
library(dplyr)
library(clusterProfiler)
library(msigdbr)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== GSEA run (HALLMARK) ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

source(snakemake@input[["ranking_lib"]])

# --- Load DESeq2 results ---------------------------------------------
results_df <- read_tsv(snakemake@input[["results"]], show_col_types = FALSE)
contrast_name <- snakemake@params[["contrast_name"]]
gsea_config <- snakemake@params[["gsea_config"]]
cat("Contrast:", contrast_name, "\n")
cat("Analysis name:", snakemake@params[["gsea_name"]], "\n")
cat("Database:", gsea_config$database, "\n\n")

if (gsea_config$database != "HALLMARK") {
  stop(
    "ERROR: gsea_run_hallmark.R was invoked with database = '",
    gsea_config$database, "'; this script only handles HALLMARK. ",
    "Check the Snakefile's wildcard_constraints for ",
    "gsea_run_hallmark."
  )
}

# --- Build ranking (shared helper) -----------------------------------
ranking <- validate_and_rank(results_df, gsea_config)

# --- Run GSEA -----------------------------------------------------
cat("\nFetching MSigDB HALLMARK gene sets (Homo sapiens)...\n")
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
cat(
  "Loaded", length(unique(hallmark_sets$gs_name)),
  "HALLMARK gene set(s)\n"
)

cat("\nRunning GSEA (HALLMARK)...\n")
gsea_result <- GSEA(
  geneList = ranking,
  TERM2GENE = hallmark_sets[, c("gs_name", "gene_symbol")],
  minGSSize = gsea_config$min_set_size,
  maxGSSize = gsea_config$max_set_size,
  pvalueCutoff = gsea_config$padj_threshold,
  by = "fgsea",
  seed = TRUE
)
cat(
  "\nGSEA complete:", nrow(as.data.frame(gsea_result)),
  "significant term(s) at padj <", gsea_config$padj_threshold, "\n"
)

# --- Save gseaResult object for postprocessing ------------------------
output_rds <- snakemake@output[["gsea_result_rds"]]
dir.create(dirname(output_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(gsea_result, output_rds)
cat("Saved:", output_rds, "\n")

cat("\n=== GSEA run (HALLMARK) complete ===\n")
sink(type = "message")
sink()
close(logfile)
