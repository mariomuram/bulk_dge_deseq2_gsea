#!/usr/bin/env Rscript
# GSEA (preranked) using clusterProfiler::gseGO() against the
# configured GO ontology, backed by fgsea. Ranking statistic comes
# directly from the DESeq2 results table (no reconstruction). This
# script only runs the expensive enrichment step -- it saves the
# live gseaResult object to RDS for scripts/gsea_postprocess.R to
# pick up (simplify(), tables, plots all happen there instead, so
# tweaking a plot never reruns fgsea).
library(readr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== GSEA run (GO) ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

source(snakemake@input[["ranking_lib"]])

# --- Load DESeq2 results ---------------------------------------------
results_df <- read_tsv(snakemake@input[["results"]], show_col_types = FALSE)
contrast_name <- snakemake@params[["contrast_name"]]
gsea_config <- snakemake@params[["gsea_config"]]
cat("Contrast:", contrast_name, "\n")
cat("Analysis name:", snakemake@params[["gsea_name"]], "\n")
cat("Database:", gsea_config$database, "\n")
cat("Ontology:", gsea_config$ontology, "\n\n")

if (gsea_config$database != "GO") {
  stop(
    "ERROR: gsea_run_go.R was invoked with database = '",
    gsea_config$database, "'; this script only handles GO. Check ",
    "the Snakefile's wildcard_constraints for gsea_run_go."
  )
}

# --- Build ranking (shared helper) -----------------------------------
ranking <- validate_and_rank(results_df, gsea_config)

# --- Run GSEA -----------------------------------------------------
cat("\nRunning gseGO (ontology:", gsea_config$ontology, ")...\n")
# pvalueCutoff = 1 -- fgsea computes NES/p-value for every gene set
# within the size bounds regardless of this value; it's a post-hoc
# filter on what's returned, not a computational shortcut. Keeping
# everything here (rather than gsea_config$padj_threshold) lets
# postprocess apply the significance call explicitly, and lets
# scripts/gsea_summary.R report a gene set's real NES/padj in every
# contrast it was tested in, not just where it happened to clear
# threshold. Already-computed p.adjust values for anything
# significant are unaffected -- BH correction runs over the full
# tested universe either way.
gsea_result <- gseGO(
  geneList = ranking,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = gsea_config$ontology,
  minGSSize = gsea_config$min_set_size,
  maxGSSize = gsea_config$max_set_size,
  pvalueCutoff = 1,
  by = "fgsea",
  seed = TRUE
)
n_sig <- sum(as.data.frame(gsea_result)$p.adjust < gsea_config$padj_threshold)
cat(
  "\ngseGO complete:", nrow(as.data.frame(gsea_result)),
  "gene set(s) tested,", n_sig, "significant at padj <",
  gsea_config$padj_threshold, "(significance applied in",
  "postprocess, not here)\n"
)

# --- Save gseaResult object for postprocessing ------------------------
output_rds <- snakemake@output[["gsea_result_rds"]]
dir.create(dirname(output_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(gsea_result, output_rds)
cat("Saved:", output_rds, "\n")

cat("\n=== GSEA run (GO) complete ===\n")
sink(type = "message")
sink()
close(logfile)
