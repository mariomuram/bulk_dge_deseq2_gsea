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

contrast_name <- snakemake@params[["contrast_name"]]
gsea_name <- snakemake@params[["gsea_name"]]
database <- snakemake@params[["database"]]
gene_id_style <- snakemake@params[["gene_id_style"]]
ranking_stat <- snakemake@params[["ranking_stat"]]
min_set_size <- snakemake@params[["min_set_size"]]
max_set_size <- snakemake@params[["max_set_size"]]
padj_threshold <- snakemake@params[["padj_threshold"]]

cat("Contrast:", contrast_name, "\n")
cat("Analysis name:", gsea_name, "\n")
cat("Database:", database, "\n")
cat("Gene ID style:", gene_id_style, "\n")
cat("Ranking stat:", ranking_stat, "\n")
cat("Min set size:", min_set_size, "\n")
cat("Max set size:", max_set_size, "\n")
cat("Adjusted p-value threshold:", padj_threshold, "\n\n")

if (database != "HALLMARK") {
  stop(
    "ERROR: gsea_run_hallmark.R was invoked with database = '",
    database, "'; this script only handles HALLMARK. ",
    "Check the Snakefile's wildcard_constraints for ",
    "gsea_run_hallmark."
  )
}

# --- Load DESeq2 results ---------------------------------------------
results_df <- read_tsv(snakemake@input[["results"]], show_col_types = FALSE)

# --- Build ranking (shared helper) -----------------------------------
ranking <- validate_and_rank(results_df, ranking_stat, gene_id_style)

# --- Run GSEA -----------------------------------------------------
cat("\nFetching MSigDB HALLMARK gene sets (Homo sapiens)...\n")
hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H")
if (gene_id_style == "SYMBOL") {
  term2gene <- hallmark_sets[, c("gs_name", "gene_symbol")]
} else if (gene_id_style == "ENSEMBL") {
  library(org.Hs.eg.db)
  # Convert Ensembl IDs to symbols
  gene_conversion <- bitr(
    names(ranking),
    fromType = "ENSEMBL",
    toType = "SYMBOL",
    OrgDb = org.Hs.eg.db
  )
  # Map HALLMARK sets to union of corresponding gene symbols
  term2gene <- hallmark_sets[, c("gs_name", "gene_symbol")] %>%
    inner_join(gene_conversion, by = c("gene_symbol" = "SYMBOL"),
               relationship = "many-to-many") %>%
    distinct(gs_name, ENSEMBL)
} else {
  stop(
    "ERROR: gsea_run_hallmark.R was invoked with gene_id_style = '",
    gene_id_style, "'; this script only supports `SYMBOL` or ",
    "`ENSEMBL` style gene IDs."
  )
}
cat(
  "Loaded", length(unique(hallmark_sets$gs_name)),
  "HALLMARK gene set(s)\n"
)

cat("\nRunning GSEA (HALLMARK)...\n")
# pvalueCutoff = 1 -- see the matching comment in gsea_run_go.R;
# same reasoning applies here, database-independent.
gsea_result <- GSEA(
  geneList = ranking, TERM2GENE = term2gene, minGSSize = min_set_size,
  maxGSSize = max_set_size, pvalueCutoff = 1, by = "fgsea", seed = TRUE
)
n_sig <- sum(as.data.frame(gsea_result)$p.adjust < padj_threshold)
cat(
  "\nGSEA complete:", nrow(as.data.frame(gsea_result)),
  "gene set(s) tested,", n_sig, "significant at padj <",
  padj_threshold, "(significance applied in postprocess,",
  "not here)\n"
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
