#!/usr/bin/env Rscript

# SPDX-FileCopyrightText: © 2026 Marcel Schilling <mschilling@idibell.cat>
# SPDX-License-Identifier: AGPL-3.0-or-later


#########
# About #
#########

# File:     workflow/scripts/gsea_run_hallmark.R
# Author:   Marcel Schilling <mschilling@idibell.cat>
# Created:  2026-09-01
# Modified: 2026-10-15
# License:  GNU Affero General Public License version 3 or later
# Purpose:  HALLMARK GSEA run R script for the bulk RNA-seq DGE (DESeq2 +
#           GSEA) Snakemake workflow.


###############
# Description #
###############

# GSEA (preranked) using clusterProfiler::GSEA() against MSigDB
# HALLMARK gene sets, backed by fgsea. Ranking statistic comes
# directly from the DESeq2 results table (no reconstruction). This
# script only runs the expensive enrichment step -- it saves the
# live gseaResult object to RDS for scripts/gsea_postprocess.R to
# pick up (tables and plots happen there instead, so tweaking a
# plot never reruns fgsea).


###########
# License #
###########

# HALLMARK GSEA run R script for bulk RNA-seq DGE (DESeq2 + GSEA) Snakemake
# workflow.
# Copyright (C) 2026  Marcel Schilling
#
# This file is part of the bulk RNA-seq DGE (DESeq2 + GSEA) Snakemake
# workflow.
# The bulk RNA-seq DGE (DESeq2 + GSEA) Snakemake workflow is free software:
# you can redistribute it and/or modify it under the terms of the GNU
# Affero General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


###########
# Imports #
###########

library(readr)
library(dplyr)
library(clusterProfiler)
library(msigdbr)


###########
# Logging #
###########

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== GSEA run (HALLMARK) ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


###############
# Shared code #
###############

source(snakemake@input[["ranking_lib"]])


##############
# Parameters #
##############

contrast_name <- snakemake@params[["contrast_name"]]
gsea_name <- snakemake@params[["gsea_name"]]
database <- snakemake@params[["database"]]
gene_id_style <- snakemake@params[["gene_id_style"]]
ranking_stat <- snakemake@params[["ranking_stat"]]
min_set_size <- snakemake@params[["min_set_size"]]
max_set_size <- snakemake@params[["max_set_size"]]
deseq2_results <- snakemake@input[["results"]]
output_rds <- snakemake@output[["gsea_results_rds"]]

cat("Contrast:", contrast_name, "\n")
cat("Analysis name:", gsea_name, "\n")
cat("Database:", database, "\n")
cat("Gene ID style:", gene_id_style, "\n")
cat("Ranking stat:", ranking_stat, "\n")
cat("Min set size:", min_set_size, "\n")
cat("Max set size:", max_set_size, "\n\n")

if (database != "HALLMARK") {
  stop(
    "ERROR: gsea_run_hallmark.R was invoked with database = '",
    database, "'; this script only handles HALLMARK. ",
    "Check the Snakefile's wildcard_constraints for ",
    "gsea_run_hallmark."
  )
}


#########
# Input #
#########

# Load DESeq2 results.
results_df <- read_tsv(deseq2_results, show_col_types = FALSE)


################
# Gene ranking #
################

# Build ranking (shared helper).
ranking <- validate_and_rank(results_df, ranking_stat, gene_id_style)


#################
# HALLMARK sets #
#################

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


############
# Run GSEA #
############

cat("\nRunning GSEA (HALLMARK)...\n")
# pvalueCutoff = 1 -- see the matching comment in gsea_run_go.R;
# same reasoning applies here, database-independent.
gsea_results <- GSEA(
  geneList = ranking, TERM2GENE = term2gene, minGSSize = min_set_size,
  maxGSSize = max_set_size, pvalueCutoff = 1, by = "fgsea", seed = TRUE
)
cat(
  "\nGSEA complete:", nrow(as.data.frame(gsea_results)),
  "gene set(s) tested (significance is applied in postprocess,",
  "not here)\n"
)


##########
# Output #
##########

# Save gseaResult object for postprocessing.
dir.create(dirname(output_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(gsea_results, output_rds)
cat("Saved:", output_rds, "\n")


###########
# Cleanup #
###########

cat("\n=== GSEA run (HALLMARK) complete ===\n")
sink(type = "message")
sink()
close(logfile)
