#!/usr/bin/env Rscript

# SPDX-FileCopyrightText: © 2026 Marcel Schilling <mschilling@idibell.cat>
# SPDX-License-Identifier: AGPL-3.0-or-later


#########
# About #
#########

# File:     workflow/scripts/gsea_run_go.R
# Author:   Marcel Schilling <mschilling@idibell.cat>
# Created:  2026-09-01
# Modified: 2026-10-15
# License:  GNU Affero General Public License version 3 or later
# Purpose:  GO GSEA run R script for the bulk RNA-seq DGE (DESeq2 + GSEA)
#           Snakemake workflow.


###############
# Description #
###############

# GSEA (preranked) using clusterProfiler::gseGO() against the
# configured GO ontology, backed by fgsea. Ranking statistic comes
# directly from the DESeq2 results table (no reconstruction). This
# script only runs the expensive enrichment step -- it saves the
# live gseaResult object to RDS for scripts/gsea_postprocess.R to
# pick up (simplify(), tables, plots all happen there instead, so
# tweaking a plot never reruns fgsea).


###########
# License #
###########

# GO GSEA run R script for bulk RNA-seq DGE (DESeq2 + GSEA) Snakemake
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
library(org.Hs.eg.db)
library(GOSemSim)


###########
# Logging #
###########

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== GSEA run (GO) ===\n")
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
ontology <- snakemake@params[["ontology"]]
gene_id_style <- snakemake@params[["gene_id_style"]]
ranking_stat <- snakemake@params[["ranking_stat"]]
min_set_size <- snakemake@params[["min_set_size"]]
max_set_size <- snakemake@params[["max_set_size"]]
deseq2_results <- snakemake@input[["results"]]
output_rds <- snakemake@output[["gsea_results_rds"]]

cat("Contrast:", contrast_name, "\n")
cat("Analysis name:", gsea_name, "\n")
cat("Database:", database, "\n")
cat("Ontology:", ontology, "\n")
cat("Gene ID style:", gene_id_style, "\n")
cat("Ranking stat:", ranking_stat, "\n")
cat("Min set size:", min_set_size, "\n")
cat("Max set size:", max_set_size, "\n\n")

if (database != "GO") {
  stop(
    "ERROR: gsea_run_go.R was invoked with database = '",
    database, "'; this script only handles GO. Check ",
    "the Snakefile's wildcard_constraints for gsea_run_go."
  )
}


#########
# Input #
#########

# Load DESeq2 results 
results_df <- read_tsv(deseq2_results, show_col_types = FALSE)


################
# Gene ranking #
################

# Build ranking (shared helper)
ranking <- validate_and_rank(results_df, ranking_stat, gene_id_style)


############
# Run GSEA #
############

cat("\nRunning gseGO (ontology:", ontology, ")...\n")
# pvalueCutoff = 1 -- fgsea computes NES/p-value for every gene set
# within the size bounds regardless of this value; it's a post-hoc
# filter on what's returned, not a computational shortcut. Keeping
# everything here (rather than padj_threshold) lets postprocess apply
# the significance call explicitly, and lets scripts/gsea_summary.R
# report a gene set's real NES/padj in every contrast it was tested
# in, not just where it happened to clear threshold. Already-computed
# p.adjust values for anything significant are unaffected -- BH
# correction runs over the full tested universe either way.
gsea_results <- gseGO(
  geneList = ranking,
  OrgDb = org.Hs.eg.db,
  keyType = gene_id_style,
  ont = ontology,
  minGSSize = min_set_size,
  maxGSSize = max_set_size,
  pvalueCutoff = 1,
  by = "fgsea",
  seed = TRUE
)
cat(
  "\ngseGO complete:", nrow(as.data.frame(gsea_results)),
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

cat("\n=== GSEA run (GO) complete ===\n")
sink(type = "message")
sink()
close(logfile)
