#!/usr/bin/env Rscript

# SPDX-FileCopyrightText: © 2026 Marcel Schilling <mschilling@idibell.cat>
# SPDX-License-Identifier: AGPL-3.0-or-later


#########
# About #
#########

# File:     workflow/scripts/analysis_pca.R
# Author:   Marcel Schilling <mschilling@idibell.cat>
# Created:  2026-10-01
# Modified: 2026-10-29
# License:  GNU Affero General Public License version 3 or later
# Purpose:  PCA analysis R script for the bulk RNA-seq DGE (DESeq2 + GSEA)
#           Snakemake workflow.


###############
# Description #
###############

# QC: PCA across all samples (not contrast-specific), for one
# correction stage (`"uncorrected"` = blind VST; otherwise the
# batch-corrected VST). The `prcomp()` object is stashed so
# `pc_field_association` doesn't recompute it.


###########
# License #
###########

# PCA R script for bulk RNA-seq DGE (DESeq2 + GSEA) Snakemake
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

library(DESeq2)


###########
# Logging #
###########

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== QC: PCA Analysis ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


#############
# Constants #
#############

# Same default as `DESeq2::plotPCA()` (which this script used before):
# top 500 genes (see below).
ntop <- 500L


##############
# Parameters #
##############

# Extract parameterized variables
dds <- snakemake@input[["dds_rds"]]
vst_mat <- snakemake@input[["assay"]]
scores_file <- snakemake@output[["scores"]]


#########
# Input #
#########

dds <- readRDS(dds)
cat("Loaded dds:", ncol(dds), "samples,", nrow(dds), "genes\n")
vst_mat <- readRDS(vst_mat)
stopifnot(identical(colnames(vst_mat), colnames(dds)))


#############
# PCA genes #
#############

# Same computation as `DESeq2::plotPCA()` (which this script used before)
# with its defaults: top 500 genes (see above) by variance, centered,
# unscaled. Done by hand (see below) because `plotPCA()` returns only
# PC1/PC2 and not the `prcomp` object.
gene_var <- apply(vst_mat, 1, var)
top_genes <- order(gene_var, decreasing = TRUE)[
  seq_len(min(ntop, length(gene_var)))
]


#######
# PCA #
#######

pca <- prcomp(t(vst_mat[top_genes, ]))


##########
# Output #
##########

saveRDS(pca, scores_file)
cat("Saved:", scores_file, "\n")


###########
# Cleanup #
###########

cat("\n=== QC: PCA Analysis Complete ===\n")
sink(type = "message")
sink()
close(logfile)
