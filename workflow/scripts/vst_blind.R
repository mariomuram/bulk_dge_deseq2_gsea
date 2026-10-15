#!/usr/bin/env Rscript

# SPDX-FileCopyrightText: © 2026 Marcel Schilling <mschilling@idibell.cat>
# SPDX-License-Identifier: AGPL-3.0-or-later


#########
# About #
#########

# File:     workflow/scripts/vst_blind.R
# Author:   Marcel Schilling <mschilling@idibell.cat>
# Created:  2026-10-10
# Modified: 2026-10-15
# License:  GNU Affero General Public License version 3 or later
# Purpose:  (Blind) VST R script for the bulk RNA-seq DGE (DESeq2 + GSEA)
#           Snakemake workflow.


###############
# Description #
###############

# Blind variance-stabilizing transform of the fitted dds, stored once as a
# plain assay matrix. Every consumer (plot_pca, batch_correct) reads this file
# instead of recomputing the transform, so they all start from identical values.


###########
# License #
###########

# (Blind) VST R script for bulk RNA-seq DGE (DESeq2 + GSEA) Snakemake
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

cat("=== VST (blind) ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


##############
# Parameters #
##############

dds <- snakemake@input[["dds_rds"]]
output_file <- snakemake@output[["vst"]]


#########
# Input #
#########

dds <- readRDS(dds)
cat("Loaded dds:", ncol(dds), "samples,", nrow(dds), "genes\n")


#######
# VST #
#######

# blind = TRUE: QC should not use the fitted design, so sample
# structure isn't partially "explained away" before we look at
# it. (blind = FALSE is appropriate later, for design-aware
# downstream use, not for QC.)
cat("\nRunning variance-stabilizing transform (blind = TRUE)...\n")
vst_data <- vst(dds, blind = TRUE)


##########
# Output #
##########

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
saveRDS(assay(vst_data), output_file)
cat("Saved:", output_file, "\n")


###########
# Cleanup #
###########

cat("\n=== VST (blind) Complete ===\n")
sink(type = "message")
sink()
close(logfile)
