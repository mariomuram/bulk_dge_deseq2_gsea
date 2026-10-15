#!/usr/bin/env Rscript

# SPDX-FileCopyrightText: © 2026 Marcel Schilling <mschilling@idibell.cat>
# SPDX-License-Identifier: AGPL-3.0-or-later


#########
# About #
#########

# File:     workflow/scripts/deseq2_relevel.R
# Author:   Marcel Schilling <mschilling@idibell.cat>
# Created:  2026-08-27
# Modified: 2026-10-15
# License:  GNU Affero General Public License version 3 or later
# Purpose:  DESeq2 column releveling R script for the bulk RNA-seq DGE
#           (DESeq2 + GSEA) Snakemake workflow.


###############
# Description #
###############

# Relevel + Wald retest: take the primary (interaction-design)
# DESeqDataSet, relevel `condition` to a new reference, and
# rerun only nbinomWaldTest(). Dispersions and size factors are
# reused as-is: the primary design is saturated (as many free
# parameters as distinct condition groups), so reparametrizing
# to ~condition changes which coefficient represents which
# comparison, but not the fitted group means or dispersions
# themselves. This makes every pairwise contrast against this
# reference extractable -- and apeglm-shrinkable -- as an
# ordinary named coefficient.


###########
# License #
###########

# DESeq2 column releveling R script for bulk RNA-seq DGE (DESeq2 + GSEA)
# Snakemake workflow.
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

cat("=== DESeq2 Relevel + Wald Retest ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


###############
# Shared code #
###############

source(snakemake@input[["confounders_lib"]])


##############
# Parameters #
##############

dds_file <- snakemake@input[["dds_rds"]]
reference_factor <- snakemake@params[["reference_factor"]]
confounder_field_names <- unlist(snakemake@params[["confounder_fields"]])
confounder_field_types <- snakemake@params[["confounder_field_types"]]


#########
# Input #
#########

cat("Loading primary DESeqDataSet:", dds_file, "\n")
dds <- readRDS(dds_file)

available_conditions <- levels(factor(colData(dds)$condition))
if (!(reference_factor %in% available_conditions)) {
  stop(
    "ERROR: reference_factor '", reference_factor, "' is not an ",
    "available condition. Available: ",
    paste(available_conditions, collapse = ", ")
  )
}


###############
# Relevelling #
###############

cat("Relevelling 'condition' to reference:", reference_factor, "\n")
colData(dds)$condition <- relevel(
  factor(colData(dds)$condition), ref = reference_factor
)

# The saved colData holds every field as a factor (see deseq2_fit.R); a
# "continuous" confounder has to be numeric again before it enters the
# design, exactly as in the primary fit (see lib/confounders.R).
colData(dds) <- cast_continuous_confounders(
  colData(dds), confounder_field_names, confounder_field_types
)

relevel_formula <- as.formula(
  paste("~", paste(c(confounder_field_names, "condition"), collapse = " + "))
)
cat("Releveled design formula:", deparse(relevel_formula), "\n")
design(dds) <- relevel_formula


#############
# Refitting #
#############

cat("\nRerunning nbinomWaldTest() only (dispersions and size",
    "factors reused from the primary fit)...\n")
dds <- tryCatch(
  nbinomWaldTest(dds),
  error = function(e) {
    stop("ERROR: nbinomWaldTest() failed: ", e$message)
  }
)

cat("\nModel coefficients (resultsNames) after relevel:\n")
print(resultsNames(dds))


##########
# Output #
##########

output_dds <- snakemake@output[["dds_rds"]]
dir.create(dirname(output_dds), showWarnings = FALSE, recursive = TRUE)
saveRDS(dds, output_dds)
cat("\nSaved:", output_dds, "\n")


###########
# Cleanup #
###########

cat("\n=== DESeq2 Relevel + Wald Retest Complete ===\n")
sink(type = "message")
sink()
close(logfile)
