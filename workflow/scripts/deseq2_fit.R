#!/usr/bin/env Rscript

# SPDX-FileCopyrightText: © 2026 Marcel Schilling <mschilling@idibell.cat>
# SPDX-License-Identifier: AGPL-3.0-or-later


#########
# About #
#########

# File:     workflow/scripts/deseq2_fit.R
# Author:   Marcel Schilling <mschilling@idibell.cat>
# Created:  2026-08-24
# Modified: 2026-10-15
# License:  GNU Affero General Public License version 3 or later
# Purpose:  DESeq2 model fitting R script for the bulk RNA-seq DGE
#           (DESeq2 + GSEA) Snakemake workflow.


###############
# Description #
###############

# DESeq2 fit: read the merged, pre-filtered count matrix and sample sheet
# (see merge_counts.R), run normalization, dispersion estimation and model
# fitting.
#
# The design is the full crossing of every parsed sample field
# (all main effects and interactions), built and fit as an
# explicit model matrix rather than a formula string. This is
# necessary, not just tidier: if the fields don't form a
# complete factorial (e.g. a background x dox-status cell with
# zero samples, as will happen once 'clean' samples with no
# GATA2R396Q counterpart are added), the naive formula produces
# an all-zero model-matrix column and DESeq2 refuses to fit
# (checkFullRank error). Building the matrix directly lets us
# detect and drop exactly those structurally-empty columns --
# something a formula string cannot express -- following
# DESeq2's own documented workaround for this situation.
#
# The cleaned matrix is stashed in metadata(dds) so
# deseq2_contrast.R can reuse the exact matrix that was fit,
# rather than recomputing (and risking drift from) it.


###########
# License #
###########

# DESeq2 model fitting R script for bulk RNA-seq DGE (DESeq2 + GSEA)
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
library(readr)


###########
# Logging #
###########

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== DESeq2 Fit ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


###############
# Shared code #
###############

source(snakemake@input[["confounders_lib"]])


##############
# Parameters #
##############

counts_file <- snakemake@input[["counts"]]
coldata_file <- snakemake@input[["coldata"]]
reference_levels <- snakemake@params[["reference_levels"]]
condition_field_names <- unlist(snakemake@params[["condition_fields"]])
confounder_field_names <- unlist(snakemake@params[["confounder_fields"]])
confounder_field_types <- snakemake@params[["confounder_field_types"]]


#########
# Input #
#########

# Read the merged, pre-filtered counts and the sample sheet (both written
# by `merge_counts.R` - see there for the file layout).
cat("Reading merged counts:", counts_file, "\n")
counts_df <- read_tsv(
  counts_file,
  col_types = cols(gene_id = col_character(), .default = col_double())
)
counts_matrix <- as.matrix(counts_df[, -1, drop = FALSE])
rownames(counts_matrix) <- counts_df$gene_id
if (any(counts_matrix != round(counts_matrix))) {
  stop("ERROR: ", counts_file, " contains non-integer counts.")
}
storage.mode(counts_matrix) <- "integer"
cat("Counts:", nrow(counts_matrix), "genes x", ncol(counts_matrix),
    "samples\n")

# Every column is read back as plain strings (na = character(): a field value
# that happens to read "NA" is a value, not a missing one); the factor levels
# (and thereby the reference levels) are set below from the `reference_levels`
# param.
cat("Reading sample sheet:", coldata_file, "\n")
coldata <- as.data.frame(
  read_tsv(coldata_file, col_types = cols(.default = col_character()),
           na = character(), trim_ws = FALSE)
)
rownames(coldata) <- coldata$sample
stopifnot(identical(colnames(counts_matrix), coldata$sample))

# Snakefile: REFERENCE_LEVELS has one entry per parsed sample field.
field_names <- names(reference_levels)
stopifnot(all(field_names %in% colnames(coldata)))
cat("Parsed sample fields:", paste(field_names, collapse = ", "), "\n")

cat("Condition fields (crossed):", paste(condition_field_names,
                                         collapse = ", "), "\n")
confounder_mat <- NULL
if (length(confounder_field_names) > 0) {
  cat("Confounder fields (additive only):",
      paste(confounder_field_names, collapse = ", "), "\n")
} else {
  # Either none are configured, or ComBat-Seq has removed all of them from
  # the counts beforehand.
  cat("No confounder fields in the design.\n")
}
stopifnot(all(union(condition_field_names, confounder_field_names)
              %in% field_names))

cat("\nSample sheet:\n")
print(coldata)


####################
# Reference levels #
####################

# Affects resultsNames()/coefficient signs and which
# combinations end up as the "fully reference" corner of the
# full-crossing matrix below; does not affect fitted group means
# or dispersions.

cat("\nReference levels:\n")
for (f in names(reference_levels)) {
  ref <- reference_levels[[f]]
  cat(" -", f, ":", ref, "\n")
  coldata[[f]] <- relevel(factor(coldata[[f]]), ref = ref)
}


##########
# Design #
##########

# Build the full-crossing design matrix, dropping any structurally-empty
# (all-zero) columns.

condition_formula <- as.formula(
  paste("~", paste(condition_field_names, collapse = "*"))
)
cat("\nCondition-crossing formula (before rank check):",
    deparse(condition_formula), "\n")

condition_mat <- model.matrix(condition_formula, data = coldata)
empty_cols <- colSums(condition_mat != 0) == 0

if (any(empty_cols)) {
  cat(
    "\nNOTE: the parsed condition fields do not form a complete",
    "factorial -- the following coefficient(s) have no",
    "observations and are being dropped before fitting (expected,",
    "not an error, e.g. a construct-free condition existing for",
    "only one background):\n"
  )
  print(colnames(condition_mat)[empty_cols])
  condition_mat <- condition_mat[, !empty_cols, drop = FALSE]
}

n_distinct_conditions <- length(unique(coldata$condition))
if (ncol(condition_mat) != n_distinct_conditions) {
  stop(
    "ERROR: cleaned condition-crossing matrix has ", ncol(condition_mat),
    " columns but there are ", n_distinct_conditions,
    " distinct observed conditions. The condition part of the model ",
    "is not saturated relative to the observed groups -- this would ",
    "break the relevel-based contrast extraction used elsewhere in ",
    "this pipeline. This should not happen for a matrix built purely ",
    "from crossing the parsed condition fields; investigate before ",
    "trusting downstream results."
  )
}

# Confounders enter additively only -- never crossed with condition
# or with each other -- so they shift the fitted baseline without
# changing which linear combination of condition columns identifies
# which comparison. This is what keeps every contrast apeglm-
# shrinkable exactly as before confounders were introduced.
if (length(confounder_field_names) > 0) {
  confounder_formula <- as.formula(
    paste("~", paste(confounder_field_names, collapse = " + "))
  )
  cat("\nConfounder formula (additive, before rank check):",
      deparse(confounder_formula), "\n")
  # Every field was turned into a factor by relevel() above; "continuous"
  # confounders must be numeric again (see lib/confounders.R).
  confounder_coldata <- cast_continuous_confounders(
    coldata, confounder_field_names, confounder_field_types
  )
  confounder_mat <- model.matrix(confounder_formula, data = confounder_coldata)
  confounder_mat <- confounder_mat[
    , colnames(confounder_mat) != "(Intercept)", drop = FALSE
  ]
  empty_confounder_cols <- colSums(confounder_mat != 0) == 0
  if (any(empty_confounder_cols)) {
    # Not expected: every confounder level was observed in at least
    # one sample by construction (it came from parsing an actual
    # sample label). An all-zero column here means a coding bug,
    # not a legitimate missing combination.
    stop(
      "ERROR: confounder design matrix has unexpected all-zero ",
      "column(s): ",
      paste(colnames(confounder_mat)[empty_confounder_cols], collapse = ", "),
      ". Investigate before trusting downstream results."
    )
  }
  full_mat <- cbind(condition_mat, confounder_mat)   # condition (incl. Intercept) first
} else {
  full_mat <- condition_mat
}
stopifnot(colnames(full_mat)[1] == "(Intercept)")

cat("\nFinal design matrix columns (", ncol(full_mat), "):\n", sep = "")
print(colnames(full_mat))


#################
# DESeq dataset #
#################

# `design` is a placeholder (`~condition`, always full rank by
# construction) at object-creation time: `DESeqDataSetFromMatrix`
# would itself reject the full, uncleaned crossing formula via
# its own rank check. The actual (cleaned) design matrix is
# supplied to `DESeq()` below via `full =`.

dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = coldata,
  design = ~condition
)

# Low-count genes were already removed by merge_counts.R.
cat("DESeqDataSet created with", nrow(dds), "genes and",
    ncol(dds), "samples\n")


##############
# Run DESeq2 #
##############

cat("\nRunning DESeq2 (normalize, estimate dispersions, fit)",
    "with the cleaned full-crossing design matrix...\n")
dds <- tryCatch(
  DESeq(dds, full = full_mat, betaPrior = FALSE),
  error = function(e) {
    stop("ERROR: DESeq() failed to fit the model: ", e$message)
  }
)

cat("Normalization, dispersion estimation and model fitting",
    "complete.\n")
cat("\nModel coefficients (resultsNames):\n")
print(resultsNames(dds))

if ("betaConv" %in% colnames(mcols(dds))) {
  n_not_converged <- sum(!mcols(dds)$betaConv, na.rm = TRUE)
  if (n_not_converged > 0) {
    cat(
      "WARNING:", n_not_converged,
      "gene(s) did not converge during beta (coefficient)",
      "fitting. Their results should be treated with caution;",
      "see mcols(dds)$betaConv for the affected genes.\n"
    )
  }
}

if ("dispOutlier" %in% colnames(mcols(dds))) {
  n_disp_outliers <- sum(mcols(dds)$dispOutlier, na.rm = TRUE)
  if (n_disp_outliers > 0) {
    cat(
      "NOTE:", n_disp_outliers,
      "gene(s) had dispersion estimates flagged as outliers",
      "relative to the fitted mean-dispersion trend (not",
      "shrunk). Expected for some genes; not necessarily a",
      "problem, but worth a look for genes of interest.\n"
    )
  }
}


########################
# Column names syncing #
########################

# Sync column names to `resultsNames(dds)` before stashing.

# DESeq2 relabels raw model.matrix() coefficient names for its
# own resultsNames()/results(name=) interface -- notably ":" in
# interaction terms becomes "." (e.g. "backgroundGATA2R396Q:doxdox"
# -> "backgroundGATA2R396Q.doxdox"). This is a deliberate DESeq2
# display choice, not an R syntax requirement, and it means the
# raw model.matrix() colnames on full_mat do NOT match
# resultsNames(dds) for the identical coefficients. The pipeline
# itself resolves coefficients by position, not by name, so this
# mismatch was never a functional bug -- but syncing here means
# anyone loading this RDS for a manual, ad hoc check sees the
# same names in metadata(dds)$full_model_matrix as resultsNames()
# reports, and isn't caught out by it the way I was when I first
# encountered this.

if (length(resultsNames(dds)) != ncol(full_mat)) {
  stop(
    "ERROR: resultsNames(dds) has ", length(resultsNames(dds)),
    " entries but full_mat has ", ncol(full_mat), " columns; ",
    "cannot safely synchronize column names. Investigate before ",
    "trusting the stashed model matrix."
  )
}
colnames(full_mat) <- resultsNames(dds)

cat("\nStashed model matrix columns (synced to resultsNames):\n")
print(colnames(full_mat))


#################################
# Matrix & field names stashing #
#################################

metadata(dds)$full_model_matrix <- full_mat
metadata(dds)$condition_model_matrix <- condition_mat
# NULL (no key) when there are no confounders
metadata(dds)$confounder_model_matrix <- confounder_mat
metadata(dds)$condition_field_names <- condition_field_names
metadata(dds)$confounder_field_names <- confounder_field_names


##########
# Output #
##########

output_dds <- snakemake@output[["dds_rds"]]
dir.create(dirname(output_dds), showWarnings = FALSE, recursive = TRUE)
saveRDS(dds, output_dds)
cat("\nSaved DESeqDataSet:", output_dds, "\n")


###########
# Cleanup #
###########

cat("\n=== DESeq2 Fit Complete ===\n")
sink(type = "message")
sink()
close(logfile)
