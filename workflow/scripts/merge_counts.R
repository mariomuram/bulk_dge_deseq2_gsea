#!/usr/bin/env Rscript

# SPDX-FileCopyrightText: © 2026 Marcel Schilling <mschilling@idibell.cat>
# SPDX-License-Identifier: AGPL-3.0-or-later


#########
# About #
#########

# File:     workflow/scripts/merge_counts.R
# Author:   Marcel Schilling <mschilling@idibell.cat>
# Created:  2026-10-03
# Modified: 2026-10-15
# License:  GNU Affero General Public License version 3 or later
# Purpose:  Count matrix merging R script for the bulk RNA-seq DGE (DESeq2
#           + GSEA) Snakemake workflow.


###############
# Description #
###############

# Merge counts: read the count matrices of all configured files into one
# gene x sample matrix, build the sample sheet (coldata) from the labels
# parsed in the Snakefile, and drop low-count genes.
#
# Everything here is independent of how confounders are handled downstream
# (modelled in the DESeq2 design, or removed from the counts beforehand), so
# every analysis path can start from this one matrix and sample sheet.
#
# Outputs are plain gzipped TSV, not RDS (same convention as the pipeline's
# input and results files):
#   * counts:  first column `gene_id`, then one column per sample (integers),
#              already pre-filtered (genes with rowSums < 10 are dropped);
#   * coldata: first column `sample`, then `condition` and one column per
#              parsed field, all plain strings. No reference levels are
#              applied here; deseq2_fit.R relevels from the Snakefile-supplied
#              `reference_levels` param, never from the stored file.


###########
# License #
###########

# Count matrix merging R script for bulk RNA-seq DGE (DESeq2 + GSEA)
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

library(readr)


###########
# Logging #
###########

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== Merge Counts ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


##############
# Parameters #
##############

count_files <- snakemake@input[["count_files"]]  # vector of file paths
file_to_samples <- snakemake@params[["file_to_samples"]]  # file path -> samples mapping
samples_info <- snakemake@params[["samples_info"]]    # sample -> label
sample_fields <- snakemake@params[["sample_fields"]]   # sample -> list
condition_field_names <- unlist(snakemake@params[["condition_fields"]])
confounder_field_names <- unlist(snakemake@params[["confounder_fields"]])
counts_file <- snakemake@output[["counts"]]
coldata_file <- snakemake@output[["coldata"]]


#########
# Input #
#########

# Read and merge count matrices (from multiple files).


if (length(count_files) == 1L) {
  cat("Reading count matrix...\n")
} else {
  cat("Reading count matrices from multiple files...\n")
  cat("Expected files:", paste(count_files, collapse = ", "), "\n")
}

# Extract all samples from all files
all_samples <- do.call(c, lapply(file_to_samples,
                                  function(file) { names(file$samples) }))

# Check for duplicates within each file
for (i in seq_along(file_to_samples)) {
  file_samples <- names(file_to_samples[[i]]$samples)
  duplicated_in_file <- file_samples[duplicated(file_samples)]
  if (length(duplicated_in_file) > 0) {
    stop(
      "ERROR: the following sample(s) appear multiple times in file ",
      count_files[i], ": ",
      paste(duplicated_in_file, collapse = ", ")
    )
  }
}

# Check for duplicates across files
sample_counts <- table(all_samples)
duplicated_samples <- names(sample_counts)[sample_counts > 1]
if (length(duplicated_samples) > 0) {
  stop(
    "ERROR: the following sample(s) appear in more than one file: ",
    paste(duplicated_samples, collapse = ", ")
  )
}

cat("Total samples to load:", length(sample_counts), "\n\n")

# Read all files and merge on gene_id
file_to_counts <- list()
seen_genes <- NULL

for (file_path in count_files) {
  cat("Reading:", file_path, "\n")

  file_idx <- which(sapply(file_to_samples,
                           function(file) { if ("path" %in% names(file))
                                              return(file$path)
                                            return(NULL) })
                      == file_path)
  if (length(file_idx) == 0L) {
    stop(
      "ERROR: the following file path is not listed in the file-to-samples ",
      "mapping: ", file_path
    )
  } else {
    if (length(file_idx) > 1L) {
      stop(
        "ERROR: the following file path is listed several times in the",
        "file-to-samples mapping:",
        paste(names(file_to_samples)[file_idx], collapse = ", ")
      )
    }
  }
  stopifnot(length(file_idx) == 1L)
  file_key <- names(file_to_samples)[[file_idx]]
  file_samples <- names(file_to_samples[[file_idx]]$samples)
  n_annotation_cols <- file_to_samples[[file_idx]]$annotation_columns

  df <- read_tsv(file_path, show_col_types = FALSE)

  # First column is gene_id; the next `n_annotation_cols` columns are
  # unused annotation (e.g. RSEM's transcript_id(s)); the rest are samples.
  gene_col <- df[[1]]
  sample_cols <- colnames(df)[-seq_len(1 + n_annotation_cols)]
  file_idx <- which(sapply(file_to_samples,
                           function(file) { if ("path" %in% names(file))
                                              return(file$path)
                                            return(NULL) })
                      == file_path)
  if (length(file_idx) == 0L) {
    stop(
      "ERROR: the following file path is not listed in the file-to-samples ",
      "mapping: ", file_path
    )
  } else {
    if (length(file_idx) > 1L) {
      stop(
        "ERROR: the following file path is listed several times in the",
        "file-to-samples mapping:",
        paste(names(file_to_samples)[file_idx], collapse = ", ")
      )
    }
  }
  stopifnot(length(file_idx) == 1L)
  file_key <- names(file_to_samples)[[file_idx]]
  file_samples <- names(file_to_samples[[file_idx]]$samples)
  if (length(file_samples) == 0L) {
    stop(
      "ERROR: the following file path is not associated to any sample(s): ",
      file_path
    )
  }
  missing_samples <- setdiff(file_samples, sample_cols)
  if (length(missing_samples) > 0L) {
    stop(
      "ERROR: the following sample(s) is/are missing from `", file_path, "`:",
      paste(missing_samples, collapse = ", ")
    )
  }
  extra_samples <- setdiff(sample_cols, file_samples)
  if (length(extra_samples) > 0L) {
    cat(
      "WARNING: `", file_path, "` contains the following extra sample(s):\n  ",
      paste(extra_samples, collapse = ", "), "\n...dropping\n", sep = ""
    )
  }

  cat("  Genes:", nrow(df), ", Samples:", length(sample_cols), "(kept:",
      length(file_samples), ")\n")

  # Check gene consistency
  if (is.null(seen_genes)) {
    seen_genes <- gene_col
  } else {
    if (!identical(gene_col, seen_genes)) {
      stop(
        "ERROR: gene set in ", file_path, " does not match ",
        "the first file. Genes must be identical across all files ",
        "(same set, same order). Use an explicit gene_id join if ",
        "gene sets differ."
      )
    }
  }

  file_to_samples[[file_key]] <- file_samples
  file_to_counts[[file_key]] <- df[, file_samples, drop = FALSE]
}


###########
# Merging #
###########

gene_ids <- seen_genes

# Merge all count matrices side-by-side
# (same gene order across all files means cbind is safe)
counts_matrix <- do.call(cbind, unname(file_to_counts))
counts_matrix <- as.matrix(counts_matrix)
rownames(counts_matrix) <- gene_ids

# Double-check that columns match file_to_samples keys (parsed sample IDs)
expected_samples <- names(sample_counts)
actual_samples <- colnames(counts_matrix)
stopifnot(setequal(expected_samples, actual_samples))

cat("\nMerged dimensions:", nrow(counts_matrix), "genes x",
    ncol(counts_matrix), "samples (total across all files)\n")



############
# Metadata #
############

# Build sample metadata (fields parsed in the Snakefile).

missing_samples <- setdiff(names(samples_info), colnames(counts_matrix))
if (length(missing_samples) > 0) {
  stop(
    "ERROR: sample(s) listed in config.yaml 'samples' are not ",
    "present in the count matrix: ",
    paste(missing_samples, collapse = ", ")
  )
}

cat("\nUsing", length(samples_info), "of",
    ncol(counts_matrix), "samples found in count matrix.\n")

counts_matrix <- counts_matrix[, names(samples_info), drop = FALSE]

field_names <- names(sample_fields[[1]])
cat("Parsed sample fields:", paste(field_names, collapse = ", "), "\n")

cat("Condition fields (crossed):", paste(condition_field_names,
                                         collapse = ", "), "\n")
if (length(confounder_field_names) > 0) {
  cat("Confounder fields (additive only):",
      paste(confounder_field_names, collapse = ", "), "\n")
} else {
  cat("No confounder fields configured.\n")
}
stopifnot(all(union(condition_field_names, confounder_field_names)
              %in% field_names))

coldata <- data.frame(
  sample = colnames(counts_matrix),
  condition = unlist(samples_info)[colnames(counts_matrix)],
  stringsAsFactors = FALSE
)
for (f in field_names) {
  coldata[[f]] <- vapply(
    sample_fields[colnames(counts_matrix)],
    function(x) x[[f]], character(1)
  )
}
rownames(coldata) <- coldata$sample

cat("\nSample sheet:\n")
print(coldata)

field_has_na <- vapply(coldata[field_names], anyNA, logical(1))
if (anyNA(coldata$condition) || any(field_has_na)) {
  stop(
    "ERROR: could not map all samples to their condition/",
    "fields. Check config.yaml 'samples' and ",
    "'sample_label_pattern'."
  )
}

cat("\nReplicates per condition:\n")
print(table(coldata$condition))

n_per_condition <- table(coldata$condition)
if (any(n_per_condition < 2)) {
  cat(
    "WARNING: at least one condition has fewer than 2",
    "replicates. DESeq2 requires n >= 2 per group; results",
    "for affected contrasts will likely fail or be unreliable.\n"
  )
} else if (any(n_per_condition == 2)) {
  cat(
    "NOTE: minimum replicate number is 2. Dispersion",
    "estimates from n = 2 are noisy; treat results as",
    "exploratory and interpret with caution.\n"
  )
}


###################
# Counts rounding #
###################

# DESeq2 needs integer counts, so fractional (e.g. RSEM expected) counts are
# rounded; the filter then runs on the rounded values, exactly as it did when
# it lived in deseq2_fit.R.
if (anyNA(counts_matrix)) {
  stop("ERROR: the merged count matrix contains NA value(s).")
}
counts_matrix <- round(counts_matrix)


#################
# Pre-filtering #
#################

# Remove genes with very low counts.
cat("\nPre-filtering low-count genes (rowSums >= 10)...\n")
keep <- rowSums(counts_matrix) >= 10
cat("Genes removed:", sum(!keep), "\n")
counts_matrix <- counts_matrix[keep, , drop = FALSE]
cat("Genes retained after filtering:", nrow(counts_matrix), "\n")
if (nrow(counts_matrix) == 0L) {
  stop("ERROR: no gene passed the low-count pre-filter.")
}


##########
# Output #
##########

# Stored as integers: readr would otherwise write large doubles in scientific
# notation (e.g. 1e+05).
if (max(counts_matrix) > .Machine$integer.max) {
  stop("ERROR: a count exceeds the integer range; cannot store as integer.")
}
storage.mode(counts_matrix) <- "integer"


# write_tsv() gzip-compresses on its own, going by the .gz file extension.
dir.create(dirname(counts_file), showWarnings = FALSE, recursive = TRUE)
write_tsv(
  data.frame(gene_id = rownames(counts_matrix), counts_matrix,
             check.names = FALSE, row.names = NULL,
             stringsAsFactors = FALSE),
  counts_file
)
cat("\nSaved counts:", counts_file, "\n")

dir.create(dirname(coldata_file), showWarnings = FALSE, recursive = TRUE)
write_tsv(coldata, coldata_file)
cat("Saved coldata:", coldata_file, "\n")


###########
# Cleanup #
###########

cat("\n=== Merge Counts Complete ===\n")
sink(type = "message")
sink()
close(logfile)
