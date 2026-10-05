#!/usr/bin/env Rscript
# ComBat-Seq correction of the merged count matrix: removes the effect of the
# configured confounders from the counts themselves, BEFORE the DESeq2 fit.
# (The alternative -- modelling them in the DESeq2 design -- is the
# `confounders_modelled` analysis; see the Snakefile's `{confounder_handling}`.)
#
# The biological signal is protected in every call: `group` is always
# coldata$condition, whichever strategy is active.
#
# Constraint, all strategies: ComBat_seq()'s `batch` is strictly discrete. A
# confounder declared "continuous" can therefore never be the batch; it can
# only ever be preserved (as a covariate), never corrected for. Whatever
# ComBat-Seq does not remove from the counts stays a term of the DESeq2 design
# of the corrected fit (see COMBAT_SEQ_RESIDUAL_CONFOUNDERS in the Snakefile).
#
# `multi_confounder_strategy` decides what happens with several confounders
# (with a single discrete one and nothing else, all three are the same plain
# ComBat_seq() call):
#   sequential (default)
#     One ComBat_seq() call per discrete confounder, in declaration order. At
#     each step the not-yet-processed discrete confounders and all continuous
#     ones are covariates, so they are not stripped out as noise before their
#     own turn (discrete) or at all (continuous). Order-dependent, best-effort
#     approximation, not a validated joint method.
#   preserve_extra
#     One call; batch = the first discrete confounder; every other confounder
#     (discrete or continuous) is a covariate, i.e. preserved and NOT
#     corrected for. The log says which.
#   composite
#     One call; batch = the interaction of all discrete confounders' levels.
#     Requires all confounders to be discrete. Interaction cells can get very
#     small (n ~ 10-35 samples in total), so their sizes are logged.
#
# Output: same genes and samples, in the same order, as the input (later
# steps compare column names with the DESeq2 objects), as integers.

library(readr)
library(sva)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== ComBat-Seq Correction ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

source(snakemake@input[["confounders_lib"]])

# Batches (and composite cells) with fewer samples than this are flagged.
MIN_BATCH_SIZE <- 3L

# =====================================================================
# Functions
# =====================================================================

# Integer count matrix (genes x samples) from a TSV as written by
# merge_counts.R: first column `gene_id`, then one column per sample.
read_count_matrix <- function(path) {
  df <- read_tsv(
    path,
    col_types = cols(gene_id = col_character(), .default = col_double())
  )
  counts <- as.matrix(df[, -1, drop = FALSE])
  rownames(counts) <- df$gene_id
  if (any(counts != round(counts))) {
    stop("ERROR: ", path, " contains non-integer counts.")
  }
  storage.mode(counts) <- "integer"
  counts
}

# Model matrix of `fields` without the intercept (ComBat_seq() wants none),
# or NULL if there are no fields. `typed_coldata` must already have its
# continuous fields cast to numeric (cast_continuous_confounders()).
build_covar_mod <- function(typed_coldata, fields) {
  if (length(fields) == 0) return(NULL)
  covar_formula <- as.formula(
    paste("~", paste0("`", fields, "`", collapse = " + "))
  )
  covar_mod <- model.matrix(covar_formula, data = typed_coldata)
  covar_mod[, colnames(covar_mod) != "(Intercept)", drop = FALSE]
}

# The thing that would silently corrupt every downstream DESeq2 fit.
validate_corrected_counts <- function(corrected, original) {
  if (!identical(dim(corrected), dim(original))) {
    stop(
      "ERROR: ComBat-Seq output has dimensions ",
      paste(dim(corrected), collapse = " x "), " but its input has ",
      paste(dim(original), collapse = " x "), "."
    )
  }
  if (!identical(dimnames(corrected), dimnames(original))) {
    stop("ERROR: ComBat-Seq output has different gene or sample names (or ",
         "a different order) than its input.")
  }
  if (anyNA(corrected)) {
    stop("ERROR: ComBat-Seq output contains NA value(s).")
  }
  if (any(corrected < 0)) {
    stop("ERROR: ComBat-Seq output contains negative value(s).")
  }
  if (any(corrected != round(corrected))) {
    stop("ERROR: ComBat-Seq output contains non-integer value(s).")
  }
  invisible(TRUE)
}

# One validated ComBat_seq() call. `label` says what the call is for.
run_combat_seq <- function(counts, batch, group, covar_mod, label) {
  cat("\n--- ComBat-Seq:", label, "---\n")
  batch <- factor(batch)
  if (nlevels(batch) < 2) {
    stop("ERROR: the batch has a single level (", levels(batch), "); ",
         "ComBat-Seq needs at least two (", label, ").")
  }
  batch_sizes <- table(batch)
  cat("Batch sizes:\n")
  print(batch_sizes)
  if (any(batch_sizes < MIN_BATCH_SIZE)) {
    cat(
      "WARNING: batch(es) with fewer than", MIN_BATCH_SIZE, "samples:",
      paste(names(batch_sizes)[batch_sizes < MIN_BATCH_SIZE], collapse = ", "),
      "-- the correction for them is poorly determined.\n"
    )
  }

  corrected <- tryCatch(
    ComBat_seq(counts, batch = batch, group = group, covar_mod = covar_mod),
    error = function(e) {
      stop("ERROR: ComBat_seq() failed (", label, "): ", conditionMessage(e))
    }
  )
  if (is.null(dimnames(corrected))) dimnames(corrected) <- dimnames(counts)
  validate_corrected_counts(corrected, counts)
  cat("Entries changed by this step:", sum(corrected != counts), "of",
      length(counts), "\n")
  corrected
}

correct_sequential <- function(counts, typed_coldata, group, discrete_fields,
                               continuous_fields) {
  if (length(discrete_fields) > 1) {
    cat(
      "WARNING: strategy 'sequential' with several discrete confounders is an",
      "order-dependent, best-effort approximation (declaration order:",
      paste(discrete_fields, collapse = " -> "),
      "), not a validated joint correction.\n"
    )
  }
  if (length(continuous_fields) > 0) {
    cat(
      "NOTE: continuous confounder(s) preserved, NOT corrected for:",
      paste(continuous_fields, collapse = ", "), "\n"
    )
  }
  for (i in seq_along(discrete_fields)) {
    batch_field <- discrete_fields[i]
    pending <- discrete_fields[-seq_len(i)]
    protected <- c(pending, continuous_fields)
    counts <- run_combat_seq(
      counts, batch = typed_coldata[[batch_field]], group = group,
      covar_mod = build_covar_mod(typed_coldata, protected),
      label = sprintf(
        "step %d of %d, batch = %s, protected: %s", i, length(discrete_fields),
        batch_field,
        if (length(protected) > 0) paste(protected, collapse = ", ")
        else "condition only"
      )
    )
  }
  counts
}

correct_preserve_extra <- function(counts, typed_coldata, group,
                                   confounder_fields, discrete_fields) {
  batch_field <- discrete_fields[1]
  preserved <- setdiff(confounder_fields, batch_field)
  cat("Corrected for:", batch_field, "\n")
  if (length(preserved) > 0) {
    cat(
      "WARNING: strategy 'preserve_extra' corrects for the first discrete",
      "confounder only. Preserved, NOT corrected for:",
      paste(preserved, collapse = ", "), "\n"
    )
  }
  run_combat_seq(
    counts, batch = typed_coldata[[batch_field]], group = group,
    covar_mod = build_covar_mod(typed_coldata, preserved),
    label = sprintf(
      "batch = %s, protected: %s", batch_field,
      if (length(preserved) > 0) paste(preserved, collapse = ", ")
      else "condition only"
    )
  )
}

correct_composite <- function(counts, typed_coldata, group, discrete_fields,
                              continuous_fields) {
  if (length(continuous_fields) > 0) {
    stop(
      "ERROR: strategy 'composite' needs every confounder to be discrete, ",
      "but these are continuous: ", paste(continuous_fields, collapse = ", "),
      ". Use 'sequential' or 'preserve_extra' instead."
    )
  }
  batch <- interaction(
    lapply(discrete_fields, function(f) typed_coldata[[f]]), drop = TRUE
  )
  cat("Composite batch of:", paste(discrete_fields, collapse = " x "), "\n")
  cat("Cell sizes (", nlevels(batch), " cells):\n", sep = "")
  print(table(batch))
  run_combat_seq(
    counts, batch = batch, group = group, covar_mod = NULL,
    label = paste("composite batch of", paste(discrete_fields, collapse = ", "))
  )
}

apply_combat_seq <- function(counts, coldata, confounder_fields,
                             confounder_field_types, discrete_fields,
                             strategy = "sequential") {
  stopifnot(identical(colnames(counts), rownames(coldata)))
  if (length(confounder_fields) == 0) {
    stop(
      "ERROR: no confounders configured (design.confounders is empty); there ",
      "is nothing to correct for and this rule should not have been requested."
    )
  }
  stopifnot(all(discrete_fields %in% confounder_fields))
  continuous_fields <- setdiff(confounder_fields, discrete_fields)
  cat("Strategy:", strategy, "\n")
  cat("Discrete confounders:",
      if (length(discrete_fields) > 0) paste(discrete_fields, collapse = ", ")
      else "(none)", "\n")
  cat("Continuous confounders:",
      if (length(continuous_fields) > 0) paste(continuous_fields, collapse = ", ")
      else "(none)", "\n")

  if (strategy != "composite" && length(discrete_fields) == 0) {
    stop(
      "ERROR: all confounders are continuous, so ComBat-Seq has nothing to ",
      "correct for (its `batch` must be discrete). Model them in the DESeq2 ",
      "design instead, or use the limma diagnostic (PCA after limma ",
      "correction) to look at their effect."
    )
  }

  # Continuous fields numeric, everything else untouched (strings).
  typed_coldata <- cast_continuous_confounders(
    coldata, confounder_fields, confounder_field_types
  )
  group <- factor(coldata$condition)

  switch(
    strategy,
    sequential = correct_sequential(
      counts, typed_coldata, group, discrete_fields, continuous_fields
    ),
    preserve_extra = correct_preserve_extra(
      counts, typed_coldata, group, confounder_fields, discrete_fields
    ),
    composite = correct_composite(
      counts, typed_coldata, group, discrete_fields, continuous_fields
    ),
    stop("ERROR: unknown multi_confounder_strategy: ", strategy)
  )
}

# =====================================================================
# Main
# =====================================================================

counts_file <- snakemake@input[["counts"]]
coldata_file <- snakemake@input[["coldata"]]
confounder_fields <- unlist(snakemake@params[["confounder_fields"]])
confounder_field_types <- snakemake@params[["confounder_field_types"]]
discrete_confounder_fields <-
  unlist(snakemake@params[["discrete_confounder_fields"]])
strategy <- snakemake@params[["multi_confounder_strategy"]]

cat("Reading merged counts:", counts_file, "\n")
counts <- read_count_matrix(counts_file)
cat("Counts:", nrow(counts), "genes x", ncol(counts), "samples\n")

# Plain strings, as in deseq2_fit.R ("NA" is a value, not a missing one).
cat("Reading sample sheet:", coldata_file, "\n")
coldata <- as.data.frame(
  read_tsv(coldata_file, col_types = cols(.default = col_character()),
           na = character(), trim_ws = FALSE)
)
rownames(coldata) <- coldata$sample
if (!identical(colnames(counts), coldata$sample)) {
  stop("ERROR: the sample columns of the counts do not match the samples of ",
       "the sample sheet (same samples, same order expected).")
}
stopifnot(all(confounder_fields %in% colnames(coldata)))
stopifnot("condition" %in% colnames(coldata))

corrected <- apply_combat_seq(
  counts, coldata, confounder_fields, confounder_field_types,
  discrete_confounder_fields, strategy
)

# Final check of the combined result (every step was checked already).
validate_corrected_counts(corrected, counts)
if (max(corrected) > .Machine$integer.max) {
  stop("ERROR: a corrected count exceeds the integer range; cannot store ",
       "as integer.")
}
storage.mode(corrected) <- "integer"

cat("\nTotal entries changed:", sum(corrected != counts), "of",
    length(counts), "\n")
# Not re-filtered: the confounders_modelled and confounders_corrected analyses
# deliberately share one gene universe (the filter ran in merge_counts.R).
cat("Genes with rowSums < 10 after correction:",
    sum(rowSums(corrected) < 10),
    "(kept, so that all analyses use the same genes)\n")

# write_tsv() gzip-compresses on its own, going by the .gz file extension.
output_file <- snakemake@output[["counts"]]
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write_tsv(
  data.frame(gene_id = rownames(corrected), corrected,
             check.names = FALSE, row.names = NULL,
             stringsAsFactors = FALSE),
  output_file
)
cat("\nSaved corrected counts:", output_file, "\n")

cat("\n=== ComBat-Seq Correction Complete ===\n")
sink(type = "message")
sink()
close(logfile)
