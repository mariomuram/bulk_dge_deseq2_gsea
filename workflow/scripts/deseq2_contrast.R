#!/usr/bin/env Rscript
# DESeq2 contrast extraction, with apeglm LFC shrinkage
# attempted for every contrast.
#
# Pairwise contrasts (numerator/denominator) run against a dds
# already releveled to `denominator` under ~condition (see
# deseq2_relevel.R) -- both results(contrast = c("condition",
# ...)) and apeglm's coef= shrinkage are available via the
# ordinary, well-supported path there.
#
# groups-style contrasts (the interaction term, or any future
# multi-term combination) run against the primary
# interaction-design dds via a group-averaged model-matrix
# contrast. Shrinkage for these is only possible when the
# resolved weights happen to exactly match a single existing
# model coefficient (true for the interaction term today);
# otherwise log2FoldChange_shrunk/lfcSE_shrunk are left NA with
# a clear warning -- the unshrunk columns are unaffected.

library(DESeq2)
library(apeglm)
library(readr)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== DESeq2 Contrast Extraction ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

dds_file <- snakemake@input[["dds_rds"]]
cat("Loading DESeqDataSet:", dds_file, "\n")
dds <- readRDS(dds_file)

spec <- snakemake@params[["contrast_spec"]]
contrast_type <- spec[["type"]]
padj_threshold <- snakemake@params[["padj_threshold"]]
log2fc_threshold <- snakemake@params[["log2fc_threshold"]]

cat("Contrast type:", contrast_type, "\n")

check_replicates <- function(dds, conditions) {
  n_per_condition <- table(colData(dds)$condition)[conditions]
  cat("\nReplicates per group:\n")
  print(n_per_condition)
  if (any(n_per_condition < 2)) {
    cat(
      "WARNING: fewer than 2 replicates for a condition used in",
      "this contrast. Results may be unreliable.\n"
    )
  }
}

# --- Extraction + coefficient resolution ----------------------------

coef_name <- NULL  # resolved iff this contrast is apeglm-shrinkable

if (contrast_type == "pairwise") {

  numerator <- spec[["numerator"]]
  denominator <- spec[["denominator"]]
  cat("Numerator:  ", numerator, "\n")
  cat("Denominator:", denominator, "\n")

  levels_condition <- levels(factor(colData(dds)$condition))
  missing <- setdiff(c(numerator, denominator), levels_condition)
  if (length(missing) > 0) {
    stop(
      "ERROR: condition(s) not found in the fitted model: ",
      paste(missing, collapse = ", "),
      ". Available: ", paste(levels_condition, collapse = ", ")
    )
  }
  if (levels_condition[1] != denominator) {
    stop(
      "ERROR: expected dds to be releveled with reference = '",
      denominator, "' but found reference = '", levels_condition[1],
      "'. Check that the Snakefile mapped this contrast to the ",
      "correct dds_object_relevel_*.rds file."
    )
  }

  check_replicates(dds, c(numerator, denominator))

  cat("\nExtracting results (character-triple contrast)...\n")
  res <- tryCatch(
    results(
      dds,
      contrast = c("condition", numerator, denominator),
      independentFiltering = TRUE,
      alpha = padj_threshold
    ),
    error = function(e) {
      stop("ERROR: results() failed for this contrast: ", e$message)
    }
  )

  n_condition_levels <- length(levels_condition)
  n_total_coefs <- length(resultsNames(dds))
  if (n_total_coefs < n_condition_levels) {
    stop(
      "ERROR: resultsNames(dds) has ", n_total_coefs, " entries, fewer ",
      "than the ", n_condition_levels, " condition levels. Positional ",
      "coefficient resolution is unsafe here; check that ",
      "deseq2_relevel.R fit '~ <confounders> + condition' as expected."
    )
  }
  if (resultsNames(dds)[1] != "Intercept") {
    stop(
      "ERROR: expected resultsNames(dds)[1] == 'Intercept', got '",
      resultsNames(dds)[1], "'. Positional coefficient resolution is ",
      "unsafe here; investigate before trusting downstream results."
    )
  }
  n_confounder_coefs <- n_total_coefs - n_condition_levels
  coef_index <- n_confounder_coefs + which(levels_condition == numerator)
  coef_name <- resultsNames(dds)[coef_index]

} else if (contrast_type == "groups") {

  groups <- spec[["groups"]]
  cat("Groups (condition: weight):\n")
  for (g in names(groups)) cat(" -", g, ":", groups[[g]], "\n")

  available_conditions <- levels(factor(colData(dds)$condition))
  missing <- setdiff(names(groups), available_conditions)
  if (length(missing) > 0) {
    stop(
      "ERROR: condition(s) not found in the fitted model: ",
      paste(missing, collapse = ", "),
      ". Available: ", paste(available_conditions, collapse = ", ")
    )
  }

  check_replicates(dds, names(groups))

  # design(dds) is only a placeholder (~condition) on the
  # primary fit -- the matrix actually used by DESeq() was a
  # cleaned full-crossing matrix built in deseq2_fit.R and
  # stashed here to avoid recomputing (and risking drift from)
  # it.
  mod_mat <- metadata(dds)[["full_model_matrix"]]
  if (is.null(mod_mat)) {
    stop(
      "ERROR: metadata(dds)$full_model_matrix is missing. This ",
      "dds does not look like it was produced by the current ",
      "deseq2_fit.R -- was it built with an older version of ",
      "the pipeline?"
    )
  }
  contrast_vec <- rep(0, ncol(mod_mat))
  names(contrast_vec) <- colnames(mod_mat)
  for (cond in names(groups)) {
    rows <- colData(dds)$condition == cond
    group_mean_row <- colMeans(mod_mat[rows, , drop = FALSE])
    contrast_vec <- contrast_vec + groups[[cond]] * group_mean_row
  }
  cat("\nResolved coefficient weights:\n")
  print(contrast_vec)

  # Determine single-coefficient reducibility BEFORE extracting
  # results, so that -- when reducible -- results() is called
  # via name= rather than the generic numeric contrast vector.
  # This matters beyond cosmetics: lfcShrink()'s own consistency
  # check rejects a `res` built from a raw numeric contrast even
  # when that vector is numerically identical to a single
  # coefficient's unit vector ("'coef' should specify same
  # coefficient as in results 'res'") -- it checks how `res` was
  # built, not just its values. name= and a matching unit-vector
  # contrast= produce numerically identical results() output, so
  # this only changes res's internal provenance tag, not any
  # extracted number.
  #
  # Intercept is excluded by POSITION (always column 1 for a
  # standard, non-"~0+..." formula), not by name -- full_mat's
  # column names now come from resultsNames(dds) (synced in
  # deseq2_fit.R), which calls it "Intercept" with no
  # parentheses, so a name-based check here would need to track
  # that convention by hand for no real benefit over position.
  if (length(resultsNames(dds)) != ncol(mod_mat)) {
    stop(
      "ERROR: resultsNames(dds) has ", length(resultsNames(dds)),
      " entries but the stored full_model_matrix has ",
      ncol(mod_mat), " columns. Positional coefficient ",
      "resolution is unsafe here; investigate before trusting ",
      "downstream results."
    )
  }
  tol <- 1e-8
  nonzero <- which(abs(contrast_vec) > tol)
  if (length(nonzero) == 1 && abs(contrast_vec[nonzero] - 1) < tol &&
        nonzero != 1) {
    coef_name <- resultsNames(dds)[nonzero]
  }

  if (!is.null(coef_name)) {
    cat(
      "\nGroups reduce to a single coefficient (", coef_name,
      "); extracting via name= for apeglm compatibility.\n",
      sep = ""
    )
    res <- tryCatch(
      results(
        dds,
        name = coef_name,
        independentFiltering = TRUE,
        alpha = padj_threshold
      ),
      error = function(e) {
        stop("ERROR: results() failed for this contrast: ", e$message)
      }
    )
  } else {
    cat("\nExtracting results (numeric contrast vector)...\n")
    res <- tryCatch(
      results(
        dds,
        contrast = contrast_vec,
        independentFiltering = TRUE,
        alpha = padj_threshold
      ),
      error = function(e) {
        stop("ERROR: results() failed for this contrast: ", e$message)
      }
    )
  }

} else {
  stop("ERROR: unknown contrast type: ", contrast_type)
}

# --- Build the results table -----------------------------------------

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res)
res_df <- res_df[, c(
  "gene_id", "baseMean", "log2FoldChange", "lfcSE",
  "stat", "pvalue", "padj"
)]

n_outlier_na <- sum(is.na(res_df$pvalue))
n_filtered_na <- sum(is.na(res_df$padj) & !is.na(res_df$pvalue))
if (n_outlier_na > 0) {
  cat(
    "\nNOTE:", n_outlier_na,
    "gene(s) have NA pvalue/padj: flagged as Cook's distance",
    "outliers.\n"
  )
}
if (n_filtered_na > 0) {
  cat(
    "NOTE:", n_filtered_na,
    "gene(s) have NA padj (but a pvalue): removed by",
    "independent filtering of low mean-count genes.\n"
  )
}

# --- Attempt apeglm shrinkage -----------------------------------------

expected_shrinkable <- isTRUE(snakemake@params[["expected_shrinkable"]])

if (!is.null(coef_name)) {
  cat("\nAttempting apeglm shrinkage using coef:", coef_name, "\n")
  res_shrunk <- tryCatch(
    lfcShrink(dds, coef = coef_name, res = res, type = "apeglm"),
    error = function(e) {
      msg <- conditionMessage(e)
      if (expected_shrinkable) {
        stop(
          "ERROR: apeglm shrinkage failed for a contrast the ",
          "Snakemake-level pre-check (SHRINKABLE_CONTRASTS) ",
          "predicted was shrinkable (coef = '", coef_name, "'). ",
          "This is a genuine disagreement between the Python ",
          "pre-check and this runtime call, not an expected ",
          "unshrinkable case -- investigate before trusting this ",
          "run. Underlying lfcShrink() error: ", msg
        )
      }
      cat("WARNING: apeglm shrinkage failed:", msg, "\n")
      NULL
    }
  )
  if (!is.null(res_shrunk)) {
    idx <- match(res_df$gene_id, rownames(res_shrunk))
    res_df$log2FoldChange_shrunk <- res_shrunk$log2FoldChange[idx]
    res_df$lfcSE_shrunk <- res_shrunk$lfcSE[idx]
    cat("Shrinkage complete.\n")
  } else {
    res_df$log2FoldChange_shrunk <- NA_real_
    res_df$lfcSE_shrunk <- NA_real_
  }
} else {
  if (expected_shrinkable) {
    stop(
      "ERROR: the Snakemake-level pre-check (SHRINKABLE_CONTRASTS) ",
      "predicted this contrast's groups reduce to a single model ",
      "coefficient, but this runtime check found no matching ",
      "unit-weight coefficient in the resolved contrast vector ",
      "(printed above). This is a genuine disagreement between ",
      "the two implementations -- investigate before trusting ",
      "this run. Compare the printed contrast_vec names against ",
      "resultsNames(dds) from the primary fit's log."
    )
  }
  cat(
    "\nWARNING: this contrast's groups are not reducible to a",
    "single model coefficient, so apeglm shrinkage is not",
    "available for it. log2FoldChange_shrunk/lfcSE_shrunk will",
    "be NA for every gene. The unshrunk log2FoldChange/stat/",
    "pvalue/padj columns are unaffected and remain fully valid.\n"
  )
  res_df$log2FoldChange_shrunk <- NA_real_
  res_df$lfcSE_shrunk <- NA_real_
}

res_df <- res_df[order(res_df$padj), ]

n_sig <- sum(
  res_df$padj < padj_threshold &
    abs(res_df$log2FoldChange) >= log2fc_threshold,
  na.rm = TRUE
)
cat(
  "\nSignificant genes (padj <", padj_threshold,
  "& |log2FC| >=", log2fc_threshold, "):", n_sig, "\n"
)

# --- Save (gzip --best compressed TSV) -----------------------------

output_file <- snakemake@output[["results"]]
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
con <- gzfile(output_file, open = "wb", compression = 9)
write_tsv(res_df, con)
close(con)
cat("\nSaved:", output_file, "\n")

cat("\n=== DESeq2 Contrast Extraction Complete ===\n")
sink(type = "message")
sink()
close(logfile)
