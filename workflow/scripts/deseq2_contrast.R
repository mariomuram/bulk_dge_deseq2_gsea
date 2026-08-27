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
      "WARNING: fewer than 2 replicates for a condition used in ",
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

  coef_index <- which(levels_condition == numerator)
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

  mod_mat <- model.matrix(design(dds), data = as.data.frame(colData(dds)))
  contrast_vec <- rep(0, ncol(mod_mat))
  names(contrast_vec) <- colnames(mod_mat)
  for (cond in names(groups)) {
    rows <- colData(dds)$condition == cond
    group_mean_row <- colMeans(mod_mat[rows, , drop = FALSE])
    contrast_vec <- contrast_vec + groups[[cond]] * group_mean_row
  }
  cat("\nResolved coefficient weights:\n")
  print(contrast_vec)

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

  # apeglm shrinkage only works when these weights exactly match
  # a single existing coefficient (weight 1 on one term, 0
  # elsewhere) -- true for the interaction term today.
  tol <- 1e-8
  nonzero <- which(abs(contrast_vec) > tol)
  if (length(nonzero) == 1 && abs(contrast_vec[nonzero] - 1) < tol &&
        names(contrast_vec)[nonzero] != "(Intercept)") {
    coef_name <- names(contrast_vec)[nonzero]
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
    "gene(s) have NA pvalue/padj: flagged as Cook's distance ",
    "outliers.\n"
  )
}
if (n_filtered_na > 0) {
  cat(
    "NOTE:", n_filtered_na,
    "gene(s) have NA padj (but a pvalue): removed by ",
    "independent filtering of low mean-count genes.\n"
  )
}

# --- Attempt apeglm shrinkage -----------------------------------------

if (!is.null(coef_name)) {
  cat("\nAttempting apeglm shrinkage using coef:", coef_name, "\n")
  res_shrunk <- tryCatch(
    lfcShrink(dds, coef = coef_name, res = res, type = "apeglm"),
    error = function(e) {
      cat("WARNING: apeglm shrinkage failed:", e$message, "\n")
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
  cat(
    "\nWARNING: this contrast's groups are not reducible to a ",
    "single model coefficient, so apeglm shrinkage is not ",
    "available for it. log2FoldChange_shrunk/lfcSE_shrunk will ",
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
