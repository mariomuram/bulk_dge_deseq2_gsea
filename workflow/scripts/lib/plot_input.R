GENES_OF_INTEREST_FALLBACK_N <- 10L

contrast_name <- snakemake@params[["contrast_name"]]
lfc_type <- snakemake@params[["lfc_type"]]
padj_threshold <- snakemake@params[["padj_threshold"]]
genes_of_interest <- snakemake@params[["genes_of_interest"]]
cat("Contrast:", contrast_name, "\n")
cat("LFC type:", lfc_type, "\n")
cat("padj threshold:", padj_threshold, "\n")
cat("Genes of interest:", genes_of_interest, "\n\n")

results_df <- read_tsv(snakemake@input[["results"]], show_col_types = FALSE)

lfc_col <- if (lfc_type == "shrunk") "log2FoldChange_shrunk" else "log2FoldChange"
if (!(lfc_col %in% colnames(results_df))) {
  stop("ERROR: column '", lfc_col, "' not found in results table.")
}
if (lfc_type == "shrunk" && all(is.na(results_df[[lfc_col]]))) {
  stop(
    "ERROR: log2FoldChange_shrunk is entirely NA for this ",
    "contrast -- it should not have been requested as a ",
    "shrinkable target. Check the Snakefile's ",
    "SHRINKABLE_CONTRASTS derivation."
  )
}

lfc_label <- if (lfc_type == "shrunk") {
  "Log2 Fold Change (apeglm-shrunk)"
} else {
  "Log2 Fold Change (unshrunk MLE)"
}

cat("\nResolving genes of interest...\n")
if (!is.null(genes_of_interest) && length(genes_of_interest) > 0) {
  tested_ids <- results_df$gene_id
  missing_ids <- setdiff(genes_of_interest, tested_ids)
  if (length(missing_ids) > 0) {
    cat(
      "WARNING:", length(missing_ids), "configured genes not found",
      "in this contrast's tested universe; omitting:",
      paste(missing_ids, collapse = ", "), "\n"
    )
  }
  genes_of_interest <- intersect(genes_of_interest, tested_ids)
  if (length(genes_of_interest) == 0) {
    stop(
      "ERROR: none of the configured `genes_of_interest` were found",
      " among the ", nrow(results_df), " gene(s) tested for",
      " contrast '", contrast_name, "'. ",
      "Check the ID(s) in `config[\"genes_of_interest\"]` against",
      " the id column of `results` for this contrast."
    )
  }
  cat(
    "Using", length(genes_of_interest), "configured gene(s) of",
    "interest.\n"
  )
} else {
  cat(
    "No genes of interest configured; falling back to the top",
    GENES_OF_INTEREST_FALLBACK_N, "by padj...\n"
  )
  genes_of_interest <- results_df %>%
    slice_min(padj, n = GENES_OF_INTEREST_FALLBACK_N) %>%
    pull(gene_id)
}

results_df <- results_df %>%
  # Note: empty string rather than `NA` to ensure unlabeled points are
  #  considered for `ggrepel`'s repulsion.
  mutate(label = ifelse(gene_id %in% genes_of_interest, gene_id, ""))
