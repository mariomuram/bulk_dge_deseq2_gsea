contrast_name <- snakemake@params[["contrast_name"]]
lfc_type <- snakemake@params[["lfc_type"]]
padj_threshold <- snakemake@params[["padj_threshold"]]
cat("Contrast:", contrast_name, "\n")
cat("LFC type:", lfc_type, "\n")
cat("padj threshold:", padj_threshold, "\n")

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
