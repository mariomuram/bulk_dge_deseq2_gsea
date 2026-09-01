#!/usr/bin/env Rscript
# Cross-contrast GSEA summary, run once per database/analysis
# (i.e. once per `gsea` wildcard -- gobp, hallmark, etc). Reads
# every contrast's already-significance-thresholded, unsimplified
# GSEA table (see scripts/gsea_postprocess.R for why it must be
# the unsimplified one), unions the pathways significant in ANY
# contrast, and builds:
#   - a pathway x contrast matrix of NES (padj alongside)
#   - a heatmap of that matrix (capped for legibility)
#   - a pairwise Jaccard-overlap table of each contrast pair's
#     significant-pathway sets
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== GSEA cross-contrast summary ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

gsea_name <- snakemake@params[["gsea_name"]]
contrast_names <- snakemake@params[["contrast_names"]]
padj_threshold <- snakemake@params[["padj_threshold"]]
heatmap_top_n <- snakemake@params[["heatmap_top_n"]]

cat("Analysis name:", gsea_name, "\n")
cat("Contrasts:", paste(contrast_names, collapse = ", "), "\n")
cat(
  "padj threshold (inherited from per-contrast runs):",
  padj_threshold, "\n\n"
)

input_paths <- unlist(snakemake@input)
if (length(input_paths) != length(contrast_names)) {
  stop(
    "ERROR: got ", length(input_paths), " input table(s) but ",
    length(contrast_names), " contrast name(s) in params. These ",
    "must line up positionally (both come from the same ",
    "CONTRAST_NAMES-ordered expand() in the Snakefile) -- check ",
    "the aggregate_gsea rule's input and params."
  )
}

# --- Load every contrast's significant-pathway table ------------------

cat("Loading per-contrast tables...\n")
per_contrast <- vector("list", length(contrast_names))
for (i in seq_along(contrast_names)) {
  contrast <- contrast_names[i]
  path <- input_paths[i]
  cat(" -", contrast, ":", path, "\n")
  df <- read_tsv(path, show_col_types = FALSE)
  if (!all(c("pathway", "NES", "padj") %in% colnames(df))) {
    stop(
      "ERROR: ", path, " is missing one of the expected columns ",
      "'pathway', 'NES', 'padj'. Got: ",
      paste(colnames(df), collapse = ", ")
    )
  }
  df$contrast <- contrast
  per_contrast[[i]] <- df %>% select(pathway, contrast, NES, padj)
}

long_df <- bind_rows(per_contrast)
n_sig_pathways <- n_distinct(long_df$pathway)
cat(
  "\nUnion of significant pathways across all contrasts:",
  n_sig_pathways, "\n"
)

if (nrow(long_df) == 0) {
  cat(
    "WARNING: no contrast had any significant pathway for '",
    gsea_name, "' (padj < ", padj_threshold, "). Writing empty ",
    "summary table, jaccard table, and placeholder heatmap.\n",
    sep = ""
  )
}

# --- Wide summary table: pathway rows, per-contrast NES/padj cols -----

cat("\nBuilding wide summary table...\n")
if (nrow(long_df) > 0) {
  nes_wide <- long_df %>%
    select(pathway, contrast, NES) %>%
    pivot_wider(
      names_from = contrast, values_from = NES,
      names_prefix = "NES_"
    )
  padj_wide <- long_df %>%
    select(pathway, contrast, padj) %>%
    pivot_wider(
      names_from = contrast, values_from = padj,
      names_prefix = "padj_"
    )
  summary_df <- nes_wide %>%
    left_join(padj_wide, by = "pathway") %>%
    arrange(pathway)
  cat(
    "NOTE: a pathway with NA for a given contrast means it did",
    "not clear that contrast's own padj <", padj_threshold,
    "threshold -- not necessarily that it was tested and found",
    "non-significant vs. never reaching significance at all;",
    "the per-contrast tables only ever contain",
    "already-thresholded results.\n"
  )
} else {
  summary_df <- data.frame(pathway = character())
  for (contrast in contrast_names) {
    summary_df[[paste0("NES_", contrast)]] <- numeric()
    summary_df[[paste0("padj_", contrast)]] <- numeric()
  }
}

cat("\nSaving summary table...\n")
output_table <- snakemake@output[["summary_table"]]
dir.create(dirname(output_table), showWarnings = FALSE, recursive = TRUE)
con <- gzfile(output_table, open = "wb", compression = 9)
write_tsv(summary_df, con)
close(con)
cat("Saved:", output_table, "\n")

# --- Heatmap (capped for legibility) -----------------------------

cat("\nGenerating summary heatmap...\n")
if (nrow(long_df) > 0) {
  pathway_rank <- long_df %>%
    group_by(pathway) %>%
    summarise(
      n_contrasts = n_distinct(contrast),
      max_abs_nes = max(abs(NES)),
      .groups = "drop"
    ) %>%
    arrange(desc(n_contrasts), desc(max_abs_nes))

  n_total_pathways <- nrow(pathway_rank)
  top_pathways <- pathway_rank %>% slice_head(n = heatmap_top_n)
  if (n_total_pathways > heatmap_top_n) {
    cat(
      "NOTE:", n_total_pathways, "significant pathway(s) found;",
      "heatmap limited to the top", heatmap_top_n, "(ranked by",
      "number of contrasts significant in, then by max |NES|).",
      "The full set is in the summary table.\n"
    )
  }

  heatmap_df <- long_df %>%
    filter(pathway %in% top_pathways$pathway) %>%
    mutate(
      pathway = factor(pathway, levels = rev(top_pathways$pathway)),
      contrast = factor(contrast, levels = contrast_names)
    )

  n_pathways_shown <- nrow(top_pathways)
  heatmap_plot <- ggplot(
    heatmap_df,
    aes(x = contrast, y = pathway, fill = NES)
  ) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0, na.value = "grey90"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = paste("Cross-contrast GSEA summary:", gsea_name),
      x = "", y = "", fill = "NES"
    )
  ggsave(
    snakemake@output[["heatmap"]], heatmap_plot,
    width = 10, height = max(6, 0.18 * n_pathways_shown + 2),
    dpi = 300, limitsize = FALSE
  )
  cat("Heatmap saved:", snakemake@output[["heatmap"]], "\n")
} else {
  cat("No pathways to plot.\n")
  empty_heatmap <- ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5,
      label = "No significant pathways found in any contrast",
      size = 5, color = "gray50"
    ) +
    theme_void()
  ggsave(
    snakemake@output[["heatmap"]], empty_heatmap,
    width = 8, height = 6, dpi = 300
  )
}

# --- Pairwise Jaccard overlap of significant-pathway sets --------

cat("\nComputing pairwise Jaccard overlap...\n")
sig_sets <- split(long_df$pathway, long_df$contrast)
# Contrasts with zero significant pathways still get an empty set,
# not a missing one, so every pair is still reported.
for (contrast in contrast_names) {
  if (is.null(sig_sets[[contrast]])) {
    sig_sets[[contrast]] <- character()
  }
}

pairs <- combn(contrast_names, 2, simplify = FALSE)
jaccard_rows <- lapply(pairs, function(p) {
  a <- unique(sig_sets[[p[1]]])
  b <- unique(sig_sets[[p[2]]])
  n_shared <- length(intersect(a, b))
  n_union <- length(union(a, b))
  jaccard <- if (n_union == 0) NA_real_ else n_shared / n_union
  data.frame(
    contrast_a = p[1], contrast_b = p[2],
    jaccard = jaccard, n_shared = n_shared, n_union = n_union
  )
})
jaccard_df <- bind_rows(jaccard_rows) %>% arrange(desc(jaccard))

cat("Pairwise Jaccard summary:\n")
print(jaccard_df)

cat("\nSaving Jaccard table...\n")
output_jaccard <- snakemake@output[["jaccard"]]
dir.create(dirname(output_jaccard), showWarnings = FALSE, recursive = TRUE)
con <- gzfile(output_jaccard, open = "wb", compression = 9)
write_tsv(jaccard_df, con)
close(con)
cat("Saved:", output_jaccard, "\n")

cat("\n=== GSEA cross-contrast summary complete ===\n")
sink(type = "message")
sink()
close(logfile)
