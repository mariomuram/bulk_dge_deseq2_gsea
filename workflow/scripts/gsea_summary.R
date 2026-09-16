#!/usr/bin/env Rscript
# Cross-contrast GSEA summary, run once per database/analysis
# (i.e. once per `gsea` wildcard -- gobp, hallmark, etc). Reads
# every contrast's gsea_table_all (ALL tested gene sets, not just
# this contrast's own significant ones -- see gsea_postprocess.R).
# Unions the pathways significant in ANY contrast, then reports
# each unioned pathway's real NES/padj from every contrast it was
# tested in -- including contrasts where that pathway didn't
# independently clear padj_threshold. Filtering per-contrast before
# unioning would make a pathway with near-identical NES look
# contrast-specific whenever it happened to sit just either side of
# the threshold in different contrasts; this avoids that. Builds:
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
gene_sets_of_interest <- snakemake@params[["gene_sets_of_interest"]]
heatmap_top_n <- snakemake@params[["heatmap_top_n"]]

cat("Analysis name:", gsea_name, "\n")
cat("Contrasts:", paste(contrast_names, collapse = ", "), "\n")
cat(
  "padj threshold (applied here, per contrast, to decide union",
  "membership):", padj_threshold, "\n\n"
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

# --- Load every contrast's full (all tested gene sets) table ----------

cat("Loading per-contrast tables...\n")
per_contrast <- vector("list", length(contrast_names))
for (i in seq_along(contrast_names)) {
  contrast <- contrast_names[i]
  path <- input_paths[i]
  cat(" -", contrast, ":", path, "\n")
  df <- read_tsv(path, show_col_types = FALSE)
    if (!all(c("id", "description", "NES", "padj") %in% colnames(df))) {
    stop(
      "ERROR: ", path, " is missing one of the expected columns ",
      "'id', 'description', 'NES', 'padj'. Got: ",
      paste(colnames(df), collapse = ", ")
    )
  }
  df$contrast <- contrast
  per_contrast[[i]] <- df %>% select(id, description, contrast, NES, padj)
}

long_df <- bind_rows(per_contrast)
cat(
  "\nTotal (pathway, contrast) rows loaded across all tested gene",
  "sets:", nrow(long_df), "\n"
)

# --- Per-contrast significant sets, and their union --------------------

sig_by_contrast <- long_df %>% filter(padj < padj_threshold)
sig_sets <- split(sig_by_contrast$id, sig_by_contrast$contrast)
for (contrast in contrast_names) {
  if (is.null(sig_sets[[contrast]])) {
    sig_sets[[contrast]] <- character()
  }
  sig_sets[[contrast]] <- unique(sig_sets[[contrast]])
}

union_sig_ids <- unique(unlist(sig_sets, use.names = FALSE))
cat(
  "\nUnion of pathways significant in at least one contrast:",
  length(union_sig_ids), "\n"
)

if (length(union_sig_ids) == 0) {
  cat(
    "WARNING: no contrast had any significant pathway for '",
    gsea_name, "' (padj < ", padj_threshold, "). Writing empty ",
    "summary table, jaccard table, and placeholder heatmap.\n",
    sep = ""
  )
}

# --- Restrict to that union, but keep every contrast's real NES/padj ---

# Deliberately NOT filtering each row by that contrast's own
# significance: a pathway significant in contrast A still shows
# its actual (possibly just-short-of-threshold) NES/padj from
# contrast B here, rather than NA.
display_df <- long_df %>% filter(id %in% union_sig_ids)

# One description per ID, for display (heatmap axis, summary
# table). An ID should map to a single description reliably
# within a run; this defensively takes the first encountered
# per ID rather than assuming it, per pipeline convention.
id_description <- display_df %>%
  distinct(id, description) %>%
  group_by(id) %>%
  summarise(description = first(description), .groups = "drop")
n_multi_desc <- sum(duplicated(
  (display_df %>% distinct(id, description))$id
))
if (n_multi_desc > 0) {
  cat(
    "WARNING:", n_multi_desc, "gene set ID(s) mapped to more",
    "than one description across contrasts; keeping the first",
    "description encountered per ID.\n"
  )
}

# --- Wide summary table: pathway rows, per-contrast NES/padj cols -----

cat("\nBuilding wide summary table...\n")
if (nrow(display_df) > 0) {
  nes_wide <- display_df %>%
    select(id, contrast, NES) %>%
    pivot_wider(
      names_from = contrast, values_from = NES,
      names_prefix = "NES_"
    )
  padj_wide <- display_df %>%
    select(id, contrast, padj) %>%
    pivot_wider(
      names_from = contrast, values_from = padj,
      names_prefix = "padj_"
    )
  summary_df <- id_description %>%
    left_join(nes_wide, by = "id") %>%
    left_join(padj_wide, by = "id") %>%
    arrange(id)
  cat(
    "NOTE: a pathway with NA for a given contrast means it was",
    "not tested there at all (e.g. too few of its member genes",
    "survived that contrast's own ranking-stat filtering to",
    "clear min_set_size) -- NOT that it fell short of",
    "significance. Every contrast's real NES/padj is shown for",
    "any pathway significant in at least one contrast, even",
    "where that particular contrast didn't independently clear",
    "padj <", padj_threshold, ".\n"
  )
} else {
  summary_df <- data.frame(id = character(), description = character())
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
if (nrow(display_df) > 0) {
  # Ranked by how many contrasts each pathway is SIGNIFICANT in
  # (not merely reported in -- with display_df now carrying every
  # unioned pathway's NES/padj across nearly every contrast, "rows
  # present" is no longer a useful ranking signal on its own).
  sig_membership <- bind_rows(lapply(contrast_names, function(ct) {
    data.frame(id = sig_sets[[ct]], contrast = ct)
  }))
  pathway_rank <- sig_membership %>%
    count(id, name = "n_contrasts_significant") %>%
    left_join(
      display_df %>%
        group_by(id) %>%
        summarise(
          max_abs_nes = max(abs(NES), na.rm = TRUE),
          .groups = "drop"
        ),
      by = "id"
    ) %>%
    left_join(id_description, by = "id") %>%
    arrange(desc(n_contrasts_significant), desc(max_abs_nes))

  if (!is.null(gene_sets_of_interest) && length(gene_sets_of_interest) > 0) {
    missing_ids <- setdiff(gene_sets_of_interest, union_sig_ids)
    if (length(missing_ids) > 0) {
      cat(
        "WARNING:", length(missing_ids), "configured gene set ID(s)",
        "not significant in any contrast; omitting:",
        paste(missing_ids, collapse = ", "), "\n"
      )
    }
    gene_sets_of_interest <- intersect(gene_sets_of_interest, union_sig_ids)
    if (length(gene_sets_of_interest) == 0) {
      stop(
        "ERROR: none of the configured `gene_sets_of_interest` for '",
        gsea_name, "' were found among the ", length(union_sig_ids),
        " gene set(s) significant across contrasts. ",
        "Check the ID(s) in `config[\"gsea\"][\"", gsea_name, "\"]",
        "[\"gene_sets_of_interest\"]` against the id column of ",
        "`gsea_table_signif` for all contrasts."
      )
    }
    cat(
      "Using", length(gene_sets_of_interest), "configured gene set(s)",
      "of interest.\n"
    )
  } else {
    cat(
      "No gene sets of interest configured; falling back to the top",
      heatmap_top_n, "by padj...\n"
    )
    gene_sets_of_interest <- pathway_rank %>%
      slice_head(n = heatmap_top_n) %>%
      pull(id)

    n_total_pathways <- nrow(pathway_rank)
    if (n_total_pathways > heatmap_top_n) {
      cat(
        "NOTE:", n_total_pathways, "significant pathway(s) found;",
        "heatmap limited to the top", heatmap_top_n, "(ranked by",
        "number of contrasts significant in, then by max |NES|).",
        "The full set is in the summary table.\n"
      )
    }
  }


  gene_sets_of_interest <- pathway_rank %>%
    filter(id %in% gene_sets_of_interest)

  heatmap_df <- display_df %>%
    filter(id %in% gene_sets_of_interest$id) %>%
    select(-description) %>%
    left_join(id_description, by = "id") %>%
    mutate(
      description = factor(description,
                           levels = rev(gene_sets_of_interest$description)),
      contrast = factor(contrast, levels = contrast_names),
      significant = factor(
        padj < padj_threshold,
        levels = c(TRUE, FALSE)
      )
    )

  n_pathways_shown <- nrow(gene_sets_of_interest)
  heatmap_plot <- ggplot(
    heatmap_df,
    aes(x = contrast, y = description, fill = NES, alpha = significant)
  ) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0, na.value = "grey90"
    ) +
    scale_alpha_manual(
      values = c(`TRUE` = 1, `FALSE` = 0.35),
      labels = c(
        `TRUE` = paste0("padj < ", padj_threshold),
        `FALSE` = paste0("padj \u2265 ", padj_threshold)
      ),
      name = "Significant",
      drop = FALSE
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
pairs <- combn(contrast_names, 2, simplify = FALSE)
jaccard_rows <- lapply(pairs, function(p) {
  a <- sig_sets[[p[1]]]
  b <- sig_sets[[p[2]]]
  n_shared <- length(intersect(a, b))
  n_union <- length(union(a, b))
  jaccard <- if (n_union == 0) NA_real_ else n_shared / n_union
  data.frame(
    contrast_a = p[1], contrast_b = p[2],
    jaccard = jaccard, n_shared = n_shared, n_union = n_union
  )
})
jaccard_df <- bind_rows(jaccard_rows) %>% arrange(desc(jaccard))

cat("Pairwise Jaccard summary (on significant-pathway sets):\n")
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
