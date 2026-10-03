#!/usr/bin/env Rscript
# Cross-group GSEA summary, run once per database/analysis (i.e. once per
# `gsea` wildcard -- gobp, hallmark, etc). "Group" is whatever
# `grouping_var` names: "contrast" today, potentially "strategy" later
# (reuse via `use rule aggregate_gsea as ... with: params:`).
#
# Reads every group's gsea_table_all (ALL tested gene sets, not just that
# group's own significant ones -- see gsea_postprocess.R). Unions the
# pathways significant in ANY group, then reports each unioned pathway's
# real NES/padj from every group it was tested in -- including groups where
# it didn't independently clear padj_threshold. Filtering per group before
# unioning would make a pathway with near-identical NES look group-specific
# whenever it sat just either side of the threshold in different groups.
# NA in the output means "not tested", never "not significant".
# Builds:
#   - a pathway x group matrix of NES (padj alongside)
#   - a heatmap of that matrix (capped for legibility)
#   - a pairwise Jaccard-overlap table of each group pair's
#     significant-pathway sets (unrestricted by the heatmap cap)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# =====================================================================
# Functions
# =====================================================================

# Read every group's table, validate columns, row-bind with a column named
# by the VALUE of grouping_var holding each table's group label.
load_group_tables <- function(input_paths, group_names, grouping_var) {
  required <- c("id", "description", "NES", "padj")
  if (grouping_var %in% required) {
    stop("ERROR: grouping_var '", grouping_var, "' collides with a ",
         "required column name (", paste(required, collapse = ", "), ").")
  }
  if (length(input_paths) != length(group_names)) {
    stop(
      "ERROR: got ", length(input_paths), " input table(s) but ",
      length(group_names), " group name(s) in params. These must line ",
      "up positionally (both come from the same ordered expand() in the ",
      "Snakefile) -- check the aggregate_gsea rule's input and params. ",
      "Input paths seen: ", paste(input_paths, collapse = ", ")
    )
  }
  per_group <- vector("list", length(group_names))
  for (i in seq_along(group_names)) {
    path <- input_paths[i]
    cat(" -", group_names[i], ":", path, "\n")
    df <- read_tsv(path, show_col_types = FALSE)
    if (!all(required %in% colnames(df))) {
      stop(
        "ERROR: ", path, " is missing one of the expected columns ",
        "'id', 'description', 'NES', 'padj'. Got: ",
        paste(colnames(df), collapse = ", ")
      )
    }
    df[[grouping_var]] <- group_names[i]
    per_group[[i]] <- df %>%
      select(all_of(c("id", "description", grouping_var, "NES", "padj")))
  }
  bind_rows(per_group)
}

# Per-group significant sets, their union, and the union-restricted table.
# display_df deliberately keeps every group's real NES/padj for each
# unioned pathway (not filtered by that group's own significance).
compute_group_union <- function(long_df, grouping_var, group_names,
                                padj_threshold) {
  sig_by_group <- long_df %>% filter(padj < padj_threshold)
  sig_sets <- split(sig_by_group$id, sig_by_group[[grouping_var]])
  for (g in group_names) {
    if (is.null(sig_sets[[g]])) {
      sig_sets[[g]] <- character()
    }
    sig_sets[[g]] <- unique(sig_sets[[g]])
  }
  sig_sets <- sig_sets[group_names]

  union_sig_ids <- unique(unlist(sig_sets, use.names = FALSE))
  display_df <- long_df %>% filter(id %in% union_sig_ids)

  # One description per ID (first encountered), defensively.
  id_description <- display_df %>%
    distinct(id, description) %>%
    group_by(id) %>%
    summarise(description = first(description), .groups = "drop")
  n_multi_desc <- sum(duplicated(
    (display_df %>% distinct(id, description))$id
  ))

  list(
    sig_sets = sig_sets, union_sig_ids = union_sig_ids,
    display_df = display_df, id_description = id_description,
    n_multi_desc = n_multi_desc
  )
}

# Wide table: pathway rows, NES_{group}/padj_{group} columns.
build_summary_table <- function(display_df, id_description, group_names,
                                grouping_var) {
  if (nrow(display_df) > 0) {
    nes_wide <- display_df %>%
      select(id, all_of(grouping_var), NES) %>%
      pivot_wider(
        names_from = all_of(grouping_var), values_from = NES,
        names_prefix = "NES_"
      )
    padj_wide <- display_df %>%
      select(id, all_of(grouping_var), padj) %>%
      pivot_wider(
        names_from = all_of(grouping_var), values_from = padj,
        names_prefix = "padj_"
      )
    id_description %>%
      left_join(nes_wide, by = "id") %>%
      left_join(padj_wide, by = "id") %>%
      arrange(id)
  } else {
    summary_df <- data.frame(id = character(), description = character())
    for (g in group_names) {
      summary_df[[paste0("NES_", g)]] <- numeric()
      summary_df[[paste0("padj_", g)]] <- numeric()
    }
    summary_df
  }
}

# --- Ranking heuristics ----------------------------------------------
# All rank_by_* share one signature (display_df, sig_sets, group_names) so
# rank_pathways() can dispatch uniformly; each returns a tibble with `id`
# plus its metric columns, already in rank order.

# Default: n groups SIGNIFICANT in (desc), then max |NES| (desc). Ranked by
# significance count, not "rows present": display_df carries every unioned
# pathway across nearly every group, so presence is not informative.
rank_by_breadth <- function(display_df, sig_sets, group_names) {
  tibble(id = unlist(sig_sets[group_names], use.names = FALSE)) %>%
    count(id, name = "n_groups_significant") %>%
    left_join(
      display_df %>%
        group_by(id) %>%
        summarise(
          max_abs_nes = max(abs(NES), na.rm = TRUE),
          .groups = "drop"
        ),
      by = "id"
    ) %>%
    arrange(desc(n_groups_significant), desc(max_abs_nes))
}

# NES range across groups (desc): which pathways' calls are most sensitive
# to the group/strategy used. Pathways tested in <2 groups get NA (nothing
# to diverge from) and sort last. NOTE: untested as a *choice* until a real
# strategy comparison exists to look at -- inclusion is not validation.
rank_by_divergence <- function(display_df, sig_sets, group_names) {
  display_df %>%
    group_by(id) %>%
    summarise(
      n_groups_tested = sum(!is.na(NES)),
      nes_range = if (sum(!is.na(NES)) < 2) {
        NA_real_
      } else {
        diff(range(NES, na.rm = TRUE))
      },
      .groups = "drop"
    ) %>%
    arrange(desc(nes_range), id)
}

rank_pathways <- function(display_df, sig_sets, id_description, group_names,
                          rank_fn = rank_by_breadth) {
  rank_fn(display_df, sig_sets, group_names) %>%
    left_join(id_description, by = "id")
}

# gene_sets_of_interest override vs. top-N. Returns the ranked rows
# (id, description, ...) to draw, in rank order.
select_heatmap_pathways <- function(pathway_rank, union_sig_ids,
                                    gene_sets_of_interest, heatmap_top_n,
                                    gsea_name = "") {
  if (!is.null(gene_sets_of_interest) && length(gene_sets_of_interest) > 0) {
    missing_ids <- setdiff(gene_sets_of_interest, union_sig_ids)
    if (length(missing_ids) > 0) {
      cat(
        "WARNING:", length(missing_ids), "configured gene set ID(s)",
        "not significant in any group; omitting:",
        paste(missing_ids, collapse = ", "), "\n"
      )
    }
    chosen <- intersect(gene_sets_of_interest, union_sig_ids)
    if (length(chosen) == 0) {
      stop(
        "ERROR: none of the configured `gene_sets_of_interest` for '",
        gsea_name, "' were found among the ", length(union_sig_ids),
        " gene set(s) significant across groups. ",
        "Check the ID(s) in `config[\"gsea\"][\"", gsea_name, "\"]",
        "[\"gene_sets_of_interest\"]` against the id column of ",
        "`gsea_table_signif` for all groups."
      )
    }
    cat("Using", length(chosen), "configured gene set(s) of interest.\n")
  } else {
    cat(
      "No gene sets of interest configured; falling back to the top",
      heatmap_top_n, "by the selected ranking heuristic...\n"
    )
    chosen <- pathway_rank %>% slice_head(n = heatmap_top_n) %>% pull(id)
    if (nrow(pathway_rank) > heatmap_top_n) {
      cat(
        "NOTE:", nrow(pathway_rank), "significant pathway(s) found;",
        "heatmap limited to the top", heatmap_top_n,
        "(rank order per the heuristic logged above).",
        "The full set is in the summary table.\n"
      )
    }
  }
  pathway_rank %>% filter(id %in% chosen)
}

# Layout for the GSEA heatmap. Row pitch shrinks with row count (0.28 in up
# to 30 rows, floored at 0.2 in); y-label size follows the pitch (capped at
# the theme's base size); width grows with the longest pathway label.
# Best-guess constants -- tweak after eyeballing real output.
gsea_heatmap_layout <- function(selected_pathways, group_names) {
  n_rows <- nrow(selected_pathways)
  row_in <- max(0.2, min(0.28, 8.4 / n_rows))
  y_text_pt <- min(base_size, 0.7 * row_in * 72)  # base_size: plot_theme.R
  max_chars <- max(0, nchar(as.character(selected_pathways$description)),
                   na.rm = TRUE)
  list(
    y_text_pt = y_text_pt,
    height = max(5, row_in * n_rows + 3),
    width = max(8, 3 + 1.2 * length(group_names) +
                  max_chars * y_text_pt * 0.55 / 72)
  )
}

plot_gsea_heatmap <- function(display_df, id_description, selected_pathways,
                              group_names, grouping_var, padj_threshold,
                              color_up, color_down, color_none, color_na,
                              gsea_name, layout, confounder_handling = NULL) {
  heatmap_df <- display_df %>%
    filter(id %in% selected_pathways$id) %>%
    select(-description) %>%
    left_join(id_description, by = "id") %>%
    mutate(
      description = factor(description,
                           levels = rev(selected_pathways$description))
    )
  heatmap_df[[grouping_var]] <- factor(heatmap_df[[grouping_var]],
                                       levels = group_names)

  ggplot(
    heatmap_df,
    aes(x = .data[[grouping_var]], y = description, fill = NES)
  ) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(
      aes(label = significance_asterisks(padj, padj_threshold)),
      vjust = 0.75, size = default_text_size
    ) +
    scale_fill_gradient2(
      low = color_down, mid = color_none, high = color_up,
      midpoint = 0, na.value = color_na
    ) +
    theme_heatmap_no_grid() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = layout$y_text_pt)
    ) +
    labs(
      title = paste0(
        "Cross-", grouping_var, " GSEA summary: ",
        format_plot_title(gsea_name, confounder_handling)
      ),
      x = "", y = "", fill = "NES",
      caption = paste0("* = padj < ", padj_threshold)
    )
}

# Output columns are {grouping_var}_a / {grouping_var}_b
# (contrast_a / contrast_b for today's grouping_var = "contrast").
compute_pairwise_jaccard <- function(sig_sets, group_names, grouping_var) {
  pairs <- combn(group_names, 2, simplify = FALSE)
  rows <- lapply(pairs, function(p) {
    a <- sig_sets[[p[1]]]
    b <- sig_sets[[p[2]]]
    n_shared <- length(intersect(a, b))
    n_union <- length(union(a, b))
    data.frame(
      group_a = p[1], group_b = p[2],
      jaccard = if (n_union == 0) NA_real_ else n_shared / n_union,
      n_shared = n_shared, n_union = n_union
    )
  })
  out <- bind_rows(rows) %>% arrange(desc(jaccard))
  colnames(out)[1:2] <- paste0(grouping_var, c("_a", "_b"))
  out
}

# =====================================================================
# Main
# =====================================================================

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== GSEA cross-group summary ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

require_param <- function(name) {
  v <- snakemake@params[[name]]
  if (is.null(v)) {
    stop("ERROR: required param '", name, "' missing from the rule.")
  }
  v
}
get_input <- function(name) {
  v <- snakemake@input[[name]]
  if (is.null(v)) {
    stop("ERROR: required input '", name, "' missing from the rule.")
  }
  v
}

source(get_input("plot_theme_lib"))
source(get_input("heatmap_theme_lib"))
source(get_input("plot_title_lib"))

gsea_name <- require_param("gsea_name")
grouping_var <- require_param("grouping_var")
group_names <- require_param("group_names")
padj_threshold <- require_param("padj_threshold")
heatmap_top_n <- require_param("heatmap_top_n")
color_up <- require_param("color_up")
color_down <- require_param("color_down")
color_none <- require_param("color_none")
color_na <- require_param("color_na")
gene_sets_of_interest <- snakemake@params[["gene_sets_of_interest"]]
rank_heuristic <- snakemake@params[["rank_heuristic"]]
if (is.null(rank_heuristic)) rank_heuristic <- "breadth"
# NULL if the rule has no `confounder_handling` wildcard (a summary that
# compares confounder handlings has them as groups instead).
confounder_handling <- snakemake@wildcards[["confounder_handling"]]

rank_fns <- list(breadth = rank_by_breadth, divergence = rank_by_divergence)
if (!rank_heuristic %in% names(rank_fns)) {
  stop(
    "ERROR: unknown rank_heuristic '", rank_heuristic, "'. Expected one ",
    "of: ", paste(names(rank_fns), collapse = ", ")
  )
}

cat("Analysis name:", gsea_name, "\n")
cat("Confounder handling:", confounder_handling, "\n")
cat("Grouping variable:", grouping_var, "\n")
cat("Groups:", paste(group_names, collapse = ", "), "\n")
cat("Ranking heuristic:", rank_heuristic, "\n")
cat(
  "padj threshold (applied here, per group, to decide union",
  "membership):", padj_threshold, "\n\n"
)

# The per-group tables are the named `gsea_tables` input, ordered as
# group_names (both come from CONTRAST_NAMES in the Snakefile). Looked up by
# name so the lib-file inputs never enter the picture.
input_paths <- unlist(get_input("gsea_tables"), use.names = FALSE)

cat("Loading per-group tables...\n")
long_df <- load_group_tables(input_paths, group_names, grouping_var)
cat(
  "\nTotal (pathway, group) rows loaded across all tested gene",
  "sets:", nrow(long_df), "\n"
)

grp <- compute_group_union(long_df, grouping_var, group_names, padj_threshold)
sig_sets <- grp$sig_sets
union_sig_ids <- grp$union_sig_ids
display_df <- grp$display_df
id_description <- grp$id_description
cat(
  "\nUnion of pathways significant in at least one group:",
  length(union_sig_ids), "\n"
)
if (length(union_sig_ids) == 0) {
  cat(
    "WARNING: no group had any significant pathway for '",
    gsea_name, "' (padj < ", padj_threshold, "). Writing empty ",
    "summary table, jaccard table, and placeholder heatmap.\n",
    sep = ""
  )
}
if (grp$n_multi_desc > 0) {
  cat(
    "WARNING:", grp$n_multi_desc, "gene set ID(s) mapped to more",
    "than one description across groups; keeping the first",
    "description encountered per ID.\n"
  )
}

# --- Wide summary table -------------------------------------------------

cat("\nBuilding wide summary table...\n")
summary_df <- build_summary_table(
  display_df, id_description, group_names, grouping_var
)
if (nrow(display_df) > 0) {
  cat(
    "NOTE: a pathway with NA for a given group means it was",
    "not tested there at all (e.g. too few of its member genes",
    "survived that group's own ranking-stat filtering to",
    "clear min_set_size) -- NOT that it fell short of",
    "significance. Every group's real NES/padj is shown for",
    "any pathway significant in at least one group, even",
    "where that particular group didn't independently clear",
    "padj <", padj_threshold, ".\n"
  )
}
cat("\nSaving summary table...\n")
summary_table_file <- snakemake@output[["summary_table"]]
dir.create(dirname(summary_table_file), showWarnings = FALSE, recursive = TRUE)
write_tsv(summary_df, summary_table_file)
cat("Saved:", summary_table_file, "\n")

# --- Heatmap (capped for legibility) -----------------------------------

cat("\nGenerating summary heatmap...\n")
if (nrow(display_df) > 0) {
  pathway_rank <- rank_pathways(
    display_df, sig_sets, id_description, group_names,
    rank_fn = rank_fns[[rank_heuristic]]
  )
  selected_pathways <- select_heatmap_pathways(
    pathway_rank, union_sig_ids, gene_sets_of_interest, heatmap_top_n,
    gsea_name
  )
  layout <- gsea_heatmap_layout(selected_pathways, group_names)
  heatmap_plot <- plot_gsea_heatmap(
    display_df, id_description, selected_pathways, group_names,
    grouping_var, padj_threshold, color_up, color_down, color_none,
    color_na, gsea_name, layout, confounder_handling
  )
  # save_plot(filename, plot, width, height, ...): `...` -> ggsave
  save_plot(
    snakemake@output[["heatmap"]], heatmap_plot,
    width = layout$width, height = layout$height, limitsize = FALSE
  )
  cat("Heatmap saved:", snakemake@output[["heatmap"]], "\n")
} else {
  cat("No pathways to plot.\n")
  empty_heatmap <- ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5,
      label = paste0("No significant pathways found in any ", grouping_var),
      size = 5, color = "gray50"
    ) +
    theme_void()
  ggsave(
    snakemake@output[["heatmap"]], empty_heatmap,
    width = 8, height = 6, dpi = 300
  )
}

# --- Pairwise Jaccard overlap of significant-pathway sets --------------

cat("\nComputing pairwise Jaccard overlap...\n")
jaccard_df <- compute_pairwise_jaccard(sig_sets, group_names, grouping_var)
cat("Pairwise Jaccard summary (on significant-pathway sets):\n")
print(jaccard_df)
cat("\nSaving Jaccard table...\n")
jaccard_file <- snakemake@output[["jaccard"]]
dir.create(dirname(jaccard_file), showWarnings = FALSE, recursive = TRUE)
write_tsv(jaccard_df, jaccard_file)
cat("Saved:", jaccard_file, "\n")

cat("\n=== GSEA cross-group summary complete ===\n")
sink(type = "message")
sink()
close(logfile)
