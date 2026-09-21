# Shared helper: validate the configured ranking statistic and
# build a de-duplicated, sorted ranking vector from a DESeq2
# results table. Sourced (not library()'d) by gsea_run_go.R and
# gsea_run_hallmark.R -- kept identical across both since ranking
# logic doesn't depend on which database the enrichment step uses.
library(dplyr)

validate_and_rank <- function(results_df, ranking_stat,
                              gene_id_style="SYMBOL") {
  stat_col <- ranking_stat
  cat("Ranking statistic:", stat_col, "\n")
  if (!(stat_col %in% colnames(results_df))) {
    stop(
      "ERROR: gsea.ranking_stat = '", stat_col, "' is not a ",
      "column in the DESeq2 results table. Available columns: ",
      paste(colnames(results_df), collapse = ", ")
    )
  }

  # Filter on the ranking column itself, not padj: independent
  # filtering leaves padj NA but stat/log2FoldChange intact, and
  # preranked GSEA should use the full ranked transcriptome, not
  # only genes that cleared the DE significance test.
  n_before <- nrow(results_df)
  results_df <- results_df %>%
    filter(
      !is.na(.data[[stat_col]]),
      !is.infinite(.data[[stat_col]])
    )
  cat(
    "Dropped", n_before - nrow(results_df),
    "gene(s) with NA/infinite", stat_col, "\n"
  )

  if (any(duplicated(results_df$gene_id))) {
    n_dup <- sum(duplicated(results_df$gene_id))
    cat(
      "WARNING:", n_dup, "duplicate gene_id(s) found; keeping",
      "the entry with the largest |", stat_col, "| for each.",
      "fgsea requires unique gene names.\n"
    )
    results_df <- results_df %>%
      mutate(.abs_stat = abs(.data[[stat_col]])) %>%
      arrange(desc(.abs_stat)) %>%
      distinct(gene_id, .keep_all = TRUE) %>%
      select(-.abs_stat)
  }

  ranking <- results_df[[stat_col]]
  names(ranking) <- results_df$gene_id
  if (gene_id_style == "ENSEMBL") {
    # Clean version suffixes
    names(ranking) <- sub("\\.\\d+$", "", names(ranking))
  }
  ranking <- sort(ranking, decreasing = TRUE)

  cat("\nGenes ranked:", length(ranking), "\n")
  cat("Top 5 genes:\n")
  print(head(ranking, 5))
  cat("Bottom 5 genes:\n")
  print(tail(ranking, 5))

  ranking
}
