#!/usr/bin/env Rscript
# QC: field x field collinearity matrix (not correction-specific: metadata
# doesn't change across correction stages).
#
# Covers every field in params$field_types (condition + confounder +
# unclassified), minus constant/identifier types, so an undeclared field
# aliased with a declared one is still surfaced.
#
# HOW TO USE IT
#   A high association between two fields (e.g. Cramer's V near 1 for two
#   categorical fields) means correcting for either should look similar in the
#   after-correction PC-field matrix, so the two-experiment comparison (try
#   confounder = A, then confounder = B, compare the pc_field_association_limma
#   heatmaps) is expected to show little difference. A moderate value means
#   the comparison is more likely to show a real, actionable difference.
#   The two diagnostics are complements: this matrix tells you what to
#   EXPECT from the comparison; the comparison itself tells you what actually
#   happens for gene expression.
#
# CAVEATS
#   * Fill is |estimate|. Signed correlations (Spearman, point-biserial, phi,
#     Kendall) and unsigned associations (Cramer's V) are not on a common
#     scale; the `method` column of the table says which was used per pair.
#   * Cramer's V is biased upward for sparse tables (few samples, many
#     levels): two categorical fields can show a high V by chance at small n.

library(DESeq2)
library(ggplot2)
library(smartcor)
library(readr)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== QC: Field Collinearity Matrix ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

source(snakemake@input[["plot_theme_lib"]])
source(snakemake@input[["association_lib"]])

condition_field_names <- unlist(snakemake@params[["condition_fields"]])
confounder_field_names <- unlist(snakemake@params[["confounder_fields"]])

compute_field_collinearity_matrix <- function(coldata, field_types) {
  prepared <- prepare_association_data(coldata, field_types)
  compute_associations(prepared$data, prepared$types)
}

dds <- readRDS(snakemake@input[["dds_rds"]])
coldata <- as.data.frame(colData(dds))
field_types <- diagnostic_field_types(coldata, snakemake@params[["field_types"]])
cat("Fields compared:", paste(names(field_types), collapse = ", "), "\n")
if (length(field_types) < 2) {
  stop("ERROR: fewer than two fields left after excluding constant/identifier ",
       "types; nothing to compare.")
}

collinearity <- compute_field_collinearity_matrix(coldata, field_types)

n_valid <- sum(!is.na(collinearity$estimate))
cat("Valid pairs:", n_valid, "of", nrow(collinearity), "\n")
if (n_valid == 0) {
  # TODO: replace with a placeholder output + warning() once a dummy-output
  # strategy exists (Snakemake needs a real output file, not a skipped one).
  stop(
    "ERROR: no valid field-field associations were computed; see the ",
    "smartcor warnings earlier in this log."
  )
}

table_file <- snakemake@output[["table"]]
dir.create(dirname(table_file), showWarnings = FALSE, recursive = TRUE)
readr::write_tsv(collinearity, table_file)
cat("Saved:", table_file, "\n")

# Order: condition, confounder, other.
role_rank <- ifelse(names(field_types) %in% condition_field_names, 1,
             ifelse(names(field_types) %in% confounder_field_names, 2, 3))
field_order <- names(field_types)[order(role_rank, names(field_types))]

# Mirror the upper triangle for readability; diagonal = 1.
heat_data <- rbind(
  data.frame(a = collinearity$var_x, b = collinearity$var_y,
             estimate = abs(collinearity$estimate), p = collinearity$p_value),
  data.frame(a = collinearity$var_y, b = collinearity$var_x,
             estimate = abs(collinearity$estimate), p = collinearity$p_value),
  data.frame(a = field_order, b = field_order, estimate = 1, p = NA_real_)
)
heat_data$a <- factor(heat_data$a, levels = field_order)
heat_data$b <- factor(heat_data$b, levels = rev(field_order))
heat_data$label <- ifelse(
  is.na(heat_data$estimate), "NA",
  sprintf("%.2f%s", heat_data$estimate,
          ifelse(!is.na(heat_data$p) & heat_data$p < 0.05, "*", ""))
)

heatmap <- ggplot(heat_data, aes(x = a, y = b, fill = estimate)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = label), size = default_text_size) +
  scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, 1),
                      na.value = "grey85", name = "|association|") +
  labs(x = NULL, y = NULL,
       caption = paste("Signed correlations shown as |r|; Cramer's V is",
                       "unsigned. * = raw p < 0.05.")) +
  theme(panel.grid = element_blank(), axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

output_file <- snakemake@output[["heatmap"]]
# Sized for the 20 pt projector theme (~1 in per field).
n_fields <- length(field_order)
save_plot(output_file, heatmap, width = 4 + n_fields, height = 3 + n_fields)
cat("Saved:", output_file, "\n")

cat("\n=== QC: Field Collinearity Matrix Complete ===\n")
sink(type = "message")
sink()
close(logfile)
