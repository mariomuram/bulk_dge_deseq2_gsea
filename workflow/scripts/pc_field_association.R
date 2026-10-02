#!/usr/bin/env Rscript
# QC: PC-field association matrix, for one correction stage.
#
# One cell = one PC vs. one field, tested univariately via smartcor (method
# selection: see lib/association.R).
#
# Heatmap fill is -log10(p), the only quantity comparable across the
# different association measures smartcor selects; estimate and method are
# in the table.
# Asterisks mark raw p < 0.05 (no multiple-testing correction). A colored
# strip along the field axis shows each field's role in the design.

library(DESeq2)
library(ggplot2)
library(patchwork)
library(smartcor)
library(readr)

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== QC: PC-Field Association Matrix ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

source(snakemake@input[["plot_theme_lib"]])
source(snakemake@input[["association_lib"]])
source(snakemake@input[["heatmap_theme_lib"]])

correction <- snakemake@wildcards[["correction"]]
condition_field_names <- unlist(snakemake@params[["condition_fields"]])
confounder_field_names <- unlist(snakemake@params[["confounder_fields"]])
n_pcs_max <- 10  # PCs shown; the rest explain little and clutter the heatmap
cat("Correction stage:", correction, "\n")

# One smart_cormat() call over [PC scores | fields]; only the PC-vs-field pairs
# are kept (the PC-PC and field-field pairs are not needed here).
compute_pc_field_associations <- function(scores, coldata, field_types) {
  pc_names <- colnames(scores)
  if (any(pc_names %in% names(field_types))) {
    stop("ERROR: a field is named like a PC (PC1, PC2, ...); rename it.")
  }
  prepared <- prepare_association_data(coldata, field_types)
  data <- cbind(as.data.frame(scores), prepared$data)
  types <- c(setNames(rep(list("continuous"), length(pc_names)), pc_names),
             prepared$types)

  pairs <- compute_associations(data, types)
  x_is_pc <- pairs$var_x %in% pc_names
  y_is_pc <- pairs$var_y %in% pc_names
  pairs <- pairs[xor(x_is_pc, y_is_pc), , drop = FALSE]
  x_is_pc <- pairs$var_x %in% pc_names

  field <- ifelse(x_is_pc, pairs$var_y, pairs$var_x)
  data.frame(
    pc = ifelse(x_is_pc, pairs$var_x, pairs$var_y),
    field = field,
    field_type = unname(field_types[field]),
    pairs[c("estimate", "p_value", "method", "method_label", "n_complete")],
    stringsAsFactors = FALSE
  )
}

plot_pc_field_association <- function(assoc, var_explained, roles) {
  pc_names <- names(var_explained)
  role_order <- c("condition", "confounder", "unused")
  field_order <- names(roles)[order(match(roles, role_order), names(roles))]
  field_levels <- rev(field_order)  # first field at the top

  assoc$field <- factor(assoc$field, levels = field_levels)
  assoc$pc <- factor(assoc$pc, levels = pc_names)
  assoc$neglogp <- pmin(-log10(assoc$p_value), 4)
  assoc$star <- significance_asterisks(assoc$p_value)

  main <- ggplot(assoc, aes(x = pc, y = field, fill = neglogp)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = star), size = default_text_size) +
    scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, 4),
                        na.value = "grey85",
                        name = expression(-log[10](p)~"(cap 4)")) +
    labs(x = NULL, y = NULL) +
    theme_heatmap_no_grid() + theme(axis.text = element_blank())

  role_data <- data.frame(
    field = factor(names(roles), levels = field_levels),
    role = factor(unname(roles), levels = role_order)
  )
  role_strip <- ggplot(role_data, aes(x = 1, y = field, fill = role)) +
    geom_tile(colour = "white") +
    scale_fill_manual(
      values = c(condition = "#D55E00", confounder = "#0072B2",
                 unused = "#BBBBBB"),
      drop = FALSE, name = "Field role"
    ) +
    labs(x = NULL, y = NULL) +
    theme_heatmap_no_grid() + theme(axis.text.x = element_blank())

  var_data <- data.frame(
    pc = factor(pc_names, levels = pc_names),
    percent = 100 * var_explained,
    label = sprintf("%s\n%.1f%%", pc_names, 100 * var_explained)
  )
  var_strip <- ggplot(var_data, aes(x = pc, y = 1, fill = percent)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = label), size = default_text_size) +
    scale_fill_gradient(low = "white", high = "steelblue",
                        name = "% variance") +
    labs(x = NULL, y = NULL) +
    theme_heatmap_no_grid() + theme(axis.text = element_blank())

  wrap_plots(R = role_strip, H = main, V = var_strip, design = "RH\n#V") +
    plot_layout(widths = c(1, length(pc_names)),
                heights = c(length(field_order), 2), guides = "collect")
}

dds <- readRDS(snakemake@input[["dds_rds"]])
coldata <- as.data.frame(colData(dds))
pca <- readRDS(snakemake@input[["pca_scores"]])

field_types <- diagnostic_field_types(coldata, snakemake@params[["field_types"]])
cat("Fields tested:", paste(names(field_types), collapse = ", "), "\n")

n_pcs <- min(n_pcs_max, ncol(pca$x))
scores <- pca$x[, seq_len(n_pcs), drop = FALSE]
stopifnot(setequal(rownames(scores), rownames(coldata)))
coldata <- coldata[rownames(scores), , drop = FALSE]
var_explained <- (pca$sdev^2 / sum(pca$sdev^2))[seq_len(n_pcs)]
names(var_explained) <- colnames(scores)

assoc <- compute_pc_field_associations(scores, coldata, field_types)

n_valid <- sum(!is.na(assoc$p_value))
cat("Valid cells:", n_valid, "of", nrow(assoc), "\n")
if (n_valid == 0) {
  # TODO: replace with a placeholder output + warning() once a dummy-output
  # strategy exists (Snakemake needs a real output file, not a skipped one).
  stop(
    "ERROR: no valid PC-field associations were computed; see the smartcor ",
    "warnings earlier in this log."
  )
}

roles <- ifelse(names(field_types) %in% condition_field_names, "condition",
         ifelse(names(field_types) %in% confounder_field_names, "confounder",
                "unused"))
names(roles) <- names(field_types)

table_file <- snakemake@output[["table"]]
dir.create(dirname(table_file), showWarnings = FALSE, recursive = TRUE)
write_tsv(assoc, table_file)
cat("Saved:", table_file, "\n")

output_file <- snakemake@output[["heatmap"]]
# Sized for the 20 pt projector theme: ~1 in per PC column, ~0.6 in per field
# row, plus room for the field labels, legends and the % variance strip.
save_plot(
  output_file, plot_pc_field_association(assoc, var_explained, roles),
  width = 5 + n_pcs, height = 3 + 0.6 * length(field_types)
)
cat("Saved:", output_file, "\n")

cat("\n=== QC: PC-Field Association Matrix Complete ===\n")
sink(type = "message")
sink()
close(logfile)
