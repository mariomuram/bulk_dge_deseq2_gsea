#!/usr/bin/env Rscript

# SPDX-FileCopyrightText: © 2026 Marcel Schilling <mschilling@idibell.cat>
# SPDX-License-Identifier: AGPL-3.0-or-later


#########
# About #
#########

# File:     workflow/scripts/plot_volcano.R
# Author:   Marcel Schilling <mschilling@idibell.cat>
# Created:  2026-09-21
# Modified: 2026-10-20
# License:  GNU Affero General Public License version 3 or later
# Purpose:  Volcano plot R script for the bulk RNA-seq DGE (DESeq2 + GSEA)
#           Snakemake workflow.


###############
# Description #
###############

# Plot: Volcano plot for a single contrast, unshrunk or shrunk LFC.


###########
# License #
###########

# Volcano plot R script for bulk RNA-seq DGE (DESeq2 + GSEA) Snakemake
# workflow.
# Copyright (C) 2026  Marcel Schilling
#
# This file is part of the bulk RNA-seq DGE (DESeq2 + GSEA) Snakemake
# workflow.
# The bulk RNA-seq DGE (DESeq2 + GSEA) Snakemake workflow is free software:
# you can redistribute it and/or modify it under the terms of the GNU
# Affero General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


###########
# Imports #
###########

library(ggplot2)
library(ggrepel)


###########
# Logging #
###########

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== Plot: Volcano Plot ===\n")


###############
# Shared code #
###############

source(snakemake@input[["plot_theme_lib"]])
source(snakemake@input[["plot_input_lib"]])
source(snakemake@input[["plot_title_lib"]])


##############
# Parameters #
##############

confounder_handling <- snakemake@wildcards[["confounder_handling"]]
log2fc_threshold <- snakemake@params[["log2fc_threshold"]]
plot_file <- snakemake@output[["volcano"]]

cat("Confounder handling:", confounder_handling, "\n")
cat("log2FC threshold:", log2fc_threshold, "\n\n")


#########
# Input #
#########

n_before <- nrow(results_df)
results_df <- results_df %>%
  filter(!is.na(padj)) %>%
  mutate(
    lfc_value = .data[[lfc_col]],
    neg_log10_padj = -log10(padj),
    # Overwrites the sign-only `category` from plot_input_lib with the
    # stricter (padj + |LFC|) classification volcano plots conventionally use.
    category = case_when(
      category %in% "Up" & lfc_value >= log2fc_threshold ~ "Up",
      category %in% "Down" & lfc_value <= -log2fc_threshold ~ "Down",
      TRUE ~ "Not Significant"
    ) %>% factor(levels = c("Up", "Down", "Not Significant"))
  )
cat(
  "Excluded", n_before - nrow(results_df),
  "gene(s) with NA padj (cannot plot -log10(NA)).\n"
)


################
# Volcano plot #
################

lfc_max <- max(abs(results_df$lfc_value), na.rm = TRUE)
volcano_plot <- ggplot(
  order_for_plotting(results_df),
  aes(x = lfc_value, y = neg_log10_padj, color = category, label = label)
) +
  geom_point(alpha = 0.6, size = default_point_size) +
  geom_vline(
    xintercept = c(-log2fc_threshold, log2fc_threshold),
    linetype = "dashed", color = "gray50", linewidth = default_line_width
  ) +
  geom_hline(
    yintercept = -log10(padj_threshold),
    linetype = "dashed", color = "gray50", linewidth = default_line_width
  ) +
  geom_label_repel(fill="white", alpha=.5,
    size = default_text_size, fontface = "bold", box.padding = .5,
    point.padding = 0, point.size = default_point_size, force = 3,
    force_pull = .8, max.iter = 1e6, max.time = 60,
    max.overlaps = Inf, min.segment.length = 0, segment.color = "grey40",
    segment.size = 0.4, seed = 42L, show.legend = FALSE
  ) +
  # +10% breathing room past the extremes
  coord_cartesian(xlim = c(-lfc_max, lfc_max) * 1.05) +
  scale_color_manual(values = category_colors,
                     labels = get_category_labels(results_df)) +
  labs(
    title = format_plot_title(contrast_name, confounder_handling, lfc_type),
    x = lfc_label,
    y = "-Log10(Adjusted p-value)",
    color = ""
  )


##########
# Output #
##########

dir.create(dirname(plot_file), showWarnings = FALSE, recursive = TRUE)
save_plot(plot_file, volcano_plot)
cat("Saved:", plot_file, "\n")


###########
# Cleanup #
###########

cat("\n=== Plot: Volcano Plot Complete ===\n")
sink(type = "message")
sink()
close(logfile)
