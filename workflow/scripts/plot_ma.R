#!/usr/bin/env Rscript

# SPDX-FileCopyrightText: © 2026 Marcel Schilling <mschilling@idibell.cat>
# SPDX-License-Identifier: AGPL-3.0-or-later


#########
# About #
#########

# File:     workflow/scripts/plot_ma.R
# Author:   Marcel Schilling <mschilling@idibell.cat>
# Created:  2026-09-21
# Modified: 2026-10-15
# License:  GNU Affero General Public License version 3 or later
# Purpose:  MA plot R script for the bulk RNA-seq DGE (DESeq2 + GSEA)
#           Snakemake workflow.


###############
# Description #
###############

# Plot: MA plot for a single contrast, unshrunk or shrunk LFC.


###########
# License #
###########

# MA plot R script for bulk RNA-seq DGE (DESeq2 + GSEA) Snakemake workflow.
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
library(dplyr)
library(readr)
library(ggrepel)


###########
# Logging #
###########

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== Plot: MA Plot ===\n")


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
output_file <- snakemake@output[["ma"]]

cat("Confounder handling:", confounder_handling, "\n")


#########
# Input #
#########

n_na <- sum(is.na(results_df$padj))
if (n_na > 0) {
  cat(
    "NOTE:", n_na, "gene(s) have NA padj (independent",
    "filtering and/or Cook's outliers); shown as",
    "'Not Significant' rather than dropped, so the low-count",
    "filtering funnel stays visible on the plot.\n"
  )
}


###########
# MA plot #
###########

lfc_max <- max(abs(results_df$lfc_value), na.rm = TRUE)
ma_plot <- ggplot(
  order_for_plotting(results_df),
  aes(x = baseMean, y = lfc_value, color = category, label = label)
) +
  geom_point(alpha = 0.5, size = default_point_size) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_label_repel(fill="white", alpha=.5,
    size = default_text_size, fontface = "bold", box.padding = .5,
    point.padding = 0, point.size = default_point_size, force = 3,
    force_pull = .8, max.iter = 1e6, max.time = 60,
    max.overlaps = Inf, min.segment.length = 0, segment.color = "grey40",
    segment.size = 0.4, seed = 42L, show.legend = FALSE
  ) +
  scale_x_log10() +
  # +10% breathing room past the extremes
  coord_cartesian(ylim = c(-lfc_max, lfc_max) * 1.05) +
  scale_color_manual(values = category_colors,
                     labels = get_category_labels(results_df)) +
  labs(
    title = format_plot_title(contrast_name, confounder_handling, lfc_type),
    x = "Mean Expression (baseMean, log10 scale)",
    y = lfc_label,
    color = paste0("padj < ", padj_threshold)
  )


##########
# Output #
##########

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
save_plot(output_file, ma_plot)
cat("Saved:", output_file, "\n")


###########
# Cleanup #
###########

cat("\n=== Plot: MA Plot Complete ===\n")
sink(type = "message")
sink()
close(logfile)
