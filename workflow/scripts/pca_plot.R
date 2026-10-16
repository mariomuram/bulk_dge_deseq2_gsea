#!/usr/bin/env Rscript

# SPDX-FileCopyrightText: © 2026 Marcel Schilling <mschilling@idibell.cat>
# SPDX-License-Identifier: AGPL-3.0-or-later


#########
# About #
#########

# File:     workflow/scripts/pca_plot.R
# Author:   Marcel Schilling <mschilling@idibell.cat>
# Created:  2026-10-01
# Modified: 2026-10-21
# License:  GNU Affero General Public License version 3 or later
# Purpose:  PCA plot R script for the bulk RNA-seq DGE (DESeq2 + GSEA)
#           Snakemake workflow.


###############
# Description #
###############

# QC: PCA plot across all samples (not contrast-specific), for one
# correction stage (`"uncorrected"` = blind VST; otherwise the
# batch-corrected VST). 


###########
# License #
###########

# PCA plot R script for bulk RNA-seq DGE (DESeq2 + GSEA) Snakemake
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

library(DESeq2)
library(ggplot2)


###########
# Logging #
###########

logfile <- snakemake@log[[1]]
dir.create(dirname(logfile), showWarnings = FALSE, recursive = TRUE)
logfile <- file(logfile, open = "w")
sink(logfile, type = "output")
sink(logfile, type = "message", append = TRUE)

cat("=== QC: PCA Plot ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")


###############
# Shared code #
###############

source(snakemake@input[["plot_theme_lib"]])


##############
# Parameters #
##############

dds <- snakemake@input[["dds_rds"]]
scores_file <- snakemake@input[["scores"]]
correction <- snakemake@wildcards[["correction"]]
color_var <- snakemake@params[["color"]]
color_label <- unname(color_var[[1]])
color_var <- names(color_var)[1]
shape_var <- snakemake@params[["shape"]]
shape_label <- unname(shape_var[[1]])
shape_var <- names(shape_var)[1]
label_var <- snakemake@params[["label"]]
palette <- snakemake@params[["palette"]]
title_suffix <- snakemake@params[["title_suffix"]]
plot_file <- snakemake@output[["pca"]]
pca_x <- snakemake@params[["pca_x"]]
pca_y <- snakemake@params[["pca_y"]]
cat("Correction stage:", correction, "\n")
cat("Plotting with: color =", color_var, ", shape =", shape_var, "\n\n")


#########
# Input #
#########

pca <- readRDS(scores_file)
dds <- readRDS(dds)
cat("Loaded dds:", ncol(dds), "samples,", nrow(dds), "genes\n")


######################
# Variance explained #
######################

percent_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2))


####################
# PCA scatter plot #
####################

# TODO: Allow configuring which PCs are plotted (default: PC1/PC2).
# Consider splitting PCA computation (above) into separate script to allow
# plotting severale PC pairs from the same PCA without re-computing the PCA
# by calling the PCA script once and this plotting script several times
# with different parameters from Snakemake.
fields <- unique(c(color_var, shape_var, label_var))

pca_data <- cbind(
  setNames(
    as.data.frame(pca$x[, c(pca_x, pca_y), drop=FALSE]),
    paste0("PC", c(pca_x, pca_y))
  ),
  as.data.frame(colData(dds))[
    rownames(pca$x), fields, drop = FALSE])

pca_plot <- ggplot(pca_data,
                   aes(x = .data[[paste0("PC", pca_x)]], y = .data[[paste0("PC", pca_y)]],
                   color = .data[[color_var]], 
                   shape = .data[[shape_var]])) +
  geom_point(size = default_point_size, alpha = 0.8) +
  scale_color_brewer(palette = palette) +
  labs(
    title = paste0("PCA", title_suffix),
    x = paste0("PC", pca_x, ": ", percent_var[pca_x], "% variance"),
    y = paste0("PC", pca_y, ": ", percent_var[pca_y], "% variance"),
    color = color_label, shape = shape_label
  )
if (!is.null(label_var)) {
  library(ggrepel)
  pca_plot <- pca_plot +
    geom_label_repel(
      aes(label = .data[[label_var]]), fill = "white", alpha = .5,
      size = default_text_size, fontface = "bold", box.padding = .5,
      point.padding = 0, point.size = default_point_size, force = 3,
      force_pull = .8, max.iter = 1e6, max.time = 60, max.overlaps = Inf,
      min.segment.length = 0, segment.color = "grey40", segment.size = 0.4,
      seed = 42L, show.legend = FALSE
    )
}


##########
# Output #
##########

dir.create(dirname(plot_file), showWarnings = FALSE, recursive = TRUE)
save_plot(plot_file, pca_plot)
cat("Saved:", plot_file, "\n")


###########
# Cleanup #
###########

cat("\n=== QC: PCA Plot Complete ===\n")
sink(type = "message")
sink()
close(logfile)
