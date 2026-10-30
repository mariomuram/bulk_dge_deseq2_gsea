#!/usr/bin/env Rscript

# SPDX-FileCopyrightText: © 2026 Marcel Schilling <mschilling@idibell.cat>
# SPDX-License-Identifier: AGPL-3.0-or-later


#########
# About #
#########

# File:     workflow/scripts/plot_pca.R
# Author:   Marcel Schilling <mschilling@idibell.cat>
# Created:  2026-10-01
# Modified: 2026-10-29
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


#############
# Constants #
#############

# Default prefix for principal component columns in PCA scores data frame.
pc_prefix <- "PC"


###############
# Shared code #
###############

source(snakemake@input[["plot_theme_lib"]])


##############
# Parameters #
##############

pcs <- strsplit(snakemake@wildcards[["pcs"]], "_")[[1]]
pca_x <- as.integer(sub("^PC", "", pcs[1]))
pca_y <- as.integer(sub("^PC", "", pcs[2]))
dds <- snakemake@input[["dds_rds"]]
scores_file <- snakemake@input[["scores"]]
correction <- snakemake@wildcards[["correction"]]
colour_config <- snakemake@params[["colour_setting"]]
colour_var <- snakemake@params[["colour_field"]]
colour_label <- snakemake@params[["colour_label"]]
shape_var <- snakemake@params[["shape"]][["field"]]
shape_label <- snakemake@params[["shape"]][["label"]]
label_var <- snakemake@params[["label"]]
title_suffix <- snakemake@params[["title_suffix"]]
plot_file <- snakemake@output[["pca"]]
cat("Correction stage:", correction, "\n")
cat("Plotting with: color =", colour_var, ", shape =", shape_var, "\n\n")


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


################
# Checking PCs #
################

#Checking that PCs requested are within the range of PCs calculated.
n_pcs <- ncol(pca$x)
if (max(pca_x, pca_y) > n_pcs) {
  stop("Requested PCs (", pca_x, ", ", pca_y, ") exceed number of computed PCs (",
       n_pcs, ").")
} 


####################
# PCA scatter plot #
####################

fields <- unique(c(colour_var, shape_var, label_var))

pca_data <- cbind(
  setNames(
    as.data.frame(pca$x[, c(pca_x, pca_y), drop = FALSE]),
    paste0(pc_prefix, c(pca_x, pca_y))
  ),
  as.data.frame(colData(dds))[
    rownames(pca$x), fields, drop = FALSE])

###For the future this should be checked and changed
# Try to convert the selected colour variable to numeric
colour_values <- pca_data[[colour_var]]

numeric_values <- suppressWarnings(
  as.numeric(as.character(colour_values))
)

# Determine whether the variable can be treated as continuous
is_continuous <- (
  any(!is.na(colour_values)) &&
  all(!is.na(numeric_values) | is.na(colour_values))
)

# Convert the variable if it can be treated as continuous
if (is_continuous) {
  pca_data[[colour_var]] <- numeric_values
}

pca_plot <- ggplot(pca_data,
                   aes(x = .data[[paste0("PC", pca_x)]], y = .data[[paste0("PC", pca_y)]],
                   color = .data[[colour_var]], 
                   shape = .data[[shape_var]])) +
  geom_point(size = default_point_size, alpha = 0.8) +
  labs(
    title = paste0("PCA", title_suffix),
    x = paste0(pc_prefix, pca_x, ": ", percent_var[pca_x], "% variance"),
    y = paste0(pc_prefix, pca_y, ": ", percent_var[pca_y], "% variance"),
    color = colour_label, shape = shape_label
  )

# Select the color scale depending on the variable type
if (is_continuous) {
  if ("palette" %in% names(colour_config)) {
    pca_plot <- pca_plot + scale_color_distiller(palette = colour_config[["palette"]])
  } else {
    if ("mid" %in% names(colour_config)) {
      pca_plot <- pca_plot + scale_color_gradient2(
        low = colour_config[["low"]],
        mid = colour_config[["mid"]],
        high = colour_config[["high"]],
        midpoint = median(pca_data[[colour_var]], na.rm = TRUE)
      )
    } else {
      pca_plot <- pca_plot + scale_color_gradient(
        low = colour_config[["low"]],
        high = colour_config[["high"]]
      )
    }
  }
} else {
  if ("palette" %in% names(colour_config)) {
    pca_plot <- pca_plot +
      scale_color_brewer(palette = colour_config[["palette"]])
  } else {
    pca_plot <- pca_plot +
      scale_color_brewer(palette = "Dark2")  # Default palette for categorical variables
  }
}

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
