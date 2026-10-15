# SPDX-FileCopyrightText: © 2026 Marcel Schilling <mschilling@idibell.cat>
# SPDX-License-Identifier: AGPL-3.0-or-later


#########
# About #
#########

# File:     workflow/scripts/lib/plot_theme.R
# Author:   Marcel Schilling <mschilling@idibell.cat>
# Created:  2026-09-19
# Modified: 2026-10-15
# License:  GNU Affero General Public License version 3 or later
# Purpose:  Plot theme helper R script for the bulk RNA-seq DGE (DESeq2 +
#           GSEA) Snakemake workflow.


###############
# Description #
###############

# Shared plot theme helpers to unifiy plot styling in a central place.


###########
# License #
###########

# Plot theme helper R script for bulk RNA-seq DGE (DESeq2 + GSEA) Snakemake
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


#############
# Constants #
#############

base_size <- 20
default_text_size <- (base_size * 0.75) / .pt
default_point_size <- 4
default_line_width <- 1.2


#############
# Functions #
#############

# TODO: Parameterize / default to constants!
# Global theme for projector presentations
theme_projector <- function(base_size = base_size, base_family = "sans") {
  theme_minimal(
    base_size = base_size,
    base_family = base_family,
    base_line_size = 0.6,
    base_rect_size = 0.6
  ) +
  theme(
    plot.title = element_text(face = "bold", size = rel(1.3), margin = margin(b = 10)),
    plot.subtitle = element_text(size = rel(1.1), margin = margin(b = 8), color = "#666666"),
    axis.line = element_line(color = "#333333", linewidth = .75 * default_line_width),
    axis.title = element_text(size = rel(1.1), face = "bold"),
    axis.text = element_text(size = rel(1.0), color = "#333333"),  # Inherit 20 pt
    legend.text = element_text(size = rel(0.95)),
    legend.title = element_text(size = rel(1.05), face = "bold"),
    panel.grid.major = element_line(color = "#6e6e6e", linewidth = .35 * default_line_width),
    panel.grid.minor = element_blank()
  )
}

# Export function
save_plot <- function(filename, plot, width = 8, height = 6, ...) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    ...,
    device = function(filename, ...) {
      pdf(file = filename, useDingbats = FALSE, ...)
    }
  )
}


#########
# Setup #
#########

# Set global theme
theme_set(theme_projector(base_size))
