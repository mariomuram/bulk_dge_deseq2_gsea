# SPDX-FileCopyrightText: © 2026 Marcel Schilling <mschilling@idibell.cat>
# SPDX-License-Identifier: AGPL-3.0-or-later


#########
# About #
#########

# File:     workflow/scripts/lib/heatmap_theme.R
# Author:   Marcel Schilling <mschilling@idibell.cat>
# Created:  2026-10-02
# Modified: 2026-10-15
# License:  GNU Affero General Public License version 3 or later
# Purpose:  Heatmap theme helpers R script for the bulk RNA-seq DGE (DESeq2
#           + GSEA) Snakemake workflow.


###############
# Description #
###############

# Shared heatmap helpers for the GSEA summary (gsea_summary.R) and the
# Phase 1 PC/field heatmaps (pc_field_association.R, field_collinearity.R).
# Layer theme_heatmap_no_grid() on top of theme_projector() (from
# plot_theme_lib).


###########
# License #
###########

# Heatmap theme helpers R script for bulk RNA-seq DGE (DESeq2 + GSEA)
# Snakemake workflow.
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


#############
# Functions #
#############

# "*" where p_value < threshold, "" otherwise. NA p-values get "" (not NA):
# an NA padj means "not tested", and must neither draw an asterisk nor be
# silently dropped by geom_text as a missing label.
significance_asterisks <- function(p_value, threshold = 0.05) {
  ifelse(!is.na(p_value) & p_value < threshold, "*", "")
}

theme_heatmap_no_grid <- function() {
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank()
  )
}
