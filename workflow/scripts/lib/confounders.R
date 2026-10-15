# SPDX-FileCopyrightText: © 2026 Marcel Schilling <mschilling@idibell.cat>
# SPDX-License-Identifier: AGPL-3.0-or-later


#########
# About #
#########

# File:     workflow/scripts/lib/confounders.R
# Author:   Marcel Schilling <mschilling@idibell.cat>
# Created:  2026-10-03
# Modified: 2026-10-15
# License:  GNU Affero General Public License version 3 or later
# Purpose:  Confounder field helpers R script for the bulk RNA-seq DGE
#           (DESeq2 + GSEA) Snakemake workflow.


###############
# Description #
###############

# Shared helper for scripts that put confounder fields into a design.
#
# Every coldata column is a string parsed from the sample label (and
# deseq2_fit.R turns each field into a factor via relevel()), while
# model.matrix() decides factor-vs-numeric from the R type. A confounder
# declared "continuous" must therefore be cast explicitly, or it gets one
# dummy column per distinct value.


###########
# License #
###########

# Confounder field helpers R script for bulk RNA-seq DGE (DESeq2 + GSEA)
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

# Returns `coldata` with the columns of all confounders declared
# "continuous" converted to numeric; every other column is left untouched.
# `coldata` can be a data.frame or a DataFrame (e.g. colData(dds)); the
# return value has the same class.
cast_continuous_confounders <- function(coldata, confounder_field_names,
                                        confounder_field_types) {
  for (f in confounder_field_names) {
    if (identical(confounder_field_types[[f]], "continuous")) {
      # as.character() first: as.numeric() on a factor returns the level
      # codes, not the values.
      values <- as.character(coldata[[f]])
      numeric_values <- suppressWarnings(as.numeric(values))
      if (anyNA(numeric_values)) {
        stop(
          "ERROR: confounder '", f, "' is declared continuous but has ",
          "non-numeric value(s): ",
          paste(unique(values[is.na(numeric_values)]), collapse = ", ")
        )
      }
      coldata[[f]] <- numeric_values
    }
  }
  coldata
}
