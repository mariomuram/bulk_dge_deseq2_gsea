# SPDX-FileCopyrightText: © 2026 Marcel Schilling <mschilling@idibell.cat>
# SPDX-License-Identifier: AGPL-3.0-or-later


#########
# About #
#########

# File:     workflow/scripts/lib/plot_title.R
# Author:   Marcel Schilling <mschilling@idibell.cat>
# Created:  2026-10-03
# Modified: 2026-10-15
# License:  GNU Affero General Public License version 3 or later
# Purpose:  Plot title helper R script for the bulk RNA-seq DGE (DESeq2 +
#           GSEA) Snakemake workflow.


###############
# Description #
###############

# Shared plot title helpers to unifiy plot labeling in a central place.

# Plot titles that say how confounders were handled in the analysis behind
# the plot, and (for LFC plots) whether the LFCs are shrunk.


###########
# License #
###########

# Plot title helper R script for bulk RNA-seq DGE (DESeq2 + GSEA) Snakemake
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


#############
# Functions #
#############

# `confounder_handling` is the value of the Snakefile's `{confounder_handling}`
# wildcard, i.e. snakemake@wildcards[["confounder_handling"]]:
#   "no_confounders"        nothing to say; omitted from the title
#   NULL                    a rule without that wildcard (e.g. a cross-
#                           strategy summary); omitted as well
#   "confounders_<how>"     labelled "conf. <how>", e.g. "conf. modelled"
# Any other value is an error rather than a silently misleading title: if a
# new handling does not follow the `confounders_<how>` convention, teach this
# function about it.
#
# `shrinkage` is the `lfc_type` ("shrunk" / "unshrunk") or NULL for plots
# without LFCs. Only "shrunk" is spelled out; unshrunk is the default.
#
# The parts are joined with "; " inside one pair of parentheses:
#   format_plot_title("X", "confounders_modelled", "shrunk")   "X (conf. modelled; shrunk)"
#   format_plot_title("X", "confounders_modelled", "unshrunk") "X (conf. modelled)"
#   format_plot_title("X", "no_confounders", "shrunk")         "X (shrunk)"
#   format_plot_title("X", "no_confounders", "unshrunk")       "X"
format_plot_title <- function(base, confounder_handling = NULL,
                              shrinkage = NULL) {
  parts <- character()

  if (!is.null(confounder_handling) &&
        !identical(confounder_handling, "no_confounders")) {
    if (!startsWith(confounder_handling, "confounders_")) {
      stop(
        "ERROR: cannot label confounder handling '", confounder_handling,
        "' in a plot title. If it is a new value of the Snakefile's ",
        "`{confounder_handling}` wildcard, extend format_plot_title() in ",
        "scripts/lib/plot_title.R."
      )
    }
    parts <- c(parts, sub("^confounders_", "conf. ", confounder_handling))
  }

  if (!is.null(shrinkage)) {
    if (!shrinkage %in% c("shrunk", "unshrunk")) {
      stop("ERROR: unknown LFC type '", shrinkage, "' in a plot title.")
    }
    if (shrinkage == "shrunk") parts <- c(parts, "shrunk")
  }

  if (length(parts) == 0) return(base)
  paste0(base, " (", paste(parts, collapse = "; "), ")")
}
