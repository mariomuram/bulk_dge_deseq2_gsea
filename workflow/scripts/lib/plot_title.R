# Plot titles that say how confounders were handled in the analysis behind
# the plot, and (for LFC plots) whether the LFCs are shrunk.
#
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
