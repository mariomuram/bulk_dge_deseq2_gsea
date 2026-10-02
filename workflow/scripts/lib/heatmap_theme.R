# Shared heatmap helpers for the GSEA summary (gsea_summary.R) and the
# Phase 1 PC/field heatmaps (pc_field_association.R, field_collinearity.R).
# Layer theme_heatmap_no_grid() on top of theme_projector() (from
# plot_theme_lib).

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
