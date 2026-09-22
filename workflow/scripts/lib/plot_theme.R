library(ggplot2)

base_size <- 20
default_text_size <- (base_size * 0.75) / .pt
default_point_size <- 4
default_line_width <- 1.2

# Global theme for projector presentations
theme_projector <- function(base_size = 20, base_family = "sans") {
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
save_plot <- function(filename, plot, width = 8, height = 6) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = function(filename, ...) {
      pdf(file = filename, useDingbats = FALSE, ...)
    }
  )
}


# Set global theme
theme_set(theme_projector(base_size))
