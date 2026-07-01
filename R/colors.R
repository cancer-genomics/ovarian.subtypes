#' Use consistent color scheme for cancer subtypes
#'
#' @export
tumor_colors <- function() {
  dx.colors <- c(
    "Uterine endometrioid" = "#DDCC7F",
    "Uterine endometrial" = "#DDCC7F",
    "Ovarian endometrioid" = "#0F7554",
    "Ovarian mucinous" = "#44AA99",
    "Colorectal mucinous" = "#882255",
    "Pancreatic mucinous" = "#AA4499",
    "Stomach mucinous" = "#D695D0"
  )
}

tile_theme <- function() {
  theme(
    axis.title = element_blank(),
    strip.placement = "outside",
    panel.grid = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white"),
    axis.line = element_blank(),
    strip.text = element_text(
      size = 9,
      hjust = 0.5,
      vjust = 0.5
    ),
    strip.text.y.left = element_text(
      angle = 0,
      size = 11
    ),
    strip.background = element_rect(
      fill = "grey88",
      color = "grey88"
    ),
    panel.spacing.x = unit(1, "lines")
  )
}
