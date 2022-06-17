#' Area proportional venn diagrams
#'
#' This functions uses eulerr::euler to plot area proportional venn diagramms
#' but plots it using ggplot2
#'
#' @param combinations set relationships as a named numeric vector, matrix, or data.frame(See `eulerr::euler`)
#' @param quantities whether to show number of intersecting elements
#' @param ... further arguments passed to eulerr::euler
ggeulerr <- function(combinations, show_quantities = TRUE, show_labels = TRUE, ...) {
  data <-
    eulerr::euler(combinations = combinations, ...) %>%
    plot(quantities = show_quantities) %>%
    pluck("data")
  
  tibble() %>%
    ggplot() +
    ggforce::geom_ellipse(
      data = data$ellipses %>% as_tibble(rownames = "Set"),
      mapping = aes(x0 = h, y0 = k, a = a, b = b, angle = 0, fill = Set),
      alpha = 0.5
    ) +
    geom_text(
      data = {
        data$centers %>%
          mutate(
            label = labels %>% map2(quantities, ~ {
              if (!is.na(.x) && !is.na(.y) && show_labels) {
                paste0(.x, "\n", sprintf(.y, fmt = "%.2g"))
              } else if (!is.na(.x) && show_labels) {
                .x
              } else if (!is.na(.y)) {
                .y
              } else {
                ""
              }
            })
          )
      },
      mapping = aes(x = x, y = y, label = label),
      size = 2.5
    ) +
    theme_void() +
    coord_fixed() +
    guides(fill = guide_legend(override.aes = list(alpha = 1))) +
    scale_fill_hue()
}
