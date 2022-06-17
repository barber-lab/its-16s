#' Facet samples and taxa grouping into patchwork
#'
#' @param data tibble with columns taxa_group, samples_group and plt
facet_groupings_plots <- function(data, title) {
  plot_title <- function(label, size = 9, ...) {
    ggplot() +
      annotate(geom = "text", x = 1, y = 1, label = label, size = size, ...) +
      theme_void()
  }

  rows <- data$samples_group %>% unique()
  cols <- data$taxa_group %>% unique()

  nrow <- length(rows)
  ncol <- length(cols)

  res <-
    data %>%
    group_by(samples_group) %>%
    nest() %>%
    mutate(plts = samples_group %>% map2(data, ~ {
      c(
        list(plot_title(.x, angle = 90)),
        .y$plt
      )
    }))

  c(
    list(plot_title("")),
    data$taxa_group %>% unique() %>% map(plot_title),
    res %>% unnest(plts) %>% pull(plts)
  ) %>%
    wrap_plots(ncol = ncol + 1, nrow = nrow + 1) +
    plot_layout(
      guides = "collect",
      widths = c(0.15, rep(1, nrow)),
      heights = c(0.15, rep(1, ncol))
    )
}
