#' Wrap plots on top of each other to a stack with interleaved labels
#' @param data data frame with columns provided by `plot_col` and `label_col`
#' @param plot_col name of column in `data` containing plots
#' @param label_col name of character column in `data` containing labels for the plots
wrap_plots_vertically <- function(data, plot_col, label_col) {
  data %>%
    {
      map2(
        .x = .[[label_col]] %>% map(~ ggpubr::text_grob(.x, size = 30)),
        .y = .[[plot_col]],
        .f = ~ list(.x, .y)
      ) %>% purrr::flatten()
    } %>%
    wrap_plots(ncol = 1, heights = {
      c(0.05, 1) %>% rep(length(.))
    })
}
