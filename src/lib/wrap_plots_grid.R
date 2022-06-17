#' Arrange ggplots into a grid and adds titles
#' Like facet_grid but for whole plot objects
#' 
#' @param data data frame containing columns mentioned in `formula` and `plot_column`
#' @param formula
wrap_plots_grid <- function(data, formula, plot_column = NULL, title_size = 13, label_size = 10) {
  if(is.null(plot_column)) stop("Must provide plot_column")
  
  cols_name <-  formula %>% as.list() %>% pluck(2) 
  cols <- data[[cols_name]] %>% unique() %>% sort()
  
  rows_name <- formula %>% as.list() %>% pluck(3)
  rows <- data[[rows_name]] %>% unique() %>% sort()
  
  empty_plt <- ggplot() + theme_void()
  
  data <-
    data %>%
    ungroup() %>%
    complete(
      !!sym(cols_name), !!sym(rows_name),
      fill = list() %>% inset2(plot_column, list(empty_plt))
    )
  
  cols_labels_plt <-
    cols %>%
    map(~ text_grob(.x, size = label_size)) %>%
    wrap_plots(nrow = 1) %>%
    list(text_grob(cols_name, size = title_size), .) %>%
    wrap_plots(ncol = 1)
  
  rows_labels_plt <-
    rows %>%
    map(~ text_grob(.x, rot = 90, size = label_size)) %>%
    wrap_plots(ncol = 1) %>%
    list(text_grob(rows_name, rot = 90, size = title_size), .) %>%
    patchwork::wrap_plots(nrow = 1)
  
  plots_plt <-
    data %>%
    arrange(!!sym(rows_name), !!sym(cols_name)) %>%
    pluck(plot_column) %>%
    patchwork::wrap_plots(ncol = length(cols), nrow = length(rows)) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  list(
    list(empty_plt, cols_labels_plt) %>% wrap_plots(widths = c(0.05, 1)),
    list(rows_labels_plt, plots_plt) %>% wrap_plots(nrow = 1, widths = c(0.05, 1))
  )  %>%
    patchwork::wrap_plots(ncol = 1, heights = c(0.05, 1))
}