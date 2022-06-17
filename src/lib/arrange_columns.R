#' Arrange columns using provided ordering
#'
#' @param data: data farme to be arranged
#' @param ordering: character vector of columns names in the desired order. Superset allowed.
#' @return `data` with arranged columns. Columns present in data but not in `ordering` will be attached at the right end
arrange_columns <- function(
                            data,
                            ordering = c(
                              "samples_grouping", "samples_group", "taxa_grouping", "taxa_group",
                              "norm_method", "dist_method", "dist_type", "ord_method"
                            )) {
  data_cols <- colnames(data)
  ordered_cols <- ordering %>% intersect(data_cols)
  other_cols <- data_cols %>% setdiff(ordering)

  data %>%
    dplyr::select(
      tidyselect::all_of(ordered_cols), other_cols
    )
}
