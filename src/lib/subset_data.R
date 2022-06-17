subset_data <- function(data, samples_grouping, samples_group, taxa_grouping, taxa_group) {
  res <- data

  if (samples_group != "all") {
    res %<>%
      dplyr::filter_at("samples_grouping", ~ .x == samples_grouping) %>%
      dplyr::filter_at("samples_group", ~ .x == samples_group)
  }

  if (taxa_group != "all") {
    res %<>% filter(!!sym(taxa_grouping) == taxa_group)
  }

  res
}
