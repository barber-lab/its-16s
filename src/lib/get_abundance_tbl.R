get_abundance_tbl <- function(taxrank = "genus", normalization_method = "tss", filter_zeros = FALSE) {
  res <- list.files(
    path = "raw/abundances",
    pattern = glue("{taxrank}\\.{normalization_method}"),
    full.names = TRUE
  ) %>%
    map({
      ~ .x %>%
        read_csv() %>%
        pivot_longer(-sample_id, names_to = taxrank, values_to = "abundance")
    }) %>%
    bind_rows() %>%
    rename(taxon = !!taxrank)

  if (filter_zeros) res %<>% dplyr::filter(abundance != 0)

  res
}
