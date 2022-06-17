#' Read a biom file as a tibble with columns run_id, otu, and abundance
read_biom_tbl <- function(biom_path) {
  biom <- biomformat::read_biom(biom_path)

  biom %>%
    purrr::pluck("data") %>%
    as.matrix() %>%
    tibble::as_tibble() %>%
    dplyr::rename(abundance = V1) %>%
    tidyr::unnest(abundance) %>%
    dplyr::mutate(
      otu = biom$rows %>% purrr::map_chr(~ pluck(.x, "id")),
      run_id = biom$columns[[1]]$id
    ) %>%
    dplyr::select(run_id, otu, abundance)
}
