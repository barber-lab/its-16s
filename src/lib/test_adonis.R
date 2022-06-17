#' Test using adonis PERMANOVA
#' @param formula formula with LHS abundance. All variables used in RHS must be columns of `data`
#' @param data tibble with columns taxon, sample_id, taxon, abundance and all covariates used in RHS of `formula`
#' @param parallel number of threads for `vegan::adonis`
#' @param type Features to test. Either 'taxa' or 'samples'
#' @param ... additional arguments forwarded to `vegan::adonis`
#' @return tidy tibble of adonis test results
test_adonis <- function(
                        formula = formula(abundance ~ environment_group),
                        data,
                        parallel = 3,
                        type = "samples",
                        ...) {
  covariates <-
    formula %>%
    as.list() %>%
    magrittr::extract2(3) %>%
    as.character() %>%
    stringr::str_split("\\+") %>%
    purrr::simplify() %>%
    purrr::discard(~ .x == "")

  data %<>% dplyr::filter_at(covariates, ~ !is.na(.x))

  abundance <-
    data %>%
    dplyr::select(taxon, sample_id, abundance) %>%
    tidyr::pivot_wider(
      names_from = "taxon", values_from = "abundance",
      values_fill = list(abundance = 0)
    ) %>%
    magrittr::set_rownames(.$sample_id) %>%
    dplyr::select(-sample_id) %>%
    as.matrix()

  if (type == "taxa") abundance %<>% t()

  covariates_tbl <-
    switch(type,
      "samples" = {
        data %>%
          distinct_at(c("sample_id", covariates)) %>%
          set_rownames(.$sample_id)
      },
      "taxa" = {
        data %>%
          distinct_at(c("taxon", covariates)) %>%
          set_rownames(.$taxon)
      }
    )

  adonis_res <-
    vegan::adonis(
      formula = formula,
      data = covariates_tbl,
      parallel = parallel,
      ...
    )

  adonis_res %>%
    purrr::pluck("aov.tab") %>%
    tibble::as_tibble(rownames = "var") %>%
    dplyr::mutate(method = "Adonis PERMANOVA") %>%
    dplyr::rename(p.value = `Pr(>F)`) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(q.value = p.value %>% p.adjust(method = "fdr"))
}
