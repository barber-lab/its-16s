parse_ancombc_result <- function(x, name) {
  x %>%
    pluck("res", name) %>%
    as_tibble(rownames = "taxon") %>%
    pivot_longer(-taxon, names_to = "samples_covariate_level", values_to = name)
}