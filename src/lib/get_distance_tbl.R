#' Get distance tibble from phyloseq object
get_distance_tbl <- function(phy, method = "bray", ...) {
  phy %>%
    phyloseq::distance(method = method, ...) %>%
    broom::tidy() %>%
    dplyr::rename(from = item1, to = item2)
}
