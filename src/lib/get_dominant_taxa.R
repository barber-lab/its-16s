# Get dominant taxa
#
# @param data: tibble with columns abundance
# @param quantile: dominance threshold by quantile: e.g. 0.1 to get the 10% most dominant taxa
get_dominant_taxa <- function(data, quantile = 0.1) {
  data %>%
    group_by(taxon) %>%
    dplyr::summarise(abundance = mean(abundance)) %>%
    dplyr::arrange(-abundance) %>%
    purrr::pluck("taxon") %>%
    unique() %>%
    head(length(.) * quantile)
}
