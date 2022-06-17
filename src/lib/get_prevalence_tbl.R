# Get prevalence tibble
#
# @param data tibble with columns taxon, sample_id, and abundance
# @param min_abundance a taxon must have a greater abundance than this value to be considered as prevalent
get_prevalence_tbl <- function(data, min_abundance = 0) {
  data %>%
    filter(abundance > min_abundance) %>%
    select(sample_id, taxon) %>%
    distinct() %>%
    group_by(taxon) %>%
    count() %>%
    rename(prevalence_n = n) %>%
    mutate(prevalence_perc = prevalence_n / length(data$sample_id %>% unique()) * 100) %>%
    arrange(-prevalence_n) %>%
    ungroup()
}
