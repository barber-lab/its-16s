#' Get samples meta data with abundances for constrained ordination e.g. dbRDA
get_samples_abundance_meta <- function(sub_abundance, other_kingdom) {
  sub_abundance %>%

    # pool counts
    group_by(sample_id, genus) %>%
    summarise(abundance = sum(abundance)) %>%
    select(sample_id, genus, abundance) %>%
    filter(sample_id %in% features) %>%
    pivot_wider(names_from = genus, values_from = abundance, values_fill = list(abundance = 0))
}
