get_dysbalances <- function(abundances, samples, lineages,
                            pooling_col, dist_method = "bray", groupings = c("bioproject_id", "kingdom"), max_clusters, max_balanced_quantile) {
  abundances %>%
    # account for library size when pooling
    filter(norm_method == "tss") %>%
    pull(data) %>%
    first() %>%
    # only join required columns for speed and space
    left_join(samples %>% select(sample_id, bioproject_id, environment_group), by = "sample_id") %>%
    left_join(lineages %>% select(genus, kingdom), by = c("taxon" = pooling_col)) %>%
    group_by_at(groupings) %>%
    nest() %>%
    mutate(
      dysbalance_tbl = data %>% map(possibly(~ {
        .x %>%
          tbl_to_phy(samples_tbl = samples, lineages_tbl = lineages) %>%
          phyloseq::distance(method = dist_method) %>%
          get_dysbalance(max_clusters = max_clusters, max_balanced_quantile = max_balanced_quantile) %>%
          pluck("dysbalance_tbl")
      }, NULL))
    ) %>%
    select(dysbalance_tbl) %>%
    unnest(dysbalance_tbl)
}
