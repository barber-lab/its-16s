get_dominances <- function(abundances, samples, lineages, taxa_grouping, samples_grouping, pooling_col, dominance_quantile) {
  abundances %>%
    # account for library size when pooling
    filter(norm_method == "tss") %>%
    pull(data) %>%
    first() %>%
    left_join(
      samples %>%
        select_at(c("sample_id", samples_grouping)),
      by = "sample_id"
    ) %>%
    left_join(
      lineages %>%
        select_at(c("kingdom", pooling_col)),
      by = c("taxon" = pooling_col)
    ) %>%
    group_by_at(c(samples_grouping, taxa_grouping)) %>%
    nest() %>%
    mutate(
      dominant = data %>% map(~ get_dominant_taxa(.x, quantile = dominance_quantile)),
      taxon = data %>% map(~ .x$taxon %>% unique())
    ) %>%
    unnest(taxon) %>%
    mutate(
      dominance = map2_chr(taxon, dominant, ~ {
        .x %in% .y %>% ifelse("dominant", "other")
      })
    ) %>%
    select_at(c(samples_grouping, taxa_grouping, "taxon", "dominance")) %>%
    arrange_columns()
}
