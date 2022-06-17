stat_abundances <- function(sub_abundances) {
  # Either samples or taxa group must be all
  if (
    list(
      sub_abundances$samples_group[[1]],
      sub_abundances$taxa_group[[1]]
    ) %>%
      flatten_chr() %>%
      keep(~ .x == "all") %>%
      length() != 1
  ) {
    return(tibble())
  }

  sub_abundances %>%
    mutate(data = data %>% map(~ .x %>% select(colnames(.x) %>% setdiff(colnames(sub_abundances))))) %>%
    unnest(data) %>%
    group_by(norm_method, samples_grouping, taxa_grouping, taxon) %>%
    nest() %>%
    mutate(
      unique_group = data %>% map_lgl(~ .x$samples_group %>%
        unique() %>%
        length() == 1)
    ) %>%
    filter(!unique_group) %>%
    mutate(
      test = data %>% map(~ possibly(test_aov, NA)(
        formula = abundance ~ samples_group,
        data = .x
      ))
    )
}
