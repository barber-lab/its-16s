get_samples_count_tbl <- function(samples_grouping) {
  mapped_samples <- samples_tbl %>%
    pull(sample_id) %>%
    unique()
  abundant_samples <- abundance_tbl %>%
    pull(sample_id) %>%
    unique()

  count_samples <- function(samples_grouping, selected_samples) {
    samples_tbl %>%
      filter(sample_id %in% selected_samples) %>%
      group_by_at(samples_grouping) %>%
      count()
  }

  samples_count_tbl <-
    samples_tbl %>%
    left_join(count_samples(samples_grouping, mapped_samples) %>%
      rename(mapped = n)) %>%
    left_join(count_samples(samples_grouping, abundant_samples) %>%
      rename(abundant = n)) %>%
    collect() %>%
    mutate_at(c("mapped", "abundant"), ~ replace_na(.x, 0))

  samples_count_tbl
}
