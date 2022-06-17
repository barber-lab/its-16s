plot_samples_by_samples_grouping_and_stage <- function(samples_grouping) {
  samples_count_tbl <- get_samples_count_tbl(samples_grouping)

  samples_count_tbl %>%
    mutate(max_samples = mapped) %>%
    pivot_longer(cols = c("mapped", "abundant"), values_to = "n", names_to = "pipeline_stage") %>%
    mutate_at(samples_grouping, ~ .x %>% fct_reorder(max_samples)) %>%
    filter(n >= 10) %>%
    ggplot(aes_string(samples_grouping, "n")) +
    geom_bar(aes(color = pipeline_stage), stat = "identity", position = "identity") +
    coord_flip()
}
