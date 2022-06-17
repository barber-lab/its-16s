summarise_alphadiv <- function(alphadiv, alphadiv_metrics) {
  alphadiv %>%
    unnest(data) %>%
    pivot_longer(
      cols = alphadiv_metrics,
      names_to = "metric",
      values_to = "alphadiv"
    ) %>%
    group_by(metric, samples_grouping, taxa_grouping, taxa_group) %>%
    nest() %>%
    mutate(
      test = data %>% map(~
      possibly(test_aov, NULL)(alphadiv ~ samples_group, .x) %>%
        mutate(q.value = p.adjust(p.value))),
      box_plt = list(data, metric, samples_grouping) %>% pmap(~ {
        plt <-
          ..1 %>%
          ggplot2::ggplot(aes(samples_group, alphadiv)) +
          ggplot2::geom_boxplot(aes(fill = samples_group)) +
          ggplot2::labs(
            y = ..2,
            x = "",
            fill = samples_grouping
          )

        switch(..3,
          "environment_group" = {
            plt +
              stat_compare_means_environment_group() +
              scale_fill_environment_group()
          },
          "dysbalance" = {
            plt +
              stat_compare_means_dysbalance() +
              scale_fill_dysbalance()
          },
          {
            plt +
              scale_fill_npg() +
              stat_compare_means()
          }
        )
      })
    ) %>%
    select(-data) %>%
    arrange_columns()
}
