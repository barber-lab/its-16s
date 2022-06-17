#' Plot alpha diversity between sample paires of two disjunct taxa groups
get_alphadiv_plot_by_taxa_group <- function(alphadiv, taxa_grouping, taxa_group_x, taxa_group_y, metric = "Chao1", cor_method = "spearman") {
  alphadiv %>%
    unnest(data) %>%
    filter(
      taxa_grouping == taxa_grouping &
        samples_grouping == "environment_group" &
        taxa_group %in% c(taxa_group_x, taxa_group_y)
    ) %>%
    group_by(samples_grouping, samples_group, taxa_grouping) %>%
    nest() %>%
    mutate(
      alphadiv_plt = data %>% map2(samples_group, ~ {
        .x %>%
          select_at(c("taxa_group", "sample_id", metric)) %>%
          pivot_wider(names_from = taxa_group, values_from = metric) %>%
          ggplot(aes_string(taxa_group_x, taxa_group_y)) +
          geom_point(aes(color = .y), alpha = 0.3) +
          stat_smooth() +
          scale_color_environment_group() +
          stat_cor(method = cor_method) +
          guides(color = FALSE) +
          labs(
            title = .y,
            x = glue("{taxa_group_x} ({metric})"),
            y = glue("{taxa_group_y} ({metric})")
          )
      })
    ) %>%
    arrange_columns()
}
