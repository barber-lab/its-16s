#' Plot dissimilarities between sample paires of two disjunct taxa groups
get_distances_plot_by_taxa_group <- function(distances, taxa_grouping, taxa_group_x, taxa_group_y) {
  distances %>%
    filter(
      norm_method == "raw" &
        taxa_grouping == taxa_grouping &
        dist_type == "samples" &
        dist_method == "bray"
    ) %>%
    group_by(samples_grouping, samples_group, dist_method) %>%
    nest() %>%
    mutate(
      dists = data %>% map(~ {
        .x %>%
          select(taxa_group, dist) %>%
          mutate(dist = dist %>% map(tidy)) %>%
          unnest(dist) %>%
          pivot_wider(names_from = "taxa_group", values_from = "distance")
      }),
      dists_plt = list(dists, samples_group, dist_method) %>% pmap(~ {
        ..1 %>%
          sample_n(20e3) %>%
          ggplot(aes_string(taxa_group_x, taxa_group_y)) +
          geom_point(aes(color = ..2), alpha = 0.1) +
          geom_smooth(method = "lm") +
          ggpubr::stat_cor(method = "pearson") +
          scale_color_environment_group() +
          coord_fixed() +
          guides(color = FALSE) +
          labs(
            title = ..2,
            x = glue("{taxa_group_x} ({..3})"),
            y = glue("{taxa_group_y} ({..3})")
          )
      })
    ) %>%
    select(samples_grouping, samples_group, dist_method, dists_plt) %>%
    mutate(taxa_grouping = taxa_grouping) %>%
    arrange_columns()
}
