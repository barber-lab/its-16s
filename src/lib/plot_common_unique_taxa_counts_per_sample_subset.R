plot_common_unique_taxa_counts_per_sample_subset <- function(abundances, lineages, samples, kingdoms_groups) {
  prevalent_taxa <-
    abundances$data[[1]] %>%
    left_join(lineages) %>%
    left_join(samples) %>%
    filter(abundance > 0 & kingdom %in% kingdoms_groups) %>%
    distinct(environment_group, kingdom, phylum, taxon) %>%
    mutate(
      is_prevalent = TRUE,
      phylum_color = phyla_colors[phylum] %>% replace_na("grey")
    ) %>%
    pivot_wider(
      names_from = environment_group,
      values_from = is_prevalent,
      values_fill = list(is_prevalent = FALSE)
    ) %>%
    distinct(taxon, .keep_all = TRUE)
  
  kingdoms_groups %>%
    map(~ {
      data <- prevalent_taxa %>% filter(kingdom == .x)
      
      tibble() %>%
        ggplot(aes(y = ..count.. / nrow(data) * 100, fill = phylum_color)) +
        geom_bar(
          data = data %>% filter(aquatic & !soil & !host),
          mapping = aes(x = "unique to aquatic")
        ) +
        geom_bar(
          data = data %>% filter(!aquatic & !soil & host),
          mapping = aes(x = "unique to host")
        ) +
        geom_bar(
          data = data %>% filter(!aquatic & soil & !host),
          mapping = aes(x = "unique to soil")
        ) +
        geom_bar(
          data = data %>% filter(aquatic & host),
          mapping = aes(x = "common in aquatic and host")
        ) +
        geom_bar(
          data = data %>% filter(aquatic & soil),
          mapping = aes(x = "common in aquatic and soil")
        ) +
        geom_bar(
          data = data %>% filter(soil & host),
          mapping = aes(x = "common in soil and host")
        ) +
        geom_bar(
          data = data %>% filter(aquatic & host & soil),
          mapping = aes(x = "common in all")
        ) +
        labs(
          x = "",
          y = str_glue("Genera (% of all {.x})")
        ) +
        theme(
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()
        ) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 60)) +
        facet_wrap(~kingdom, ncol = 1, scales = "free") +
        coord_flip()
    }) %>%
    wrap_plots(ncol = 1)
}
