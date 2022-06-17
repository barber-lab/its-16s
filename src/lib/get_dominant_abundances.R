get_dominant_abundances <- function(sub_abundances, dominant_taxa) {
  environment_groups_dominances <-
    dominant_taxa %>%
    filter(
      norm_method == "tss" &
        samples_grouping == "environment_group" &
        taxa_grouping == "kingdom" &
        taxa_group %in% c("Bacteria", "Fungi")
    ) %>%
    unnest(cols = data) %>%
    mutate(is_dominant = TRUE) %>%
    select(samples_group, taxon = data, is_dominant)

  dominant_abundances <-
    sub_abundances %>%
    filter(
      norm_method == "tss" &
        samples_grouping == "environment_group" &
        taxa_grouping == "kingdom" &
        taxa_group %in% c("Bacteria", "Fungi")
    ) %>%
    unnest(data) %>%
    left_join(environment_groups_dominances) %>%
    replace_na(list(is_dominant = FALSE)) %>%
    group_by(taxa_group, samples_group, is_dominant) %>%
    dplyr::summarize(abundance = sum(abundance)) %>%
    nest() %>%
    mutate(
      pie_plt = data %>% map2(samples_group, ~ {
        .x %>%
          ggpie("abundance", fill = "is_dominant") +
          labs(title = .y)
      })
    )

  dominant_abundances
}
