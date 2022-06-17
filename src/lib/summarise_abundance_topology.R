summarise_abundance_topology <- function(abundance_topology) {
  abundance_topology %>%
    mutate(plt = plt %>% map(~ {
      .x +
        theme(
          axis.title = element_text(size = 13)
        )
    })) %>%
    group_by(topology_name, cor_method, samples_grouping, taxa_grouping) %>%
    nest() %>%
    mutate(plt = data %>% map2(topology_name, facet_groupings_plots)) %>%
    select(-data)
}
