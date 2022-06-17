get_topologies <- function(graphs, node_topology_metrics) {
  graphs %>%
    filter(!is.na(cor_graph)) %>%
    mutate(
      nodes_tbl = cor_graph %>% map(~ {
        .x %>%
          activate(nodes) %>%
          as_tibble()
      })
    ) %>%
    unnest(nodes_tbl) %>%
    pivot_longer(
      cols = node_topology_metrics,
      names_to = "topology_metric",
      values_to = "topology_value"
    ) %>%
    group_by(samples_grouping, taxa_grouping, taxa_group, cor_method, topology_metric) %>%
    nest() %>%
    mutate(topol_plt = list(data, topology_metric, samples_grouping) %>% pmap(~ {
      plt <-
        ..1 %>%
        ggplot(aes(samples_group, topology_value)) +
        labs(
          y = ..2,
          x = ""
        )

      switch(..3,
        "environment_group" = {
          plt +
            geom_boxplot(aes(fill = samples_group)) +
            stat_compare_means_environment_group() +
            scale_fill_environment_group()
        },
        "dysbalance" = {
          plt +
            geom_boxplot(aes(fill = samples_group)) +
            stat_compare_means_dysbalance() +
            scale_fill_dysbalance()
        },
        {
          plt +
            geom_boxplot() +
            stat_compare_means()
        }
      )
    }))
}
