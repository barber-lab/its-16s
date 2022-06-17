#' Get relationship between abundance and topology
get_abundance_topology <- function(sub_abundances, graphs, node_topology_metrics, keystones, foundations) {
  graphs %>%
    inner_join(sub_abundances) %>%
    rename(abundances = data, coabundances = cor_graph) %>%
    filter(!is.na(coabundances)) %>%
    mutate(
      abundances = abundances %>% map(~ {
        .x %>%
          group_by(taxon) %>%
          dplyr::summarise(abundance = mean(abundance))
      }),
      coabundances = coabundances %>% map(~ {
        .x %>%
          activate(nodes) %>%
          as_tibble() %>%
          select(taxon, kingdom, node_topology_metrics) %>%
          pivot_longer(cols = node_topology_metrics, names_to = "topology_name", values_to = "topology_value")
      }),
      both = coabundances %>% map2(abundances, ~ inner_join(.x, .y, by = "taxon"))
    ) %>%
    # one group per topology metric
    select(-coabundances, -abundances) %>%
    unnest(both) %>%
    group_by(samples_grouping, samples_group, taxa_grouping, taxa_group, cor_method, topology_name) %>%
    nest() %>%
    # annotate keystones and foundations
    left_join(keystones %>% rename(keystones = data)) %>%
    left_join(foundations %>% rename(foundations = data)) %>%
    mutate(data = list(data, keystones, foundations) %>% pmap(~ ..1 %>%
      left_join(..2, by = "taxon") %>%
      left_join(..3, by = "taxon"))) %>%
    mutate(
      plt = data %>% map2(topology_name, ~ {
        .x %>%
          ggplot(aes(abundance, topology_value)) +
          stat_smooth(method = "loess", color = "grey", fill = "lightgrey") +
          stat_density_2d(aes(color = kingdom)) +
          scale_color_kingdom(drop = FALSE) +
          labs(color = "Kingdom") +
          ggnewscale::new_scale_color() +
          geom_point(
            data = .x %>% arrange(-keystoneness) %>% head(10),
            mapping = aes(color = "Keystone taxon")
          ) +
          geom_point(
            data = .x %>% arrange(-foundationess) %>% head(10),
            mapping = aes(color = "Foundation taxon")
          ) +
          scale_color_d3() +
          stat_cor(method = "spearman") +
          labs(color = "Importance") +
          scale_x_log10() +
          scale_y_log10() +
          annotation_logticks() +
          labs(
            x = "Mean abundance",
            y = .y
          )
      })
    ) %>%
    select(samples_grouping, samples_group, taxa_grouping, taxa_group, cor_method, plt)
}
