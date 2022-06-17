#' Summarizes graphs with same cor method and groupings for sampels and taxa
summarise_graphs <- function(graphs) {
  graphs %>%
    group_by(cor_method, samples_grouping, taxa_grouping) %>%
    nest() %>%
    mutate(
      cor_graph = data %>% map(~ {
        .x %>%
          filter(!is.na(cor_graph)) %>%
          mutate(cor_graph = list(cor_graph, samples_group, taxa_group) %>% pmap(~ {
            ..1 %>%
              activate(edges) %>%
              mutate(
                samples_group = ..2,
                taxa_group = ..3
              )
          })) %>%
          pull(cor_graph) %>%
          reduce(graph_join)
      }),
      cor_plt = cor_graph %>% map(~ {
        .x %>%
          tidygraph::activate(edges) %>%
          ggraph::ggraph("circle") +
          ggraph::geom_edge_link(aes(color = estimate), size = 6) +
          ggraph::geom_node_point(aes(color = kingdom)) +
          ggraph::facet_edges(samples_group ~ taxa_group) +
          ggplot2::theme_void() +
          scale_color_kingdom() +
          scale_edge_color_gradient2(
            high = "#770000",
            low = "#000077",
            midpoint = 0
          ) +
          ggplot2::coord_fixed()
      })
    ) %>%
    select(-data, -cor_graph)
}
