#' Calculate keystoneness
#'
#' Keystone taxa have a high mean degree and closeness but low betweeness centrality (Nat Rev Microbiol 16, 567–576 (2018))
#' Furthermore, they are very prevalent. This method treats all properties with equal weight and calculates a keystoness in range of 0 to 1
#' @param graph coabundance tibble graph
#' @param prevalences tibble with columns taxon and prevalence_perc
#' @return tibble graph with node columns including keystoneness attached
add_keystoneness <- function(graph, prevalences) {
  graph %>%
    activate(nodes) %>%
    left_join(prevalences %>% mutate_at("taxon", as.character(), by = "taxon")) %>%
    mutate(
      graph_order = graph_order(),
      betweeness = centrality_betweenness(),
      closeness = centrality_closeness(),
      weighted_degree = centrality_degree(weights = estimate),
      # High high mean degree, high closeness centrality, low betweeness centrality and prevalence
      keystoneness = (
        rank(weighted_degree) +
          rank(closeness) -
          rank(betweeness) -
          rank(prevalence_perc)
      ) %>% scales::rescale()
    )
}
