#' Calculate foundation taxa score
#'
#' @param graph coabundance tibble graph
#' @param abundances tibble with columns sample_id, taxon abundance
#' @return tibble graph with node columns including keystoneness attached
add_foundationess <- function(graph, abundances) {
  abundances_summary <-
    abundances %>%
    dplyr::group_by(taxon) %>%
    dplyr::summarize(abundance = mean(abundance)) %>%
    dplyr::arrange(-abundance)

  pos_weighted_degrees <-
    graph %>%
    tidygraph::activate(edges) %>%
    dplyr::filter(estimate > 0) %>%
    tidygraph::activate(nodes) %>%
    mutate(
      pos_weighted_degree = centrality_degree(weights = estimate)
    ) %>%
    as_tibble() %>%
    select(taxon, pos_weighted_degree)

  graph %>%
    activate(nodes) %>%
    left_join(abundances_summary, by = "taxon") %>%
    left_join(pos_weighted_degrees, by = "taxon") %>%
    mutate(
      graph_order = graph_order(),
      betweeness = centrality_betweenness(),
      closeness = centrality_closeness(),
      weighted_degree = centrality_degree(weights = estimate),
      foundationess = (
        rank(abundance) + rank(pos_weighted_degree)
      ) %>% scales::rescale()
    )
}
