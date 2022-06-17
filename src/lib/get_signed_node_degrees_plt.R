get_signed_node_degrees_plt <- function(graph) {
  pos_weighted_degrees <-
    graph %>%
    tidygraph::activate(edges) %>%
    dplyr::filter(estimate > 0) %>%
    tidygraph::activate(nodes) %>%
    mutate(
      pos_weighted_degree = centrality_degree(weights = estimate)
    ) %>%
    as_tibble() %>%
    select(name, pos_weighted_degree)

  neg_weighted_degrees <-
    graph %>%
    tidygraph::activate(edges) %>%
    dplyr::filter(estimate < 0) %>%
    tidygraph::activate(nodes) %>%
    mutate(
      neg_weighted_degree = centrality_degree(weights = abs(estimate))
    ) %>%
    as_tibble() %>%
    select(name, neg_weighted_degree)

  weighted_degrees <-
    graph %>%
    tidygraph::activate(nodes) %>%
    mutate(
      weighted_degree = centrality_degree(weights = estimate)
    ) %>%
    as_tibble() %>%
    select(name, weighted_degree)

  degrees <-
    list(
      pos_weighted_degrees,
      neg_weighted_degrees,
      weighted_degrees
    ) %>%
    reduce(~ inner_join(.x, .y, by = "name")) %>%
    mutate(
      weighted_degree = pos_weighted_degree + neg_weighted_degree,
      ratio_weighted_degree = pos_weighted_degree / neg_weighted_degree,
      fraction_pos_degree = pos_weighted_degree / weighted_degree
    )

  lm(weighted_degree ~ ratio_weighted_degree, data = degrees) %>% tidy()
  lm(pos_weighted_degree ~ neg_weighted_degree, data = degrees) %>% tidy()
  lm(fraction_pos_degree ~ weighted_degree, data = degrees) %>% tidy()
  cor.test(~ pos_weighted_degree + neg_weighted_degree, data = degrees) %>% tidy()


  degrees %>%
    ggplot(aes(weighted_degree, fraction_pos_degree)) +
    geom_point(alpha = 0.2) +
    stat_smooth(
      data = degrees %>% mutate(group = "[0,1]"),
      mapping = aes(color = group, fill = group),
      method = "glm",
      method.args = list(family = "quasibinomial"),
    ) +
    stat_smooth(
      data = degrees %>% filter(fraction_pos_degree > 0 & fraction_pos_degree < 1) %>% mutate(group = "]0,1["),
      mapping = aes(color = group, fill = group),
      method = "glm",
      method.args = list(family = "quasibinomial"),
    ) +
    scale_y_continuous(labels = percent) +
    scale_color_npg() +
    scale_fill_npg() +
    labs(
      x = "Weighted degree",
      y = "Positive weighted degree",
      color = "y interval",
      fill = "y interval"
    )
}
