get_subsets <- function() {
  environment_groups <- c("host", "aquatic", "soil")
  kingdoms_groups <- c("Bacteria", "Fungi")
  dominance_groups <- c("dominant", "rare")
  dysbalance_groups <- c("dysbalanced", "balanced")

  environment_group_kingdom_subsets <-
    tibble(environment_group = environment_groups) %>%
    expand_grid(kingdom = kingdoms_groups) %>%
    transmute(
      subset_name = "environment_group,kingdom",
      subset_value = str_glue("{environment_group},{kingdom}"),
      subset_query = str_glue("environment_group == '{environment_group}' & kingdom == '{kingdom}'")
    )

  dysbalanced_kingdom_subsets <-
    expand_grid(dysbalance = dysbalance_groups) %>%
    expand_grid(kingdom = kingdoms_groups) %>%
    transmute(
      subset_name = "dysbalance,kingdom",
      subset_value = str_glue("{dysbalance},{kingdom}"),
      subset_query = str_glue("kingdom == '{kingdom}' & dysbalance == '{dysbalance}'")
    )

  dysbalanced_environment_group_kingdom_subsets <-
    tibble(environment_group = environment_groups) %>%
    expand_grid(dysbalance = dysbalance_groups) %>%
    expand_grid(kingdom = kingdoms_groups) %>%
    transmute(
      subset_name = "dysbalance,environment_group,kingdom",
      subset_value = str_glue("{dysbalance},{environment_group},{kingdom}"),
      subset_query = str_glue("environment_group == '{environment_group}' & kingdom == '{kingdom}' & dysbalance == '{dysbalance}'")
    )

  dominance_environment_group_kingdom_subsets <-
    tibble(environment_group = environment_groups) %>%
    expand_grid(dominance = dominance_groups) %>%
    expand_grid(kingdom = kingdoms_groups) %>%
    transmute(
      subset_name = "dominance,environment_group,kingdom",
      subset_value = str_glue("{dominance},{environment_group},{kingdom}"),
      subset_query = str_glue("dominance == '{dominance}' & environment_group == '{environment_group}' & kingdom == '{kingdom}'")
    )

  list(
    environment_group_kingdom_subsets,
    dysbalanced_kingdom_subsets,
    dysbalanced_environment_group_kingdom_subsets,
    dominance_environment_group_kingdom_subsets
  ) %>%
    # ignore type glue
    map(~ .x %>% mutate_all(as.character)) %>%
    bind_rows() %>%
    view()
}
