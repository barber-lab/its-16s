annotate_name_in_edges <- function(graph, from_names = "name", to_names_suffix = "taxon") {
  state <- graph %>% active()
  
  names <-
    graph %>% 
    activate(nodes) %>%
    mutate(id := row_number()) %>%
    as_tibble() %>%
    select(!!from_names, id)
  
  graph %>%
    activate(edges) %>%
    left_join(names %>% rename(!!paste0("from_", to_names_suffix) := name), by = c("from" = "id")) %>%
    left_join(names %>% rename(!!paste0("to_", to_names_suffix) := name), by = c("to" = "id")) %>%
    activate(!!state)
}
