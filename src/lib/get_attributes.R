#' Get list of parameters from a drake target name
#' @param target drake target name produced by drake using static brancing
#' @param  attributes vector of group names. These are the values to search in the target name
get_attributes <- function(target, attributes) {
  attributes %>%
    enframe(name = "attribute", value = "group") %>%
    rowwise() %>%
    mutate(
      group = group %>% paste(collapse = "|"),
      value = target %>% str_extract(group)
    ) %>%
    add_row(attribute = "step", value = target %>% str_extract("^[^_]+")) %>%
    select(attribute, value) %>%
    deframe()
}
