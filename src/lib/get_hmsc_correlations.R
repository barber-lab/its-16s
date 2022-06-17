get_hmsc_correlations <- function(model, ...) {
  model %>%
    getPostEstimate("OmegaCor", ...) %>%
    # for mean, support, and supportNeg
    map(function(m) {
      rownames(m) <- model$spNames
      colnames(m) <- model$spNames
      
      m %>%
        as_tibble(rownames = "from_clean_taxon") %>%
        pivot_longer(-from_clean_taxon, names_to = "to_clean_taxon")
    }) %>%
    enframe() %>%
    unnest(value) %>%
    pivot_wider(names_from = name, values_from = value)
}
