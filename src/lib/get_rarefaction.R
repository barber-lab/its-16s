# outside of pmap to prevent nesting issues
get_rarefaction <- function(sample_grouping, rarefactions_trails, samples, alphadiv_metrics, kingdoms_groups, abundances, lineages) {
  taxa_deck <-
    abundances$data[[1]] %>%
    left_join(lineages)
  
  samples %>%
    count(.data[[sample_grouping]]) %>%
    mutate(n_samples = n %>% map(~ seq(1, .x, length.out = 10) %>% round())) %>%
    select(.data[[sample_grouping]], n_samples) %>%
    unnest(n_samples) %>%
    expand_grid(rarefactions_trail = rarefactions_trails) %>%
    mutate(
      samples = list(n_samples, .data[[sample_grouping]], rarefactions_trail) %>% pmap(~ {
        set.seed(..3)

        samples %>%
          filter(.data[[sample_grouping]] == .y) %>%
          pull(sample_id) %>%
          sample(size = .x, replace = FALSE)
      })
    ) %>%
    # Use the same random sample subset for both kingdoms and all alphadiv metrics
    expand_grid(
      alphadiv_metric = c(alphadiv_metrics, "Observed"),
      kingdom = kingdoms_groups
    ) %>%
    mutate(
      alpha_diversity = list(samples, kingdom, alphadiv_metric) %>% pmap_dbl(possibly(~ {
        taxa_deck %>%
          filter(sample_id %in% .x & kingdom == .y) %>%
          # pool every sample into one to asses
          # alpha diversity of all samples together
          group_by(taxon) %>%
          summarise(
            sample_id = "pooled_sample",
            taxon = first(taxon),
            abundance = sum(abundance)
            #.groups = "taxon"
          ) %>%
          pivot_wider(names_from = taxon, values_from = abundance, values_fill = list(abundance = 0)) %>%
          column_to_rownames("sample_id") %>%
          as.matrix() %>%
          otu_table(taxa_are_rows = FALSE) %>%
          estimate_richness() %>%
          pluck(..3)
      }, NA))
    )
}
