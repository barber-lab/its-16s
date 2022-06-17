#' Long tibble to phyloseq object
#' @param abundance_tbl tibble with columns sample_id, taxon, and abundance
tbl_to_phy <- function(abundance_tbl, samples_tbl, lineages_tbl, pooling_col = "genus", taxon_name = "taxon") {
  otu_phy <-
    abundance_tbl %>%
    dplyr::group_by_at(c("sample_id", taxon_name)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select_at(c("sample_id", taxon_name, "abundance")) %>%
    tidyr::pivot_wider(
      names_from = taxon_name, values_from = "abundance",
      values_fill = list(abundance = 0)
    ) %>%
    as_matrix("sample_id") %>%
    phyloseq::otu_table(taxa_are_rows = FALSE)

  samples_phy <-
    samples_tbl %>%
    dplyr::collect() %>%
    dplyr::filter(sample_id %in% abundance_tbl$sample_id) %>%
    dplyr::distinct() %>%
    phyloseq::sample_data()
  phyloseq::sample_names(samples_phy) <- samples_phy$sample_id

  tax_tbl <-
    lineages_tbl %>%
    dplyr::filter_at(pooling_col, ~ .x %in% phyloseq::taxa_names(otu_phy)) %>%
    dplyr::select(get_taxranks(pooling_col, direction = "with_upstream")) %>%
    group_by_at(pooling_col) %>%
    slice(1) %>%
    ungroup()

  tax_phy <- phyloseq::tax_table(tax_tbl)
  phyloseq::taxa_names(tax_phy) <- tax_tbl[[pooling_col]]

  phy <-
    phyloseq::phyloseq(otu_phy, tax_phy, samples_phy) %>%
    phyloseq::prune_samples(phyloseq::sample_sums(.) > 0, .) %>%
    phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0, .)

  phy
}
