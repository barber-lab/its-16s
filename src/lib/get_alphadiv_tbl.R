#' @param tbl with raw counts and columns sample_id, abundance and pooling_col
#' @param pooling_col col name of tbl used to pool counts
get_alpha_div_tbl <- function(tbl, pooling_col, measures = c("Shannon", "Chao1", "Simpson")) {
  tbl %>%
    group_by_at(c("sample_id", pooling_col)) %>%
    dplyr::summarize(abundance = sum(abundance)) %>%
    pivot_wider(names_from = pooling_col, values_from = "abundance", values_fill = list(abundance = 0)) %>%
    set_rownames(.$sample_id) %>%
    ungroup() %>%
    column_to_rownames("sample_id") %>%
    phyloseq::otu_table(taxa_are_rows = FALSE) %>%
    phyloseq::estimate_richness(measures = measures) %>%
    as_tibble(rownames = "sample_id") %>%
    select(-se.chao1) %>%
    # phyloseq::estimate_richness reformats sample_id
    # e.g. with prefix X if sample_id is numeric
    # Undo this behavior
    mutate(sample_id = sample_id %>% str_remove("^X"))
}
