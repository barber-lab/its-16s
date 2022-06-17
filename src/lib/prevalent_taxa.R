#' Get prevalent taxa
#' @param min_prevalence min fraction of samples required to count as prevalent
#'   in any sample group defined by samples_grouping
#' @param min_abundance min fraction (TSS scaled value)
#'   required to count as prevalent in a sample
#' @param pooling_col taxonomic rank e.g. family to get prevalent families
#' @param samples_grouping taxon must be prevalent if it is prevalent
#'   in any of the groups defined in column samples_grouping
prevalent_taxa <- function(
                           min_prevalence = 0.2,
                           min_abundance = 0.1e-2,
                           pooling_col = "genus",
                           samples_grouping = "bioproject_id") {
  con <- DBI::dbConnect(RSQLite::SQLite(), "etc/db.sqlite")
  abundances_tbl <-
    dplyr::tbl(con, "abundance") %>%
    filter(normalization_method == "tss" & pooling_col == !!pooling_col)
  samples_tbl <- dplyr::tbl(con, "samples")

  n_samples_tbl <-
    samples_tbl %>%
    group_by_at(samples_grouping) %>%
    count() %>%
    rename(n_samples = n)

  prevalent_taxa <-
    abundances_tbl %>%
    left_join(samples_tbl) %>%
    filter(abundance >= min_abundance) %>%
    group_by_at(vars("taxon", samples_grouping)) %>%
    count() %>%
    left_join(n_samples_tbl) %>%
    filter(n >= min_prevalence * n_samples) %>%
    collect() %>%
    pull(taxon) %>%
    unique()

  prevalent_taxa
}
