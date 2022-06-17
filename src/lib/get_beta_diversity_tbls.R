#' Get beta diversity tibbles
#'
#' @param nested_tbl tibble with columns samples_dist and taxa_dist
#' @param type to ordinate. One of taxa or samples
#' @return list of tibbles
get_beta_diversity_tbls <- function(nested_tbl, type = "taxa") {
  dist_col <- glue("{type}_dist")

  beta_diversity_tbl <-
    nested_tbl %>%
    select(params$samples$grouping, one_of(dist_col)) %>%
    unnest(dist_col) %>%
    ungroup()

  beta_diversity_tbls <-
    beta_diversity_tbl %>%
    mutate(
      grouping = (environment_group == "all" & !is.na(environment_group)) %>%
        ifelse("between", "within")
    ) %>%
    group_by(grouping) %>%
    nest() %>%
    deframe()

  beta_diversity_tbls
}
