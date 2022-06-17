#!/usr/bin/env Rscript

library(purrr)
library(broom)
library(Hmisc)
library(RSQLite)
library(DBI)
library(SpiecEasi)
library(tidygraph)
library(ggraph)
library(optparse)
library(tidyverse)

#' Get abundance matrix
#' @param per_kingdom return a list of matricies. One per kingdom. Needed for mb.
get_abundance_matrix <- function(
                                 pooling_col, samples_filter_query,
                                 min_prevalence_perc, per_kingdom) {
  selected_downstream_cols <-
    get_taxranks(pooling_col, "downstream")

  con <- dbConnect(SQLite(), "etc/db.sqlite")
  abundances_tbl <- con %>% tbl(glue("abundance_{pooling_col}_raw"))
  samples_tbl <- con %>% tbl("samples")
  lineages_tbl <- tbl(con, "lineages")

  tbl <-
    abundances_tbl %>%
    inner_join(samples_tbl, by = "sample_id") %>%
    # filter(sample_id %in% samples) %>%
    group_by(kingdom) %>%
    collect() %>%
    filter(eval(parse(text = samples_filter_query)))

  # samples must have both fungal and bacterial data available
  valid_samples <-
    tbl %>%
    select(sample_id, kingdom) %>%
    distinct() %>%
    mutate(available = TRUE) %>%
    pivot_wider(names_from = kingdom, values_from = available) %>%
    filter(bacteria & fungi) %>%
    pull(sample_id)

  get_mat <- function(tbl) {
    tbl %>%
      select(sample_id, taxon, abundance) %>%
      pivot_wider(names_from = "taxon", values_from = "abundance") %>%
      distinct() %>%
      arrange(sample_id) %>%
      select(-sample_id) %>%
      replace(., is.na(.), 0) %>%
      as.matrix()
  }

  tbl <-
    tbl %>%
    filter(sample_id %in% valid_samples) %>%
    arrange(sample_id)

  if (per_kingdom) {
    tbl %>%
      group_by(kingdom) %>%
      nest() %>%
      mutate(mats = data %>% map(get_mat)) %>%
      pull(mats)
  } else {
    tbl %>%
      ungroup() %>%
      get_mat()
  }
}
