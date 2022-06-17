#!/usr/bin/env R

library(ggraph)
library(igraph)
library(SpiecEasi)
library(broom)
library(tidygraph)

get_graph_pulsar <- function(cor_res, min_abs_estimate = 0) {
  used_taxa <- cor_res$est$data %>% colnames()

  switch(cor_res$est$method,
    "mb" = {
      cor_mb_edge_estimates_tbl <-
        cor_res %>%
        getOptBeta() %>%
        as.matrix() %>%
        as_tibble(rownames = "from") %>%
        pivot_longer(-one_of("from"),
          names_to = "to", values_to = "estimate"
        ) %>%
        dplyr::filter(abs(estimate) >= min_abs_estimate) %>%
        mutate(to = to %>% str_remove("^V") %>% as.integer()) %>%
        mutate_at(c("from", "to"), as.integer)

      cor_mb_graph <-
        cor_res %>%
        getRefit() %>%
        adj2igraph() %>%
        as_tbl_graph() %>%
        mutate(taxon = used_taxa) %>%
        activate(edges) %>%
        left_join(cor_mb_edge_estimates_tbl) %>%
        activate(nodes)

      cor_mb_graph
    },
    "glasso" = {
      cor_glasso_edge_estimates_tbl <-
        cor_res %>%
        getOptiCov() %>%
        as.matrix() %>%
        as_tibble(rownames = "from") %>%
        pivot_longer(-one_of("from"),
          names_to = "to", values_to = "estimate"
        ) %>%
        dplyr::filter(abs(estimate) >= min_estimate) %>%
        mutate(to = to %>% str_remove("^V") %>% as.integer()) %>%
        mutate_at(c("from", "to"), as.integer)

      cor_glasso_graph <-
        cor_res %>%
        getRefit() %>%
        adj2igraph() %>%
        as_tbl_graph() %>%
        mutate(taxon = used_taxa) %>%
        activate(edges) %>%
        left_join(cor_glasso_edge_estimates_tbl) %>%
        activate(nodes)

      cor_glasso_graph
    },
    error = function(e) stop("not implemented")
  )
}

filter_and_topologize_graph <- function(tbl_graph, min_abs_estimate = 0.1) {
  state <- tbl_graph %>% active()

  tbl_graph %>%
    activate(nodes) %>%
    # sanity filter nodes
    filter(kingdom %in% c("Bacteria", "Fungi")) %>%
    # annotate edges
    activate(edges) %>%
    # sanity filter edges
    filter(abs(estimate) >= min_abs_estimate) %>%
    # filter({
    #   if ("p.value" %in% names(.)) p.value else NULL
    # } <= 0.05) %>%
    activate(nodes) %>%
    # calculate node topology features
    mutate(
      degree = centrality_degree(),
      component = group_components(),
      closeness = centrality_closeness(),
      betweeness = centrality_betweenness()
    ) %>%
    activate(!!state)
}

#' Rectify vector representing triangle of pairwise values to a matrix
#' @param up_tri_vec vector containing values of upper triangle
#' @param elm_names names of the elements
rectify_pairs_vector <- function(up_tri_vec, elm_names, diag_val = 0) {
  # transform vector to matrix
  # see: https://github.com/zdk123/SpiecEasi/issues/17

  n <- length(elm_names)
  p <- rep(diag_val, n)

  X <- diag(p)
  X[upper.tri(X, diag = FALSE)] <- up_tri_vec
  X <- X + t(X)
  colnames(X) <- elm_names
  rownames(X) <- elm_names

  X
}

annotate_node_attributes_in_edges <- function(graph) {
  state <- graph %>% active()

  graph %>%
    activate(edges) %>%
    # ensure idempotency
    select(-matches("^(from|to)_")) %>%
    left_join(
      graph %>%
        activate(nodes) %>%
        as_tibble() %>%
        mutate(name = row_number()) %>%
        rename_all(~ .x %>% paste0("from_", .)),
      by = c("from" = "from_name")
    ) %>%
    left_join(
      graph %>%
        activate(nodes) %>%
        as_tibble() %>%
        mutate(name = row_number()) %>%
        rename_all(~ .x %>% paste0("to_", .)),
      by = c("to" = "to_name")
    ) %>%
    activate(!!state)
}

#' Tidy symmetric correlation matrix
#' @param mat symmetric square matrix with row names
tidy_matrix <- function(mat) {
  if (!isSymmetric(mat)) stop("Matrix must be symmetric")

  mat %>%
    as_tibble(rownames = "column1") %>%
    pivot_longer(cols = -column1, names_to = "column2") %>%
    rowwise() %>%
    mutate(
      comparison = c(column1, column2) %>% sort() %>% paste0(collapse = "-")
    ) %>%
    group_by(comparison) %>%
    slice(1) %>%
    ungroup() %>%
    select(-comparison)
}

get_graph <- function(cor_res, pooling_col, lineages_tbl, min_abs_estimate = 0) {
  taxa <- c(
    "kingdom" = 8, "phylum" = 7, "class" = 6, "order" = 5, "family" = 4,
    "genus" = 3, "species" = 2, "otu" = 1
  )
  selected_downstream_cols <- taxa[taxa < taxa[pooling_col]]

  graph <- switch(class(cor_res),
    "pulsar.refit" = get_graph_pulsar(
      cor_res = cor_res,
      min_abs_estimate = min_abs_estimate
    ),
    "rcorr" = {
      cor_res %>%
        tidy() %>%
        dplyr::rename(from = column1, to = column2) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(q.value = p.value %>% p.adjust(method = "fdr")) %>%
        dplyr::filter(q.value < 0.05 & abs(estimate) >= min_abs_estimate) %>%
        tidygraph::as_tbl_graph(directed = FALSE) %>%
        dplyr::mutate(taxon = name, name = row_number())
    },
    "spiec_easi_sparcc_res" = {
      cor_tbl <-
        cor_res$pval$cors %>%
        rectify_pairs_vector(taxa, diag_val = 1) %>%
        tidy_matrix() %>%
        dplyr::filter(abs(estimate) >= min_estimate) %>%
        dplyr::rename(estimate = value)

      pval_tbl <-
        cor_res$pval$pvals %>%
        rectify_pairs_vector(taxa, diag_val = 0) %>%
        tidy_matrix() %>%
        rename(p.value = value)

      cor_tbl %>%
        inner_join(pval_tbl, by = c("column1", "column2")) %>%
        rename(from = column1, to = column2) %>%
        ungroup() %>%
        mutate(q.value = p.value %>% p.adjust(method = "fdr")) %>%
        filter(p.value < 0.05) %>%
        filter(abs(estimate) >= min_abs_estimate) %>%
        as_tbl_graph(directed = FALSE) %>%
        mutate(taxon = name, name = row_number())
    },
    error = function(e) stop("not implemented")
  )

  graph %>%
    activate(nodes) %>%
    rename(!!pooling_col := taxon) %>%
    left_join(lineages_tbl, by = pooling_col) %>%
    annoate_node_attributes_in_edges() %>%
    activate(nodes) %>%
    select(-name) %>%
    filter_and_topologize_graph(min_abs_estimate = min_abs_estimate)
}

annotate_edges <- function(graph, ...) {
  graph %>%
    activate(edges) %>%
    mutate(...)
}

plot_graph <- function(tbl_graph) {
  tbl_graph %>%
    # activate(nodes) %>%
    # mutate(betweeness = betweeness %>% as.double()) %>%
    ggraph(layout = "igraph", algorithm = "circle") +
    geom_edge_link(aes(color = estimate), size = 6) +
    geom_node_point(aes(color = kingdom, alpha = betweeness, size = closeness)) +
    #  scale_edge_color_distiller(palette = "vikO") +
    scale_edge_color_gradient2(
      high = "#770000", low = "#000077",
      midpoint = 0
    ) +
    theme_void() +
    coord_fixed()
}
