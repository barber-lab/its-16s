#' @param data data frame of covariates used to control correlation estiamate
#' @param estimate non partial correlation effect size
#' @param tol tolerance to calculate the pseudoinverse as tol in cor2pcor
partial_correlation <- function(from_taxon, to_taxon, estimate, covariates = c("lon", "lat", "collection_datetime"), data) {
  dummy_cols <-
    data %>%
    select_at(covariates) %>%
    keep(is.character) %>%
    names()
  
  covariates_tbl <-
    data %>%
    select(from_taxon = !!from_taxon, to_taxon = !!to_taxon, !!covariates) %>%
    fastDummies::dummy_cols(dummy_cols) %>%
    select(-dummy_cols) %>%
    # do not mutate taxa abundances
    group_by(from_taxon, to_taxon) %>%
    # correlation needs numeric variables
    mutate_at(vars(-group_cols()), as.double) %>%
    ungroup()
  
  cor_mat <-
    covariates_tbl %>%
    colnames() %>%
    combn(m = 2) %>%
    t() %>%
    as_tibble() %>%
    mutate(estimate = V1 %>%
             map2_dbl(V2, ~ cor.test(covariates_tbl[[.x]], covariates_tbl[[.y]], method = "pearson")$estimate)) %>%
    # assume no correlation e.g. if meta data is missing to do pearson correlation
    replace_na(list(estimate = 0)) %>%
    { x <- .; x %>% bind_rows(x %>% rename(V1 = V2, V2 = V1))} %>%
    pivot_wider(names_from = V2, values_from = estimate) %>%
    set_rownames(.$V1) %>%
    select(-V1) %>%
    as.matrix()
  
  cor_mat <- cor_mat[rownames(cor_mat) %>% sort(), colnames(cor_mat) %>% sort()] # sort by names
  diag(cor_mat) <- 1 # autocorrelation
  cor_mat[["from_taxon", "to_taxon"]] <- estimate # non partial correlation
  cor_mat[["to_taxon", "from_taxon"]] <- estimate # non partial correlation
  
  # try to get partial correlation using minimal tolerance possible
  tol = max(dim(cor_mat))*max(cor_mat)*.Machine$double.eps
  options(warn=-1)
  while(TRUE) {
    partial_cor_mat <- cor_mat %>% corpcor::cor2pcor(tol = tol)

    if(! partial_cor_mat %>% is.na() %>% any()) break
    tol <- tol * 2
  }
  options(warn=0)
  
  partial_cor_mat <- cor_mat %>% corpcor::cor2pcor(tol = tol)
  rownames(partial_cor_mat) <- rownames(cor_mat)
  colnames(partial_cor_mat) <- colnames(cor_mat)
  
  list(
    estimate = partial_cor_mat[["from_taxon", "to_taxon"]],
    tol = tol
  )
}
