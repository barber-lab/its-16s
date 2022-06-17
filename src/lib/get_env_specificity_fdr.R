#' Calculate p value in enviornment group coabudnance
#'
#' Permutation: Enviornment group labels were shuffeld across samples
#' Coabundances were calculated for each shuffeled environment
#' Two-sided FDR P value is calculated based on the number of permutations with a estimate of at leas a given threshold
#' @see: https://www.nature.com/articles/s41467-020-17840-y#Sec2
get_env_specificity_fdr <- function(permuted_coabundance_res, from, to, estimate) {
  fdr_counts <-
    permuted_coabundance_res %>%
    unnest(cor_res) %>%
    rename(permuted_from = from, permuted_to = to, permuted_estimate = estimate) %>%
    filter(permuted_from == from & permuted_to == to) %>%
    # tow sided test
    mutate(false_discovery = abs(permuted_estimate) >= abs(estimate)) %>%
    group_by(false_discovery) %>%
    count() %>%
    ungroup() %>%
    complete(false_discovery = c(TRUE, FALSE), fill = list(n = 0)) %>%
    deframe()

  case_when(
    # no false discovery found, so p value is zero
    sum(fdr_counts) == 0 ~ 0,

    # default case
    TRUE ~ fdr_counts[["TRUE"]] / sum(fdr_counts)
  )
}
