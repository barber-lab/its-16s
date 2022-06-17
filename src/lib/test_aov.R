#' Test significance using Analysis of Variance
#'
#' pre: Kruskal Wallis test
#' post: Dunn's test
#' @param formula formula describing test of type response ~ category
#' @param data data frame containing all columns used by formula
#' @param reject null hypothesis if p value of pre test
#'   is larger than this value
#' @param ... arguments parsed to `dunn.test`
#' @return tidy tibble of test results
test_aov <- function(formula, data, ...) {
  formula_parts <-
    formula %>%
    as.character() %>%
    stringr::str_split("~") %>%
    purrr::simplify() %>%
    purrr::discard(~ .x == "")

  response <- formula_parts[1]
  covariate <- formula_parts[2]
  covariates <- data %>%
    dplyr::pull(covariate) %>%
    unique()

  if (length(covariates) == 1) stop("Covariate must have at least two levels")

  if (length(covariates) == 2) {
    res <-
      wilcox.test(formula, data, alternative = "two.sided") %>%
      broom::tidy() %>%
      mutate(comparison = paste0(covariates, collapse = "-"))

    return(res)
  }

  pre_res <-
    kruskal.test(formula, data) %>%
    broom::tidy()

  post_res <-
    hush(
      dunn.test::dunn.test(
        x = data %>% purrr::pluck(response),
        g = data %>% purrr::pluck(covariate),
        ...
      )
    )

  post_res <-
    post_res %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(
      method = "Dunn's test",
      comparison = comparisons %>% stringr::str_replace_all(" - ", "-")
    ) %>%
    dplyr::rename(
      p.value = P,
      q.value = P.adjusted,
      statistic = Z
    ) %>%
    dplyr::select(-chi2) # remove info from pre test

  dplyr::bind_rows(pre_res, post_res)
}
