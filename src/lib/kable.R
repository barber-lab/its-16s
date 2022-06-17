#' @param remove_uniform_columns remove columns which have the same value in every row
kable <- function(data, remove_uniform_columns = TRUE, ...) {
  if (remove_uniform_columns) {
    data %<>% discard(~ .x %>%
      unique() %>%
      length() <= 1)
  }

  data %>%
    dplyr::select_if(~ !is.list(.x)) %>%
    kableExtra::kable(...) %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover"))
}
