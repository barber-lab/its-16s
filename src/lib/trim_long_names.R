trim_long_names <- function(x) {
  x %>%
  str_split("-") %>%
  purrr::simplify() %>%
  head(2) %>%
  paste0(collapse = "-") %>%
  
  str_split(" ") %>%
  purrr::simplify() %>%
  head(2) %>%
  paste0(collapse = " ")
}