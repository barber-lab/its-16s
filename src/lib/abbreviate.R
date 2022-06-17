abbreviate <- function(s, n = 2) {
  s %>% str_split("[ -]") %>% simplify() %>% str_sub(1, n) %>% paste0(sep = ".", collapse = " ")
}
