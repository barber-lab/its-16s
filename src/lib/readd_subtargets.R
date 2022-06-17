readd_subtargets <- function(pattern, n = Inf) {
  pattern <- pattern %>% paste0("_[a-f0-9]+$")
  
  drake::cached() %>%
    keep(~ .x %>% str_detect(pattern = pattern)) %>%
    discard(is.na) %>%
    head(n) %>%
    map(~ possibly(readd,tibble())(.x, character_only = TRUE)) %>%
    bind_rows()
}

