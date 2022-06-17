get_failed_traces <- function() {
  traces <-
    list.files("raw/etc/traces", pattern = "trace.tsv", full.names = TRUE) %>%
    map(read_tsv) %>%
    bind_rows()

  failed_traces <-
    traces %>%
    filter(status == "FAILED") %>%
    rowwise() %>%
    mutate(tag = tag %>% parse_json() %>% list()) %>%
    ungroup() %>%
    unnest_wider(tag)

  return(failed_traces)
}
