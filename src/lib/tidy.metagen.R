tidy.metagen <- function(x) {
  tibble(
    method = "Generic inverse variance meta-analysis",
    p.value = x$pval.Q,
    df = x$df.Q,
    statistic = x$Q
  )
}