scale_min_max <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}
