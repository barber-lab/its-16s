#' Adopts missing attributs from y to x
#' @param which specifiy which attributes of y to add. Default: all
#' @return x with attributes added from y
add_missing_attributes <- function(x, y, which = NA) {
  if (is.na(x)) {
    return(x)
  }

  attrs_x <- base::attributes(x) %>% names()
  attrs_y <- base::attributes(y) %>% names()
  attr_to_add <- attrs_y %>% setdiff(attrs_x)

  if (!is.na(which)) attr_to_add %<>% keep(~ .x %in% which)

  for (a in attr_to_add) {
    base::attr(x, a) <- base::attributes(y)[[a]]
  }
  x
}
