#' @param s string like 29.55149167 N 118.9075972 E
parse_lon_lat <- function(s) {
  parsed <- s %>%
    str_split(" ") %>%
    purrr::simplify()
  lat <- parsed[1] %>%
    as.numeric() %>%
    ifelse(parsed[2] == "N", ., -.)
  lon <- parsed[3] %>%
    as.numeric() %>%
    ifelse(parsed[4] == "E", ., -.)

  list(lat = lat, lon = lon)
}
