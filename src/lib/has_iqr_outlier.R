#' IQR method to see if there are outlieres in a vector
#' @seealso https://www.nature.com/articles/s41467-020-17840-y#Sec10
has_iqr_outlier <- function(x) {
  quantiles <- quantile(x)
  b1 <- quantiles["0%"]
  b4 <- quantiles["100%"]
  iqr2 <- (quantiles["75%"] - quantiles["25%"]) * 2
  
  if (quantiles["25%"] - b1 > iqr2) {
    return(TRUE)
  }
  if (b4 - quantiles["75%"] > iqr2) {
    return(TRUE)
  }
  return(FALSE)
}