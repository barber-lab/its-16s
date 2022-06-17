#' Supress any output of an expression
hush <- function(code) {
  sink("NUL") # use /dev/null in UNIX
  tmp <- code
  sink()
  return(tmp)
}
