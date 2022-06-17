#' Adverb to run function in parallel by setting mc.cores
#' @param f function to call
#' @param cores number of cores used to set option mc.cores
parallelly <- function(.f, cores = 1, ...) {
  .f <- as_mapper(.f)
  force(cores)

  function(...) {
    old.mc.cores <- getOption("mc.cores")
    options(mc.cores = cores)
    .f(...)
    options(mc.cores = old.mc.cores)
  }
}
