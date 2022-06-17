silently <-function (f) {
  f <- as_mapper(f)
  function(...) {
    file <- file("/dev/null", open = "wt")
    sink(file, type = "output")
    sink(file, type = "message")
    f(...)
    sink()
    sink()
  }
}
