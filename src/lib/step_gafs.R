#' Feature selection using Genetic Algorithms
step_select_gafs <- function(
  recipe,
  ...,
  outcome = NULL,
  role = "predictor",
  trained = FALSE,
  exclude = NULL,
  options = list(),
  res = NULL,
  skip = FALSE,
  id = recipes::rand_id("select_gafs")) {
  
  recipes::recipes_pkg_check("gafs")
  
  recipes::add_step(
    recipe,
    step_select_gafs_new(
      terms = recipes::ellipse_check(...),
      trained = trained,
      outcome = outcome,
      role = role,
      exclude = exclude,
      options = options,
      res = res,
      skip = skip,
      id = id
    )
  )
}

# wrapper around 'step' function that sets the class of new step objects
#' @importFrom recipes step
step_select_gafs_new <- function(terms, role, trained, outcome, exclude,
                                   options, res, skip, id) {
  recipes::step(
    subclass = "select_gafs",
    terms = terms,
    role = role,
    trained = trained,
    outcome = outcome,
    exclude = exclude,
    options = options,
    res = res,
    skip = skip,
    id = id
  )
}

#' @export
prep.step_select_gafs <- function(x, training, info = NULL, ...) {
  
  # translate the terms arguments
  x_names <- recipes::terms_select(terms = x$terms, info = info)
  y_name <- recipes::terms_select(x$outcome, info = info)
  y_name <- y_name[1]
  
  if (length(x_names) > 0) {
    
    call <- rlang::call2(
      .fn = "gafs",
      .ns = "gafs",
      x = rlang::quo(training[, x_names]),
      y = rlang::quo(training[[y_name]]),
      !!!x$options
    )
    
    res <- rlang::eval_tidy(call)
    
    exclude <- names(res$finalDecision[res$finalDecision == "Rejected"])
    
  } else {
    exclude <- character()
  }
  
  step_select_gafs_new(
    terms = x$terms,
    trained = TRUE,
    role = x$role,
    outcome = y_name,
    exclude = exclude,
    options = x$options,
    res = res,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_select_gafs <- function(object, new_data, ...) {
  if (length(object$exclude) > 0) {
    new_data <- new_data[, !colnames(new_data) %in% object$exclude]
  }
  as_tibble(new_data)
}

#' @export
print.step_select_gafs <- function(x, width = max(20, options()$width - 30), ...) {
  cat("gafs feature selection")
  
  if(recipes::is_trained(x)) {
    n <- length(x$exclude)
    cat(paste0(" (", n, " excluded)"))
  }
  cat("\n")
  
  invisible(x)
}

#' @rdname step_select_gafs
#' @param x A `step_select_gafs` object.
#' @export
tidy.step_select_gafs <- function(x, ...) {
  if (recipes::is_trained(x)) {
    res <- tibble(terms = x$exclude)
  } else {
    term_names <- recipes::sel2char(x$terms)
    res <- tibble(terms = rlang::na_chr)
  }
  res$id <- x$id
  res
}
