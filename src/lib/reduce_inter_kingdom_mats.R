#'@param mats list of matricies which should be combined into one matrix
reduce_inter_kingdom_mats <- function(mats) {
  if(class(mats)[[1]] == "matrix") {
    mats
  } else {
    mats %>% map(as.data.frame) %>% reduce(bind_cols) %>% as.matrix()
  }
}
