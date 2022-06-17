as_matrix <- function(data, rownames_col = NULL) {
  if(! is.null(rownames_col)) {
    res <-
      data %>% 
      dplyr::select(-!!rownames_col) %>%
      as.matrix()
    
    rownames(res) <- data[[rownames_col]]
  } else {
    res <- data %>% as.matrix()
  }
  
  res
}