#' Filter samples with highly correlated abundance profile
#'
#' @param max_abs_estimate max aboslute correlation coefficient allowed between two samples
#' @param data data frame with columns sample_id, taxon and abundance
filter_samples_correlation <- function(data, max_abs_estimate = 0.95, method = "spearman") {
    removed_samples <- 
      data %>%
      select(sample_id, taxon, abundance) %>%
      pivot_wider(
        id_cols = taxon,
        names_from = sample_id,
        values_from = abundance,
        values_fill = list(abundance = 0)
      ) %>%
      select(-taxon) %>%
      as.matrix() %>%
      cor(method = method) %>%
      as_tibble(rownames = "from") %>%
      pivot_longer(-from, names_to = "to", values_to = "estimate") %>%
      # get too high correlations
      filter(from != to & abs(estimate) > max_abs_estimate) %>%
      # remove one sample of the pair
      pull(from) %>%
      unique()
    
     data %>% filter(! sample_id %in% removed_samples)
}
