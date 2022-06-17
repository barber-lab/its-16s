#' Median alpha diverity in samples lacking specific taxa
#' 
#' Do not need to calculate alpha diversity of samples having these specific taxa,
#' because the sum of alpahdiversities of samples with and without those taxa is constant
#'
#' @param taxa selected samples must lack these taxa
#' @param permute remove random taxa (as many as given in the argument taxa) instead
get_median_alpha_div_of_samples_lacking_taxa <- function(permute = FALSE, seed = 1, taxa, sub_abundances, alphadiv) {
  if (permute) {
    set.seed(seed)
    
    taxa <-
      sub_abundances$data[[2]]$taxon %>%
      unique() %>%
      sample(length(taxa))
  }
  
  samples_having_any_generalists <-
    sub_abundances %>%
    filter(subset_name == "kingdom" & norm_method == "tss") %>%
    transmute(kingdom = subset_value, data) %>%
    select(-kingdom) %>%
    unnest(data) %>%
    mutate(has_generalist = taxon %in% taxa) %>%
    filter(abundance > 0) %>%
    group_by(sample_id, kingdom) %>%
    summarise(has_generalist = any(has_generalist))
  
  samples_having_any_generalists %>%
    inner_join(
      alphadiv %>%
        filter(subset_name == "kingdom") %>%
        select(kingdom = subset_value, data) %>%
        unnest(data) %>%
        select(-Simpson)
    ) %>%
    pivot_longer(cols = c(Chao1, Shannon), names_to = "alphadiv_metric") %>%
    group_by(kingdom, has_generalist, alphadiv_metric) %>%
    summarise(median_alpha_div = median(value))
}
