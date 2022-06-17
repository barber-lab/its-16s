#' Count abundant samples
#' @param  abundance data frame with columns sample_id, taxon and abundance
#' @param samples data frame with columns sample_id and bioproject_id and the one provided in `samples_grouping`
#' @param bioprojects data frame with columns bioproject_id
#' @param samples_grouping column name to group samples
count_abundant_samples <- function(abundance, samples, bioprojects, samples_grouping = "bioproject_id") {
  abundance %>%
    left_join(samples %>% select_at(c("sample_id", "bioproject_id", samples_grouping))) %>%
    left_join(bioprojects) %>%
    select_at(c("sample_id", samples_grouping)) %>%
    distinct() %>%
    group_by_at(samples_grouping) %>%
    count() %>%
    rename(abundant = n)
}
