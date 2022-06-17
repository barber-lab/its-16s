#' Get list to translate flatten bioproject_ids to merged ones
#' Some studies have multiple bioprojects
get_merged_bioprojects_id_dict <- function() {
  # some studies have separate bioprojects for fungi and bacteria
  merged_bioprojects <-
    read_csv("metadata/samples.csv") %>%
    pull(bioproject_id) %>%
    unique() %>%
    keep(~ .x %>% str_detect("_"))

  merged_bioprojects_flatten <-
    merged_bioprojects %>%
    map(~ .x %>%
      str_split("_") %>%
      pluck(1)) %>%
    purrr::simplify()

  get_merged_bioproject_id <- function(flat_bioproject_id) {
    merged_bioprojects[merged_bioprojects %>% str_detect(flat_bioproject_id)]
  }

  merged_bioprojects_id_dict <-
    merged_bioprojects_flatten %>%
    set_names(merged_bioprojects_flatten) %>%
    map_chr(get_merged_bioproject_id)

  merged_bioprojects_id_dict
}
