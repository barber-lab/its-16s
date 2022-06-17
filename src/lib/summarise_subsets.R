#' @param  data tibble with columns subset_name, subset_value and data
#' @param sep character used to split multiple groups of subset_name and subset_value
#' @param grouping one value of subset_name describing grouping used to summarise
summarise_subsets <- function(data, grouping = "environment_group", sep = ",") {
  get_group_pos <- function(x, grouping) {
    x %>%
      str_split(sep) %>%
      purrr::simplify() %>%
      map_lgl(~ .x == grouping) %>%
      which(TRUE)
  }

  clean_commas <- function(x) {
    x %>%
      str_replace_all(",+", ",") %>%
      str_remove_all("(^,)|(,$)")
  }

  start_groups <- data %>%
    group_vars() %>%
    setdiff(c("subset_name", "subset_value"))

  data %>%
    mutate(
      group = subset_name %>% map2_chr(
        subset_value, ~ {
          group_pos <- get_group_pos(.x, grouping)
          .y %>%
            str_split(sep) %>%
            purrr::simplify() %>%
            possibly(extract2, NA)(group_pos)
        }
      ),
      subset_name = subset_name %>% str_remove(grouping) %>% clean_commas(),
      subset_value = subset_value %>% str_remove(group) %>% clean_commas()
    ) %>%
    unnest(data) %>%
    group_by_at(c(start_groups, "subset_name", "subset_value")) %>%
    nest() %>%
    arrange_columns()
}
