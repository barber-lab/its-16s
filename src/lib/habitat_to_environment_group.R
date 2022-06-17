habitat_to_environment_group <- function(habitat) {
  env_habitats <- list(
    soil = c(
      "soil",
      "temperate forests", 
      "conifer forests", 
      "mediterranean forests", 
      "boreal forests", 
      "tropical forests"
    ),
    aquatic = c(
      "aquatic",
      "marine",
      "other freshwater",
      "other aquatic",
      "large lakes"
    ),
    host = c(
      "host",
      "other hosts",
      "insects",
      "lung",
      "gut",
      "plants",
      "skin",
      "mouth"
    )
  )
  
  env_habitats %>% simplify() %>% keep(~ .x == habitat) %>% names() %>% str_remove("[0-9]+$")
}
