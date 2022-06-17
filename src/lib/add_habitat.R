add_habitat <- function(samples, min_samples_per_habitat = 10) {
  ecoregions <-
    c(
      Freshwater = "raw/databases/DX4-ecoregions/Freshwater_Ecoregions.shp",
      Terrestrial = "raw/databases/DX4-ecoregions/Terrestrial_Ecoregions.shp"
    ) %>%
    map(read_ecoregion) %>%
    map(~ .x %>% mutate(eco_region = eco_region %>% str_to_lower()))

  res_tbl <-
    samples %>%
    bind_rows(
      # host samples
      samples %>%
        filter(environment_group == "host") %>%
        mutate(
          habitat = case_when(
            biosample_scientific_name == "endophyte metagenome" ~ "endophytes",

            biosample_scientific_name %in% c("ant metagenome", "insect metagenome") ~ "insects",

            biosample_scientific_name %>% str_detect("gut") ~ "gut",
            biosample_scientific_name %>% str_detect("lung") ~ "lung",
            biosample_scientific_name %>% str_detect("oral") ~ "mouth",
            bioproject_id == "PRJNA473079" ~ "skin", # scalp


            biosample_scientific_name %in% c("phyllosphere metagenome") ~ "plants",
            bioproject_id == "PRJEB23282" ~ "plants", # pollen

            TRUE ~ "other hosts"
          )
        ),

      # spatial samples
      samples %>%
        filter(environment_group != "host") %>%
        filter(!is.na(lon) & !is.na(lat)) %>%
        group_by(environment_group) %>%
        nest() %>%
        mutate(
          eco_type = environment_group %>% recode(aquatic = "Freshwater", soil = "Terrestrial"),
          data = data %>% map2(eco_type, ~ {
            .x %>%
              # kepp copy of columns lat and lon
              mutate(lon2 = lon, lat2 = lat) %>%
              st_as_sf(coords = c("lon2", "lat2"), crs = st_crs(4326)) %>%
              mutate(eco_type = .y) %>%
              st_join(ecoregions[[.y]], by = "eco_type") %>%
              mutate(habitat = eco_region)
          })
        ) %>%
        unnest(data),

      # other samples
      samples
    ) %>%
    # Remove duplicate sanmples
    group_by(sample_id) %>%
    arrange(habitat) %>%
    slice(1) %>%
    ungroup()

  too_small_habitats <-
    res_tbl %>%
    group_by(habitat) %>%
    count() %>%
    filter(n < min_samples_per_habitat) %>%
    pull(habitat)

  res_tbl %>%
    mutate(
      habitat = case_when(
        !is.na(habitat) ~ habitat,
        
        environment_group == "aquatic" & biosample_scientific_name == "freshwater metagenome" ~ "other freshwater",
        environment_group == "aquatic" & biosample_scientific_name == "riverine metagenome" ~ "river",
        environment_group == "aquatic" & biosample_scientific_name == "marine metagenome" ~ "marine",
        
        TRUE ~ habitat
      ),
      
      habitat = habitat %>% recode(
        "temperate broadleaf and mixed forests" = "temperate forests",
        "temperate conifer forests" = "conifer forests",
        "mediterranean forests, woodlands and scrub" = "mediterranean forests",
        "boreal forests/taiga" = "boreal forests",
        "tropical and subtropical moist broadleaf forests" = "tropical forests"
      ),
      
      habitat = case_when(
        habitat %in% too_small_habitats ~ str_glue("other {environment_group}"),
        is.na(habitat) ~ str_glue("other {environment_group}"),
        TRUE ~ habitat
      ),
      habitat = habitat %>% str_remove("[ ]?metagenome$")
    ) %>%
    # filter habitats with too few samples
    group_by(habitat) %>%
    filter(n() >= min_samples_per_habitat) %>%
    ungroup()
}
