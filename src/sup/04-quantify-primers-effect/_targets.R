#!/usr/bin/env R


work_dir <- "04-quantify-primers-effect"
source(here::here(paste("src", work_dir, "defaults.R", sep = "/")))


tar_option_set(packages = c("tidyverse",
                            "tarchetypes",
                            "reshape",
                            "patchwork"))

list(
  tar_target(
    valid_samples,
    read_csv("raw/valid_samples/samples_used1580.csv")
  ),
  tar_target(
    bioproject_primers,
    read_csv("raw/primers/primers_comparison_16s_its.csv")
  ),
  tar_target(
    runs,
    read_csv("raw/runs.csv")
  ),
  tar_target(
    raw_abundance_bac,
    read_csv("raw/abundances/bacteria.genus.raw.csv")
  ),
  tar_target(
    raw_abundance_fun,
    read_csv("raw/abundances/fungi.genus.raw.csv")
  ),
  tar_target(
    biosamples,
    read_tsv("raw/biosamples.tsv") %>% 
      select(bioproject_id, environment_group) %>% 
      distinct() %>% 
      drop_na(environment_group)
  ),
  tar_target(
    metadata,
    annotate_valid_samples(valid_samples, bioproject_primers, biosamples)
  ),
  tar_target(
    bacteria_by_primers,
    group_primers(metadata, raw_abundance_bac, "bac"),
    deployment = "main"
  ),
  tar_target(
    bacterial_commonalities,
    intersection_between_primers(bacteria_by_primers)
  ),
  tar_target(
    correlate_16s_primers,
    correlate_primers(bacterial_commonalities),
    pattern = map(bacterial_commonalities)
  ),
  tar_target(
    bacterial_primers_sim,
    calculate_similarity(correlate_16s_primers),
    pattern = map(correlate_16s_primers)
  ),
  tar_target(
    bacterial_primers_heat,
    primers_heatmap(bacterial_primers_sim, "bacterial-primers-group-corr-env", "bacteria")
  ),
  tar_target(
    fungi_by_primers,
    group_primers(metadata, raw_abundance_fun, "fun"),
    deployment = "main"
  ),
  tar_target(
    fungal_commonalities,
    intersection_between_primers(fungi_by_primers)
  ),
  tar_target(
    correlate_18s_primers,
    correlate_primers(fungal_commonalities),
    pattern = map(fungal_commonalities)
  ),
  tar_target(
    fungal_primers_sim,
    calculate_similarity(correlate_18s_primers),
    pattern = map(correlate_18s_primers)
  ),
  tar_target(
    fungal_primers_heat,
    primers_heatmap(fungal_primers_sim, "fungal-primers-group-corr-env", "fungi")
  ),
  tar_target(
    fungal_primers_heat_rds,
    primers_heatmap_rds(fungal_primers_sim, "fungi")
  ),
  tar_target(
    bacterial_primers_heat_rds,
    primers_heatmap_rds(bacterial_primers_sim, "bacteria")
  )
)