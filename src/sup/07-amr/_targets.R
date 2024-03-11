#!/usr/bin/env R

work_dir <- "07-amr"
Sys.setenv(TAR_PROJECT = "amr_detection")
source(here::here(paste("src", work_dir, "defaults.R", sep = "/")))

# Need to install qs package for efficient storage
tar_option_set(packages = c("tidyverse","tarchetypes","jsonlite", "patchwork", "ggpubr"),
               format = "qs", 
               memory = "transient", 
               garbage_collection = TRUE, 
               storage = "worker", 
               retrieval = "worker")

### Targets results dir
gen_spe_tar_dir <- "cache/01-r-generalist_specialist_genome_metadata"

### Import results from other targets that were already analised
genomes2find_BGCs <- tar_read(name = genomes2find_BGCs, store = gen_spe_tar_dir)


library(future)
library(future.callr)

plan(callr)

list(
  tar_target(
    amr_results_files,
    list_amr_files("cache/16-amr-gff")
  ),
  tar_target(
    get_amr,
    import_amr_results(amr_results_files),
    pattern = map(amr_results_files)
  ),
  tar_target(
    valid_amr,
    tidy_amr(get_amr)
  ),
  tar_target(
    amr_comparisons_plot,
    amr_plots(valid_amr, genomes2find_BGCs)
  ),
  tar_target(
    amr2sup_file,
    export_amr(valid_amr)
  )
)

