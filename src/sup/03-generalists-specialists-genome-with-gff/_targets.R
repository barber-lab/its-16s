#!/usr/bin/env R

# This pipeline analise the antismash results with gff3 file suplied together with the genome.

source(here::here("src/03-generalists-specialists-genome-with-gff/defaults.R"))

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
#bgcs_by_genome    <- tar_read(name = bgcs_by_genome, store = gen_spe_tar_dir)
genomes2find_BGCs <- tar_read(name = genomes2find_BGCs, store = gen_spe_tar_dir)
generalist_genus  <- tar_read(name = generalist_genus, store = gen_spe_tar_dir)
specialist_genus  <- tar_read(name = specialist_genus, store = gen_spe_tar_dir)
genome_gff        <- tar_read(name = genome_gff, store = gen_spe_tar_dir)
cds_by_genome     <- tar_read(name = cds_by_genome, store = gen_spe_tar_dir)
genome_length     <- tar_read(name = genome_length, store = gen_spe_tar_dir)

library(future)
library(future.callr)

plan(callr)

list(
  tar_target(
    list_bgcs_files,
    antismash_files("cache/14-antismash-results-gff")
  ),
  tar_target(
    import_bgcs,
    parse_antismash(list_bgcs_files),
    pattern = map(list_bgcs_files)
    # ,
    # deployment = "main"
  ),
  tar_target(
    bgcs_by_genome,
    group_bgcs(import_bgcs)
  ),
  tar_target(
    generalist_specialists_bgcs,
    generate_annotation(bgcs_by_genome, genomes2find_BGCs, generalist_genus, specialist_genus)
  ),
  tar_target(
    figure1,
    plot_bgcs_gen_size_cds(generalist_specialists_bgcs, genomes2find_BGCs, cds_by_genome)
  ),
  tar_target(
    figure2,
    ratio_bgc_cds(generalist_specialists_bgcs, cds_by_genome)
  ),
  tar_target(
    figure3,
    plot_cds(generalist_specialists_bgcs, cds_by_genome, genome_length)
  ),
  tar_target(
    bcgs_groups,
    read_tsv("raw/02-antismash-groups/antismash-bcgs-products.tsv")
  ),
  tar_target(
    figure4,
    plot_bcgs_gen_spec(generalist_specialists_bgcs, bcgs_groups)
  ),
  tar_target(
    export_bgcs_list,
    bgcs2tsv(generalist_specialists_bgcs)
  ),
  tar_target(
    bcgs_compare_means,
    bcgs_statistics(generalist_specialists_bgcs, bcgs_groups)
  )
)

