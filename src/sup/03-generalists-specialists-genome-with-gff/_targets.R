#!/usr/bin/env R

# This pipeline analyse the antismash results with gff3 file suplied together with the genome.

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
amr_dir         <- "cache/01-r-07-amr"

### Import results from other targets that were already analised
genomes2find_BGCs <- tar_read(name = genomes2find_BGCs, store = gen_spe_tar_dir)
generalist_genus  <- tar_read(name = generalist_genus, store = gen_spe_tar_dir)
specialist_genus  <- tar_read(name = specialist_genus, store = gen_spe_tar_dir)
genome_gff        <- tar_read(name = genome_gff, store = gen_spe_tar_dir)
cds_by_genome     <- tar_read(name = cds_by_genome, store = gen_spe_tar_dir)
genome_length     <- tar_read(name = genome_length, store = gen_spe_tar_dir)
valid_amr         <- tar_read(name = valid_amr, store = amr_dir)

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
  ),
  tar_target(
    genome_length2,
    import_genome_length("cache/17-genome-lenght")
  ),
  tar_target(
    paper_main_fig2,
    plot_main_fig2_paper(generalist_specialists_bgcs, genomes2find_BGCs, cds_by_genome, genome_length, valid_amr)
  ),
  tar_target(
    cds2permutate,
    cds_to_plot(genomes2find_BGCs, cds_by_genome, genome_length)
  ),
  tar_target(
    permutate_cds,
    cds2permutate %>% 
      filter(kingdom == "Bacteria") %>% 
      permutation_test(y = "cds_n", x = "Group", n_iter = 10000)
  ),
  tar_target(
    permutate_cds_norm,
    cds2permutate %>% 
      filter(kingdom == "Bacteria") %>% 
      permutation_test(y = "cds_norm", x = "Group", n_iter = 10000)
  ),
  tar_target(
    bcgs2permutate,
    bgcs_df_paper(generalist_specialists_bgcs, genome_length)
  ),
  tar_target(
    permutate_bcgs,
    bcgs2permutate %>% 
      filter(kingdom == "Bacteria") %>% 
      permutation_test(y = "value", x = "Groups", n_iter = 10000)
  ),
  tar_target(
    amr2permutate,
    amr_df_paper(valid_amr, genomes2find_BGCs, genome_length)
  ),
  tar_target(
    permutate_amr,
    amr2permutate %>% 
      filter(kingdom == "Bacteria") %>% 
      permutation_test(y = "norm_amr", x = "Groups", n_iter = 10000)
  ),
  tar_target(
    compare_bcgs_products,
    compare_bcgs_products_bac_gen(generalist_specialists_bgcs, genome_length, bcgs_groups)
  )
)