#!/usr/bin/env R

# This pipeline works in three stages, it means that some results need to be
# generated outside R and then used.

# 1st stage: Collect genomes metadata
# - Ends in the target genomes2find_BGCs

# 2nd stage: Get BGCs
# - Ends in the target bgcs_plot_comparison

# 3rd stage: Get CDS and genome information 
# - Ends in the target export_cds_genome_plots

source(here::here("src/03-generalists-specialists-genome/defaults.R"))

# Need to install qs package for efficient storage
tar_option_set(packages = c("tidyverse","tarchetypes","jsonlite","ggpubr", "patchwork"),
               format = "qs", 
               memory = "transient", 
               garbage_collection = TRUE, 
               storage = "worker", 
               retrieval = "worker")

list(
  tar_target(
    generalist_genus,
    read_tsv("raw/generalists/generalist.txt")
  ),
  tar_target(
    specialist_genus,
    read_tsv("raw/specialists/specialists.txt")
  ),
  tar_target(
    lineages,
    read_csv("raw/lineages.csv")
  ),
  tar_target(
    metadata,
    get_annotation(generalist_genus, specialist_genus, lineages)
  ),
  tar_target(
    generalist_list,
    get_files_list("cache/11-generalists-genome-metadata", ".jsonl$")
  ),
  tar_target(
    tidy_generalists,
    parse_jsonl(generalist_list)
  ),
  tar_target(
    generalists_pattern,
    generalists_no_genome(generalist_genus)
  ),
  tar_target(
    generalists_with_genome,
    filter_generalists(tidy_generalists, generalists_pattern)
  ),
  tar_target(
    random_generalist_genomes,
    shuffle_genome(generalists_with_genome, metadata)
  ),
  tar_target(
    specialist_list,
    get_files_list("cache/12-specialists-genome-metadata", ".jsonl$")
  ),
  tar_target(
    tidy_specialists,
    parse_jsonl(specialist_list)
  ),
  tar_target(
    specialists_with_genome,
    filter_specialists(tidy_specialists)
  ),
  tar_target(
    random_specialist_genomes,
    shuffle_genome(specialists_with_genome, metadata)
  ),
  tar_target(
    failed_genomes,
    remove_genome_id()
  ),
  tar_target(
    genomes2find_BGCs,
    export_genome_list(failed_genomes, random_generalist_genomes, random_specialist_genomes, lineages)
  ),
  tar_target(
    list_bgcs_files,
    antismash_files("cache/14-antismash-results")
  ),
  tar_target(
    import_bgcs,
    parse_antismash(list_bgcs_files),
    pattern = map(list_bgcs_files),
    deployment = "main"
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
    bgcs_plot_comparison,
    plot_bgcs(generalist_specialists_bgcs)
  ),
  tar_target(
    genome_gff,
    get_files_list("cache/13-generalists-specialists-gff3", ".gff$")
  ),
  tar_target(
    cds_by_genome,
    count_cds(genome_gff),
    pattern = map(genome_gff)
    # ,
    # deployment = "main"
  ),
  tar_target(
    genome_fasta,
    get_files_list("cache/13-generalists-specialists-genome", ".fna$")
  ),
  tar_target(
    genome_length,
    get_genome_size(genome_fasta),
    pattern = map(genome_fasta)
    # ,
    # deployment = "main"
  ),
  tar_target(
    export_cds_genome_plots,
    plot_cds_gen_metrics(genomes2find_BGCs, cds_by_genome, genome_length)
  ),
  tar_target(
    plot_ratio,
    ratio_bgc_cds(generalist_specialists_bgcs, cds_by_genome)
  )
)

