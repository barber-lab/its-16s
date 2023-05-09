#!/usr/bin/env R

source(here::here("src/00-aln-metrics/defaults.R"))

tar_option_set(packages = c("tidyverse",
                            "tarchetypes"))

list(
  tar_target(
    bed_files,
    list_files("cache/10-bam2bed")
  ),
  tar_target(
   import_bed,
   mapped_reads_from_bed(bed_files),
   deployment = "main"
   ),
  tar_target(
    reads_length_files,
    list_files("cache/09-reads-length")
  ),
  tar_target(
    get_reads_length,
    import_reads_length(reads_length_files),
    deployment = "main"
  ),
  tar_target(
    coverage_files,
    list_files("cache/08-coverage")
  ),
  tar_target(
    import_aln_coverage,
    aln_covered_bases(coverage_files),
    deployment = "main"
  ),
  tar_target(
    runs,
    read_csv("raw/runs.csv")
  ),
  tar_target(
    reads_coverage,
    merge_metrics(import_bed, get_reads_length, import_aln_coverage, runs),
    deployment = "main"
  ),
  tar_target(
    coverage_by_reads_16s,
    plot_reads_coverage(reads_coverage),
    deployment = "main"
  ),
  tar_target(
    coverage_by_regerence_16s,
    plot_reference_coverage(import_aln_coverage, runs),
    deployment = "main"
  ),
  tar_target(
    plot_its_primers,
    extracted_its("cache/00-tmp/extracted_its_seqs_NAs_primers.txt", runs),
    deployment = "main"
  )
)
