#!/usr/bin/env Rscript

options(
  Ncpus = parallel::detectCores()
)

install.packages(c("devtools", "tidyverse"))

library(tidyverse)

get_requirements <- function(path) {
  path %>%
    readr::read_lines() %>%
    purrr::map_chr(~ .x %>%
      stringr:::str_remove("#.*$|^$") %>%
      stringr::str_trim()) %>%
    purrr::discard(~ .x %>% stringr::str_detect("^$"))
}

# fix dependency of jaccard: Must be installed first
devtools::install_bioc("3.12/qvalue")

# install from CRAN
get_requirements("src/install/requirements/cran.txt") %>%
  devtools::install_cran(upgrade_dependencies = FALSE)

# install from Bioconductor
get_requirements("src/install/requirements/bioconductor.txt") %>%
  devtools::install_bioc(upgrade_dependencies = FALSE)

# install from GitHub
get_requirements("src/install/requirements/git.txt") %>%
  devtools::install_git(upgrade_dependencies = FALSE)

# install versions
devtools::install_version("vctrs", "0.3.6", repos = "https://cloud.r-project.org", upgrade = "never")
devtools::install_version("ggplot2", "3.3.2", upgrade = "never")

# Fix: WGCNA dependencies not available using devtools::install_bioc
install.packages("BiocManager")
BiocManager::install("WGCNA", update = FALSE)
install.packages("DGCA")

# Install complexHeatmap afterwards
devtools::install_bioc("3.11/ComplexHeatmap")