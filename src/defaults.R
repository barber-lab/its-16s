#!/usr/bin/env Rscript

library(latex2exp)
library(Hmsc)
library(themis)
library(recipes)
library(workflows)
library(tidymodels)
library(rsample)
#library(recipeselectors)
library(yardstick)
library(fastDummies)
library(tidygraph)
library(optparse)
library(meta)
library(ggsci)
library(igraph)
library(ggpubr)
library(pbmcapply)
library(here)
library(memoise)
library(yaml)
library(drake)
library(patchwork)
library(qs)
library(storr)
library(fpc)
library(lme4)
library(lubridate)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(phyloseq)
library(jsonlite)
library(scales)
library(clustermq)
library(dplyr)
library(magrittr)
library(purrr)
library(glue)
library(ggnewscale)
library(ComplexHeatmap)
library(vegan)
#library(FEAST)
library(ggalluvial)
library(ggfittext)
library(ggvenn)
library(ggfortify)
# vip needs predict functions in global NS
library(kknn)
library(glmnet)
library(ranger)
library(parsnip)
library(kernlab)
library(nnet)
library(readxl)
library(writexl)
library(ggtext)
# call tidyverse last to overwrite default namespace
library(tidyverse)

select <- dplyr::select
filter <- dplyr::filter
simplify <- purrr::simplify

results_path <- function(x) {
  paste0("/sbidata/dloos/prj/AX3-its-16s-stat/results/", x)
}

set.seed(1337)

# set working directory to project directory
# especially if R shell started in a sub directory
setwd(here::here())

Sys.setenv(PATH = paste0("/miniconda3/envs/fastspar/bin/:", Sys.getenv("PATH")))

c("src/lib", "src/plans/") %>%
  list.files(pattern = "\\.R$", full.names = TRUE) %>%
  walk(source)

jobs <- min(60, parallel::detectCores() - 1)

# Resitrict Linear Algebra solver threads to minimize overhead
#' @see https://github.com/hmsc-r/HMSC/issues/23
RhpcBLASctl::blas_set_num_threads(4)

# internal parallelization
# jobs %>%
#   parallel::makePSOCKcluster() %>%
#   doParallel::registerDoParallel()

weightsRef <- NULL # object must be provided for WGCNA::modulePreservation

options(
  mc.cores = jobs,
  knitr.kable.NA = "",
  future.globals.maxSize = 10 * 1024^3,
  clustermq.scheduler = "multicore"
)

plan <-
  list(
    "base" = get_base_plan(),
    "feast" = get_feast_plan(),
    "coabundance" = get_coabundance_plan(),
    "graph" = get_graph_plan(),
    "ordinations" = get_ordinations_plan(),
    "plots" = get_plots_plan(),
    "generalists_plots" = get_generalists_plots_plan(),
    "ml" = get_ml_plan(),
    "selbal" = get_selbal_plan(),
    "wgcna" = get_wgcna_plan(),
    "enrichment" = get_enrichment_plan(),
    "lmm" = get_lmm_plan(),
    "hmsc" = get_hmsc_plan(),
    "dgca" = get_dgca_plan(),
    "rarefaction" = get_rarefaction_plan(),
    "generalists" = get_generalists_plan(),
    "files" = get_files_plan()
  ) %>%
  map2(names(.), ~ .x %>% mutate(part = .y)) %>%
  bind_plans()

#
# Theme
#

scale_fill_discrete <- function(...) {
  ggsci::scale_color_npg(..., na.value = "grey")
}

scale_fill_continuous <- function(...) {
  ggplot2::scale_fill_viridis_c(..., na.value = "grey")
}

scale_color_discrete <- function(...) {
  ggsci:.scale_color_npg(..., na.value = "grey")
}

scale_color_continuous <- function(...) {
  ggplot2::scale_color_viridis_c(..., na.value = "grey")
}

scale_colour_discrete <- function(...) {
  ggsci::scale_colour_npg(..., na.value = "grey")
}

scale_colour_continuous <- function(...) {
  ggplot2::scale_colour_viridis_c(..., na.value = "grey")
}

environment_group_colors <- c(
  "soil" = "#8b786d",
  "aquatic" = "#78a1bb",
  "host" = "#e0aa5a",
  "all" = "black"
)

environment_group_sets_colors <- c(
  "host-aquatic-soil" = "#000000",
  "aquatic-soil" = "#155A82",
  "host-aquatic" = "#83155F",
  "host-soil" = "#826516"
) %>%
  append(environment_group_colors)

kingdoms_colors <- c(
  "Bacteria" = "#3f78a9",
  "Fungi" = "#A93F55",
  "all" = "black"
)

dominance_colors <- c(
  "dominant" = "#03045e",
  "other" = "#00b4d8"
)

dysbalance_colors <- c(
  "balanced" = "#83c5be",
  "dysbalanced" = "#ffddd2"
)

modeltype_colors <- c(
  "univariate" = "#adb5bd",
  "multivariate" = "#343a40"
)

binary_colors <- c(
  "FALSE" = "#adb5bd",
  "TRUE" = "#343a40"
)

prevalence_group_colors <- c(
  "generalist" = "#3a86ff",
  "specialist" = "#ffbe0b"
)

phyla_colors <- c(
  # https://coolors.co/03045e-023e8a-0077b6-0096c7-00b4d8-48cae4-90e0ef-ade8f4-caf0f8
  Proteobacteria = "#03045e",
  Firmicutes = "#0077b6",
  Actinobacteria = "#00b4d8",
  Bacteroidetes = "#90e0ef",
  `other Bacteria` = "#8baec1",

  # https://coolors.co/03071e-370617-6a040f-9d0208-d00000-dc2f02-e85d04-f48c06-faa307-ffba08
  Ascomycota = "#9D0208",
  Basidiomycota = "#E85D04",
  `other Fungi` = "#c18b8d"
)

habitat_colors <- c(
  # aquatic
  `marine` = "#6fc6e5",
  `other freshwater` = "#316088",
  `other aquatic` = "#5f9fdd",
  `large lakes` = "#4f8dac",

  # soil
  `conifer forests` = "#5e9a43",
  `mediterranean forests` = "#cccb48",
  `temperate forests` = "#5f682f",
  `boreal forests` = "#7adb53",
  `tropical forests` = "#bfcf8f",

  # host
  `endophytes` = "#8b3554",
  `other hosts` = "#d57d35",
  `insects` = "#d83f8a",
  `lung` = "#7a482f",
  `gut` = "#e1453a",
  `plants` = "#d9849c",
  `skin` = "#a53d32",
  `mouth` = "#d6977c"
)

scale_color_environment_group <- function(drop = FALSE, ...) {
  scale_color_manual(
    values = environment_group_colors,
    na.value = "#d9d9d9",
    drop = drop,
    ...
  )
}

scale_fill_environment_group <- function(drop = FALSE, na.value = "#d9d9d9", ...) {
  scale_fill_manual(
    values = environment_group_colors,
    na.value = na.value,
    drop = drop,
    ...
  )
}

scale_fill_dominance <- function(drop = FALSE, ...) {
  scale_fill_manual(
    values = dominance_colors,
    na.value = "#d9d9d9",
    drop = drop,
    ...
  )
}

scale_color_dominance <- function(drop = FALSE, ...) {
  scale_color_manual(
    values = dominance_colors,
    na.value = "#d9d9d9",
    drop = drop,
    ...
  )
}

scale_fill_dysbalance <- function(drop = FALSE, ...) {
  scale_fill_manual(
    values = dysbalance_colors,
    na.value = "#d9d9d9",
    drop = drop,
    ...
  )
}

scale_fill_binary <- function(drop = FALSE, ...) {
  scale_fill_manual(
    values = binary_colors,
    na.value = "#d9d9d9",
    drop = drop,
    ...
  )
}

scale_color_dysbalance <- function(drop = FALSE, ...) {
  scale_color_manual(
    values = dysbalance_colors,
    na.value = "#d9d9d9",
    drop = drop,
    ...
  )
}

scale_color_kingdom <- function(drop = FALSE, na.value = "#d9d9d9", ...) {
  scale_color_manual(
    values =  purrr::simplify(kingdoms_colors),
    na.value = na.value,
    drop = drop,
    ...
  )
}

scale_fill_kingdom <- function(drop = FALSE, na.value = "#d9d9d9", ...) {
  scale_fill_manual(
    values =  purrr::simplify(kingdoms_colors),
    na.value = na.value,
    drop = drop,
    ...
  )
}

scale_fill_phyla <- function(drop = FALSE, ...) {
  scale_fill_manual(
    values =  purrr::simplify(phyla_colors),
    na.value = "#d9d9d9",
    drop = drop,
    ...
  )
}

scale_color_phyla <- function(drop = FALSE, ...) {
  scale_color_manual(
    values =  purrr::simplify(phyla_colors),
    na.value = "#d9d9d9",
    drop = drop,
    ...
  )
}

scale_color_binary <- function(drop = FALSE, ...) {
  scale_color_manual(
    values =  purrr::simplify(binary_colors),
    na.value = "#d9d9d9",
    drop = drop,
    ...
  )
}

scale_color_habitat <- function(drop = FALSE, ...) {
  scale_color_manual(
    values =  purrr::simplify(habitat_colors),
    na.value = "#d9d9d9",
    drop = drop,
    ...
  )
}

scale_fill_habitat <- function(drop = FALSE, ...) {
  scale_fill_manual(
    values =  purrr::simplify(habitat_colors),
    na.value = "#d9d9d9",
    drop = drop,
    ...
  )
}

scale_fill_modeltype <- function(drop = FALSE, ...) {
  scale_fill_manual(
    values =  purrr::simplify(modeltype_colors),
    na.value = "#d9d9d9",
    drop = drop,
    ...
  )
}


scale_color_modeltype <- function(drop = FALSE, ...) {
  scale_color_manual(
    values =  purrr::simplify(modeltype_colors),
    na.value = "#d9d9d9",
    drop = drop,
    ...
  )
}

scale_fill_prevalence_group <- function(drop = FALSE, ...) {
  scale_fill_manual(
    values =  purrr::simplify(prevalence_group_colors),
    na.value = "#000000",
    drop = drop,
    ...
  )
}


scale_color_prevalence_group <- function(drop = FALSE, ...) {
  scale_color_manual(
    values =  purrr::simplify(prevalence_group_colors),
    na.value = "#000000",
    drop = drop,
    ...
  )
}

stat_compare_means <- function(size = 2.5, symnum.args = list(
                                 cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                 symbols = c("\u2605\u2605\u2605", "\u2605\u2605", "\u2605", "ns")
                               ), ...) {
  ggpubr::stat_compare_means(
    size = size,
    symnum.args = symnum.args,
    ...
  )
}

stat_compare_means_environment_group <- function(...) {
  stat_compare_means(
    comparisons = list(
      c("aquatic", "host"),
      c("aquatic", "soil"),
      c("host", "soil")
    ),
    ...
  )
}

stat_compare_means_dysbalance <- function(method = "wilcox", ...) {
  stat_compare_means(
    method = method,
    comparisons = list(
      c("balanced", "dysbalanced")
    ),
    ...
  )
}

stat_compare_means_kingdom <- function(...) {
  stat_compare_means(
    comparisons = list(c("Fungi", "Bacteria")),
    ...
  )
}

theme_my <-
  ggplot2::theme_minimal(base_size = 20) +
  ggplot2::theme(
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_line(colour = "black", size = 0.8),
    axis.text = ggplot2::element_text(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(size = 0.8, fill = "transparent")
  )

theme_pub <- function(show_axes = TRUE) {
  line_width <- 0.5
  small_line_width <- line_width / 2

  theme_base <- function() {
    ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        axis.line.x = ggplot2::element_line(size = small_line_width),
        axis.line.y = ggplot2::element_line(size = small_line_width),
        axis.ticks = ggplot2::element_line(colour = "black", size = small_line_width),

        axis.title = element_text(size = 9),
        strip.text = element_text(size = 9),
        title = element_text(size = 9),

        legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(),
        legend.box.margin = margin(),

        plot.margin = margin(2, 2, 2, 2),
        panel.border = ggplot2::element_rect(fill = "transparent", size = line_width)
      )
  }

  if (show_axes) {
    theme_base()
  } else {
    theme_base() +
      ggplot2::theme(
        axis.title = ggplot2::element_blank()
      )
  }
}

ggplot2::theme_set(theme_pub())
ggplot2::update_geom_defaults("bar", base::list(fill = "lightgrey"))
ggplot2::update_geom_defaults("point", list(shape = 20, size = 1))

#
# Other functions
#

read_rdata <- function(file) {
  res <- new.env()
  load(file, envir = res)
  res
}

write_rds <- function(x, path, compress = "none", ...) {
  # enable compression if not set and file is > 1MB
  if (missing("compress") && object.size(x) > 10^6) {
    compress <- "gz"
  }

  readr::write_rds(x = x, path = path, compress = compress, ...)
}

filter <- function(...) dplyr::filter(...)
