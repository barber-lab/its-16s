#!/usr/bin/env Rscript

source(here::here("src/defaults.R"))

cache_path <- "/cache/drake/"
cache <- drake::new_cache(cache_path)
drake::drake_cache(cache_path)$unlock()

system("rm log/hmsc*.log")
system("rm *.log")

make_args <-
  list(
    plan = plan,
    parallelism = "clustermq",

    # allow concurrent workers accessing the env
    lock_envir = FALSE,
    recoverable = TRUE,
    recover = TRUE,

    # caching
    cache = cache,
    garbage_collection = TRUE,
    memory_strategy = "preclean",

    # config
    seed = 1337,
    keep_going = TRUE,
    jobs = jobs,

    # ouput
    verbose = 1,
    format = "qs",
    log_make = "log/drake.log",
    log_worker = TRUE
  )

# some targets e.g. file serialization can not be parallelized
# resulting in heavy overhead
default_parallel_targets <- c(
  "selected_generalists_specialists",
  "stepwise_constrained_ordinations_contribution",
  "stepwise_constrained_ordinations_anova",
  "generalists_alphadiv_test",
  "preselected_constrained_ordinations",
  "selected_specialists_specialness",
  "generalists_specialists_abundance_test",
  "generalists_specialists_overall_prevalences",
  "generalists_common_coabundance_graph_plt",
  "correlation_common_taxa_across_environment_group_topology_plt"
)
default_serial_targets <- c(
  "generalists_alphadiv_test_plt",
  "all_files"
)

if (Sys.getenv("TARGETS") != "") {
  # Targets can be overwritten by
  my_targets <- Sys.getenv("TARGETS") %>%
    str_split(",") %>%
    simplify()
} else {
  my_targets <- default_parallel_targets
}

message("TARGETS: ", my_targets %>% paste0(collapse = ", "))

make_args %>%
  inset2("jobs", 24) %>%
  inset2("targets", my_targets) %>%
  inset2("recover", TRUE) %>%
  inset2("recoverable", FALSE) %>%
  do.call(make, args = .)

make_args %>%
  inset2("parallelism", "loop") %>%
  inset2("targets", default_serial_targets) %>%
  inset2("recover", TRUE) %>%
  inset2("recoverable", FALSE) %>%
  do.call(make, args = .)