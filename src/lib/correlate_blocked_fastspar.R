#' Coabundance analysis using SparCC as implemented in fastspar
#'
#' This implementation is way faster than `correlate_spiec_easi_sparcc` but requires linux and the external shell command `fastspar`.
#' @param data integer matrix of abundance count data. One sample per row and one taxon per column
correlate_blocked_fastspar <- function(data, iterations = 50, exclude_iterations = 10, permutations = 200, strata = NULL, threads = getOption("mc.cores"), ...) {
  system <- function(...) base::system(ignore.stdout = TRUE, ignore.stderr = TRUE, ...)
  
  threads <- min(parallel::detectCores(), threads)
  
  # sanity checks
  if (! "matrix" %in% class(data)) stop("data must be of type matrix")
  if (str_glue("fastspar --version") %>% system() != 0) {
    stop("Command fastspar not found")
  }
  
  dir <- tempfile(pattern = "fastspar")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  
  data_path <- str_glue("{dir}/data.tsv")
  
  data %>%
    t() %>%
    as_tibble(rownames = "#OTU ID") %>%
    write_tsv(data_path)
  
  paste(
    "fastspar",
    "--yes",
    "--iterations", iterations,
    "--exclude_iterations", exclude_iterations,
    "--otu_table", data_path,
    "--correlation", paste0(dir, "/cor.tsv"),
    "--covariance", paste0(dir, "median_covariance.tsv"),
    "--threads", threads,
    sep = " "
  ) %>%
    system()
  
  
  paste0(dir, "/bootstraps_counts") %>% dir.create()
  paste0(dir, "/bootstraps_cor") %>% dir.create()
  
  
  permutations %>%
    seq() %>%
    walk(
      ~ data %>%
        as_tibble(rownames = "sample_id") %>%
        mutate(strata = strata) %>%
        group_by(strata) %>%
        # sample group_by does not work, nesting needed
        nest() %>%
        mutate(
          data = data %>% map(
            ~ .x %>%
              mutate(sample_id = sample_id %>% sample(replace = FALSE)) %>%
              arrange(sample_id)
          )
        ) %>%
        unnest(data) %>%
        ungroup() %>%
        pivot_longer(cols = matches("^V"), names_to = "#OTU ID", values_to = "abundance") %>%
        select(sample_id, `#OTU ID`, abundance) %>%
        arrange(sample_id, `#OTU ID`) %>%
        mutate(sample_id = sample_id %>% paste0("bs_", .)) %>%
        pivot_wider(names_from = sample_id, values_from = abundance) %>%
        write_tsv(
          str_glue("{dir}/bootstraps_counts/data_{.x}.tsv")
        )
    )
  
  paste(
    "parallel",
    "--jobs", threads,
    
    "fastspar",
    "--yes",
    "--otu_table {}",
    "--correlation", paste0(dir, "/bootstraps_cor/cor_{/}"),
    "--covariance", paste0(dir, "/bootstraps_cor/cov_{/}"),
    "--iterations", iterations,
    "--exclude_iterations", exclude_iterations,
    "--threads", 1,
    ":::",
    paste0(dir, "/bootstraps_counts/*"),
    sep = " "
  ) %>%
    system()
  
  paste(
    "fastspar_pvalues",
    "--otu_table", data_path,
    "--correlation", paste0(dir, "/cor.tsv"),
    "--prefix", paste0(dir, "/bootstraps_cor/cor_data_"),
    "--permutations", bootstraps,
    "--outfile", paste0(dir, "/pvals.tsv"),
    sep = " "
  ) %>%
    system()
  
  pval_tbl <-
    paste0(dir, "/pvals.tsv") %>%
    read_tsv(col_types = cols(.default = "c")) %>%
    dplyr::rename(from = `#OTU ID`) %>%
    tidyr::pivot_longer(-from, names_to = "to", values_to = "p.value")
  
  cor_tbl <-
    paste0(dir, "/cor.tsv") %>%
    readr::read_tsv(col_types = cols(.default = "c")) %>%
    dplyr::rename(from = `#OTU ID`) %>%
    tidyr::pivot_longer(-from, names_to = "to", values_to = "estimate")
  
  paste0("rm -rf ", dir) %>% system()
  
  pval_tbl %>%
    dplyr::inner_join(cor_tbl, by = c("from", "to")) %>%
    # only keep triangle
    dplyr::mutate(comp = from %>% map2_chr(to, ~ c(.x, .y) %>%
                                             sort() %>%
                                             paste0(collapse = "-"))) %>%
    dplyr::group_by(comp) %>%
    dplyr::slice(1) %>%
    dplyr::select(-comp) %>%
    readr::type_convert() %>%
    dplyr::ungroup() %>%
    dplyr::filter(from != to) %>%
    dplyr::mutate(q.value = p.adjust(p.value, method = "fdr")) %>%
    dplyr::select(from, to, p.value, q.value, estimate) %>%
    as_tbl_graph(...) %>%
    `attr<-`(., "method", "sparcc")
}
