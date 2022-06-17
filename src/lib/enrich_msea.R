#' Enrichment analysis using MSEA
#' 
#' @seealso  https://www.nature.com/articles/s41598-020-78511-y
#' @param term2gene data frame with columns term (name of set) and gene (name of element)
enrich_msea <- function(term2gene, test_taxa, universe = 10e3) {
  dir <- tempfile(pattern = "msea")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  gmt_file <- paste0(dir, "/db.gmt")
  test_taxa_file <- paste0(dir, "/test-taxa.txt")
  out_csv_file <- paste0(dir, "/out.csv")

  test_taxa %>% write_lines(test_taxa_file)

  # convert to GMT format
  term2gene %>%
    group_by(term) %>%
    summarise(text = term[[1]] %>% c(gene) %>% paste0(collapse = "\t"), .groups = 'drop') %>%
    pull(text) %>%
    write_lines(gmt_file)
  
  str_glue("python src/lib/do_msea.py {gmt_file} {test_taxa_file} {out_csv_file} {universe}") %>%
    system()
  
  res <-
    read_csv(out_csv_file, col_types = cols(
      term = col_character(),
      oddsratio = col_double(),
      pvalue = col_double(),
      qvalue = col_double(),
      shared = col_character(),
      n_shared = col_double()
    )) %>%
    rename(p.value = pvalue) %>%
    mutate(
      shared = shared %>% map(~ {
        .x %>%
          str_split("', '") %>%
          simplify() %>%
          map_chr(~ .x %>% str_remove("'\\]") %>% str_remove("\\['")) 
      }),
      
      q.value = p.value %>% p.adjust("fdr")
    )
    
  unlink(dir, recursive = TRUE)
  res
}