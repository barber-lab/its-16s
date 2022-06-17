get_taxranks <- function(current_rank, direction = "downstream") {
  taxranks <- c(
    "kingdom" = 8, "phylum" = 7, "class" = 6, "order" = 5, "family" = 4,
    "genus" = 3, "species" = 2, "otu" = 1
  )

  switch(direction,
    "with_downstream" = c(
      current_rank, taxranks[taxranks < taxranks[current_rank]] %>% names()
    ),
    "downstream" = taxranks[taxranks < taxranks[current_rank]] %>% names(),
    "with_upstream" = c(
      current_rank,
      taxranks[taxranks > taxranks[current_rank]] %>% names()
    ),
    "upstream" = taxranks[taxranks > taxranks[current_rank]] %>% names(),
    stop("direction unknown")
  )
}
