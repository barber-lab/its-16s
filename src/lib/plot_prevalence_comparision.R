#' @param type one of "Intersection" or "Union"
plot_prevalence_comparision <- function(prevalence_summary, lineages, environment_groups, type = "Intersection") {
  phyla_counts_tbl <-
    prevalence_summary %>%
    filter(subset_name == "kingdom") %>%
    select(subset_name, subset_value, taxa) %>%
    unnest_longer(taxa) %>%
    rename(taxon = taxa, group = taxa_id) %>%
    unnest(taxon) %>%
    inner_join(lineages, by = "taxon") %>%
    group_by(subset_name, subset_value, group, phylum)

  count_phyla <- function(from, to, kingdom, phyla_counts_tbl, phyla_colors, type) {
    query <- switch(type, "Intersection" = "from & to", "Union" = "from | to")

    phyla_counts_tbl %>%
      ungroup() %>%
      dplyr::filter(group %in% c(from, to) & kingdom == !!kingdom) %>%
      mutate(existent = TRUE) %>%
      pivot_wider(names_from = group, values_from = existent) %>%
      rename_at(from, ~"from") %>%
      rename_at(to, ~"to") %>%
      # get intersection
      filter(eval(parse(text = query))) %>%
      # count phyla
      group_by(phylum) %>%
      count() %>%
      arrange(-n) %>%
      # summarise other
      ungroup() %>%
      # mutate(phylum = ifelse(row_number() < 5, phylum, "other")) %>%
      mutate(
        phylum = phylum %>% map2_chr(kingdom, phylum_to_color_name) %>% factor(level = names(phyla_colors))
      ) %>%
      group_by(phylum) %>%
      summarise(n = sum(n)) %>%
      arrange(-n)
  }

  plot_pie <- function(counts_tbl) {
    plt <-
      counts_tbl %>%
      ggplot(aes(x = "", y = n, fill = phylum)) +
      geom_bar(width = 1, stat = "identity") +
      coord_polar("y", start = 0) +
      scale_fill_phyla() +
      theme_void() +
      theme(legend.position = "none")
  }

  # # square layoutr with two triangles
  # tribble(
  #   ~from, ~to, ~kingdom,
  #   "aquatic", "soil", "Bacteria",
  #   "soil", "aquatic", "Fungi",
  # 
  #   "host", "aquatic", "Bacteria",
  #   "aquatic", "host", "Fungi",
  # 
  #   "host", "soil", "Bacteria",
  #   "soil", "host", "Fungi"
  # ) %>%
  
  environments <- c("host", "aquatic", "soil")
  kingdoms <- c("Bacteria", "Fungi")

  tribble(
    ~from, ~to, ~kingdom,
    
    "host", "aquatic", "Bacteria",
    "host", "aquatic", "Fungi",
    
    "host", "soil", "Bacteria",
    "host", "soil", "Fungi",
    
    "aquatic", "soil", "Bacteria",
    "aquatic", "soil", "Fungi"
    ) %>%  
    mutate(
      counts_tbl = list(from, to, kingdom) %>% pmap(~
      count_phyla(from = ..1, to = ..2, kingdom = ..3, phyla_counts_tbl = phyla_counts_tbl, phyla_colors = phyla_colors, type = type)),
      from = from %>% factor(levels = environment_groups),
      to = to %>% factor(levels = environment_groups %>% rev())
    ) %>%
    unnest(counts_tbl) %>%
    mutate(
      phylum = phylum %>% factor(levels = names(phyla_colors))
    ) %>%
    rename(Kingdom = kingdom, Phylum = phylum) %>%
    mutate(comparison = from %>% map2_chr(to, ~ paste0(.x, "∩", .y, collapse = ""))) %>%
  
    ggplot(aes(x = comparison, n)) +
    #geom_rect(aes(fill = Kingdom), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1) +
    #scale_fill_kingdom() +
    new_scale_fill() +
    geom_bar(aes(fill = Phylum), stat = "identity", position = "stack") +
    scale_fill_phyla() +
    #scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank()
    ) +
    labs(x = "Intersection", y = "Genera")
}
