#!/usr/bin/env R

# Start with valid samples and add metadata
annotate_valid_samples <- function(valid_samples, bioproject_primers, biosamples) {
  
  valid_samples %>% 
    distinct(sample_id, bioproject_id) %>% 
    left_join(bioproject_primers) %>% 
    left_join(
      biosamples %>% 
        mutate(environment_group = case_when(
          bioproject_id == "PRJNA550037" & environment_group == "host" ~ "soil",
          bioproject_id == "PRJNA418896" & environment_group == "host" ~ "soil",
          TRUE ~ environment_group)
        ) %>% 
        distinct()
      )

}

# Put df in a longer format, merging metadata and abundance
# Remove not detected taxa
# .data[[primer_bac]] is replaced by a variable of interest
group_primers <- function(metadata, raw_abundance_bac, king_abrv) {
  
  abundant_primers <- metadata %>% 
    left_join(
      raw_abundance_bac %>% 
        pivot_longer(cols = -sample_id, 
                     names_to = "taxon", 
                     values_to = "abundance")
    ) %>%
    filter(abundance > 0) %>%
    distinct() %>% 
    dplyr::select(contains("set"), sample_id, taxon, abundance, environment_group)
  
  # Expression to filter primers
  exp_ <- c('V1-V4' = "set1bac == 'V1-V4'",
            'V3-V5' = "set2bac == 'V3-V5'",
            'ITS1' = "set1fun == 'ITS1'",
            'ITS2' = "set2fun == 'ITS2'")
  
 var2 <- exp_ %>% 
    map(rlang::parse_expr) %>% 
    map(
      function(pick_these) 
      abundant_primers %>% 
        filter(!!pick_these) %>% 
        select(-contains("set"))
      ) %>% 
      tibble(data = .) %>% 
      mutate(primers = names(exp_)
      ) %>% 
    left_join(
      colnames(abundant_primers) %>%
        as_tibble() %>%
        filter(str_detect(value, "set")) %>% 
        mutate(primers = c("V1-V4", "V3-V5", "ITS1", "ITS2"))
    ) %>%
   unnest() %>% 
   group_by(primers, environment_group, value) %>% 
   nest() %>% 
   filter(str_detect(value, king_abrv))

}

# Detect common taxa between set of primers 
common_taxa <- function(data1, data2) {
  
  data1 %>% 
    pull(taxon) %>%
    unique() %>% 
    intersect(
      data2 %>% 
        pull(taxon) %>%
        unique()) %>% 
    sort()
}

intersection_between_primers <- function(data) {
  
  # Create sets to perform correlation
  set1 <- data %>%
    unite(col = set1,
          c("primers", "environment_group"), 
          sep = "__") %>% 
    ungroup() %>% 
    select(-value) %>% 
    dplyr::rename(data1 = data)
  
  set2 <- data %>%
    unite(col = set2,
          c("primers", "environment_group"), 
          sep = "__") %>% 
    ungroup() %>% 
    select(-value) %>% 
    dplyr::rename(data2 = data)
  
  # function filter taxa
  common <- set1 %>%
    expand_grid(set2) %>%
    mutate(intersect_ = map2(
      .x = data1,
      .y = data2,
      .f = ~ common_taxa(.x, .y)))

  # replace data by intersection
  fltrd <- common %>%
    mutate(data1 = map2(
      .x = data1,
      .y = intersect_,
      .f = ~ .x %>%
        filter(abundance > 0) %>%
        filter(taxon  %in% .y)
    )) %>%
    mutate(data2 = map2(
      .x = data2,
      .y = intersect_,
      .f = ~ .x %>%
        filter(abundance > 0) %>%
        filter(taxon  %in% .y)
    ))
}

# Calculate correlation
# Get upper triangle
# Convert matrix to vector
# Remove NA
# EXTREMELY IMPORTANT TO SORT. ORDER MATTERS!
do_correlation <- function(df) {
  
  df %>% 
    distinct() %>% 
    #mutate(abundance = log10(abundance)) %>% 
    pivot_wider(names_from = taxon, 
                values_from = abundance, 
                values_fill = 0, 
                names_sort = T) %>%
    column_to_rownames("sample_id") %>%
    cor(method = "pearson") %>%
    get_upper_tri() %>%
    as.vector()
  # %>%
  #   .[!is.na(.)]
        
  # %>% 
  #   # Remove columns to make object smaller
  #   select(-c(data1, data2, intersect_))
}

correlate_primers <- function(df) {
  
  df %>% 
    mutate(cor1 = map(.x = data1, .f = ~ do_correlation(.x))) %>% 
    mutate(cor2 = map(.x = data2, .f = ~ do_correlation(.x)))
    # Remove columns to make object smaller
    # select(-c(matrix))
  
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

# https://towardsdatascience.com/how-to-measure-similarity-between-two-correlation-matrices-ce2ea13d8231
calculate_similarity <- function(df) {
  
  a <- df %>% 
    mutate(sim = map2(.x = cor1, 
                      .y = cor2,
                      .f = ~ cor.test(x = .x, y = .y, method = "pearson") %>% 
                        broom::tidy())) %>% 
    unnest(sim) %>%
    select(-contains("data"), -contains("cor"), -intersect_)
  
  # p-value adjustment
  a %>% 
    mutate(p.adj = a %>% 
             select(p.value) %>% 
             as.matrix %>% 
             as.vector %>% 
             p.adjust(method='fdr'))
    
  # %>% 
  #   mutate(p.value = round(p.value, digits = 3))
}

tidy_corr <- function(df) {
  
  df %>% 
    select(set1, set2, sim) %>% 
    mutate(df = map(sim, ~ .x %>% 
                      unclass() %>%
                      as_tibble() %>%
                      slice(1) %>%
                      select(method, estimate, p.value))) %>% 
    unnest(df)
  
}

  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
      
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
primers_heatmap <- function(df, file_name, rm_aqua) {
  
  # We don't have samples from aquatic with its2, then we decided to remove aquatic from fungi.
  if (rm_aqua == "fungi") {
    
    df <- df %>% 
      filter(!str_detect(set1, "aquatic")) %>% 
      filter(!str_detect(set2, "aquatic"))
  }
  
  mat <- df %>% 
    select(contains("set"), estimate) %>%
    # Some correlations were NA
    mutate(estimate = replace_na(estimate, 0)) %>% 
    mutate(estimate = round(estimate, digits = 2)) %>% 
    pivot_wider(names_from = set2, values_from = estimate) %>% 
    column_to_rownames("set1") %>% 
    as.matrix()
  
  # Reorder the correlation matrix
  cormat <- reorder_cormat(mat)
  
  lower_tri <- get_lower_tri(cormat)
  
  # Melt the correlation matrix
  melted_cormat <- melt(lower_tri, na.rm = TRUE)
  
  h <- melted_cormat %>% 
    ggplot(mapping = aes(X1, X2, fill = value), 
           na.rm = T) +
    geom_tile() +
    scale_fill_distiller(
      palette = "PuOr",
      na.value = "white",
      direction = 1,
      limits = c(-1, 1),
      name = "Pearson\nCorrelation:") +
    geom_text(mapping = aes(X1, X2, label = value), 
              color = "black", 
              size = 5) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_text(size=15),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  mat <- df %>% 
    select(contains("set"), p.adj) %>%
    pivot_wider(names_from = set2, values_from = p.adj) %>% 
    column_to_rownames("set1") %>% 
    as.matrix()
  
  # Reorder the correlation matrix
  cormat <- reorder_cormat(mat)
  
  lower_tri <- get_lower_tri(cormat)
  
  # Melt the correlation matrix
  melted_cormat <- melt(lower_tri, na.rm = TRUE)
  
  p <- melted_cormat %>% 
    ggplot(mapping = aes(X1, X2, fill = value), 
           na.rm = T) +
    geom_tile() +
    scale_fill_distiller(
      palette = "PuOr",
      na.value = "white",
      direction = 1,
      limits = c(-1, 1),
      name = "Adjusted\np-values:") +
    geom_text(mapping = aes(X1, X2, label = value), 
              color = "black", 
              size = 5) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_text(size=15),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  w <- wrap_plots(h,p)
  
  ext <- ".png"
  heat_name <- paste0("plots/03-primers-effect/", file_name, ext)

  ggsave(filename = heat_name,
         plot = w,
         device = "png",
         width = 16,
         height = 9,
         dpi = 300)
  
  ext <- ".pdf"
  heat_name <- paste0("plots/03-primers-effect/", file_name, ext)
  
  ggsave(filename = heat_name,
         plot = w,
         device = "pdf",
         width = 16,
         height = 9,
         dpi = 300)
  
}

primers_heatmap_rds <- function(df, rm_aqua) {
  
  # We don't have samples from aquatic with its2, then we decided to remove aquatic from fungi.
  if (rm_aqua == "fungi") {
    
    df <- df %>% 
      filter(!str_detect(set1, "aquatic")) %>% 
      filter(!str_detect(set2, "aquatic"))
  }
  
  mat <- df %>% 
    select(contains("set"), estimate) %>%
    # Some correlations were NA
    mutate(estimate = replace_na(estimate, 0)) %>% 
    mutate(estimate = round(estimate, digits = 2)) %>% 
    pivot_wider(names_from = set2, values_from = estimate) %>% 
    column_to_rownames("set1") %>% 
    as.matrix()
  
  # Reorder the correlation matrix
  cormat <- reorder_cormat(mat)
  
  lower_tri <- get_lower_tri(cormat)
  
  # Melt the correlation matrix
  melted_cormat <- melt(lower_tri, na.rm = TRUE)
  
  h <- melted_cormat %>% 
    ggplot(mapping = aes(X1, X2, fill = value), 
           na.rm = T) +
    geom_tile() +
    scale_fill_distiller(
      palette = "PuOr",
      na.value = "white",
      direction = 1,
      limits = c(-1, 1),
      name = "Pearson\nCorrelation:") +
    geom_text(mapping = aes(X1, X2, label = value), 
              color = "black", 
              size = 5) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_text(size=15),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  
  
}
