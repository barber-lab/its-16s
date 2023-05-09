#!/usr/bin/env R

# Import generalists metadata retrieved from NCBI in order to identify
# generalist species for downstream analysis ... BGCs identification.

# Taxons with no metadata from ncbi were manually removed.

get_annotation <- function(generalist_genus, specialist_genus, lineages) {
  
  king2complete <- tibble(genus = c(
    "Burkholderia-Caballeronia-Paraburkholderia",
    "Escherichia-Shigella",
    "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium")
  ) %>% 
    mutate(genus = str_split(genus, "-")) %>% 
    unnest(genus) %>% 
    mutate(kingdom = "Bacteria")
  
  a <- generalist_genus %>% 
    rename(genus = generalist) %>% 
    mutate(type = "Generalist") %>% 
    bind_rows(
      specialist_genus %>% 
        rename(genus = specialists) %>% 
        mutate(type = "Specialist")
    ) %>% 
    left_join(lineages %>% 
                select(kingdom, genus) %>% 
                distinct() %>% 
                bind_rows(king2complete)) %>% 
    add_row(genus = c("Lachnospiraceae", "Erythrobacter"),
            kingdom = c("Bacteria", "Bacteria"),
            type = c("Generalist", "Specialist"))
}

# List jsonl files
get_files_list <- function(path, pttrn) {
  
  list.files(path = path, full.names = T, pattern = pttrn) %>% 
    as_tibble()
  
}

# Read jsonl, unnest desired levels and select desired fields
import_jsonl <- function(path) {
  
  fromJSON(path) %>% 
    as_tibble() %>%
    unnest(reports) %>%
    unnest(assembly_info) %>%
    unnest(organism) %>%
    unnest(assembly_stats) %>%
    select(accession, 
           assembly_level, 
           refseq_category, 
           organism_name,
           contains("n50"))
  
}

# Tidy metadata in a tibble format and remove genus not available in NCBI
# During import if the jsonl is empty (no genome available), NA avoid crashing
parse_jsonl <- function(jsonl) {
  
  jsonl %>%
    mutate(metadata = map(
      .x = value, 
      .f = possibly(~ import_jsonl(.x),NA))
    ) %>%
    unnest(metadata) %>%
    drop_na(organism_name)
  
}

generalists_no_genome <- function(genus) {
  
  # Subsetting genus with only one word to create filtering pattern
  # Some genus has weird annotation, like Blattella germanica (2 words)
  
  patt1 <- genus %>% 
    separate(col = generalist, 
             into = paste("n", rep(1:2), sep = ""), 
             sep = " ", 
             remove = F) %>% 
    filter(is.na(n2)) %>% 
    select(-n2) %>% 
    pull(n1)
  
  # After subsetting pieces of the generalist, it is been removed those not found.
  # Manual curation done beforehand
  
  patt2 <- genus %>% 
    mutate(generalist = str_remove(generalist, " group")) %>% 
    separate(col = generalist, 
             into = paste("n", rep(1:2), sep = ""), 
             sep = " ", 
             remove = F) %>% 
    drop_na(n2) %>% 
    filter(!n2 == "Solibacter") %>% 
    filter(!n2 == "R-7") %>% 
    filter(!n2 == "sensu") %>% 
    filter(!n1 == "Prevotella") %>% 
    filter(!n2 == "UCG-014") %>% 
    pull(n2)
  
  # There were identified 53 taxons in the NCBI genome out of 58
  #2023-01-11
  patt <- c(patt1, patt2) %>% paste(collapse = "|")
  
}

filter_generalists <- function(generalists, patt) {
  
  generalists %>% 
    filter(str_detect(organism_name, patt)) %>%
    mutate(genus = str_replace(organism_name, "(?s) .*", "")) %>%
    filter(!str_detect(genus, "\\[")) %>%
    filter(!str_detect(genus, "^'")) %>% 
    filter(!str_detect(genus, "synthetic")) %>% 
    filter(!str_detect(genus, "Candidatus")) %>% 
    filter(!str_detect(genus, "uncultured"))
  # %>% 
  #   group_by(genus, refseq_category) %>% 
  #   count()
  
}

# Remove organisms/genus that doesn't match with specialists genus list
filter_specialists <- function(tidy_specialists, patt) {
  
  tidy_specialists %>% 
    mutate(genus = str_replace(organism_name, "(?s) .*", "")) %>% 
    filter(!str_detect(genus, "Candidatus")) %>% 
    filter(!str_detect(genus, "uncultured")) %>% 
    filter(!str_detect(genus, "Blattella")) %>% 
    filter(!str_detect(genus, "Blattella")) %>% 
    filter(!str_detect(genus, "\\[")) %>% 
    filter(!str_detect(genus, "synthetic"))
  # %>% 
  #   group_by(genus, refseq_category) %>% 
  #   count()
  
}

shuffle_genome <- function(genome_df, metadata) {
  
    # set.seed(123) needs to be called each time before sample_n is performed. 
  set.seed(123)
  genome_df %>% 
      left_join(metadata) %>% 
      group_by(kingdom, genus, type) %>%
      sample_n(size = 60, replace = T) %>% 
      distinct()

}

remove_genome_id <- function() {
  # When trying to download the genome, the zip file was empty
  # Or failed at antismash
  c(
    "GCF_008567345.1",
    "GCF_008565675.1",
    "GCF_001276935.1",
    "GCA_003029985.1"
  )
}
  
export_genome_list <- function(failed_genomes, random_generalist_genomes, random_specialist_genomes, lineages) {
  
  # Taxa that doesn't match the lineages table anymore needs to be fixed.
  king2complete <- tibble(genus = c(
    "Burkholderia-Caballeronia-Paraburkholderia",
    "Escherichia-Shigella",
    "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
    "Lachnospiraceae")
  ) %>% 
    mutate(genus = str_split(genus, "-")) %>% 
    unnest(genus) %>% 
    mutate(kingdom = "Bacteria")
  
  list <- random_specialist_genomes %>% 
    bind_rows(random_generalist_genomes) %>% 
    #select(accession, genus) %>% 
    left_join(lineages %>% 
                select(kingdom, genus) %>%
                distinct() %>% 
                bind_rows(king2complete)) %>% 
    ungroup() %>% 
    filter(!accession  %in% failed_genomes)
  
  list %>% 
    filter(str_detect(kingdom, "Bacteria")) %>% 
    select(accession) %>% 
    write_tsv("raw/generalist_specialist_genome/bacteria.tsv")
  
  list %>% 
    filter(str_detect(kingdom, "Fungi")) %>% 
    select(accession) %>% 
    write_tsv("raw/generalist_specialist_genome/fungi.tsv")
  
  return(list)
}

antismash_files <- function(path) {
  
  list.files(path = path, full.names = T, pattern = ".json$", recursive = T) %>% 
    as_tibble() %>% 
    mutate(genome_id = str_split(string = value, 
                                 pattern = "/", 
                                 simplify = T)[,3]) 
  
}

is_there_results <- function(file) {
  
  flat <- fromJSON(file, flatten = T) %>% 
    .$records %>% 
    flatten(recursive = T) %>%
    # Will remove empty columns
    # If there are hits, areas must not be empty
    select(where(~!all(lengths(.) == 0)))
  
  # If areas was empty, it was removed
  if("areas"  %in% colnames(flat)) {
    
    flat %>% 
      select(areas)

    } else {
      NA
    }
}

parse_antismash <- function(files) {
  
  files %>% 
    mutate(bgcs = map(
      .x = value, 
      .f = ~ is_there_results(.x))
    )
}

group_bgcs <- function(df) {
  
  df %>% 
    unnest(bgcs) %>% 
    unnest(areas) %>% 
    select(genome_id, products) %>% 
    unnest(cols = c(products)) %>% 
    group_by(genome_id, products) %>% 
    count() %>% 
    arrange(desc(n))
  
}

generate_annotation <- function(bgcs_by_genome, genomes2find_BGCs, generalist_genus, specialist_genus) {

  type2complete <- tibble(genus = c("Erythrobacter","Lachnospiraceae"),
                          type = c("Specialist","Generalist"))

  g <- generalist_genus %>% 
    mutate(type = "Generalist") %>% 
    rename(genus = generalist)
  
  s <- specialist_genus %>% 
    mutate(type = "Specialist") %>% 
    rename(genus = specialists)
  
  
  metadata <- genomes2find_BGCs %>% 
    left_join(bind_rows(g,s, type2complete))
  
  metadata %>% 
    left_join(bgcs_by_genome, by = c("accession" = "genome_id")) %>% 
    mutate(n = replace_na(n, 0))
  # %>% 
    # Remove taxa with no antismash hits
    #drop_na(products)

  
}

scale_fill_prevalence_type <- function(drop = FALSE, ...) {
  scale_fill_manual(
    values =  purrr::simplify(prevalence_type),
    na.value = "#000000",
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



plot_bgcs <- function(generalist_specialists_bgcs) {
  
  p <- generalist_specialists_bgcs %>% 
    #unite(col = Categories, c("kingdom", "type"), sep = " ") %>% 
    rename(Groups = type) %>% 
    group_by(accession, Groups, kingdom) %>% 
    summarise(value = sum(n)) %>%
    ggplot(aes(Groups, value, fill = Groups)) +
    geom_boxplot(width=0.5) +
    stat_compare_means(method = "wilcox", 
                       comparisons = list(c("Specialist", "Generalist"))) +
    scale_y_log10() +
    facet_wrap(c("kingdom")) +
    scale_fill_prevalence_type() +
    annotation_logticks(sides = "l") +
    ylab("Number of BGCs (Log10)") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 10))
    
  
  ggsave(plot = p,
         filename = "plots/02-generalist-specialist-bgcs/bgcs_boxp_no_gff.png",
         device = "png", 
         width = 9, 
         height = 7, 
         dpi = 150)
  
}

count_cds <- function(df) {
  
  df %>%
    mutate(genome_id = str_split(string = value, 
                                 pattern = "/", 
                                 simplify = T)[,3]) %>% 
    mutate(genome_id = str_extract(genome_id, "[^_]*_[^_]*")) %>% 
    mutate(cds_n = map(.x = value, 
                       .f = ~ read.delim(.x, header=F, comment.char="#") %>% 
                         filter(V3 == "CDS") %>% 
                         nrow())) %>% 
    unnest(cds_n)
                         
  
}

# Remove the header and count total bases
get_genome_size <- function(df) {
  
  df %>% 
    mutate(genome_id = str_split(string = value, 
                                 pattern = "/", 
                                 simplify = T)[,3]) %>% 
    mutate(genome_id = str_extract(genome_id, "[^_]*_[^_]*")) %>% 
    mutate(cmd = str_glue("grep -v '>' {value} | wc -c")) %>%
    mutate(g_size = map(.x = cmd, 
                        .f = ~ system(.x, intern = T)))
}

plot_cds_gen_metrics <- function(genomes2find_BGCs, cds_by_genome, genome_length) {
  
  cds <- genomes2find_BGCs %>% 
    distinct() %>% 
    select(accession, genus, type, kingdom) %>% 
    left_join(cds_by_genome, by = c("accession" = "genome_id")) %>% 
    select(-value)
  
  a <- cds %>% 
    rename(Group = type) %>% 
    ggplot(aes(Group, cds_n)) +
    geom_boxplot(aes(fill = Group)) +
    facet_wrap(c("kingdom")) +
    stat_compare_means(method = "wilcox", 
                       comparisons = list(c("Generalist", "Specialist"))) +
    ylab("CDS") +
    theme(axis.title.x = element_blank()) +
    scale_fill_prevalence_type()
  
  gen <- genomes2find_BGCs %>% 
    distinct() %>% 
    select(accession, genus, type, kingdom) %>% 
    left_join(genome_length, by = c("accession" = "genome_id")) %>% 
    select(-c(value, cmd)) %>% 
    unnest(g_size) %>% 
    mutate(g_size = as.numeric(g_size)) %>% 
    mutate(g_size = g_size/10^6)
  
  b <- gen %>% 
    rename(Group = type) %>% 
    ggplot(aes(Group, g_size)) +
    geom_boxplot(aes(fill = Group)) +
    facet_wrap(c("kingdom")) +
    stat_compare_means(method = "wilcox", 
                       comparisons = list(c("Generalist", "Specialist"))) +
    ylab("Genome size (MB)") +
    theme(axis.title.x = element_blank()) +
    scale_fill_prevalence_type()
  
  
  genVScds <- cds %>% 
    left_join(gen)
  
  c <- genVScds %>% 
    rename(Group = type) %>% 
    ggplot(aes(g_size, cds_n)) +
    geom_point(aes(color = Group)) +
    facet_wrap(c("kingdom")) +
    scale_color_prevalence_group() +
    xlab("Genome size (MB)") +
    ylab("CDS")
  
  d <- genVScds %>% 
    mutate(ratio_cds_g = cds_n/g_size) %>%
    rename(Group = type) %>% 
    ggplot(aes(Group, ratio_cds_g)) +
    geom_boxplot(aes(fill = Group)) +
    facet_wrap(c("kingdom")) +
    stat_compare_means(method = "wilcox", 
                       comparisons = list(c("Generalist", "Specialist"))) +
    ylab("Ratio CDS/Genome size (MB)") +
    scale_fill_prevalence_type()
    theme(axis.title.x = element_blank())
  
  p <- wrap_plots(a, b, guides = "collect")
  #"plots/02-generalist-specialist-bgcs/cds-gensize.png"
  ggsave(filename = "plots/02-generalist-specialist-bgcs/cds-gensize.png", 
         plot = p, 
         device = "png", 
         width = 12, 
         height = 9,
         dpi = 150)
}

ratio_bgc_cds <- function(generalist_specialists_bgcs, cds_by_genome) {
  
  p <- generalist_specialists_bgcs %>% 
    select(accession, genus, type, kingdom, n) %>% 
    rename(bgcs_n = n) %>% 
    left_join(
      cds_by_genome %>% 
        rename(accession = genome_id) %>% 
        select(-value)
    ) %>% 
    mutate(ratio_bgc_cds = bgcs_n/cds_n*100) %>% 
    rename(Groups = type) %>%
    ggplot(aes(Groups, ratio_bgc_cds)) +
    geom_boxplot(aes(fill = Groups), width=0.5) +
    facet_wrap(c("kingdom")) +
    stat_compare_means(method = "wilcox", 
                       comparisons = list(c("Specialist", "Generalist"))) +
    scale_fill_prevalence_type() +
    ylab("Ratio BGCs/CDS (%)") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 10))
  
  ggsave(filename = "plots/02-generalist-specialist-bgcs/ratio-bgcs-cds-no-gff.png", 
         plot = p, 
         device = "png", 
         width = 12, 
         height = 9,
         dpi = 150)
  
}

