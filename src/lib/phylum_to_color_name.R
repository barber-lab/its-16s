phylum_to_color_name <- function(phylum, kingdom) {
  case_when(
    phylum %in% names(phyla_colors) ~ phylum,
    kingdom == "Bacteria" ~ "other Bacteria",
    kingdom == "Fungi" ~ "other Fungi"
   )
}
