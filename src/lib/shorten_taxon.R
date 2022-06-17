shorten_taxon <- function(taxon) {
  recode(
    taxon,
    "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium" = "Rhizobium",
    "Burkholderia-Caballeronia-Paraburkholderia" = "Burkholderia",
    "Escherichia-Shigella" = "Escherichia"
  )
}
