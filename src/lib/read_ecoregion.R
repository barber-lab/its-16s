#' @param x file name of shapefile
#' @param crs output coordinate reference system: integer with the EPSG code, or character with proj4string
read_ecoregion <- function(x, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") {
  tbl <- x %>%
    read_sf() %>%
    st_transform(crs = crs)
  type <- x %>% str_extract("(Freshwater|Terrestrial)")

  switch(type,
    "Terrestrial" = tbl %>% transmute(ECO_ID_U, eco_type = "Terrestrial", eco_region = WWF_MHTNAM),
    "Freshwater" = tbl %>% transmute(ECO_ID_U, eco_type = "Freshwater", eco_region = MHT_TXT),
    tbl
  )
}
