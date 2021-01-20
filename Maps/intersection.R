library(sf)
library(lwgeom)
library(dplyr)

# Import data
ctry <- st_read("ctry_PROJ.shp")
pa <- st_read("pa_PROJ.shp")

# Find intersection
mat <- st_intersects(st_make_valid(pa), ctry, sparse=FALSE)
pa_in_ctry <- apply(mat, 1, any)
pa_int <- pa %>%
  dplyr::filter(pa_in_ctry)

# Tidying feature geometries
# https://www.r-spatial.org/r/2017/03/19/invalid.html
st_geometry_type(pa_int)
pa_AUS_QLD <- pa_int %>%
  dplyr::filter(st_is(pa_int, c("MULTIPOLYGON", "POLYGON"))) %>%
  st_cast("MULTIPOLYGON")

# Write data
st_write(pa_AUS_QLD, "results/pa_PROJ.shp")
