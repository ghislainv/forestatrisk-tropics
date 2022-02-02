#!/usr/bin/R

## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
## web             :https://ghislainv.github.io
## license         :GPLv3
## ==============================================================================


library(sf)
library(lwgeom)
library(dplyr)
library(here)

dataset <- "jrc2020"

# Import data for Australia
f_ctry <- here("Data", dataset, "Asia",
               "AUS-QLD", "data", "ctry_PROJ.shp")
f_pa <- here("Data", dataset, "Asia",
             "AUS-QLD", "data", "pa_PROJ.shp")
ctry <- st_read(f_ctry)
pa <- st_read(f_pa)

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
f_out <- here("Maps", dataset, "AUS-QLD", "pa_PROJ_intersect.shp")
st_write(pa_AUS_QLD, f_out)

# EOF
