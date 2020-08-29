#!/usr/bin/R

## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
## web             :https://ghislainv.github.io
## license         :GPLv3
## ==============================================================================

## For creating maps with R, see:
## https://geocompr.robinlovelace.net/adv-map.html
## https://cran.r-project.org/web/packages/sf/vignettes/sf5.html

## Libraries
require(dplyr)
require(sf)
require(tmap)
require(raster)
require(stars)

## Some variables
##dataset <- "gfc2020_70" 
dataset <- "jrc2020"
dir.create(file.path("Maps", dataset, "maps"), recursive=TRUE)
dir_fdb <- "/home/forestatrisk-tropics"

## Countries and continent
data("World")
cont_world <- c("Africa")

## iso3 for Africa
iso3_cont <- World %>%
	st_drop_geometry() %>%
	filter(continent %in% cont_world) %>%
	select(iso_a3) %>%
	pull() %>%
	as.character()
iso3_cont <- c(as.character(iso3_cont), "REU", "MUS")

## COD border
ctry_PROJ_shp <- file.path(dir_fdb, dataset, "Africa", "COD", "data", "ctry_PROJ.shp")
ctry_PROJ <- st_read(ctry_PROJ_shp)
proj <- st_crs(ctry_PROJ)

## Resample r at 500m resolution with gdal
in_f <- file.path(dir_fdb, dataset, "Africa", "COD", "data", "forest", "fcc123.tif")
out_f <- file.path("Maps", dataset, "fcc123_COD_500m.tif")
system(paste0('gdalwarp -r near -tr 500 500 -tap -overwrite \\
							-co "COMPRESS=LZW" -co "PREDICTOR=2" ', in_f, ' ', out_f))

## Raster of historical deforestation 2000-2010-2020
r <- read_stars(out_f)

## Deforestation color
orange <- rgb(255, 165, 0, 255, maxColorValue=255)
red <- rgb(227, 26, 28, 255, maxColorValue=255)
green <- rgb(34, 139, 34, 255, maxColorValue=255)
	
## Plot with tmap
tm <- 
	tm_shape(r) +
	  tmap_options(max.raster=c(plot=1e8, view=1e8)) +
	  tm_raster(palette=c(orange, red, green),
		  				style="cat", legend.show=FALSE) +
  tm_shape(ctry_PROJ) +
	  tm_borders(col="black")

## Save plot
tmap_save(tm, file=file.path("Maps", dataset, "maps", "fcc123_COD.png"))

