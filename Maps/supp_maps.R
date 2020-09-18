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
require(grid)

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
	dplyr::select(iso_a3) %>%
	pull() %>%
	as.character()
# iso3_cont <- c(iso3_cont, "REU", "MUS")

## COD border
ctry_PROJ_shp <- file.path(dir_fdb, dataset, "Africa", "COD", "data", "ctry_PROJ.shp")
ctry_PROJ <- st_read(ctry_PROJ_shp)
proj <- st_crs(ctry_PROJ)

## Equator
eq <- st_linestring(rbind(c(-180,0), c(180,0)))
eq_geom <- st_sfc(eq, crs=4326)
eq_sf <- st_sf(id=1, geometry=eq_geom)

## Tropics
cancer <- st_linestring(rbind(c(-180,23.43661), c(180,23.43661)))
capricorn <- st_linestring(rbind(c(-180,-23.43661), c(180,-23.43661)))
trop_geom <- st_sfc(cancer, capricorn, crs=4326)
trop_df <- data.frame(id=c(1,2),name=c("Cancer","Capricorn"))
trop_sf <- st_sf(trop_df, row.names=trop_df$id, geometry=trop_geom)

## Load GADM level0 data
gadm0 <- st_read(file.path("Maps", "GADM_data", "gadm36_level0.gpkg"))

#=========================================
# Africa with COD border
#=========================================

## Extract countries
gadm0_cont <- gadm0 %>%
	filter(GID_0 %in% iso3_cont)

## Map with tmap
tm_Afr <- 
	#tm_shape(eq_sf) +
	#tm_lines(lty=1,lwd=0.5) +
	#tm_shape(trop_sf) +
	#tm_lines(lty=2, lwd=0.5) +
	tm_shape(gadm0_cont, is.master=TRUE) +
	tm_fill(grey(0.9)) +
	tm_shape(ctry_PROJ) +
	tm_fill(col="black") +
	tm_layout(frame=FALSE)

## Aspect ratio
xy <- st_bbox(gadm0_cont)
asp <- (xy$xmax - xy$xmin)/(xy$ymax - xy$ymin)

## Viewport for inset
w <- 0.25
h <- asp * w
vp <- viewport(x=0.025, y=0.975, width=w, height=h, just=c("left", "top"))

#=========================================
# Historical deforestation map
#=========================================

## Resample r at 500m resolution with gdal
out_f <- file.path("Maps", dataset, "maps", "fcc123_COD_500m.tif")
if (!file.exists(out_f)) {
  in_f <- file.path(dir_fdb, dataset, "Africa", "COD", "data", "forest", "fcc123.tif")
  system(paste0('gdalwarp -r near -tr 500 500 -tap -overwrite \\
							  -co "COMPRESS=LZW" -co "PREDICTOR=2" ', in_f, ' ', out_f))
}

## Raster of historical deforestation 2000-2010-2020
r <- read_stars(out_f)

## Deforestation color
orange <- rgb(255, 165, 0, 255, maxColorValue=255)
red <- rgb(227, 26, 28, 255, maxColorValue=255)
green <- rgb(34, 139, 34, 255, maxColorValue=255)
	
## Plot with tmap
tm_COD_fcc <- 
	tm_shape(r) +
	  tmap_options(max.raster=c(plot=1e8, view=1e8)) +
	  tm_raster(palette=c(orange, red, green),
		  				style="cat", legend.show=FALSE) +
  tm_shape(ctry_PROJ) +
	  tm_borders(col="black")

## Save plot
tmap_save(tm_COD_fcc, file=file.path("Maps", dataset, "maps", "fcc123_COD.png"),
					insets_tm=tm_Afr, insets_vp=vp)

#=========================================
# Zoom with sample points
#=========================================

#' Function to compute zoom extent on same grid as raster
#'
#' @param rast The raster to zoom in
#' @param x_zoom_start x coordinate of the bottom left corner for the zoom
#' @param y_zoom_start y coordinate of the bottom left corner for the zoom
#' @param size_x_start size of the zoom on the x axis (in m)
#' @param size_y_start size of the zoom on the y axis (in m)
#' @return The correct zoom extent for gdal: xmin, ymin, xmax, ymax
zoom_grid <- function(rast, x_zoom_start, y_zoom_start, size_x_start, size_y_start) {
	# rast=r; x_zoom_start=3280000; y_zoom_start=55000; size_x_start=16000; size_y_start=16000
	require(raster)
	r <- raster::raster(rast)
	e <- raster::extent(r)
	res <- raster::res(r)
	r_xmin <- e@xmin
	r_ymin <- e@ymin
	x_zoom <- r_xmin + ((x_zoom_start - r_xmin) %/% res[1]) * res[1]
	y_zoom <- r_ymin + ((y_zoom_start - r_ymin) %/% res[2]) * res[2]
	size_x <- size_x_start - (size_x_start %% res[1])
	size_y <- size_y_start - (size_y_start %% res[2])
	return(paste(x_zoom, y_zoom, x_zoom+size_x, y_zoom+size_y))
}

## Zoom on region with gdal
in_f <- file.path(dir_fdb, dataset, "Africa", "COD", "data", "forest", "fcc123.tif")
out_f <- file.path("Maps", dataset, "maps", "fcc123_COD_zoom.tif")
z_ext <- zoom_grid(rast=raster(in_f), x_zoom_start=3280000, y_zoom_start=55000, 
									 size_x_start=16000, size_y_start=16000)
if (!file.exists(out_f)) {
	system(paste0('gdalwarp -te ', z_ext,' -overwrite \\
							  -co "COMPRESS=LZW" -co "PREDICTOR=2" ', in_f, ' ', out_f))
}

## Import zoom raster
r <- read_stars(out_f)

## Deforestation color
orange <- rgb(255, 165, 0, 130, maxColorValue=255)
red <- rgb(227, 26, 28, 130, maxColorValue=255)
green <- rgb(34, 139, 34, 130, maxColorValue=255)

## Sampled points for COD
in_f_sp <- file.path(dir_fdb, dataset, "Africa", "COD", "output", "sample.txt")
sp_df <- read.table(in_f_sp, header=TRUE, sep=",")
sp_COD <- st_as_sf(sp_df, coords = c("X", "Y"), crs=3395)

## Plot with tmap
tm_COD_zoom <- 
	tm_shape(r) +
	tm_raster(palette=c(orange, red, green),
						style="cat", legend.show=FALSE) +
	tm_shape(ctry_PROJ) +
	tm_borders(col="black") +
	tm_shape(sp_COD) +
	tm_dots(size=1, shape=21)

# EOF