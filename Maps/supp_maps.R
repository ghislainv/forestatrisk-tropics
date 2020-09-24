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
vp_Afr <- viewport(x=0.025, y=0.975, width=w, height=h, just=c("left", "top"))

#=========================================
# Zoom on fcc
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
	return(list(xmin=x_zoom, ymin=y_zoom, xmax=x_zoom+size_x, ymax=y_zoom+size_y))
}

## Zoom on region with gdal
in_f <- file.path(dir_fdb, dataset, "Africa", "COD", "data", "forest", "fcc123.tif")
z_ext <- zoom_grid(rast=raster(in_f), x_zoom_start=3248000, y_zoom_start=23000, 
									 size_x_start=48000, size_y_start=48000)

## Rectangle polygon for zoom
rect_ext <- extent(z_ext$xmin, z_ext$xmax, z_ext$ymin, z_ext$ymax)
rect_bbox <- st_bbox(rect_ext, crs=3395)
rect_geom <- st_as_sfc(rect_bbox)
rect <- st_sf(id=1, geometry=rect_geom)

## Zoom on region with gdal
out_f <- file.path("Maps", dataset, "maps", "fcc123_COD_zoom.tif")
z_ext_gdal <- paste(z_ext$xmin, z_ext$ymin, z_ext$xmax, z_ext$ymax)
if (!file.exists(out_f)) {
	system(paste0('gdalwarp -te ', z_ext_gdal,' -overwrite \\
							  -co "COMPRESS=LZW" -co "PREDICTOR=2" ', in_f, ' ', out_f))
}

## Deforestation color
orange <- rgb(255, 165, 0, 255, maxColorValue=255)
red <- rgb(227, 26, 28, 255, maxColorValue=255)
green <- rgb(34, 139, 34, 255, maxColorValue=255)

## Plot
r_zoom <- read_stars(out_f)
tm_COD_fcc_zoom <- 
	tm_shape(r_zoom) +
	tmap_options(max.raster=c(plot=1e8, view=1e8)) +
  tm_raster(palette=c(orange, red, green),
  					style="cat", legend.show=FALSE)

## Aspect ratio
asp <- (z_ext$xmax - z_ext$xmin)/(z_ext$ymax - z_ext$ymin)

## Viewport for inset
w <- 0.25
h <- asp * w
vp_COD_fcc_zoom <- viewport(x=0.025, y=0.025, width=w, height=h, just=c("left", "bottom"))

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

## Plot with tmap
tm_COD_fcc <- 
	tm_shape(r) +
	  tmap_options(max.raster=c(plot=1e8, view=1e8)) +
	  tm_raster(palette=c(orange, red, green),
		  				style="cat", legend.show=FALSE) +
  tm_shape(ctry_PROJ) +
	  tm_borders(col="black") +
	tm_shape(rect) +
	  tm_borders(col="black", lwd=2)

## Save plot
f <- file.path("Maps", dataset, "maps", "fcc123_COD.png")
tmap_save(tm_COD_fcc, file=f,
					insets_tm=list(tm_Afr,tm_COD_fcc_zoom), insets_vp=list(vp_Afr,vp_COD_fcc_zoom))

## Copy for manuscript
f_doc <- file.path("Manuscript", "Supplementary_Materials", "figures", "fcc123_COD.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

#=========================================
# Extra zoom on fcc
#=========================================

## Zoom on region with gdal
in_f <- file.path(dir_fdb, dataset, "Africa", "COD", "data", "forest", "fcc123.tif")
z_ext <- zoom_grid(rast=raster(in_f), x_zoom_start=3274440, y_zoom_start=60800, 
									 size_x_start=1795, size_y_start=1795)

## Rectangle polygon for zoom
rect_ext <- extent(z_ext$xmin, z_ext$xmax, z_ext$ymin, z_ext$ymax)
rect_bbox <- st_bbox(rect_ext, crs=3395)
rect_geom <- st_as_sfc(rect_bbox)
rect <- st_sf(id=1, geometry=rect_geom)

## Zoom on region with gdal
out_f <- file.path("Maps", dataset, "maps", "fcc123_COD_zoom_extra.tif")
z_ext_gdal <- paste(z_ext$xmin, z_ext$ymin, z_ext$xmax, z_ext$ymax)
if (!file.exists(out_f)) {
	system(paste0('gdalwarp -te ', z_ext_gdal,' -overwrite \\
							  -co "COMPRESS=LZW" -co "PREDICTOR=2" ', in_f, ' ', out_f))
}

## Deforestation color
orange <- rgb(255, 165, 0, 255, maxColorValue=255)
red <- rgb(227, 26, 28, 255, maxColorValue=255)
green <- rgb(34, 139, 34, 255, maxColorValue=255)

## Deforestation color with transparency
orange_transp <- rgb(255, 165, 0, 130, maxColorValue=255)
red_transp <- rgb(227, 26, 28, 130, maxColorValue=255)
green_transp <- rgb(34, 139, 34, 130, maxColorValue=255)

## Sampled points for COD
in_f_sp <- file.path(dir_fdb, dataset, "Africa", "COD", "output", "sample.txt")
sp_df <- read.table(in_f_sp, header=TRUE, sep=",")
sp_COD <- st_as_sf(sp_df, coords = c("X", "Y"), crs=3395)

## Plot
in_f <- file.path("Maps", dataset, "maps", "fcc123_COD_zoom_extra.tif")
r_zoom_extra <- read_stars(in_f)
tm_COD_fcc_zoom_extra <- 
	tm_shape(r_zoom_extra) +
	  tmap_options(max.raster=c(plot=1e8, view=1e8)) +
	  tm_raster(palette=c(orange_transp, red_transp, green_transp),
						  style="cat", legend.show=FALSE) +
	tm_shape(sp_COD) +
	  tm_dots(col="fcc23", size=0.25, shape=21, style="cat", n=2, 
					  palette=c(red, green), legend.show=FALSE)

#=========================================
# Zoom with sample points
#=========================================

## Import zoom raster
in_f <- file.path("Maps", dataset, "maps", "fcc123_COD_zoom.tif")
r <- read_stars(in_f)

## Plot with tmap
tm_COD_fcc_zoom <- 
	tm_shape(r) +
	  tmap_options(max.raster=c(plot=1e8, view=1e8)) +
	  tm_raster(palette=c(orange_transp, red_transp, green_transp),
						  style="cat", legend.show=FALSE) +
	tm_shape(rect) +
	  tm_borders(col="black", lwd=2) +
	tm_shape(sp_COD) +
	  tm_dots(col="fcc23", size=0.25, shape=21, style="cat", n=2, 
					  palette=c(red, green), legend.show=FALSE)

## Arrange plot with grid package
f <- file.path("Maps", dataset, "maps", "sample_COD.png")
png(filename=f, width=1000, height=500)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2)))
print(tm_COD_fcc_zoom, vp=viewport(layout.pos.col=1))
print(tm_COD_fcc_zoom_extra, vp=viewport(layout.pos.col=2))
dev.off()

## Copy for manuscript
f_doc <- file.path("Manuscript", "Supplementary_Materials", "figures", "sample_COD.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

# EOF