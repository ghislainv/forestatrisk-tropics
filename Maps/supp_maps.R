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
require(here)

## Declare some variables
##dataset <- "gfc2020_70" 
dataset <- "jrc2020"
ctry_iso <- "COD"

## Load country info (encoding pb for ctry names)
ctry_df <- read.csv2(here("Analysis", "ctry_run.csv"), header=TRUE, sep=";", encoding="UTF-8")

## Identify continent
cont <- as.character(ctry_df$cont_run[ctry_df$iso3==ctry_iso])

## Create directory for country maps
dir.create(here("Maps", dataset, ctry_iso), recursive=TRUE)

## Source data directory
#dir_fdb <- "/home/forestatrisk-tropics"
dir_fdb <- here("Data")

## Countries and continent
data("World")

## iso3 for continent
if (cont=="Africa") {cont_world <- c("Africa")}
if (cont=="America") {cont_world <- c("North America", "South America")}
if (cont=="Asia") {cont_world <- c("Asia", "Oceania")}
iso3_cont <- World %>%
	st_drop_geometry() %>%
	filter(continent %in% cont_world) %>%
	dplyr::select(iso_a3) %>%
	pull() %>%
	as.character()
# iso3_cont <- c(iso3_cont, "REU", "MUS")

## Country border
ctry_PROJ_shp <- file.path(dir_fdb, dataset, cont, ctry_iso, "data", "ctry_PROJ.shp")
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

## Deforestation color
orange <- rgb(255, 165, 0, 255, maxColorValue=255)
red <- rgb(227, 26, 28, 255, maxColorValue=255)
green <- rgb(34, 139, 34, 255, maxColorValue=255)

## tmap_options
tmap_opt <- function(npix=1e5, ...) {
  tmap_options_reset()
  tmap_options(max.raster=c(plot=npix, view=npix), ...)
  cat("max.raster set to: ", npix, "\n")
}

## textwidth (in cm) for figure width
textwidth <- 16.6

## Load GADM level0 data
gadm0 <- st_read(here("Maps", "GADM_data", "gadm36_level0.gpkg"))

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
w <- 0.275
h <- asp * w
vp_Afr <- viewport(x=0.01, y=0.01, width=w, height=h, just=c("left", "bottom"))

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
in_f <- file.path(dir_fdb, dataset, cont, ctry_iso, "data", "forest", "fcc123.tif")
z_ext <- zoom_grid(rast=raster(in_f), x_zoom_start=3248000, y_zoom_start=23000, 
									 size_x_start=48000, size_y_start=48000)

## Rectangle polygon for zoom
rect_ext <- extent(z_ext$xmin, z_ext$xmax, z_ext$ymin, z_ext$ymax)
rect_bbox <- st_bbox(rect_ext, crs=3395)
rect_geom <- st_as_sfc(rect_bbox)
rect <- st_sf(id=1, geometry=rect_geom)

## Zoom on region with gdal
out_f <- here("Maps", dataset, ctry_iso, "fcc123_zoom.tif")
z_ext_gdal <- paste(z_ext$xmin, z_ext$ymin, z_ext$xmax, z_ext$ymax)
if (!file.exists(out_f)) {
	system(paste0('gdalwarp -te ', z_ext_gdal,' -overwrite \\
							  -co "COMPRESS=LZW" -co "PREDICTOR=2" ', in_f, ' ', out_f))
}

## Plot
r_zoom <- read_stars(out_f)
tm_fcc_zoom <- 
	tm_shape(r_zoom) +
  tm_raster(palette=c(orange, red, green),
  					style="cat", legend.show=FALSE) +
	tm_scale_bar(breaks=c(0,5,10), text.size=1,
	             position=c(0.5,0), just=c("center", "bottom"))

## Viewport for inset
vp_fcc_zoom <- viewport(x=0.01, y=0.99, width=0.275, height=0.275, just=c("left", "top"))

#=========================================
# Historical deforestation map
#=========================================

## Resample r at 500m resolution with gdal
in_f <- file.path(dir_fdb, dataset, cont, ctry_iso, "data", "forest", "fcc123.tif")
out_f <- here("Maps", dataset, ctry_iso, "fcc123_500m.tif")
if (!file.exists(out_f)) {
  system(paste0('gdalwarp -r near -tr 500 500 -tap -overwrite \\
							  -co "COMPRESS=LZW" -co "PREDICTOR=2" ', in_f, ' ', out_f))
}

## Raster of historical deforestation 2000-2010-2020
r <- read_stars(out_f)

## Aspect
bbox_r <- st_bbox(r)
asp_fcc <- (bbox_r$xmax - bbox_r$xmin)/(bbox_r$ymax - bbox_r$ymin)

## Plot with tmap
tm_fcc <- 
	tm_shape(r) +
	  tm_raster(palette=c(orange, red, green),
		  				style="cat", legend.show=FALSE) +
  tm_shape(ctry_PROJ) +
	  tm_borders(col="black") +
	tm_shape(rect) +
	  tm_borders(col="black", lwd=2) +
	tm_scale_bar(c(0,250,500), text.size=1,
	             position=c(0.5,0), just=c("center", "bottom"))

## Save plot
tmap_opt(1e8, outer.margins=c(0,0,0,0))
f <- here("Maps", dataset, ctry_iso, "fcc123.png")
tmap_save(tm_fcc, file=f, width=textwidth, height=textwidth/asp_fcc, units="cm", dpi=300,
					insets_tm=list(tm_Afr,tm_fcc_zoom), insets_vp=list(vp_Afr,vp_fcc_zoom))

## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "figures", "fcc123.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

#=========================================
# Extra zoom on fcc
#=========================================

## Zoom on region with gdal
in_f <- file.path(dir_fdb, dataset, cont, ctry_iso, "data", "forest", "fcc123.tif")
z_ext <- zoom_grid(rast=raster(in_f), x_zoom_start=3274440, y_zoom_start=60800, 
									 size_x_start=1795, size_y_start=1795)

## Rectangle polygon for zoom
rect_ext <- extent(z_ext$xmin, z_ext$xmax, z_ext$ymin, z_ext$ymax)
rect_bbox <- st_bbox(rect_ext, crs=3395)
rect_geom <- st_as_sfc(rect_bbox)
rect_extra <- st_sf(id=1, geometry=rect_geom)

## Zoom on region with gdal
out_f <- here("Maps", dataset, ctry_iso, "fcc123_zoom_extra.tif")
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
in_f_sp <- file.path(dir_fdb, dataset, cont, ctry_iso, "output", "sample.txt")
sp_df <- read.table(in_f_sp, header=TRUE, sep=",")
sp <- st_as_sf(sp_df, coords = c("X", "Y"), crs=3395)

## Plot
in_f <- here("Maps", dataset, ctry_iso, "fcc123_zoom_extra.tif")
r_zoom_extra <- read_stars(in_f)
tm_fcc_zoom_extra <- 
	tm_shape(r_zoom_extra) +
	  tmap_options(max.raster=c(plot=1e8, view=1e8)) +
	  tm_raster(palette=c(orange_transp, red_transp, green_transp),
						  style="cat", legend.show=FALSE) +
	tm_shape(sp) +
	  tm_dots(col="fcc23", size=0.2, shape=21, style="cat", n=2, 
					  palette=c(red, green), legend.show=FALSE) +
	tm_scale_bar(c(0,0.5,1), text.size=0.8,
	             position=c(0.5,0.14), just=c("center", "top")) +
  tm_layout(outer.margins=c(0,0.015,0,0))

#=========================================
# Zoom with sample points
#=========================================

## Import zoom raster
in_f <- here("Maps", dataset, ctry_iso, "fcc123_zoom.tif")
r <- read_stars(in_f)

## Plot with tmap
tm_fcc_zoom <- 
	tm_shape(r) +
	  tmap_options(max.raster=c(plot=1e8, view=1e8)) +
	  tm_raster(palette=c(orange_transp, red_transp, green_transp),
						  style="cat", legend.show=FALSE) +
	tm_shape(sp) +
	  tm_dots(col="fcc23", size=0.2, shape=21, style="cat", n=2, 
					  palette=c(red, green), legend.show=FALSE) +
	tm_shape(rect_extra) +
	  tm_borders(col="black", lwd=2) +
	tm_scale_bar(c(0,5,10), text.size=0.8,
	             position=c(0.5,0.14), just=c("center", "top")) +
  tm_layout(outer.margins=c(0,0,0,0.015))

## Arrange plots with grid package
f <- here("Maps", dataset, ctry_iso, "sample.png")
png(filename=f, width=textwidth, height=textwidth/2, units="cm", res=300)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2)))
print(tm_fcc_zoom, vp=viewport(layout.pos.col=1))
print(tm_fcc_zoom_extra, vp=viewport(layout.pos.col=2))
dev.off()

## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "figures", "sample.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

#=========================================
# Grid for spatial random effects
#=========================================

## File paths
in_fcc_500 <- here("Maps", dataset, ctry_iso, "fcc123_500m.tif")
in_fcc_zoom <- here("Maps", dataset, ctry_iso, "fcc123_zoom.tif")
in_fcc123 <- file.path(dir_fdb, dataset, cont, ctry_iso, "data",
                       "forest", "fcc123.tif")

## Rasters
fcc_500 <- read_stars(in_fcc_500)
fcc_zoom <- read_stars(in_fcc_zoom)
fcc123 <- read_stars(in_fcc123)

## Grid
g <- st_make_grid(st_as_sfc(st_bbox(fcc123), crs=3395), cellsize=10000)

## Small grid for diagram
g_neigh <- st_sf(target=c(rep(0,3), c(0,1,0), rep(0,3)),
                 geometry=g[c(33147:33149, 33361:33363, 33575:33577)])

## Zoom
in_f <- file.path(dir_fdb, dataset, cont, ctry_iso, "data", "forest", "fcc123.tif")
z_ext <- zoom_grid(rast=raster(in_f), x_zoom_start=3248000, y_zoom_start=23000, 
									 size_x_start=48000, size_y_start=48000)

## Plot
tm_rho_grid <- 
	tm_shape(fcc_500) +
	tmap_options(max.raster=c(plot=1e8, view=1e8)) +
	tm_raster(palette=c(orange, red, green),
						style="cat", legend.show=FALSE) +
	tm_shape(ctry_PROJ) +
	tm_borders(col="black") +
	tm_shape(g) +
	tm_borders(col="black", lwd=0.1) +
	tm_shape(rect) +
	tm_borders(col="black", lwd=2) + 
	tm_scale_bar(c(0,250,500), text.size=1,
	             position=c(0.5,0), just=c("left", "bottom")) +
  tm_layout(outer.margins=c(0,0,0,0))

## Zoom
tm_rho_grid_zoom <- 
	tm_shape(fcc_zoom) +
	tmap_options(max.raster=c(plot=1e8, view=1e8)) +
	tm_raster(palette=c(orange_transp, red_transp, green_transp),
						style="cat", legend.show=FALSE) +
	tm_shape(sp) +
	tm_dots(col="fcc23", size=0.1, shape=21, style="cat", n=2, 
					palette=c(red, green), legend.show=FALSE) +
	tm_shape(g) +
	tm_borders(col="black", lwd=1.5) + 
	tm_scale_bar(c(0,5,10), text.size=1,
	             position=c(0.5,0), just=c("center", "bottom")) +
  tm_layout(outer.margins=c(0,0,0,0))

## Diagram
tm_grid_diag <- 
	tm_shape(g, bbox=st_bbox(fcc_zoom)) +
	  tm_borders(col=grey(0.7), lwd=1.5) + 
  tm_shape(g_neigh, bbox=st_bbox(fcc_zoom)) +
	  tm_fill(col="target", style="cat", n=2, title="", palette=c(grey(0.9), grey(0.7)),
	          labels=c("Neighbouring cells with rho_j'", "Target cell with rho_j")) +
    tm_borders(col="black", lwd=1.5) +
  tm_layout(outer.margins=c(0,0,0,0)) +
  tm_legend(legend.text.size=0.8,
            legend.position=c("center","bottom"),
            legend.just=c("center","bottom"), legend.width=1)

## Save plot
vp_grid_zoom <- viewport(x=0.01, y=0.99, width=0.45, height=0.45, just=c("left", "top"))
vp_grid_diag <- viewport(x=0.01, y=0.01, width=0.45, height=0.45, just=c("left", "bottom"))
f <- here("Maps", dataset, ctry_iso, "grid.png")
tmap_save(tm_rho_grid, file=f, width=textwidth, height=textwidth/asp_fcc, units="cm", dpi=300,
					insets_tm=list(tm_rho_grid_zoom, tm_grid_diag), insets_vp=list(vp_grid_zoom, vp_grid_diag))

## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "figures", "grid.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

#=======================================
# Spatial random effects
#=======================================

## File paths
in_rho_orig <- file.path(dir_fdb, dataset, cont, ctry_iso, "output", "rho_orig.tif")
in_rho <- file.path(dir_fdb, dataset, cont, ctry_iso, "output", "rho.tif")

## Rasters
rho_orig <- read_stars(in_rho_orig)
rho <- read_stars(in_rho)

## Ncells for COD
ncell <- ncell(rho_orig)
ncell

## Quantiles
q <- quantile(rho_orig[[1]], c(0.01, 0.99))
q_max <- ceiling(max(abs(q)))

## Colors
red_rho <- "#cc3300"
yellow_rho <- "#ffffbf"
green_rho <- green # "#228B22FF"

## tm_options
tmap_opt(npix=1e8)

## Plot rho_orig
tm_rho_orig <- 
	tm_shape(rho_orig) +
	tm_raster(palette=c(green_rho, yellow_rho, red_rho), title="",
						legend.reverse=TRUE,
						style="cont", midpoint=0, breaks=seq(-q_max, q_max, b=1)) +
	tm_shape(rect) +
	tm_borders(col="black", lwd=2) +
	tm_shape(ctry_PROJ) +
	tm_borders(col="black") +
	tm_scale_bar(c(0,250,500), text.size=0.8,
	             position=c(0.5,0.14), just=c("center", "top")) +
  tm_layout(outer.margins=c(0.015,0,0,0.015),
            legend.position=c("left", "bottom"),
            legend.just=c("left","bottom"))
## Zoom rho_orig
tm_rho_orig_zoom <- 
	tm_shape(rho_orig, bbox=st_bbox(rect)) +
	tm_raster(palette=c(green_rho, yellow_rho, red_rho), title="",
						legend.reverse=TRUE,
						style="cont", midpoint=0, breaks=seq(-q_max, q_max, b=1)) +
	tm_shape(ctry_PROJ) +
	tm_borders(col="black") +
	tm_scale_bar(c(0,5,10), text.size=0.8,
	             position=c(0.5,0.14), just=c("center", "top")) +
  tm_layout(outer.margins=c(0,0,0.015,0.015),
            legend.position=c("left", "bottom"),
            legend.just=c("left","bottom"))

## Plot interpolated rho
tm_rho <- 
	tm_shape(rho) +
	tm_raster(palette=c(green_rho, yellow_rho, red_rho), title="",
						legend.reverse=TRUE, legend.show=FALSE,
						style="cont", midpoint=0, breaks=seq(-q_max, q_max, b=1)) +
	tm_shape(rect) +
	tm_borders(col="black", lwd=2) +
	tm_shape(ctry_PROJ) +
	tm_borders(col="black") +
  tm_scale_bar(c(0,250,500),  text.size=0.8,
	             position=c(0.5,0.14), just=c("center", "top")) +
  tm_layout(outer.margins=c(0.015,0.015,0,0))
## Zoom interpolated rho
tm_rho_zoom <- 
	tm_shape(rho, bbox=st_bbox(rect)) +
	tm_raster(palette=c(green_rho, yellow_rho, red_rho), title="",
						legend.reverse=TRUE, legend.show=FALSE,
						style="cont", midpoint=0, breaks=seq(-q_max, q_max, b=1)) +
	tm_scale_bar(c(0,5,10),  text.size=0.8,
	             position=c(0.5,0.14), just=c("center", "top")) +
  tm_layout(outer.margins=c(0,0.015,0.015,0))

## Arrange plots with grid package
f <- here("Maps", dataset, ctry_iso, "rho.png")
png(filename=f, width=textwidth, height=textwidth, units="cm", res=300)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,2)))
print(tm_rho_orig, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(tm_rho_orig_zoom, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(tm_rho, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(tm_rho_zoom, vp=viewport(layout.pos.row=2, layout.pos.col=2))
dev.off()

## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "figures", "rho.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

#===========================================
# Spatial probability of deforestation: zoom
#===========================================

## Zoom on region with gdal
in_f <- file.path(dir_fdb, dataset, cont, ctry_iso, "output", "prob.tif")
z_ext <- zoom_grid(rast=raster(in_f), x_zoom_start=3248000, y_zoom_start=23000, 
									 size_x_start=48000, size_y_start=48000)

## Rectangle polygon for zoom
rect_ext <- extent(z_ext$xmin, z_ext$xmax, z_ext$ymin, z_ext$ymax)
rect_bbox <- st_bbox(rect_ext, crs=3395)
rect_geom <- st_as_sfc(rect_bbox)
rect <- st_sf(id=1, geometry=rect_geom)

## Zoom on region with gdal
out_f <- here("Maps", dataset, ctry_iso, "prob_zoom.tif")
z_ext_gdal <- paste(z_ext$xmin, z_ext$ymin, z_ext$xmax, z_ext$ymax)
if (!file.exists(out_f)) {
	system(paste0('gdalwarp -te ', z_ext_gdal,' -overwrite \\
							  -co "COMPRESS=LZW" -co "PREDICTOR=2" ', in_f, ' ', out_f))
}

## Plot
r_zoom <- read_stars(out_f)
tm_prob_zoom <- 
	tm_shape(r_zoom) +
	  tm_raster(style="cont", title="", legend.reverse=TRUE, legend.show=FALSE,
	            palette=c("#228b22", "#ffa500", "#e31a1c", "#000000"),
	            breaks=c(1, 39322, 54249, 65535), labels=c("0","","","1")) +
	tm_scale_bar(breaks=c(0,5,10), text.size=1,
	             position=c(0.5,0), just=c("center", "bottom"))

## Viewport for inset
vp_prob_zoom <- viewport(x=0.01, y=0.99, width=0.275, height=0.275,
                         just=c("left", "top"))

#=====================================
# Spatial probability of deforestation
#=====================================

## Resample r at 500m resolution with gdal
in_f <- file.path(dir_fdb, dataset, cont, ctry_iso, "output", "prob.tif")
out_f <- here("Maps", dataset, ctry_iso, "prob_500m.tif")
## Use near resampling for plot (including min/max), not bilinear
if (!file.exists(out_f)) {
  system(paste0('gdalwarp -r near -ot UInt16 -tr 500 500 -tap -overwrite \\
							  -co "COMPRESS=LZW" -co "PREDICTOR=2" ', in_f, ' ', out_f))
}

## Raster of historical deforestation 2000-2010-2020
r <- read_stars(out_f)

## Aspect
bbox_r <- st_bbox(r)
asp_prob <- (bbox_r$xmax - bbox_r$xmin)/(bbox_r$ymax - bbox_r$ymin)

## Plot with tmap
tm_prob <- 
	tm_shape(r) +
	  tm_raster(style="cont", title="", legend.reverse=TRUE,
	            palette=c("#228b22", "#ffa500", "#e31a1c", "#000000"),
	            breaks=c(1, 39322, 54249, 65535), labels=c("0","","","1")) +
  tm_shape(ctry_PROJ) +
	  tm_borders(col="black") +
	tm_shape(rect) +
	  tm_borders(col="black", lwd=2) +
	tm_scale_bar(c(0,250,500), text.size=1,
	             position=c(0.5,0), just=c("center", "bottom")) +
  tm_legend(position=c("left","bottom"), just=c("left","bottom")) +
  tm_layout(legend.text.size=1.2)

## Save plot
tmap_opt(1e8, outer.margins=c(0,0,0,0))
f <- here("Maps", dataset, ctry_iso, "prob.png")
tmap_save(tm_prob, file=f, width=textwidth, height=textwidth/asp_prob, units="cm", dpi=300,
					insets_tm=tm_prob_zoom, insets_vp=vp_prob_zoom)

## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "figures", "prob.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

#==================================
# Future forest cover in 2050: zoom
#==================================

## Zoom on region with gdal
in_f <- file.path(dir_fdb, dataset, cont, ctry_iso, "output", "fcc_2050.tif")
z_ext <- zoom_grid(rast=raster(in_f), x_zoom_start=3248000, y_zoom_start=23000, 
									 size_x_start=48000, size_y_start=48000)

## Rectangle polygon for zoom
rect_ext <- extent(z_ext$xmin, z_ext$xmax, z_ext$ymin, z_ext$ymax)
rect_bbox <- st_bbox(rect_ext, crs=3395)
rect_geom <- st_as_sfc(rect_bbox)
rect <- st_sf(id=1, geometry=rect_geom)

## Zoom on region with gdal
out_f <- here("Maps", dataset, ctry_iso, "fcc2050_zoom.tif")
z_ext_gdal <- paste(z_ext$xmin, z_ext$ymin, z_ext$xmax, z_ext$ymax)
if (!file.exists(out_f)) {
	system(paste0('gdalwarp -te ', z_ext_gdal,' -overwrite \\
							  -co "COMPRESS=LZW" -co "PREDICTOR=2" ', in_f, ' ', out_f))
}

## Plot
r_zoom <- read_stars(out_f)
tm_prob_zoom <- 
	tm_shape(r_zoom) +
	  tm_raster(style="cat", n=2, title="", legend.show=FALSE,
	            palette=c(red, green)) +
	tm_scale_bar(breaks=c(0,5,10), text.size=1,
	             position=c(0.5,0), just=c("center", "bottom"))

## Viewport for inset
vp_prob_zoom <- viewport(x=0.01, y=0.99, width=0.275, height=0.275,
                         just=c("left", "top"))

#============================
# Future forest cover in 2050
#============================

## Resample r at 500m resolution with gdal
in_f <- file.path(dir_fdb, dataset, cont, ctry_iso, "output", "fcc_2050.tif")
out_f <- here("Maps", dataset, ctry_iso, "fcc2050_500m.tif")
if (!file.exists(out_f)) {
  system(paste0('gdalwarp -r near -tr 500 500 -tap -overwrite \\
							  -co "COMPRESS=LZW" -co "PREDICTOR=2" ', in_f, ' ', out_f))
}

## Raster of historical deforestation 2000-2010-2020
r <- read_stars(out_f)

## Aspect
bbox_r <- st_bbox(r)
asp_prob <- (bbox_r$xmax - bbox_r$xmin)/(bbox_r$ymax - bbox_r$ymin)

## Plot with tmap
tm_fcc2050 <- 
	tm_shape(r) +
	  tm_raster(style="cat", n=2, legend.show=FALSE,
	            palette=c(red, green)) +
  tm_shape(ctry_PROJ) +
	  tm_borders(col="black") +
	tm_shape(rect) +
	  tm_borders(col="black", lwd=2) +
	tm_scale_bar(c(0,250,500), text.size=1,
	             position=c(0.5,0), just=c("center", "bottom")) +
  tm_legend(position=c("left","bottom"), just=c("left","bottom")) +
  tm_layout(legend.text.size=1.2, title="2050",
            title.position=c(0.06,0.06), title.size=1.8)

## Save plot
tmap_opt(1e8, outer.margins=c(0,0,0,0))
f <- here("Maps", dataset, ctry_iso, "fcc2050.png")
tmap_save(tm_fcc2050, file=f, width=textwidth, height=textwidth/asp_prob, units="cm", dpi=300,
					insets_tm=tm_prob_zoom, insets_vp=vp_prob_zoom)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "figures", "fcc2050.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

#==================================
# Future forest cover in 2100: zoom
#==================================

## Zoom on region with gdal
in_f <- file.path(dir_fdb, dataset, cont, ctry_iso, "output", "fcc_2100.tif")
z_ext <- zoom_grid(rast=raster(in_f), x_zoom_start=3248000, y_zoom_start=23000, 
									 size_x_start=48000, size_y_start=48000)

## Rectangle polygon for zoom
rect_ext <- extent(z_ext$xmin, z_ext$xmax, z_ext$ymin, z_ext$ymax)
rect_bbox <- st_bbox(rect_ext, crs=3395)
rect_geom <- st_as_sfc(rect_bbox)
rect <- st_sf(id=1, geometry=rect_geom)

## Zoom on region with gdal
out_f <- here("Maps", dataset, ctry_iso, "fcc2100_zoom.tif")
z_ext_gdal <- paste(z_ext$xmin, z_ext$ymin, z_ext$xmax, z_ext$ymax)
if (!file.exists(out_f)) {
	system(paste0('gdalwarp -te ', z_ext_gdal,' -overwrite \\
							  -co "COMPRESS=LZW" -co "PREDICTOR=2" ', in_f, ' ', out_f))
}

## Plot
r_zoom <- read_stars(out_f)
tm_prob_zoom <- 
	tm_shape(r_zoom) +
	  tm_raster(style="cat", n=2, title="", legend.show=FALSE,
	            palette=c(red, green)) +
	tm_scale_bar(breaks=c(0,5,10), text.size=1,
	             position=c(0.5,0), just=c("center", "bottom"))

## Viewport for inset
vp_prob_zoom <- viewport(x=0.01, y=0.99, width=0.275, height=0.275,
                         just=c("left", "top"))

#============================
# Future forest cover in 2100
#============================

## Resample r at 500m resolution with gdal
in_f <- file.path(dir_fdb, dataset, cont, ctry_iso, "output", "fcc_2100.tif")
out_f <- here("Maps", dataset, ctry_iso, "fcc2100_500m.tif")
if (!file.exists(out_f)) {
  system(paste0('gdalwarp -r near -tr 500 500 -tap -overwrite \\
							  -co "COMPRESS=LZW" -co "PREDICTOR=2" ', in_f, ' ', out_f))
}

## Raster of historical deforestation 2000-2010-2020
r <- read_stars(out_f)

## Aspect
bbox_r <- st_bbox(r)
asp_prob <- (bbox_r$xmax - bbox_r$xmin)/(bbox_r$ymax - bbox_r$ymin)

## Plot with tmap
tm_fcc2100 <- 
	tm_shape(r) +
	  tm_raster(style="cat", n=2, legend.show=FALSE,
	            palette=c(red, green)) +
  tm_shape(ctry_PROJ) +
	  tm_borders(col="black") +
	tm_shape(rect) +
	  tm_borders(col="black", lwd=2) +
	tm_scale_bar(c(0,250,500), text.size=1,
	             position=c(0.5,0), just=c("center", "bottom")) +
  tm_legend(position=c("left","bottom"), just=c("left","bottom")) +
  tm_layout(legend.text.size=1.2, title="2100",
            title.position=c(0.06,0.06), title.size=1.8)

## Save plot
tmap_opt(1e8, outer.margins=c(0,0,0,0))
f <- here("Maps", dataset, ctry_iso, "fcc2100.png")
tmap_save(tm_fcc2100, file=f, width=textwidth, height=textwidth/asp_prob, units="cm", dpi=300,
					insets_tm=tm_prob_zoom, insets_vp=vp_prob_zoom)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "figures", "fcc2100.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## Montage with ImageMagick
f1 <- here("Manuscript", "Supplementary_Materials", "figures", "fcc2050.png")
f2 <- here("Manuscript", "Supplementary_Materials", "figures", "fcc2100.png")
f_out <- here("Manuscript", "Supplementary_Materials", "figures",
              "fcc2050_2100.png")
system(paste0("montage ",f1," ",f2," -tile 2x1 -geometry +25 ",f_out))
system(paste0("convert ",f_out," -crop +25+0 +repage ",f_out))
system(paste0("convert ",f_out," -crop -25+0 +repage ",f_out))

#======================
# Explanatory variables
#======================

## Shapefiles
ctry_PROJ_shp <- file.path(dir_fdb, dataset, cont, ctry_iso, "data", "ctry_PROJ.shp")
pa_shp <- file.path(dir_fdb, dataset, cont, ctry_iso, "data", "pa_PROJ.shp")
roads_shp <- file.path(dir_fdb, dataset, cont, ctry_iso, "data", "roads_PROJ.shp")
towns_shp <- file.path(dir_fdb, dataset, cont, ctry_iso, "data", "towns_PROJ.shp")
rivers_shp <- file.path(dir_fdb, dataset, cont, ctry_iso, "data", "rivers_PROJ.shp")

## Altitude
in_f <- file.path(dir_fdb, dataset, cont, ctry_iso, "data", "altitude.tif")
out_f1 <- here("Maps", dataset, ctry_iso, "elev_500m.tif")
out_f2 <- here("Maps", dataset, ctry_iso, "elev_500m_crop.tif")
if (!file.exists(out_f2)) {
  system(paste0('gdalwarp -r near \\
                -ot Int16 -tr 500 500 -tap -overwrite \\
							  -co "COMPRESS=LZW" -co "PREDICTOR=2" ', in_f, ' ', out_f1))
  system(paste0('gdalwarp -cutline ', ctry_PROJ_shp,' -crop_to_cutline \\
                -overwrite \\
							  -co "COMPRESS=LZW" -co "PREDICTOR=2" ', out_f1, ' ', out_f2))
}
## Slope
in_f <- file.path(dir_fdb, dataset, cont, ctry_iso, "data", "slope.tif")
out_f1 <- here("Maps", dataset, ctry_iso, "slope_500m.tif")
out_f2 <- here("Maps", dataset, ctry_iso, "slope_500m_crop.tif")
if (!file.exists(out_f2)) {
  system(paste0('gdalwarp -r near \\
                -ot Int16 -tr 500 500 -tap -overwrite \\
							  -co "COMPRESS=LZW" -co "PREDICTOR=2" ', in_f, ' ', out_f1))
  system(paste0('gdalwarp -cutline ', ctry_PROJ_shp,' -crop_to_cutline \\
                -overwrite \\
							  -co "COMPRESS=LZW" -co "PREDICTOR=2" ', out_f1, ' ', out_f2))
}

# Load as sf object
ctry_PROJ <- st_read(ctry_PROJ_shp)
pa <- st_read(pa_shp)
roads <- st_read(roads_shp)
towns <- st_read(towns_shp)
rivers <- st_read(rivers_shp)

# Load as stars object
elev_tif <- here("Maps", dataset, ctry_iso, "elev_500m_crop.tif")
elev <- read_stars(elev_tif)
slope_tif <- here("Maps", dataset, ctry_iso, "slope_500m_crop.tif")
slope <- read_stars(slope_tif)

## Elevation
elev_pal <- colorRampPalette(c("darkolivegreen4","yellow","brown"))(6)
tm_elev <- 
  tm_shape(elev) +
    tm_raster(style="quantile", n=6, title="", midpoint=NA,
              palette=elev_pal) +
  tm_shape(ctry_PROJ, is.master=TRUE) +
	  tm_borders(col="black", lwd=0.5) +
  tm_legend(position=c("left", "top"),
            just=c("left","bottom"),
            text.size=0.4) +
  tm_layout(outer.margins=c(0,0,0,0.015),
            title="Elevation (m)",
            title.position=c(0.06,0.06),
            title.size=0.8)

## Slope
slope_pal <- colorRampPalette(grey((9:1)/10))(6)
tm_slope <-
  tm_shape(slope) +
    tm_raster(style="quantile", n=6, title="", midpoint=NA,
              palette=slope_pal) +
  tm_shape(ctry_PROJ, is.master=TRUE) +
	  tm_borders(col="black", lwd=0.5) +
  tm_legend(position=c("left", "top"),
            just=c("left","bottom"),
            text.size=0.4) +
  tm_layout(outer.margins=c(0,0,0,0.015),
            title="Slope (°)",
            title.position=c(0.06,0.06),
            title.size=0.8)

## Roads
tm_roads <- 
  tm_shape(roads) +
    tm_lines(lwd=0.25, col="orange4") +
  tm_shape(ctry_PROJ, is.master=TRUE) +
	  tm_borders(col="black", lwd=0.5) +
  tm_layout(outer.margins=c(0,0,0,0.015),
            title="Roads",
            title.position=c(0.06,0.06),
            title.size=0.8)

## Towns
tm_towns <- 
  tm_shape(roads) +
    tm_lines(lwd=0.25, col="orange4") +
  tm_shape(towns) +
    tm_dots(size=0.025, col=grey(0.5)) +
  tm_shape(ctry_PROJ, is.master=TRUE) +
	  tm_borders(col="black", lwd=0.5) +
  tm_layout(outer.margins=c(0,0,0,0.015),
            title="Towns",
            title.position=c(0.06,0.06),
            title.size=0.8)

## Rivers
tm_rivers <- 
  tm_shape(rivers) +
    tm_lines(lwd=0.25, col="blue") +
  tm_shape(ctry_PROJ, is.master=TRUE) +
	  tm_borders(col="black", lwd=0.5) +
  tm_scale_bar(c(0,250,500), text.size=0.5,
	             position=c(0.5,0.14), just=c("center", "top")) +
  tm_layout(outer.margins=c(0,0,0,0.015),
            title="Rivers",
            title.position=c(0.06,0.06),
            title.size=0.8)

## Protectead areas
tm_pa <- 
  tm_shape(roads) +
    tm_lines(lwd=0.25, col="orange4") +
	tm_shape(pa) +
    tm_fill(col="olivedrab3") +
  tm_shape(ctry_PROJ, is.master=TRUE) +
	  tm_borders(col="black", lwd=0.5) +
  tm_layout(outer.margins=c(0,0,0,0.015),
            title="Protected areas",
            title.position=c(0.06,0.06),
            title.size=0.8)

## Arrange plots with grid package
tmap_opt(1e8)
f <- here("Maps", dataset, ctry_iso, "var.png")
png(filename=f, width=textwidth, height=textwidth*2/3, units="cm", res=300)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,3)))
print(tm_elev, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(tm_slope, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(tm_roads, vp=viewport(layout.pos.row=1, layout.pos.col=3))
print(tm_towns, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(tm_rivers, vp=viewport(layout.pos.row=2, layout.pos.col=2))
print(tm_pa, vp=viewport(layout.pos.row=2, layout.pos.col=3))
dev.off()

## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "figures", "var.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

#======================
# Carbon map
#======================

## Crop carbon map to country extent
in_f <- file.path(dir_fdb, dataset, cont, ctry_iso, "data", "emissions", "AGB.tif")
out_f <- here("Maps", dataset, ctry_iso, "AGB.tif")
if (!file.exists(out_f)) {
  system(paste0('gdalwarp -cutline ', ctry_PROJ_shp,' -crop_to_cutline \\
                -overwrite \\
							  -co "COMPRESS=LZW" -co "PREDICTOR=2" ', in_f, ' ', out_f))
}

## Load raster
r_AGB <- read_stars(out_f)

## Aspect
bbox_r <- st_bbox(r_AGB)
asp_AGB <- (bbox_r$xmax - bbox_r$xmin)/(bbox_r$ymax - bbox_r$ymin)

## Prob raster for bbox
in_f_prob <- here("Maps", dataset, ctry_iso, "prob_500m.tif")
r_prob <- read_stars(in_f_prob)

## Palette
pal <- c("#e3d5c6", "#a59e5f", "#64701d", "#4d5c20", "#23390b", "#000900")

## Plot with tmap
tm_AGB <- 
	tm_shape(r_AGB, bbox=st_bbox(r_prob)) +
	  tm_raster(palette=pal, title="AGB (Mg/ha)",
		  				style="cont", legend.show=TRUE, legend.reverse=TRUE) +
  tm_shape(ctry_PROJ) +
	  tm_borders(col="black") +
	tm_shape(rect) +
	  tm_borders(col="black", lwd=2) +
	tm_scale_bar(c(0,250,500), text.size=1,
	             position=c(0.5,0), just=c("center", "bottom")) +
  tm_legend(position=c("left","bottom"), just=c("left","bottom")) +
  tm_layout(legend.text.size=1.2, legend.title.size=1.5)

## Zoom
tm_AGB_zoom <- 
	tm_shape(r_AGB, bbox=st_bbox(rect)) +
	  tm_raster(palette=pal, style="cont", legend.show=FALSE) +
	tm_scale_bar(c(0,5,10), text.size=1,
	             position=c(0.5,0), just=c("center", "bottom"))
## Viewport for inset
vp_AGB_zoom <- viewport(x=0.01, y=0.99, width=0.275, height=0.275, just=c("left", "top"))

## Save plot
tmap_opt(1e8, outer.margins=c(0,0,0,0))
f <- here("Maps", dataset, ctry_iso, "AGB.png")
tmap_save(tm_AGB, file=f, width=textwidth, height=textwidth/asp_AGB, units="cm", dpi=300,
					insets_tm=tm_AGB_zoom, insets_vp=vp_AGB_zoom)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "figures", "AGB.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

# EOF