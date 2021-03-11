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
require(glue)

## Declare some variables
##dataset <- "gfc2020_70" 
dataset <- "jrc2020"

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

## Reproject equator and tropic in 3395
eq_sf <- st_transform(eq_sf, crs=3395)
trop_sf <- st_transform(trop_sf, crs=3395)

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

# =============================
# Maps at the continental scale
# =============================

## Continents
continent <- c("America", "Africa", "Asia")
ncont <- length(continent)

## Load GADM level0 data
v <- here("Maps", "GADM_data", "gadm36_level0.gpkg")
gadm0 <- st_read(v)

## Simplify study area borders (keeping EPSG:3395)
for (i in 1:ncont) {
	cont <- continent[i]
	in_f <- here("Maps", dataset, cont, paste0("borders_", cont, ".gpkg"))
	out_f <- here("Maps", dataset, cont, paste0("borders_", cont, "_simp.gpkg"))
	if (!file.exists(out_f)) {
		## Simplify and reproject
		cmd <- paste0("ogr2ogr -overwrite -nlt MULTIPOLYGON ", out_f, " ", in_f, " -simplify 1000")
		system(cmd)
	}
}

## Load countries and continent data
data("World")

## Bounding boxes
## Tropics in latlong
ymin_trop <- -23.43661; ymax_trop <- 23.43661
## America
bb <- st_bbox(c(xmin=-86, xmax=-34,
                ymin=ymin_trop, ymax=ymax_trop),
              crs=4326)
bbox_Ame <- bb %>% st_as_sfc() %>% st_transform(crs=3395) %>% st_bbox()
## Africa
bb <- st_bbox(c(xmin=-10, xmax=42,
                ymin=ymin_trop, ymax=ymax_trop),
              crs=4326)
bbox_Afr <- bb %>% st_as_sfc() %>% st_transform(crs=3395) %>% st_bbox()
## Asia
bb <- st_bbox(c(xmin=94, xmax=146,
                ymin=ymin_trop, ymax=ymax_trop),
              crs=4326)
bbox_Asi <- bb %>% st_as_sfc() %>% st_transform(crs=3395) %>% st_bbox()
bbox_cont <- list(bbox_Ame, bbox_Afr, bbox_Asi)

## Points for zoom (order: America, Africa, Asia)
xmin <- c(20840, 2347879, -1636620)
ymin <- c(2548500, -2289793, 1525793)
size_m <- 100000 
coords <- data.frame(X=xmin+size_m/2, Y=ymin+size_m/2)
crs_cont <- c("ESRI:102033", "ESRI:102022", "ESRI:102028")
zoom_points <- list()
for (i in 1:3) {
  sp_proj <- st_as_sf(coords[i, ], coords = c("X", "Y"), crs=crs_cont[i])
  zoom_points[[i]] <- st_transform(sp_proj, crs=3395)
}


## List to save maps
l_prob <- list()

## Loop on continent
tmap_opt(npix=1e+07)
for (i in 1:ncont) {
	
	## Continent
	cont <- continent[i]
	if (cont=="Africa") {cont_world <- c("Africa")}
	if (cont=="America") {cont_world <- c("North America", "South America")}
	if (cont=="Asia") {cont_world <- c("Asia", "Oceania")}

	## Iso code for African coutry in World database
	## NB: World database doesn't include small countries
	iso3_cont <- World %>%
		st_drop_geometry() %>%
		filter(continent %in% cont_world) %>%
		dplyr::select(iso_a3) %>%
		pull() %>%
		as.character()
	
	## Add or remove some countries
	if (cont=="Africa") {iso3_cont <- c(as.character(iso3_cont), "REU", "MUS")}
	if (cont=="America") {
		w <- which(iso3_cont %in% c("GRL", "CAN"))
		iso3_cont <- c(as.character(iso3_cont)[-w], "GUF")
	}
	
	## Extract countries
	gadm0_cont <- gadm0 %>%
		filter(GID_0 %in% iso3_cont)
	
	## Import study area borders
	f <- here("Maps", dataset, cont, paste0("borders_", cont, "_simp.gpkg"))
	borders <- st_read(f) 
	if (dataset=="jrc2020" & cont=="Africa") {borders <- borders %>% filter(GID_0 != "STP")}
	
	## Import rasters
	r_prob <- read_stars(here("Maps", dataset, cont, "prob_500m.tif"))

	## Maps prob
	l.show <- ifelse(cont=="Africa", TRUE, FALSE)
	m_prob <- 
		tm_shape(eq_sf) +
		  tm_lines(lty=1,lwd=0.5) +
		tm_shape(trop_sf) +
		  tm_lines(lty=2, lwd=0.5) +
		tm_shape(gadm0_cont) +
		  tm_fill(grey(0.9)) +
	  tm_shape(r_prob, bbox=bbox_cont[[i]], is.master=TRUE) +
	    tm_raster(style="cont", title="", legend.reverse=TRUE, legend.show=l.show,
	              palette=c("#228b22", "#ffa500", "#e31a1c", "#000000"),
	              breaks=c(1, 39322, 54249, 65535), labels=c("0","","","1")) +
	  tm_shape(borders) +
		  tm_borders(col=grey(0.5), lwd=0.5) +
	  tm_shape(zoom_points[[i]]) +
		  tm_squares(col="black", size=0.1) +
	  tm_layout(inner.margins=c(0,0,0,0),
	            outer.margins=c(0,0,0,0))
	
	## Save in list
	l_prob[[i]] <- m_prob
}

## Add scale bar to Asia
prob_Asia <- l_prob[[3]] + 
  tm_scale_bar(breaks=c(0, 1000, 2000), text.size=0.5,
	             position=c("left", "bottom"), just=c("left", "bottom"))

## Add legend to Africa
prob_Africa <- l_prob[[2]] + 
  tm_legend(position=c("left", "bottom"), just=c("left", "bottom")) +
  tm_layout(legend.text.size=0.6)

# ===============================================
# Zooming on the deforestation probability map
# ===============================================

# ----------------------
# Function f_zoom_prob()
# ----------------------

f_zoom_prob <- function(ctr_iso, cont, continent, xmin, ymin,
                        size_m, npix, col_roads="#043353") {

  borders <- st_read(here("Data", dataset, continent, ctry_iso, "data", "ctry_PROJ.shp"))
  bb_proj <- st_bbox(c(xmin=xmin, xmax=xmin+size_m,
                       ymin=ymin, ymax=ymin+size_m),
                     crs=st_crs(borders))
  bb_ll <- bb_proj %>% st_as_sfc() %>% st_transform(crs=4326) %>% st_bbox()
  
  # Prob extract
  f_out <- here("Maps", dataset, ctry_iso, "prob_zoom.tif")
  if (!file.exists(f_out)) {
    extent <- glue("{bb_ll$xmin} {bb_ll$ymax} {bb_ll$xmax} {bb_ll$ymin}")
    proj <- "EPSG:4326"
    f_in <- glue("/vsicurl/https://forestatrisk.cirad.fr/tropics/tif/prob_{cont}.tif")
    cmd <- glue("gdal_translate -projwin {extent} -projwin_srs {proj} {f_in} {f_out}")
    system(cmd)
  }
  
  # Data
  roads <- st_read(here("Data", dataset, continent, ctry_iso, "data", "roads_PROJ.shp"))
  pa <- st_read(here("Data", dataset, continent, ctry_iso, "data", "pa_PROJ.shp"))
  r_prob <- read_stars(here("Maps", dataset, ctry_iso, "prob_zoom.tif"))
  
  # Crop to extent
  roads <- roads %>% st_crop(bb_proj)
  
  ## tmap_options
  tmap_opt(npix=npix)
  
  ## Map
  m_zoom_prob <- 
    tm_shape(borders) +
      tm_fill(col=grey(0.9)) +
    tm_shape(r_prob, bbox=bb_proj, is.master=TRUE) +
      tm_raster(style="cont", alpha=0.8, title="", legend.reverse=TRUE,
                legend.show=FALSE,
  	            palette=c("#228b22", "#ffa500", "#e31a1c", "#000000"),
  	            breaks=c(1, 39322, 54249, 65535), labels=c("0","","","1")) +
    tm_shape(roads) +
  	  tm_lines(col=col_roads, lwd=0.7) +
    tm_shape(pa) +
      tm_borders(col="black", lwd=1.2) +
    tm_shape(pa) +
      tm_fill(col=grey(0.5), alpha=0.4) +
    tm_layout(inner.margins=c(0,0,0,0),
              outer.margins=c(0,0,0,0))
  
  ## Return the map
  return(m_zoom_prob)
}

# Set npix
npix <- 1e+07

# ---------------------
# America (Mato Grosso)
# ---------------------
ctry_iso <- "BRA-MT"
cont <- "AME"
continent <- "Brazil"
xmin <- 20840; ymin <- 2548500
zoom_prob_Ame <- f_zoom_prob(ctr_iso=ctry_iso, cont=cont,
                             continent=continent, xmin=xmin, ymin=ymin,
                             size_m=100000, npix=npix, col_roads="black")
# -------------------
# Africa (Madagascar)
# -------------------
ctry_iso <- "MDG"
cont <- "AFR"
continent <- "Africa"
xmin <- 2347879; ymin <- -2289793
zoom_prob_Afr <- f_zoom_prob(ctr_iso=ctry_iso, cont=cont,
                             continent=continent, xmin=xmin, ymin=ymin,
                             size_m=100000, npix=npix, col_roads="black")
# -------------------
# Asia (Indonesia)
# -------------------
ctry_iso <- "IDN"
cont <- "ASI"
continent <- "Asia"
xmin <- -1636620; ymin <- 1525793
zoom_prob_Asi <- f_zoom_prob(ctr_iso=ctry_iso, cont=cont,
                             continent=continent, xmin=xmin, ymin=ymin,
                             size_m=100000, npix=npix, col_roads="black")
## Add legend and scale bar to Asia
zoom_prob_Asi <- zoom_prob_Asi + 
  tm_scale_bar(breaks=c(0, 10, 20), text.size=0.5,
	             position=c("left", "bottom"), just=c("left", "bottom"))

# ==============================
# Arrange plot with grid package
# ==============================

## Plot (map) size in m
height_cont_m <- bbox_Ame$ymax-bbox_Ame$ymin
width_cont_m <- bbox_Ame$xmax-bbox_Ame$xmin # Same between continents (52°)

## Figure width in m (map unit)
space_npc <- 0.0075 # White space between plots (in "npc")
width_fig_m <- (3*width_cont_m)/(1-2*space_npc)
space_m <- space_npc*width_fig_m
height_fig_m <- height_cont_m + space_m + width_cont_m
ratio <- height_fig_m/width_fig_m

## Plot (map) width and height in npc
# Continents
height_cont_npc <- height_cont_m/height_fig_m
width_cont_npc <- width_cont_m/width_fig_m
# Zooms
height_zoom_npc <- 1-space_npc-height_cont_npc
width_zoom_npc <- 1-2*space_npc-2*width_cont_npc

## Viewports in npc
# Continents
vp_Ame <- viewport(x=0, y=1, width=width_cont_npc, height=height_cont_npc, just=c(0,1))
vp_Afr <- viewport(x=0.5, y=1, width=width_cont_npc, height=height_cont_npc, just=c(0.5,1))
vp_Asi <- viewport(x=1, y=1, width=width_cont_npc, height=height_cont_npc, just=c(1,1))
# Zooms
vp_Ame_z <- viewport(x=0, y=0, width=width_zoom_npc, height=height_zoom_npc, just=c(0,0))
vp_Afr_z <- viewport(x=0.5, y=0, width=width_zoom_npc, height=height_zoom_npc, just=c(0.5,0))
vp_Asi_z <- viewport(x=1, y=0, width=width_zoom_npc, height=height_zoom_npc, just=c(1,0))

## Arrange plots with grid package
## prob
f <- here("Maps", dataset, "prob_zoom.png")
png(filename=f, width=textwidth, height=textwidth*ratio, units="cm", res=300)
grid.newpage()
# Continents
print(l_prob[[1]], vp=vp_Ame)
print(prob_Africa, vp=vp_Afr)
print(prob_Asia, vp=vp_Asi)
# Zooms
print(zoom_prob_Ame, vp=vp_Ame_z)
print(zoom_prob_Afr, vp=vp_Afr_z)
print(zoom_prob_Asi, vp=vp_Asi_z)
dev.off()
## Copy for manuscript
f_doc <- here("Manuscript", "Article", "figures", "prob_zoom.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

# EOF