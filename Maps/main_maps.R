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

## Continents
continent <- c("Africa", "America", "Asia")
ncont <- length(continent)

## Load GADM level0 data
v <- here("Maps", "GADM_data", "gadm36_level0.gpkg")
gadm0 <- st_read(v)

## Simpify study area borders (keeping EPSG:3395)
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
## Africa
bb <- st_bbox(c(xmin=-19, xmax=52,
                ymin=ymin_trop, ymax=ymax_trop),
              crs=4326)
bbox_Afr <- bb %>% st_as_sfc() %>% st_transform(crs=3395) %>% st_bbox()
## America
bb <- st_bbox(c(xmin=-94, xmax=-33,
                ymin=ymin_trop, ymax=ymax_trop),
              crs=4326)
bbox_Ame <- bb %>% st_as_sfc() %>% st_transform(crs=3395) %>% st_bbox()
## Asia
bb <- st_bbox(c(xmin=68, xmax=170,
                ymin=ymin_trop, ymax=ymax_trop),
              crs=4326)
bbox_Asi <- bb %>% st_as_sfc() %>% st_transform(crs=3395) %>% st_bbox()
bbox_cont <- list(bbox_Afr, bbox_Ame, bbox_Asi)

## List to save maps
l_fcc2100 <- list()
l_prob <- list()
l_roads <- list()
l_pa <- list()
  
## Loop on continent
tmap_opt(npix=1e7)
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
	r_fcc2100 <- read_stars(here("Maps", dataset, cont, "fcc_2100_500m.tif"))
	r_prob <- read_stars(here("Maps", dataset, cont, "prob_500m.tif"))
	r_fcc123 <- read_stars(here("Maps", dataset, cont, "fcc123_500m.tif"))
	
	## Import roads
	roads <- st_read(here("Maps", dataset, cont, "roads_simp.gpkg"))
	
	## Intersect protected areas with ctry borders and import pa
	f <- here("Maps", dataset, cont, "pa_simp_crop.gpkg")
	if (!file.exists(f)) {
	  pa <- st_read(here("Maps", dataset, cont, "pa_simp.gpkg"))
	  pa_join <- st_join(pa, borders, join=st_intersects,
	                suffix=c(".x", ".y"), left=FALSE)
	  write_sf(pa_join, f)
	  rm(pa_join)
	}
	pa <- st_read(f)
	if (cont=="Asia") {
	  pa <- pa %>% 
	    dplyr::filter(!(pa_name=="Natural Park of the Coral Sea"))
	}

	## Maps fcc2050
	m_fcc2100 <- 
		tm_shape(eq_sf, bbox=bbox_cont[[i]]) +
		  tm_lines(lty=1,lwd=0.5) +
		tm_shape(trop_sf) +
		  tm_lines(lty=2, lwd=0.5) +
		tm_shape(gadm0_cont) +
		  tm_fill(grey(0.9)) +
	  tm_shape(r_fcc2100) +
	    tm_raster(style="cat", n=2, legend.show=FALSE,
	              palette=c(red, green)) +
	  tm_shape(borders) +
		  tm_borders(col=grey(0.5), lwd=0.5) +
	  tm_layout(inner.margins=c(0,0,0,0),
	            outer.margins=c(0,0,0,0))
	
	## Maps prob
	l.show <- ifelse(cont=="Asia", TRUE, FALSE)
	m_prob <- 
		tm_shape(eq_sf, bbox=bbox_cont[[i]]) +
		  tm_lines(lty=1,lwd=0.5) +
		tm_shape(trop_sf) +
		  tm_lines(lty=2, lwd=0.5) +
		tm_shape(gadm0_cont) +
		  tm_fill(grey(0.9)) +
	  tm_shape(r_prob) +
	    tm_raster(style="cont", title="", legend.reverse=TRUE, legend.show=l.show,
	              palette=c("#228b22", "#ffa500", "#e31a1c", "#000000"),
	              breaks=c(1, 39322, 54249, 65535), labels=c("0","","","1")) +
	  tm_shape(borders) +
		  tm_borders(col=grey(0.5), lwd=0.5) +
	  tm_layout(inner.margins=c(0,0,0,0),
	            outer.margins=c(0,0,0,0))
	
	## Maps roads
	m_roads <- 
		tm_shape(eq_sf, bbox=bbox_cont[[i]]) +
		  tm_lines(lty=1,lwd=0.5) +
		tm_shape(trop_sf) +
		  tm_lines(lty=2, lwd=0.5) +
		tm_shape(gadm0_cont) +
		  tm_fill(grey(0.9)) +
	  tm_shape(roads) +
	    tm_lines(lwd=0.25, col="orange4") +
	  tm_shape(borders) +
		  tm_borders(col=grey(0.5), lwd=0.5) +
	  tm_layout(inner.margins=c(0,0,0,0),
	            outer.margins=c(0,0,0,0))
	
	## Maps pa
	m_pa <- 
		tm_shape(eq_sf, bbox=bbox_cont[[i]]) +
		  tm_lines(lty=1,lwd=0.5) +
		tm_shape(trop_sf) +
		  tm_lines(lty=2, lwd=0.5) +
		tm_shape(gadm0_cont) +
		  tm_fill(grey(0.9)) +
	  tm_shape(pa) +
	    tm_fill(col="olivedrab3") +
	  tm_shape(borders) +
		  tm_borders(col=grey(0.5), lwd=0.5) +
	  tm_layout(inner.margins=c(0,0,0,0),
	            outer.margins=c(0,0,0,0))
	
	## Save in list
	l_fcc2100[[i]] <- m_fcc2100
	l_prob[[i]] <- m_prob
	l_roads[[i]] <- m_roads
	l_pa[[i]] <- m_pa
}

## Add year and scale bar to Asia
## fcc2100
fcc2100_Asia <- l_fcc2100[[3]] + 
  tm_scale_bar(breaks=c(0, 1000, 2000), text.size=0.6,
	             position=c(0.15, 0.01), just=c("left", "bottom")) +
  tm_layout(title="2100",
            title.position=c(0.02,0.08), title.size=1.2)
## prob (add legend layout)
prob_Asia <- l_prob[[3]] + 
  tm_scale_bar(breaks=c(0, 1000, 2000), text.size=0.6,
	             position=c(0.15, 0.01), just=c("left", "bottom")) +
  tm_legend(position=c(0.02, 0.02), just=c("left","bottom")) +
  tm_layout(legend.text.size=0.8)
## roads
roads_Asia <- l_roads[[3]] + 
  tm_scale_bar(breaks=c(0, 1000, 2000), text.size=0.6,
	             position=c(0.15, 0.01), just=c("left", "bottom")) +
  tm_layout(legend.text.size=0.8)
## pa
pa_Asia <- l_pa[[3]] + 
  tm_scale_bar(breaks=c(0, 1000, 2000), text.size=0.6,
	             position=c(0.15, 0.01), just=c("left", "bottom")) +
  tm_layout(legend.text.size=0.8)

## Plot (map) size in m
height_trop_m <- bbox_Ame$ymax-bbox_Ame$ymin
width_Ame_m <- bbox_Ame$xmax-bbox_Ame$xmin
width_Afr_m <- bbox_Afr$xmax-bbox_Afr$xmin
width_Asi_m <- bbox_Asi$xmax-bbox_Asi$xmin

## Figure width in m (map unit)
space_npc <- 0.015 # White space between plots (in "npc")
width_fig_m <- (width_Ame_m+width_Afr_m)/(1-space_npc)
space_m <- space_npc*width_fig_m
height_fig_m <- 2*height_trop_m + space_m
ratio <- height_fig_m/width_fig_m

## Plot (map) width and height in npc
height_trop_npc <- 0.5*(1-space_npc)
width_Ame_npc <- width_Ame_m/width_fig_m
width_Afr_npc <- width_Afr_m/width_fig_m
width_Asi_npc <- width_Asi_m/width_fig_m

## Viewports in npc
vp_Ame <- viewport(x=0, y=1, width=width_Ame_npc, height=height_trop_npc, just=c(0,1))
vp_Afr <- viewport(x=1, y=1, width=width_Afr_npc, height=height_trop_npc, just=c(1,1))
vp_Asi <- viewport(x=0.5, y=0, width=width_Asi_npc, height=height_trop_npc, just=c(0.5,0))

## Arrange plots with grid package

## fcc2100
f <- here("Maps", dataset, "fcc2100.png")
png(filename=f, width=textwidth, height=textwidth*ratio, units="cm", res=300)
grid.newpage()
print(l_fcc2100[[2]], vp=vp_Ame)
print(l_fcc2100[[1]], vp=vp_Afr)
print(fcc2100_Asia, vp=vp_Asi)
dev.off()
## Copy for manuscript
f_doc <- here("Manuscript", "Article", "figures", "fcc2100.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## prob
f <- here("Maps", dataset, "prob.png")
png(filename=f, width=textwidth, height=textwidth*ratio, units="cm", res=300)
grid.newpage()
print(l_prob[[2]], vp=vp_Ame)
print(l_prob[[1]], vp=vp_Afr)
print(prob_Asia, vp=vp_Asi)
dev.off()
## Copy for manuscript
f_doc <- here("Manuscript", "Article", "figures", "prob.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## roads
f <- here("Maps", dataset, "roads.png")
png(filename=f, width=textwidth, height=textwidth*ratio, units="cm", res=300)
grid.newpage()
print(l_roads[[2]], vp=vp_Ame)
print(l_roads[[1]], vp=vp_Afr)
print(roads_Asia, vp=vp_Asi)
dev.off()
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "figures", "roads.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## pa
f <- here("Maps", dataset, "pa.png")
png(filename=f, width=textwidth, height=textwidth*ratio, units="cm", res=300)
grid.newpage()
print(l_pa[[2]], vp=vp_Ame)
print(l_pa[[1]], vp=vp_Afr)
print(pa_Asia, vp=vp_Asi)
dev.off()
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "figures", "pa.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

# EOF