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

## Some variables
##dataset <- "gfc2020_70" 
dataset <- "jrc2020"
dir.create(file.path("Maps", dataset, "maps"), recursive=TRUE)
dir_fdb <- "/home/forestatrisk-tropics"

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

## ===============================
## Study area borders by continent
## ===============================

## Continents including Brazil
continent <- c("Africa", "America", "Brazil", "Asia")
ncont <- length(continent)

## Loop on continent
for (i in 1:ncont) {
	cont <- continent[i]
	cmd <- paste0("find ",
								file.path(dir_fdb, dataset),
								" -regextype posix-egrep -regex '.*",
								cont,
								".*/data/ctry_PROJ.shp$' -exec ogr2ogr -update -nlt MULTIPOLYGON -append ",
								file.path("Maps", dataset, "maps", paste0("borders_", cont, ".gpkg")),
								" {} \\;")
	system(cmd)
}

## Corrections for Brazil on iso
Brazil_in <- st_read(file.path("Maps", dataset, "maps", "borders_Brazil.gpkg"))
Brazil_df <- Brazil_in %>%
	st_drop_geometry() %>%
	mutate(GID_0=substr(State_ID, 5, 6), NAME_0=as.character(State_name)) %>%
	mutate(NAME_0=ifelse(GID_0=="AM", "amazonas", NAME_0)) %>%
	mutate(NAME_0=ifelse(GID_0=="SC", "santa-catarina", NAME_0)) %>%
	select(GID_0, NAME_0)
Brazil <- st_sf(Brazil_df, geom=Brazil_in$geom)

## Combine with America and export
America_in <- st_read(file.path("Maps", dataset, "maps", "borders_America.gpkg"))
America <- rbind(America_in, Brazil)
st_write(America, file.path("Maps", dataset, "maps", "borders_America.gpkg"), delete_dsn=TRUE)

## Remove Brazil
file.remove(file.path("Maps", dataset, "maps", "borders_Brazil.gpkg"))

## Corrections for Asia on iso
Asia_in <- st_read(file.path("Maps", dataset, "maps", "borders_Asia.gpkg"))
Asia_df <- Asia_in %>%
	st_drop_geometry() %>%
	mutate(NAME_0=replace(as.character(NAME_0), c(4, 11, 13, 12, 16),
												c("Western Ghats", "Andaman and Nicobar",
													"North-East India", "Queensland", "Indonesia"))) %>%
	mutate(GID_0=replace(as.character(GID_0), c(4, 11, 13, 12, 16),
												c("WG", "AN", "NE", "QLD", "IDN")))
Asia <- st_sf(Asia_df, geom=Asia_in$geom)
st_write(Asia, file.path("Maps", dataset, "maps", "borders_Asia.gpkg"), delete_dsn=TRUE)

## Reset continents
continent <- c("Africa", "America", "Asia")
ncont <- length(continent)

## Simplify borders and reproject to 4326
## ! Simplification is done on a per geometry basis. Topology is not preserved.
## See here: https://github.com/r-spatial/sf/issues/381
## This could be corrected with grass with v.generalize or v.clean tool=snap.
## But precise enough for figures.
for (i in 1:ncont) {
	cont <- continent[i]
	f <- file.path("Maps", dataset, "maps", paste0("borders_", cont, "_simp.gpkg"))
	if (!file.exists(f)) {
		## Simplify and reproject
		in_f <- file.path("Maps", dataset, "maps", paste0("borders_", cont, ".gpkg"))
		out_f <- file.path("Maps", dataset, "maps", paste0("borders_", cont, "_simp.gpkg"))
		cmd <- paste0("ogr2ogr -overwrite -nlt MULTIPOLYGON -t_srs 'EPSG:4326' ", out_f, " ", in_f, " -simplify 1000")
		system(cmd)
	}
}

## =======================
## Figures for study areas
## =======================

## Continents
continent <- c("Africa", "America", "Asia")
ncont <- length(continent)

## Load GADM level0 data
gadm0 <- st_read(file.path("Maps", "GADM_data", "gadm36_level0.gpkg"))

## Countries and continent
data("World")

## Loop on continent
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
		select(iso_a3) %>%
		pull() %>%
		as.character()
	
	## Add some countries
	if (cont=="Africa") {iso3_cont <- c(as.character(iso3_cont), "REU", "MUS")}
	if (cont=="America") {
		w <- which(iso3_cont %in% c("GRL", "CAN"))
		iso3_cont <- as.character(iso3_cont)[-w]
	}
	
	## Extract countries
	gadm0_cont <- gadm0 %>%
		filter(GID_0 %in% iso3_cont)
	
	## Import study area borders
	f <- file.path("Maps", dataset, "maps", paste0("borders_", cont, "_simp.gpkg"))
	borders <- st_read(f) 
	if (dataset=="jrc2020" & cont=="Africa") {borders <- borders %>% filter(GID_0 != "STP")}
	
	# Compute areas
	borders$area <- st_area(borders)
	borders$iso_size <- pmax(0.4, scales::rescale(as.numeric(borders$area), to=c(0,4)))
	
	## Map with tmap
	tm <- 
		tm_shape(eq_sf) +
		  tm_lines(lty=1,lwd=0.5) +
		tm_shape(trop_sf) +
		  tm_lines(lty=2, lwd=0.5) +
		tm_shape(gadm0_cont) +
		  tm_fill(grey(0.9)) +
		tm_shape(borders, is.master=TRUE) +
		  tmap_options(max.categories=nrow(borders)+1) +
		  #tm_fill(col=grey(0.8)) +
		  tm_fill(col="GID_0", legend.show=FALSE, palette="cat") +
		  tm_borders(col="black") +
		  tm_text("GID_0", size="iso_size", auto.placement=FALSE, legend.size.show=FALSE)
	
	# ## Save maps as png
	# tmap_save(tm, file=file.path("Maps", dataset, "maps", paste0("study_areas_", cont, ".png")))
	## Save maps as svg (for modifications of label position with Inkscape)
	tmap_save(tm, file=file.path("Maps", dataset, "maps", paste0("study_areas_", cont, ".svg")))
		
}



# EOF