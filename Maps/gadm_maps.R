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
require(here)

## Some variables
##dataset <- "gfc2020_70" 
dataset <- "jrc2020"
dir.create(file.path("Maps", dataset, "maps"), recursive=TRUE)
continent <- c("Africa", "America", "Asia")
ncont <- length(continent)
dir_fdb <- "/home/forestatrisk-tropics"

## =======================
## Figures for study areas
## =======================

## Load GADM level0 data
gadm0 <- st_read(here("Maps", "GADM_data", "gadm36_level0.gpkg"))

## Countrycode
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
		w <- which(iso3_cont == "GRL")
		iso3_cont <- as.character(iso3_cont)[-w]
	}
	
	## Extract countries
	gadm0_cont <- gadm0 %>%
		filter(GID_0 %in% iso3_cont)
	
	# ## Export
	# st_write(gadm0_Afr, here("Maps", "GADM_data", "gadm36_level0_Afr.gpkg"), append=FALSE)
	# 
	# ## Simplify 
	# in_f <- here("Maps", "GADM_data", "gadm36_level0_Afr.gpkg")
	# out_f <- here("Maps", "GADM_data", "gadm36_level0_Afr_simp.gpkg")
	# cmd <- paste0("ogr2ogr -overwrite -f GPKG -nlt MULTIPOLYGON ", out_f, " ", in_f, " -simplify 1000")
	# system(cmd)
	# 
	# ## Import simplified layer
	# gadm0_Afr_simp <- st_read(out_f)
	# 

	## Object cont_regex
	cont_regex <- ifelse(cont=="America", "(America|Brazil)", cont)
	
	## Combine all borders
	cmd <- paste0("find ",
								file.path(dir_fdb, dataset),
								" -regextype posix-egrep -regex '.*",
								cont_regex,
								".*/data/ctry_PROJ.shp$' -exec ogr2ogr -update -nlt MULTIPOLYGON -append ",
								file.path("Maps", dataset, "maps", paste0("borders_", cont, ".gpkg")),
								" {} \\;")
	system(cmd)
	
	## Simplify and reproject
	in_f <- file.path("Maps", dataset, "maps", paste0("borders_", cont, ".gpkg"))
	out_f <- file.path("Maps", dataset, "maps", paste0("borders_simp_", cont, ".gpkg"))
	cmd <- paste0("ogr2ogr -overwrite -nlt MULTIPOLYGON -t_srs 'EPSG:4326' ", out_f, " ", in_f, " -simplify 1000")
	system(cmd)

	## Import study area borders
	f <- file.path("Maps", dataset, "maps", paste0("borders_simp_", cont, ".gpkg"))
	borders <- st_read(f) 
	if (dataset=="jrc2020" & cont=="Africa") {borders <- borders %>% filter(GID_0 != "STP")}
	
	## Bounding box
	bb <- st_bbox(borders)
	
	## Map with tmap
	tm <- tm_shape(gadm0_cont, bbox=bb) +
		tm_fill() +
		tm_shape(borders) +
		tm_fill() +
		tm_borders() +
		tm_text("GID_0", size=0.8)
	
}

## Save maps
tmap_save(tm, file=file.path("Maps", dataset, "maps", paste0("study_areas_", cont, ".png")))

# EOF