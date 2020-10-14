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
## Tropics
ymin_trop <- -2667917; ymax_trop <- 2667917
## Africa
xmin_Afr <- -1952868; xmax_Afr <- 5700000
bbox_Afr <- st_bbox(c(xmin=xmin_Afr, xmax=xmax_Afr,
                      ymin=ymin_trop, ymax=ymax_trop),
                    crs=st_crs(3395))
x_size_Afr <- xmax_Afr-xmin_Afr
## America
xmin_Ame <- -10500000; xmax_Ame <- -3600000
bbox_Ame <- st_bbox(c(xmin=xmin_Ame, xmax=xmin_Ame + x_size_Afr,
                      ymin=ymin_trop, ymax=ymax_trop),
                    crs=st_crs(3395))
## Asia
xmin_Asi <- 7400000; xmax_Asi <- 19100000
bbox_Asi <- st_bbox(c(xmin=xmin_Asi, xmax=xmax_Asi,
                      ymin=ymin_trop, ymax=ymax_trop),
                    crs=st_crs(3395))
x_size <- 19100000-7400000
bbox_cont <- list(bbox_Afr, bbox_Ame, bbox_Asi)

## List to save maps
l_fcc2050 <- list()
  
## Loop on continent
tmap_opt(npix=1e5)
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
		iso3_cont <- as.character(iso3_cont)[-w]
	}
	
	## Extract countries
	gadm0_cont <- gadm0 %>%
		filter(GID_0 %in% iso3_cont)
	
	## Import study area borders
	f <- here("Maps", dataset, cont, paste0("borders_", cont, "_simp.gpkg"))
	borders <- st_read(f) 
	if (dataset=="jrc2020" & cont=="Africa") {borders <- borders %>% filter(GID_0 != "STP")}
	
	## Import fcc2100 raster
	r_fcc2100 <- read_stars(here("Maps", dataset, cont, "fcc_2100_500m.tif"))
	
	## Maps
	m_fcc2050 <- 
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
	  tm_layout(outer.margins=c(0,0,0,0))
	
	## Save in list
	l_fcc2050[[i]] <- m_fcc2050
}

l_fcc2050[[1]]
l_fcc2050[[2]]
l_fcc2050[[3]]

## Arrange plots with grid package
f <- here("Maps", dataset, "fcc2050.png")
png(filename=f, width=textwidth, height=textwidth, units="cm", res=300)
grid.newpage()
pushViewport(viewport(layout=grid.layout(4,4)))
print(l_fcc2050[[2]], vp=viewport(layout.pos.row=1:2, layout.pos.col=1:2))
print(l_fcc2050[[1]], vp=viewport(layout.pos.row=1:2, layout.pos.col=3:4))
print(l_fcc2050[[3]], vp=viewport(layout.pos.row=3:4, layout.pos.col=1:4))
dev.off()

# EOF
