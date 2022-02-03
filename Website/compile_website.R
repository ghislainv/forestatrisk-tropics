#!/usr/bin/R

## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
## web             :https://ghislainv.github.io
## license         :GPLv3
## ==============================================================================

library(glue)
library(rmarkdown)
library(here)

## Compile site
render_site(here::here("Website"))

## Copy website directory
## system("cp -r Website/_site/* /home/www/forestatrisk/")
system("scp -r Website/_site/* fdb:/home/www/forestatrisk/")

## ===============================
## Make symbolic link to COG files
## ===============================

## ## Continents and abbreviations
## cont <- c("America", "Africa", "Asia")
## cab <- c("AME", "AFR", "ASI")
## year <- c("2030", "2040", "2050", "2055", "2060", "2070", "2080", "2085", "2090", "2100")

## ## loop on continent and year
## for (i in 1:length(cont)) {
## 	# fcc123
## 	cmd <- glue("ln -s /home/forestatrisk-tropics/jrc2020/Maps/{cont}/fcc123.tif \\
## 		    /home/www/forestatrisk/tropics/tif/fcc123_{cab}.tif")
## 	system(cmd)
## 	# prob
## 	cmd <- glue("ln -s /home/forestatrisk-tropics/jrc2020/Maps/{cont}/prob.tif \\
## 		    /home/www/forestatrisk/tropics/tif/prob_{cab}.tif")
## 	system(cmd)
## 	# future fcc
## 	for (j in 1:length(year)) {
## 		cont <- continent[i]
## 		cab <- continent_ab[i]
## 		yr <- year[j]
## 		cmd <- glue("ln -s /home/forestatrisk-tropics/jrc2020/Maps/{cont}/fcc_{yr}.tif \\
## 			    /home/www/forestatrisk/tropics/tif/fcc_{yr}_{cab}.tif")
## 		system(cmd)
## 	}
## }

## ============================================
## Make symbolic link to COG files in EPSG:3857
## ============================================

## Continents and abbreviations
cont <- c("America", "Africa", "Asia")
cab <- c("AME", "AFR", "ASI")
year <- c("2050", "2100")

## loop on continent and year
for (i in 1:length(cont)) {
	# fcc123
	cmd <- glue("ln -s /home/forestatrisk-tropics/jrc2020/Maps/{cont}/fcc123_epsg3857.tif \\
		    /home/www/forestatrisk/tropics/tif/fcc123_{cab}_epsg3857.tif")
	system(cmd)
	# prob
	cmd <- glue("ln -s /home/forestatrisk-tropics/jrc2020/Maps/{cont}/prob_epsg3857.tif \\
		    /home/www/forestatrisk/tropics/tif/prob_{cab}_epsg3857.tif")
	system(cmd)
	# future fcc
	for (j in 1:length(year)) {
		cont <- continent[i]
		cab <- continent_ab[i]
		yr <- year[j]
		cmd <- glue("ln -s /home/forestatrisk-tropics/jrc2020/Maps/{cont}/fcc_{yr}_epsg3857.tif \\
			    /home/www/forestatrisk/tropics/tif/fcc_{yr}_{cab}_epsg3857.tif")
		system(cmd)
	}
}

## ==============================
## Make symbolic link to raw data
## ==============================

## Note for sshfs to mbb:
## sshfs -o uid=1000 -o gid=1000 mbb:/home/gvieilledent /home/ghislain/mbb/gvieilledent
## sshfs -o uid=1000 -o gid=1000 mbb:/share/nas2-amap/gvieilledent /home/ghislain/mbb/nas

## Note: we have finally used Google Drive as data repository

## system("ln -s /home/ghislain/mbb/nas/jrc2020/Africa \\
##	  /home/www/forestatrisk/tropics/raw_data/Africa")

## EOF
