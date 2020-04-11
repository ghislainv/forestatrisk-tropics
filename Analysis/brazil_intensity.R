#!/usr/bin/R

## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
## web             :https://ghislainv.github.io
## license         :GPLv3
## ==============================================================================

## Libraries
## library(readr)
library(dplyr)
## library(ggplot2)

## Set working directory
setwd("/home/gvieilledent/Code/forestatrisk-tropics/")

## =================
## Countries info
## =================

## Load country info (encoding pb for ctry names)
ctry_df <- read.csv2("Analysis/ctry_run.csv", header=TRUE, sep=";", encoding="UTF-8")

## List of countries
iso3 <- ctry_df$iso3[ctry_df$cont_run=="Brazil"]
iso3 <- as.character(iso3[-which(iso3 %in% c("BRA-DF"))])
nctry <- length(iso3)

## ===================
## Forest cover change
## ===================

## Create table to store results
fcc_tab <- data.frame(matrix(NA, nrow=nctry, ncol=5))
names(fcc_tab) <- c("cont", "iso3", "for2000", "for2010", "for2019")

## Loop on countries
for (i in 1:nctry) {
    iso <- iso3[i]
    continent <- as.character(ctry_df$cont_run[ctry_df$iso3==iso])
    dir <- paste0("/share/nas2-amap/gvieilledent/", continent)
    ## Forest cover change
    f_name <- paste0(dir, "/", iso, "/output_jrc/forest_cover.txt")
    fcc_df <- read.table(f_name, header=FALSE, sep=",", stringsAsFactors=FALSE)
    area <- round(fcc_df[, 1])
    ## Fill in the table
    fcc_tab[i,1:2] <- c(continent, iso)
    fcc_tab[i,3:5] <- area
}

## Annual defor
fcc_tab2 <- fcc_tab %>%
    mutate(andef=round((for2010-for2019)/9)) %>%
    mutate(pdef=round(100*(1-(1-(for2010-for2019)/for2010)^(1/9)), 1))

## Save results
write.table(fcc_tab2, file="Analysis/results/fcc_brazil.csv", sep=",", row.names=FALSE)

## ===================
## Neighbors
## ===================

## ## Combine country borders
## cmd <- "find ~/nas/ -regextype posix-egrep -regex '.*/Brazil/.*/data/ctry_PROJ.shp$' \\
## -exec ogr2ogr -update -append Analysis/results/borders_brazil.shp {} \\;"
## system(cmd)

## ## Simplify
## cmd <- "ogr2ogr -overwrite Analysis/results/borders_simp_brazil.shp Analysis/results/borders_brazil.shp -simplify 1000"
## system(cmd)

## Load Brazil shapefile
library(sf)
library(geobr)
library(spdep)

## Get Brazilian state borders
brazil <- geobr::read_state(year=2018)
## Write to disk
sf::st_write(brazil,"Analysis/results/borders_brazil.shp")
## Transform as SpatialPolygonsDataFrame
brazil_sp <- as(brazil, "Spatial")

## Identify neighbors
row.names(brazil_sp) <- as.character(brazil_sp$abbrev_state)
nb <- spdep::poly2nb(brazil_sp, snap=1000*1/111111)  # Can accomodate snap (here ~1km)
mat <- spdep::nb2mat(nb, style="B")
colnames(mat) <- rownames(mat)

## Remove DF state
i <- which(colnames(mat)=="DF")
mat2 <- mat[-i, -i]

## Number of neighbors by state
n.neighbor <- apply(mat2, 2, sum)

## ========================
## Contagious deforestation
## ========================

#= Short algorithm to compute the number of deforested pixels by cell assuming
#= a constant deforestation at the ecoregion scale: defor.tot.mean = 8069130 pixels

# Year and deforestation periods
for.year <- paste("for",seq(2010,2100,by=10),sep="")
defor.period <- c("defor10.20","defor20.30","defor30.40","defor40.50",
                  "defor50.60","defor60.70","defor70.80","defor80.90",
                  "defor90.100")
n.period <- length(defor.period) # 9

# Initialization
nfor <- data.di2$nfor2010
defor <- data.di2$defor.mean
mat <- nfor

# Loop on years for simulations
for (t in 1:n.period) {
    ## Two conditions: deforestation must not exceed forest cover
    while(!all(nfor>=defor) & sum(defor)<=sum(nfor)) {
        for (i in 1:ncell.grid) {
            ## We need to have nfor > defor for each cell
            if (nfor[i] < defor[i]) { # if nfor[i]=0 and defor[i]=0: nothing
                ## Attributing excess of deforestation to neighbouring cells
                n.neigh <- n.neighbors[i]
                neigh <- neighbors.mat[neighbors.mat[,1]==i,2]
                defor.rest <- defor[i]-nfor[i]
                defor.neigh <- floor(defor.rest/n.neigh)
                defor.neigh.adjust <- defor.rest-(defor.neigh*n.neigh)
                defor[neigh] <- defor[neigh]+defor.neigh
                draw.neigh <- sample(neigh,size=1)
                defor[draw.neigh] <- defor[draw.neigh]+defor.neigh.adjust
                defor[i] <- nfor[i]
            }
        }
    }
    nfor <- nfor-defor
    mat <- cbind(mat,defor,nfor)
}
mat <- as.data.frame(mat)
names(mat) <- as.vector(rbind(for.year,c(defor.period,0)))[-c(20)]

#= Export results
data.di3 <- cbind(data.di2,mat)
head(data.di3)
write.table(data.di3,file="./results/defor.intensity.csv",sep=",",row.names=FALSE)

# End
