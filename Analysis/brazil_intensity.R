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
setwd("~/Code/forestatrisk-tropics/")

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

## Load results
fcc_tab2 <- read.table("Analysis/results/fcc_brazil.csv", sep=",", header=TRUE)

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
for.year <- paste("for",c(2019, seq(2030,2100,by=10)),sep="")
defor.period <- c("defor19.30", "defor30.40", "defor40.50",
                 "defor50.60", "defor60.70", "defor70.80",
                 "defor80.90", "defor90.100")
time.interval <- c(11, rep(10, 7))
n.period <- length(defor.period)

## Initialization
nfor <- fcc_tab2$for2019
andef <- fcc_tab2$andef
mat <- nfor
ctry_for <- rep(1, nctry)

## Loop on time periods for simulations
for (t in 1:n.period) {
    defor <- andef * time.interval[t]
    ## While a country has defor > nfor
    while(!all(nfor >= defor)) {
        excess <- 0
        for (i in 1:nctry) {
            ## We need to have nfor > defor for each ctry
            if (defor[i] > nfor[i]) { # if nfor[i]=0 and defor[i]=0: nothing
                ctry_for[i] <- 0
                excess <- excess + (defor[i] - nfor[i]) # Compute excess of deforestation
                defor[i] <- nfor[i] # Set defor to nfor for number in mat
            }
        }
        ## Number of countries with forest
        ncf <- sum(ctry_for == 1)
        ## We split the excess of deforestation among countries with forest
        ## This can make defor > nfor, thus implying the while loop
        defor[ctry_for == 1] <- defor[ctry_for == 1] + excess / ncf
    }
    nfor <- nfor - defor
    mat <- cbind(mat, defor, nfor)
    andef <- defor / time.interval[t]
}

df <- as.data.frame(mat)
names(df) <- as.vector(rbind(for.year, c(defor.period, 0)))[-c(2*length(for.year))]
df2 <- cbind(fcc_tab2[,2], df)

apply(df2[,2:8], 2, sum)

## Function diffusion
fcc_diffusion <- function(forest_t0, t0, annual_defor, t) {
    ## Variables
    nctry <- length(forest_t0)
    ctry_for <- as.numeric(forest_t0 > 0)
    ti <- t-t0  # time-interval
    defor <- annual_defor * ti
    nfor <- forest_t0
    ## While a country has defor > nfor
    while(!all(nfor >= defor)) {
        excess <- 0
        for (i in 1:nctry) {
            ## We need to have nfor > defor for each ctry
            if (defor[i] > nfor[i]) { # if nfor[i]=0 and defor[i]=0: nothing
                ctry_for[i] <- 0
                excess <- excess + (defor[i] - nfor[i]) # Compute excess of deforestation
                defor[i] <- nfor[i] # Set defor to nfor for number in mat
            }
        }
        ## Number of countries with forest
        ncf <- sum(ctry_for == 1)
        ## We split the excess of deforestation among countries with forest
        ## This can make defor > nfor, thus implying the while loop
        defor[ctry_for == 1] <- defor[ctry_for == 1] + excess / ncf
    }
    nfor <- nfor - defor
    return(list(forest_t0=forest_t0, forest_t=nfor, defor_t0_t=defor))
}

brazil_fcc_diffusion <- fcc_diffusion(forest_t0=fcc_tab2$for2019, t0=2019, annual_defor=fcc_tab2$andef, t=2050)
all((brazil_fcc_diffusion$forest_t0 - brazil_fcc_diffusion$defor_t0_t) == brazil_fcc_diffusion$forest_t)

#= Export results
data.di3 <- cbind(data.di2,mat)
head(data.di3)
write.table(data.di3,file="./results/defor.intensity.csv",sep=",",row.names=FALSE)

# End
