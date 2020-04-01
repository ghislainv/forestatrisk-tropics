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
iso3 <- ctry_df$iso3
iso3 <- as.character(iso3[-which(iso3 %in% c("BRA-DF","STP"))])
nctry <- length(iso3)

## =================
## Model performance
## =================

## Create table to store model performance results
ind_tab <- data.frame(matrix(NA, nrow=2*nctry, ncol=12))
names(ind_tab) <- c("cont", "iso3", "mod", "D", "AUC", "OA",
                    "EA", "FOM", "Sen", "Spe", "TSS", "K")

## Loop on countries
for (i in 1:nctry) {
    iso <- iso3[i]
    continent <- as.character(ctry_df$cont_run[ctry_df$iso3==iso])
    dir <- paste0("/share/nas2-amap/gvieilledent/", continent)
    ## Deviance
    f_name <- paste0(dir, "/", iso, "/output_jrc/model_deviance.csv")
    dev_df <- read.table(f_name, header=TRUE, sep=",", stringsAsFactors=FALSE)
    dev_glm <- dev_df$perc[dev_df$model=="glm"]
    dev_icar <- dev_df$perc[dev_df$model=="icar"]
    ## Performance index
    f_name <- paste0(dir, "/", iso, "/output_jrc/CV_glm.csv")
    CV_glm_df <- read.table(f_name, header=TRUE, sep=",", stringsAsFactors=FALSE)
    CV_glm <- round(CV_glm_df$mean*100)
    f_name <- paste0(dir, "/", iso, "/output_jrc/CV_icar.csv")
    CV_icar_df <- read.table(f_name, header=TRUE, sep=",", stringsAsFactors=FALSE)
    CV_icar <- round(CV_icar_df$mean*100)
    ## Sample size
    f_name <- paste0(dir, "/", iso, "/output_jrc/sample_size.csv")
    samp_size_df <- read.table(f_name, header=TRUE, sep=",")
    nfor <- samp_size_df
    ndefor <- samp_size_df
    ## Fill in the table
    ind_tab[2*(i-1)+1, 1:3] <- c(continent, iso, "icar")
    ind_tab[2*(i-1)+1, 4:12] <- c(dev_icar, CV_icar)
    ind_tab[2*(i-1)+2, 1:3] <- c(continent, iso, "glm")
    ind_tab[2*(i-1)+2, 4:12] <- c(dev_glm, CV_glm)
}

## Save results
write.table(ind_tab, file="Analysis/results/performance_index.csv", sep=",", row.names=FALSE)

## =================
## Sample size
## =================

## Create table to store results
samp_size_tab <- data.frame(matrix(NA, nrow=nctry, ncol=4))
names(samp_size_tab) <- c("cont", "iso3", "nfor", "ndef")

## Loop on countries
for (i in 1:nctry) {
    iso <- iso3[i]
    continent <- as.character(ctry_df$cont_run[ctry_df$iso3==iso])
    dir <- paste0("/share/nas2-amap/gvieilledent/", continent)
    ## Sample size
    f_name <- paste0(dir, "/", iso, "/output_jrc/sample_size.csv")
    samp_df <- read.table(f_name, header=TRUE, sep=",", stringsAsFactors=FALSE)
    nfor <- samp_df$n[samp_df$var=="nfor"]
    ndef <- samp_df$n[samp_df$var=="ndefor"]
    ## Fill in the table
    samp_size_tab[i, 1:2] <- c(continent, iso)
    samp_size_tab[i, 3:4] <- c(nfor, ndef)
}

## Save results
write.table(samp_size_tab, file="Analysis/results/samp_size.csv", sep=",", row.names=FALSE, quote=FALSE)

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
    mutate(pdef=round(100*(1-(1-(for2010-for2019)/for2010)^(1/9)), 1)) %>%
    mutate(for2035=pmax(0, for2019-16*andef), for2050=pmax(0, for2019-31*andef),
           for2055=pmax(0, for2019-36*andef), for2085=pmax(0, for2019-66*andef),
           for2100=pmax(0, for2019-81*andef)) %>%
    # Year during which forest should have disappeared
    mutate(yrdis=floor(2019 + for2019/andef))

## Save results
write.table(fcc_tab2, file="Analysis/results/forest_cover_change.csv", sep=",", row.names=FALSE, quote=FALSE)

