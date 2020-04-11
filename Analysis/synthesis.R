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

## ==========================
## Dataset with all countries
## ==========================

## Create table to store model performance results
sample_allctry_tab <- data.frame()

## Loop on countries
for (i in 1:nctry) {
    iso <- iso3[i]
    continent <- as.character(ctry_df$cont_run[ctry_df$iso3==iso])
    dir <- paste0("/share/nas2-amap/gvieilledent/", continent)
    ## Sample size
    f_name <- paste0(dir, "/", iso, "/output_jrc/sample.txt")
    sample_df <- read.table(f_name, header=TRUE, sep=",", stringsAsFactors=FALSE)
    sample_df$iso3 <- iso
    sample_df$continent <- continent
    ## Fill in the table
    sample_allctry_tab <- rbind(sample_allctry_tab,sample_df)
}

## Save results
write.table(sample_allctry_tab, file="Analysis/results/sample_allctry.csv", sep=",", row.names=FALSE)

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
write.table(samp_size_tab, file="Analysis/results/samp_size.csv", sep=",", row.names=FALSE)

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
write.table(fcc_tab2, file="Analysis/results/forest_cover_change.csv", sep=",", row.names=FALSE)

## ===================
## Model parameters
## ===================

## Create table to store results
par_tab <- data.frame(matrix(NA, nrow=nctry, ncol=12))
names(par_tab) <- c("cont", "iso3", "int", "pa", "alt",
                    "slope", "ddefor", "dedge", "driver",
                    "droad", "dtown", "Vrho")
long_var_names <- c("Intercept", "C(pa)[T.1.0]", "scale(altitude)",
                    "scale(slope)", "scale(dist_defor)",
                    "scale(dist_edge)", "scale(dist_river)",
                    "scale(dist_road)", "scale(dist_town)", "Vrho")

## Loop on countries
for (i in 1:nctry) {
    iso <- iso3[i]
    continent <- as.character(ctry_df$cont_run[ctry_df$iso3==iso])
    dir <- paste0("/share/nas2-amap/gvieilledent/", continent)
    ## Parameter estimates
    f_name <- paste0(dir, "/", iso, "/output_jrc/summary_hSDM.txt")
    par <- read.table(f_name, skip=4)
    names(par) <- c("Var", "Mean", "Sd", "CI_low", "CI_high")
    ## Fill in the table
    par_tab[i, 1:2] <- c(continent, iso)
    for (r in 1:(nrow(par)-1)) {
        j <- which(long_var_names==par$Var[r])
        par_tab[i, j+2] <- par$Mean[r]
    }
}

## Save results
write.table(par_tab, file="Analysis/results/parameter_estimates.csv", sep=",", row.names=FALSE)

## ===================
## PA effect
## ===================

## Create table to store results
parea_tab <- data.frame(matrix(NA, nrow=nctry, ncol=6))
names(parea_tab) <- c("cont", "iso3", "Mean", "Sd", "CI_low", "CI_high")

## Loop on countries
for (i in 1:nctry) {
    iso <- iso3[i]
    continent <- as.character(ctry_df$cont_run[ctry_df$iso3==iso])
    dir <- paste0("/share/nas2-amap/gvieilledent/", continent)
    ## Parameter estimates
    f_name <- paste0(dir, "/", iso, "/output_jrc/summary_hSDM.txt")
    par <- read.table(f_name, skip=4)
    names(par) <- c("Var", "Mean", "Sd", "CI_low", "CI_high")
    ## Fill in the table
    parea_tab[i, 1:2] <- c(continent, iso)
    j <- which(par$Var=="C(pa)[T.1.0]")
    if (length(j)>0) {
        parea_tab[i, 3:6] <- par[j, 2:5]
    }
}

## Save results
write.table(parea_tab, file="Analysis/results/parea_estimates.csv", sep=",", row.names=FALSE)

## ===================
## Road effect
## ===================

## Create table to store results
road_tab <- data.frame(matrix(NA, nrow=nctry, ncol=6))
names(road_tab) <- c("cont", "iso3", "Mean", "Sd", "CI_low", "CI_high")

## Loop on countries
for (i in 1:nctry) {
    iso <- iso3[i]
    continent <- as.character(ctry_df$cont_run[ctry_df$iso3==iso])
    dir <- paste0("/share/nas2-amap/gvieilledent/", continent)
    ## Parameter estimates
    f_name <- paste0(dir, "/", iso, "/output_jrc/summary_hSDM.txt")
    par <- read.table(f_name, skip=4)
    names(par) <- c("Var", "Mean", "Sd", "CI_low", "CI_high")
    ## Fill in the table
    road_tab[i, 1:2] <- c(continent, iso)
    j <- which(par$Var=="scale(dist_road)")
    if (length(j)>0) {
        road_tab[i, 3:6] <- par[j, 2:5]
    }
}

## Save results
write.table(road_tab, file="Analysis/results/road_estimates.csv", sep=",", row.names=FALSE)

## ===================
## Emissions
## ===================

## ===================
## Main results
## ===================

## Forest cover and forest loss
r1 <- "Global moist tropical forest (ha)"
r2 <- "Moist tropical forest per continent (ha)"
r3 <- "Global annual forest lost (ha)"
r4 <- "Annual forest lost per continent (ha)"
r5 <- "Total forest loss at each date in the future with percentage"
r6 <- "Total forest loss per continent at each date in the future with percentage"

## Sample
r.sa.1 <- "Number of sample points per country"
r.sa.2 <- "Total number of sample points (nfor, ndefor)"

## Model performance
r6.1 <- "Percentage of deviance explained (with SD) at the global scale"
r6.2 <- "Percentage of deviance explained (with SD) by continent"
r6.3 <- "Performance index with (with SD) at the global scale"
r6.4 <- "Performance index with (with SD) by continent"

## Variables
r7 <- "Percentage of forest inside protected areas at each date"
r8 <- "Number of countries (and %) with significant PA effect"
r9 <- "Number of countries (and %) with significant road effect"
r10 <- "Decrease in the deforestation risk with distance to road. Figure and estimates at 100, 200, 500 and 1km."
r10 <- "EDGE EFFECT: Decrease in the deforestation risk with distance to edge. Figure and estimates at 100, 200, 500 and 1km."
r10.1 <- "Variable importance"

## CO2 emissions
r11 <- "Global CO2 emissions at each date in the future"
r12 <- "CO2 emissions per continent at each date in the future"

## Fragmentation?

## End
