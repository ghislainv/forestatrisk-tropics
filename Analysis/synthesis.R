#!/usr/bin/R

## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
## web             :https://ghislainv.github.io
## license         :GPLv3
## ==============================================================================

## Libraries
## require(readr)
require(dplyr)
require(here)
## require(ggplot2)

## Set working directory
setwd(here())

## Dataset
##dataset <- "gfc2020_70" 
dataset <- "jrc2020"
dir.create(file.path("Analysis", dataset, "results"), recursive=TRUE)

## Result directory
dir_fdb <- "/home/forestatrisk-tropics"

## =================
## Countries info
## =================

## Load country info (encoding pb for ctry names)
ctry_df <- read.csv2("Analysis/ctry_run.csv", header=TRUE, sep=";", encoding="UTF-8")

## List of countries
iso3 <- ctry_df$iso3
if (dataset=="gfc2020_70") {
    iso3 <- as.character(iso3[-which(iso3 %in% c("STP", "SEN", "GMB"))])
} else if (dataset=="jrc2020") {
    iso3 <- as.character(iso3[-which(iso3 %in% c("STP"))])
}
nctry <- length(iso3)

## ## ==========================
## ## Dataset with all countries
## ## ==========================

## ## Create table to store model performance results
## sample_allctry_tab <- data.frame()

## ## Loop on countries
## for (i in 1:nctry) {
##     iso <- iso3[i]
##     continent <- as.character(ctry_df$cont_run[ctry_df$iso3==iso])
##     dir <- file.path(dir_fdb, dataset, continent)
##     ## Sample size
##     f_name <- file.path(dir, iso, "/output/sample.txt")
##     sample_df <- read.table(f_name, header=TRUE, sep=",", stringsAsFactors=FALSE)
##     sample_df$iso3 <- iso
##     sample_df$continent <- continent
##     ## Fill in the table
##     sample_allctry_tab <- rbind(sample_allctry_tab,sample_df)
## }

## ## !! NA in datasets (65535), see eg VIR and ATG 

## ## Save results
## write.table(sample_allctry_tab, file=file.path("Analysis", dataset, "results/sample_allctry.csv"), sep=",", row.names=FALSE)

## =================
## Model performance
## =================

## Create table to store model performance results
ind_tab <- data.frame(matrix(NA, nrow=3*nctry, ncol=14))
names(ind_tab) <- c("cont", "iso3", "mod", "D", "AUC", "OA",
                    "EA", "FOM", "Sen", "Spe", "TSS", "K",
                    "nfor", "ndefor")

## Loop on countries
for (i in 1:nctry) {
    iso <- iso3[i]
    continent <- as.character(ctry_df$cont_run[ctry_df$iso3==iso])
    dir <- file.path(dir_fdb, dataset, continent)
    ## Deviance
    f_name <- file.path(dir, iso, "output/model_deviance.csv")
    dev_df <- read.table(f_name, header=TRUE, sep=",", stringsAsFactors=FALSE)
    dev_glm <- dev_df$perc[dev_df$model=="glm"]
    dev_icar <- dev_df$perc[dev_df$model=="icar"]
    dev_rf <- dev_df$perc[dev_df$model=="rf"]
    ## Performance index
    ## glm
    f_name <- file.path(dir, iso, "output/CV_glm.csv")
    CV_glm_df <- read.table(f_name, header=TRUE, sep=",", stringsAsFactors=FALSE)
    CV_glm <- round(CV_glm_df$mean*100)
    ## icar
    f_name <- file.path(dir, iso, "output/CV_icar.csv")
    CV_icar_df <- read.table(f_name, header=TRUE, sep=",", stringsAsFactors=FALSE)
    CV_icar <- round(CV_icar_df$mean*100)
    ## RF
    f_name <- file.path(dir, iso, "output/CV_rf.csv")
    CV_rf_df <- read.table(f_name, header=TRUE, sep=",", stringsAsFactors=FALSE)
    CV_rf <- round(CV_rf_df$mean*100)
    ## Sample size
    f_name <- file.path(dir, iso, "output/sample_size.csv")
    samp_size_df <- read.table(f_name, header=TRUE, sep=",")
    nfor <- samp_size_df$n[samp_size_df$var=="nfor"]
    ndefor <- samp_size_df$n[samp_size_df$var=="ndefor"]
    ## Fill in the table
    ind_tab[3*(i-1)+1, 1:3] <- c(continent, iso, "icar")
    ind_tab[3*(i-1)+1, 4:12] <- c(dev_icar, CV_icar)
    ind_tab[3*(i-1)+1, 13:14] <- c(nfor, ndefor)
    ind_tab[3*(i-1)+2, 1:3] <- c(continent, iso, "glm")
    ind_tab[3*(i-1)+2, 4:12] <- c(dev_glm, CV_glm)
    ind_tab[3*(i-1)+2, 13:14] <- c(nfor, ndefor)
    ind_tab[3*(i-1)+3, 1:3] <- c(continent, iso, "rf")
    ind_tab[3*(i-1)+3, 4:12] <- c(dev_rf, CV_rf)
    ind_tab[3*(i-1)+3, 13:14] <- c(nfor, ndefor)
}

## Performance per continent
## 1. Weighted percentage with forest size in 2010
fcc_tab <- read.table(file.path("Analysis", dataset, "results/forest_cover_change.csv"), header=TRUE, sep=",")
weights <- rep(fcc_tab$for2010, each=3)
## 2. Summarize per mod
perf_mod <- ind_tab %>% 
    mutate(w=weights) %>%
    group_by(mod) %>%
    summarize(D=weighted.mean(D, w), AUC=weighted.mean(AUC,w),
              OA=weighted.mean(OA,w), FOM=weighted.mean(FOM,w),
              TSS=weighted.mean(TSS,w))
## 3. Summarize per cont and mod
perf_cont_mod <- ind_tab %>% 
    mutate(w=weights) %>%
    mutate(cont=ifelse(cont=="Brazil", "America", cont)) %>%
    group_by(cont, mod) %>%
    summarize(D=weighted.mean(D, w), AUC=weighted.mean(AUC,w),
              OA=weighted.mean(OA,w), FOM=weighted.mean(FOM,w),
              TSS=weighted.mean(TSS,w)) %>%
    mutate(cont_id=ifelse(cont=="America", 1, ifelse(cont=="Africa", 2, 3))) %>%
    arrange(cont_id) %>%
    select(-cont_id)

## Save results
write.table(ind_tab, file=file.path("Analysis", dataset, "results/performance_index.csv"), sep=",", row.names=FALSE)
write.table(perf_mod, file=file.path("Analysis", dataset, "results/perf_mod.csv"), sep=",", row.names=FALSE)
write.table(perf_cont_mod, file=file.path("Analysis", dataset, "results/perf_cont_mod.csv"), sep=",", row.names=FALSE)

## =================
## Sample size
## =================

## Create table to store results
samp_size_tab <- data.frame(matrix(NA, nrow=nctry, ncol=5))
names(samp_size_tab) <- c("cont", "ctry", "iso3", "nfor", "ndef")

## Loop on countries
for (i in 1:nctry) {
    iso <- iso3[i]
    continent <- as.character(ctry_df$cont_run[ctry_df$iso3==iso])
    country <- as.character(ctry_df$ctry_run[ctry_df$iso3==iso])
    dir <- file.path(dir_fdb, dataset, continent)
    ## Sample size
    f_name <- file.path(dir, iso, "output/sample_size.csv")
    samp_size_df <- read.table(f_name, header=TRUE, sep=",", stringsAsFactors=FALSE)
    nfor <- samp_size_df$n[samp_size_df$var=="nfor"]
    ndef <- samp_size_df$n[samp_size_df$var=="ndefor"]
    ## Fill in the table
    samp_size_tab[i, 1:3] <- c(continent, country, iso)
    samp_size_tab[i, 4:5] <- c(nfor, ndef)
}

## Equivalence in ha
samp_size_tab$nforHa <- round(samp_size_tab$nfor*30*30/10000)
samp_size_tab$ndefHa <- round(samp_size_tab$ndef*30*30/10000)

## Updating continent, country, region and code
samp_size_tab <- samp_size_tab %>%
    # Continent
    mutate(cont2=ifelse(cont=="Brazil", "America", cont)) %>%
    # Country
    mutate(ctry2=ifelse(cont=="Brazil", "Brazil", ctry)) %>%
    mutate(ctry2=ifelse(iso3=="AUS-QLD", "Australia", ctry2)) %>%
    mutate(ctry2=ifelse(iso3=="IND-AND", "India", ctry2)) %>%
    mutate(ctry2=ifelse(iso3=="IND-WEST", "India", ctry2)) %>%
    mutate(ctry2=ifelse(iso3=="IND-EAST", "India", ctry2)) %>%
    # Region
    mutate(region=ifelse(cont=="Brazil", ctry, "")) %>%
    mutate(region=ifelse(iso3=="AUS-QLD", "Queensland", region)) %>%
    mutate(region=ifelse(iso3=="IND-AND", "Andaman and N.", region)) %>%
    mutate(region=ifelse(iso3=="IND-WEST", "West. Ghats", region)) %>%
    mutate(region=ifelse(iso3=="IND-EAST", "North-East", region)) %>%
    # Code
    mutate(code=ifelse(cont=="Brazil", substr(iso3,5,6), iso3)) %>%
    mutate(code=ifelse(iso3=="AUS-QLD", "QLD", code)) %>%
    mutate(code=ifelse(iso3=="IND-AND", "AN", code)) %>%
    mutate(code=ifelse(iso3=="IND-WEST", "WG", code)) %>%
    mutate(code=ifelse(iso3=="IND-EAST", "NE", code)) %>%
    # Id
    mutate(id=ifelse(cont2=="America", 1, ifelse(cont=="Brazil", 2, ifelse(cont=="Africa", 3, 4)))) %>%
    arrange(id, ctry2) %>%
    # Select columns
    dplyr::select(cont2, ctry2, region, code, nfor, ndef, nforHa, ndefHa)

## Save results
write.table(samp_size_tab, file=file.path("Analysis", dataset, "results/samp_size.csv"), sep=",", row.names=FALSE)

## ===================
## Forest cover change
## ===================

## Create table to store results
fcc_tab <- data.frame(matrix(NA, nrow=nctry, ncol=7))
names(fcc_tab) <- c("cont", "iso3", "for2000", "for2005", "for2010", "for2015", "for2020")

## Loop on countries
for (i in 1:nctry) {
    iso <- iso3[i]
    continent <- as.character(ctry_df$cont_run[ctry_df$iso3==iso])
    dir <- file.path(dir_fdb, dataset, continent)
    ## Forest cover change
    f_name <- file.path(dir, iso, "/output/forest_cover.txt")
    fcc_df <- read.table(f_name, header=FALSE, sep=",", stringsAsFactors=FALSE)
    area <- round(fcc_df[, 1])
    ## Fill in the table
    fcc_tab[i,1:2] <- c(continent, iso)
    fcc_tab[i,3:7] <- area
}

## Annual defor
TI <- 2020-2010  ## Time-interval
fcc_tab2 <- fcc_tab %>%
    mutate(andef=round((for2010-for2020)/TI)) %>%
    mutate(pdef=round(100*(1-(1-(for2010-for2020)/for2010)^(1/TI)), 1)) %>%
    mutate(for2030=pmax(0, for2020-10*andef), for2035=pmax(0, for2020-15*andef),
           for2040=pmax(0, for2020-20*andef),
           for2050=pmax(0, for2020-30*andef), for2055=pmax(0, for2020-35*andef),
           for2060=pmax(0, for2020-40*andef),
           for2070=pmax(0, for2020-50*andef),
           for2080=pmax(0, for2020-60*andef), for2085=pmax(0, for2020-65*andef),
           for2090=pmax(0, for2020-70*andef),
           for2100=pmax(0, for2020-80*andef)) %>%
    # Year during which forest should have disappeared
    mutate(yrdis=floor(2020 + for2020/andef))

## Corrections for Brazil with deforestation diffusion
if (dataset=="gfc2020_70") {
    fname_BRA <- "Brazil/fcc_BRA_gfc.csv"
} else if (dataset=="jrc2020") {
    fname_BRA <- "Brazil/fcc_BRA_jrc.csv"
}
fcc_BRA <- read.table(file.path(dir_fdb, dataset, fname_BRA), sep=",",
                      header=TRUE, stringsAsFactors=FALSE)
if (all(fcc_BRA$iso3==fcc_tab2$iso3[fcc_tab2$cont=="Brazil"])) { # Check order
    fcc_tab2[fcc_tab2$cont=="Brazil", c(10:ncol(fcc_tab2))] <- round(fcc_BRA[, c(seq(7, 27, by=2), 30)])
}

## Save results
write.table(fcc_tab2, file=file.path("Analysis", dataset, "results/forest_cover_change.csv"), sep=",", row.names=FALSE)

## Graph showing the % decrease per continent with time compared with 2000.

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
    dir <- file.path(dir_fdb, dataset, continent)
    ## Parameter estimates
    f_name <- file.path(dir, iso, "/output/summary_hSDM.txt")
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
write.table(par_tab, file=file.path("Analysis", dataset, "results/parameter_estimates.csv"), sep=",", row.names=FALSE)

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
    dir <- file.path(dir_fdb, dataset, continent)
    ## Parameter estimates
    f_name <- file.path(dir, iso, "/output/summary_hSDM.txt")
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
write.table(parea_tab, file=file.path("Analysis", dataset,"results/parea_estimates.csv"), sep=",", row.names=FALSE)

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
    dir <- file.path(dir_fdb, dataset, continent)
    ## Parameter estimates
    f_name <- file.path(dir, iso, "/output/summary_hSDM.txt")
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
write.table(road_tab, file=file.path("Analysis", dataset, "results/road_estimates.csv"), sep=",", row.names=FALSE)

## =========================
## Significance PA and roads
## =========================

## Significance PA
parea_tab <- parea_tab %>%
    mutate(sign=ifelse(CI_low * CI_high > 0 & !is.na(Mean), 1, 0))
## Percentage of country for which the effet of protected areas is significant
perc_sign_PA <- 100*sum(parea_tab$sign==1)/nrow(parea_tab)
## Weighted percentage with forest size in 2010
fcc_tab <- read.table(file.path("Analysis", dataset, "results/forest_cover_change.csv"), header=TRUE, sep=",")
weights <- fcc_tab$for2010
perc_sign_w_PA <- 100*sum((parea_tab$sign==1)*weights)/sum(weights)

## Significance road
road_tab <- road_tab %>%
    mutate(sign=ifelse(CI_low * CI_high > 0 & !is.na(Mean), 1, 0))
## Percentage of country for which the effet of protected areas is significant
perc_sign_road <- 100*sum(road_tab$sign==1)/nrow(road_tab)
## Weighted percentage with forest size in 2010
fcc_tab <- read.table(file.path("Analysis", dataset, "results/forest_cover_change.csv"), header=TRUE, sep=",")
weights <- fcc_tab$for2010
perc_sign_w_road <- 100*sum((road_tab$sign==1)*weights)/sum(weights)

## Save results
nctry <- length(fcc_tab$iso3)
perc <- c(perc_sign_PA, perc_sign_road)
perc_w <- c(perc_sign_w_PA, perc_sign_w_road) 
sign_PA_road <- data.frame(var=c("PA","road"), nctry=nctry, perc=round(perc), perc_w=round(perc_w))

## Save results
write.table(sign_PA_road, file=file.path("Analysis", dataset, "results/sign_PA_road.csv"), sep=",", row.names=FALSE)

## =====================================
## Carbon emissions (in tonnes=10e6 g C)
## =====================================

## Create table to store results
Cem_tab <- data.frame(matrix(NA, nrow=nctry, ncol=14), stringsAsFactors=FALSE)
C_var <- paste0("C", 2020 + c(0, 10, 15, 20, 30, 35, 40, 50, 60, 65, 70, 80))
names(Cem_tab) <- c("cont", "iso3", C_var)

## Loop on countries
for (i in 1:nctry) {
    iso <- iso3[i]
    continent <- as.character(ctry_df$cont_run[ctry_df$iso3==iso])
    dir <- file.path(dir_fdb, dataset, continent)
    ## Carbon emissions
    f_name <- file.path(dir, iso, "output/C_emissions.csv")
    Cem_df <- read.table(f_name, header=TRUE, sep=",", stringsAsFactors=FALSE)
    ## Fill in the table
    Cem_tab[i, 1:2] <- c(continent, iso)
    Cem_tab[i, 3:14] <- Cem_df$C
}

## Save results
write.table(Cem_tab, file=file.path("Analysis", dataset, "results/C_emissions.csv"), sep=",", row.names=FALSE)

## ========================
## Summary carbon emissions
## ========================

## Summarize results
Cem_tab2 <- Cem_tab %>%
    mutate(cont=ifelse(cont=="Brazil", "America", cont)) %>%
    group_by(cont) %>%
    summarize_at(vars(C_var), sum) %>%
    bind_rows(data.frame(cont="TOTAL",
                         summarize_at(Cem_tab, vars(C_var), sum),
                         stringsAsFactors=FALSE)) %>%
    mutate_at(vars(C2020:C2100), function(x){x*1e-9})  # Results in PgC

## Save results
write.table(Cem_tab2, file=file.path("Analysis", dataset, "results/C_emissions_summary.csv"), sep=",", row.names=FALSE)

## In 2019-2050, 22.4 billions of tonnes of C emitted = 22.4 PgC (0.72 PgC/year on 2019-2050)
##
## References for comparison:
## Baccini et al. (2017) only 0.8 PgC/year for 2003–2014. Land use and land-cover change (LULCC) are believed to release between 0.81 and 1.14 PgC/yr.
## Baccini 2012 et al. reported 1 PgC/year on the period 2000-2010.
## See Van Der Werf, G. R. et al. CO2 emissions from forest loss. Nature Geosci. 2, 737–738 (2009).
## Friedlingstein, P. et al. Update on CO2 emissions. Nature Geosci. 3, 811–812 (2010).
## Le Quéré, C. et al. Trends in the sources and sinks of carbon dioxide. Nature Geosci. 2, 831–836 (2009).

## ======================
## Carbon emission trends
## ======================

## Compute carbon emission trends in the future
C_trend <- Cem_tab2 %>%
    mutate(T10_20=C2020/10,
           T20_30=C2030/10,
           T30_40=(C2040-C2030)/10,
           T40_50=(C2050-C2040)/10,
           T50_60=(C2060-C2050)/10,
           T60_70=(C2070-C2060)/10,
           T70_80=(C2080-C2070)/10,
           T80_90=(C2090-C2080)/10,
           T90_100=(C2100-C2090)/10) %>%
    select(cont, T10_20:T90_100)

## Carbon emissions should continue to increase: from 0.66 PgC/yr on 2019-2035 to 0.814 PgC/yr on 2050-2085.
## Deforestation of forest areas with higher carbon stocks in the future.
## Higher carbon stocks because of environmental/elevation gradient (see Asner article) and because of remote, less degraded forests.
## Will decrease in Asia because many countries with no more forests after 2085.

## Deforestation => increase in C source with time (deforestation of forest with higher carbon stocks).
## Climate change => decrease in C sink with time (higher mortality), see Hubau2020.
## The result is that forests will likely become a major C source in the future. 

## Save results
write.table(C_trend, file=file.path("Analysis", dataset, "results/C_trend.csv"), sep=",", row.names=FALSE)

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
