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
require(tidyr)
require(here)
require(ggplot2)
require(wesanderson)
require(readr)
require(glue)

## Set working directory
setwd(here())

## inv_logit function
inv_logit <- function (x, min=0, max=1) {
    p <- exp(x)/(1 + exp(x))
    p <- ifelse(is.na(p) & !is.na(x), 1, p)
    p * (max - min) + min
}

## Dataset
##dataset <- "gfc2020_70" 
dataset <- "jrc2020"
dir.create(here("Analysis", dataset), recursive=TRUE)

## Result directory
dir_fdb <- "/home/forestatrisk-tropics"

## =================
## Countries info
## =================

## Load country info (encoding pb for ctry names)
ctry_df <- read.csv2(here("Analysis", "data", "ctry_run.csv"), header=TRUE, sep=";", encoding="UTF-8")

## List of countries
iso3 <- ctry_df$iso3
if (dataset=="gfc2020_70") {
    iso3 <- as.character(iso3[-which(iso3 %in% c("STP", "SEN", "GMB"))])
} else if (dataset=="jrc2020") {
    iso3 <- as.character(iso3[-which(iso3 %in% c("STP"))])
}
nctry <- length(iso3)

## ========================================================
## Data with observations and predictions for all countries
## ========================================================

## Create table to store data
data_allctry_tab <- data.frame()

## Loop on countries
for (i in 1:nctry) {
  iso <- iso3[i]
  continent <- as.character(ctry_df$cont_run[ctry_df$iso3==iso])
  dir <- file.path(dir_fdb, dataset, continent)
  ## Area info
  area_cont <- as.character(ctry_df$area_cont[ctry_df$iso3==iso])
  area_ctry <- as.character(ctry_df$area_ctry[ctry_df$iso3==iso])
  area_name <- as.character(ctry_df$area_name[ctry_df$iso3==iso])
  area_code <- as.character(ctry_df$area_code[ctry_df$iso3==iso])
  ## Sample size
  f_name <- file.path(dir, iso, "/output/obs_pred.csv")
  obs_pred_df <- read.table(f_name, header=TRUE, sep=",", stringsAsFactors=FALSE)
  obs_pred_df$iso3 <- iso
  obs_pred_df$continent <- continent
  obs_pred_df$area_cont <- area_cont
  obs_pred_df$area_ctry <- area_ctry
  obs_pred_df$area_name <- area_name
  obs_pred_df$area_code <- area_code
  
  ## Fill in the table
  data_allctry_tab <- rbind(data_allctry_tab, obs_pred_df)
}

## Corrections for "VIR"
data_allctry_tab <- data_allctry_tab %>%
  ## No river data in VIR: replace 65535 with NA
  dplyr::mutate(dist_river=ifelse(iso3=="VIR" & dist_river==65535, NaN, dist_river))

## Save results
f <- here("Analysis", dataset, "data_allctry.csv")
write.table(data_allctry_tab, file=f, sep=",", row.names=FALSE)

## ===================
## Forest cover change
## ===================

## Annual deforestation with uncertainty
f <- here("Intensity", "output", "d_uncertainty.csv")
d_ci <- read.table(f, sep=",", header=TRUE)
d_ci <- d_ci %>%
  select(area_code, d_se, d_mean, d_min, d_max)

## Create table to store results
fcc_tab <- data.frame(matrix(NA, nrow=nctry, ncol=9))
names(fcc_tab) <- c("area_cont", "area_ctry", "area_name", "area_code", "for2000", "for2005", "for2010", "for2015", "for2020")

## Loop on countries
for (i in 1:nctry) {
  ## File path
  iso <- iso3[i]
  continent <- as.character(ctry_df$cont_run[ctry_df$iso3==iso])
  dir <- file.path(dir_fdb, dataset, continent)
  ## Area info
  area_cont <- as.character(ctry_df$area_cont[ctry_df$iso3==iso])
  area_ctry <- as.character(ctry_df$area_ctry[ctry_df$iso3==iso])
  area_name <- as.character(ctry_df$area_name[ctry_df$iso3==iso])
  area_code <- as.character(ctry_df$area_code[ctry_df$iso3==iso])
  ## Forest cover change
  f_name <- file.path(dir, iso, "/output/forest_cover.txt")
  fcc_df <- read.table(f_name, header=FALSE, sep=",", stringsAsFactors=FALSE)
  area_df <- round(fcc_df[, 1])
  ## Fill in the table
  fcc_tab[i,1:4] <- cbind(area_cont, area_ctry, area_name, area_code)
  fcc_tab[i,5:9] <- area_df
}

## Join forest and deforestation tables
fcc_tab2 <- fcc_tab %>%
  dplyr::left_join(d_ci, by="area_code")

## Simulations with uncertainty
sim <- c("mean", "min", "max")
nsim <- length(sim)

## Loop on simulations
for (i in 1:nsim) {
  
  ## Simulation id
  s <- sim[i]
  d_s <- fcc_tab2[[glue("d_{s}")]]
  
  ## Annual defor
  TI <- 2020-2010  ## Time-interval
  fcc_tab3 <- fcc_tab2 %>%
    mutate(andef=d_s) %>%
    mutate(pdef=round(100*(1-(1-(andef*TI)/for2010)^(1/TI)), 1)) %>%
    mutate(for2030=pmax(0, for2020-10*andef), for2035=pmax(0, for2020-15*andef),
           for2040=pmax(0, for2020-20*andef),
           for2050=pmax(0, for2020-30*andef), for2055=pmax(0, for2020-35*andef),
           for2060=pmax(0, for2020-40*andef),
           for2070=pmax(0, for2020-50*andef),
           for2080=pmax(0, for2020-60*andef), for2085=pmax(0, for2020-65*andef),
           for2090=pmax(0, for2020-70*andef),
           for2100=pmax(0, for2020-80*andef)) %>%
    # Year during which all the forest will have disappeared
    mutate(yrdis=floor(2020 + for2020/andef))
  
  ## Corrections for Brazil with deforestation diffusion
  if (dataset=="gfc2020_70") {
    fname_BRA <- glue("Brazil/fcc_BRA_gfc_{s}.csv")
  } else if (dataset=="jrc2020") {
    fname_BRA <- glue("Brazil/fcc_BRA_jrc_{s}.csv")
  }
  fcc_BRA <- read.table(file.path(dir_fdb, dataset, fname_BRA), sep=",",
                        header=TRUE, stringsAsFactors=FALSE)
  ## Check order and replace with correct values for Brazil
  codes_fcc <- paste0("BRA-",fcc_tab3$area_code[fcc_tab3$area_ctry=="Brazil"])
  codes_BRA <- fcc_BRA$iso3
  if (all(codes_fcc==codes_BRA)) {
    fcc_tab3[fcc_tab3$area_ctry=="Brazil", c(16:ncol(fcc_tab3))] <- round(fcc_BRA[, c(seq(7, 27, by=2), 30)])
  }
  
  ## Sort continents, and select col
  fcc_tab4 <- fcc_tab3 %>%
    # Id
    mutate(id=ifelse(area_cont=="America", 1, ifelse(area_cont=="Africa", 2, 3))) %>%
    arrange(id, area_name) %>%
    # Select columns
    dplyr::select(area_cont, area_ctry, area_name, area_code, for2000:yrdis)
  
  ## Save results
  f_name <- glue("forest_cover_change_{s}.csv")
  f <- here("Analysis", dataset, f_name)
  write.table(fcc_tab4, file=f, sep=",", row.names=FALSE)
  ## Copy for manuscript
  f_doc <- here("Manuscript", "Supplementary_Materials", "tables", f_name)
  file.copy(from=f, to=f_doc, overwrite=TRUE)
  ## Copy for data
  f_doc <- here("Manuscript", "Supplementary_Data", "tables", f_name)
  file.copy(from=f, to=f_doc, overwrite=TRUE)

}

## ====================================================
## Number of countries with no more forest in 2100
## ====================================================

## Import data
f_in <- here("Analysis", dataset, "forest_cover_change_mean.csv")
df <- read.table(f_in, header=TRUE, sep=",")

## Arrange data
df <- df %>%
  dplyr::mutate(loss21=round(100*(for2000-for2100)/for2000))

## Compute number of countries (apart from India and Brazil)
nctry_loss21 <- df %>%
  dplyr::filter(!area_ctry %in% c("Brazil", "India")) %>%
  group_by(area_cont) %>%
  summarize(n=sum(loss21==100))

nctry_loss21_large <- df %>%
  dplyr::filter(!area_ctry %in% c("Brazil", "India")) %>%
  dplyr::filter(for2020 >= 1e6) %>%
  group_by(area_cont) %>%
  summarize(n=sum(loss21==100))
  
## For Brazil
nstate_BRA_loss21 <- df %>%
  dplyr::filter(area_ctry=="Brazil") %>%
  group_by(area_ctry) %>%
  summarize(n=sum(loss21==100)) %>%
  rename(area_cont=area_ctry)

nstate_BRA_loss21_large <- df %>%
  dplyr::filter(area_ctry=="Brazil") %>%
  dplyr::filter(for2020 >= 1e6) %>%
  group_by(area_ctry) %>%
  summarize(n=sum(loss21==100)) %>%
  rename(area_cont=area_ctry)

## For India
nstate_IND_loss21 <- df %>%
  dplyr::filter(area_ctry=="India") %>%
  group_by(area_ctry) %>%
  summarize(n=sum(loss21==100)) %>%
  rename(area_cont=area_ctry)

nstate_IND_loss21_large <- df %>%
  dplyr::filter(area_ctry=="India") %>%
  dplyr::filter(for2020 >= 1e6) %>%
  group_by(area_ctry) %>%
  summarize(n=sum(loss21==100)) %>%
  rename(area_cont=area_ctry)

## All regions
nall <- df %>%
  summarize(n=sum(loss21==100)) %>%
  pull()
nctry_all <- tibble(area_cont="All", n=nall)

nall <- df %>%
  dplyr::filter(for2020 >= 1e6) %>%
  summarize(n=sum(loss21==100)) %>%
  pull()
nctry_all_large <- tibble(area_cont="All", n=nall)

## Combine results
nctry_loss21_bycont <- nctry_loss21 %>%
  dplyr::bind_rows(nstate_BRA_loss21) %>%
  dplyr::bind_rows(nstate_IND_loss21) %>%
  dplyr::bind_rows(nctry_all)

nctry_loss21_bycont_large <- nctry_loss21_large %>%
  dplyr::bind_rows(nstate_BRA_loss21_large) %>%
  dplyr::bind_rows(nstate_IND_loss21_large) %>%
  dplyr::bind_rows(nctry_all_large)

## Save results
f_out <- here("Analysis", dataset, "nctry_loss21_bycont.csv")
write.table(nctry_loss21_bycont, file=f_out, sep=",", row.names=FALSE)
f_out <- here("Analysis", dataset, "nctry_loss21_bycont_large.csv")
write.table(nctry_loss21_bycont_large, file=f_out, sep=",", row.names=FALSE)

## Number of countries with more than 1 Mha forest in 2020
nctry_large <- df %>%
  dplyr::filter(!area_ctry %in% c("Brazil", "India")) %>%
  dplyr::filter(for2020 >= 1e6) %>%
  summarize(n=n()) %>% pull()
nctry_large <- nctry_large + 2 # For Brazil and India
msg <- paste0("Number of countries with more than 1 Mha ",
              "forest in 2020: ", nctry_large)
# Save
f_out <- here("Analysis", dataset, "nctry_large.txt")
sink(f_out)
cat(msg)
sink()

## ====================================================
## Forest cover change summarized per region
## ====================================================

## Simulations with uncertainty
sim <- c("mean", "min", "max")
nsim <- length(sim)

## Loop on simulations
for (i in 1:nsim) {
  
  ## Simulation id
  s <- sim[i]

  ## Load previous fcc table
  f <- here("Analysis", dataset, glue("forest_cover_change_{s}.csv"))
  fcc_df <- read.table(f, header=TRUE, sep=",")
  fcc_df <- fcc_df %>%
    select(-d_se, -d_mean, -d_min, -d_max)
  
  ## For each continent
  fcc_cont <- fcc_df %>%
    dplyr::group_by(area_cont) %>%
    dplyr::summarise_if(is.numeric, list(sum=sum, max=max)) %>%
    dplyr::select(area_cont, for2000_sum:andef_sum, for2030_sum:for2100_sum, yrdis_max) %>%
    dplyr::mutate(id_cont=c(2, 1, 3)) %>%
    dplyr::arrange(id_cont) %>%
    dplyr::select(-id_cont) %>%
    dplyr::rename_at(.vars=vars(starts_with("for")), .funs=substr, start=1, stop=7) %>%
    dplyr::rename(andef=andef_sum, yrdis=yrdis_max)
  
  ## For all continents
  fcc_all <- fcc_df %>%
    dplyr::summarise_if(is.numeric, list(sum=sum, max=max)) %>%
    dplyr::select(for2000_sum:andef_sum, for2030_sum:for2100_sum, yrdis_max) %>%
    dplyr::mutate(area_cont="All continents") %>%
    dplyr::relocate(area_cont, .before=for2000_sum) %>%
    dplyr::rename_at(.vars=vars(starts_with("for")), .funs=substr, start=1, stop=7) %>%
    dplyr::rename(andef=andef_sum, yrdis=yrdis_max)
  
  ## For Brazil
  fcc_bra <- fcc_df %>%
    dplyr::filter(area_ctry=="Brazil") %>%
    dplyr::group_by(area_ctry) %>%
    dplyr::summarise_if(is.numeric, list(sum=sum, max=max)) %>%
    dplyr::select(area_ctry, for2000_sum:andef_sum, for2030_sum:for2100_sum, yrdis_max) %>%
    dplyr::rename(area_cont=area_ctry) %>%
    dplyr::rename_at(.vars=vars(starts_with("for")), .funs=substr, start=1, stop=7) %>%
    dplyr::rename(andef=andef_sum, yrdis=yrdis_max)
  
  ## For India
  fcc_ind <- fcc_df %>%
    dplyr::filter(area_ctry=="India") %>%
    dplyr::group_by(area_ctry) %>%
    dplyr::summarise_if(is.numeric, list(sum=sum, max=max)) %>%
    dplyr::select(area_ctry, for2000_sum:andef_sum, for2030_sum:for2100_sum, yrdis_max) %>%
    dplyr::rename(area_cont=area_ctry) %>%
    dplyr::rename_at(.vars=vars(starts_with("for")), .funs=substr, start=1, stop=7) %>%
    dplyr::rename(andef=andef_sum, yrdis=yrdis_max)
  
  ## DRC and Indonesia
  fcc_DRC_IDN <- fcc_df %>% 
    dplyr::filter(area_ctry %in% c("DRC", "Indonesia")) %>%
    dplyr::select(area_ctry, for2000:andef, for2030:for2100, yrdis) %>%
    dplyr::rename(area_cont=area_ctry)
  
  ## Combine
  TI <- 2020-2010  ## Time-interval
  fcc_comb <- fcc_ind %>%
    # Add DRC and Indonesia
    dplyr::bind_rows(fcc_DRC_IDN) %>%
    # Add Brazil
    dplyr::bind_rows(fcc_bra) %>%
    # Add continents
    dplyr::bind_rows(fcc_cont) %>%
    dplyr::bind_rows(fcc_all) %>%
    # Compute pdef
    dplyr::mutate(pdef=round(100*(1-(1-(andef*TI)/for2010)^(1/TI)), 1)) %>%
    # Compute loss21
    dplyr::mutate(loss21=100*(for2000-for2100)/for2000) %>%
    # Arrange columns
    dplyr::relocate(pdef, .after=andef) %>%
    dplyr::relocate(loss21, .before=yrdis)
  
  ## Save result
  f_name <- glue("fcc_hist_region_{s}.csv")
  f <- here("Analysis", dataset, f_name)
  write.table(fcc_comb, file=f, sep=",", row.names=FALSE)
  ## Copy for manuscript
  f_doc <- here("Manuscript", "Supplementary_Materials", "tables", f_name)
  file.copy(from=f, to=f_doc, overwrite=TRUE)
  f_doc <- here("Manuscript", "Article", "tables", f_name)
  file.copy(from=f, to=f_doc, overwrite=TRUE)

}

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
fcc_tab <- read.table(here("Analysis", dataset, "forest_cover_change_mean.csv"),
                      header=TRUE, sep=",")
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
    dplyr::select(-cont_id)

## Save results
f1 <- here("Analysis", dataset, "performance_index.csv")
f2 <- here("Analysis", dataset, "perf_mod.csv")
f3 <- here("Analysis", dataset, "perf_cont_mod.csv")
write.table(ind_tab, file=f1, sep=",", row.names=FALSE)
write.table(perf_mod, file=f2, sep=",", row.names=FALSE)
write.table(perf_cont_mod, file=f3, sep=",", row.names=FALSE)
## Copy for manuscript
f1_doc <- here("Manuscript", "Supplementary_Materials", "tables", "performance_index.csv")
f2_doc <- here("Manuscript", "Supplementary_Materials", "tables", "perf_mod.csv")
f3_doc <- here("Manuscript", "Supplementary_Materials", "tables", "perf_cont_mod.csv")
file.copy(from=f1, to=f1_doc, overwrite=TRUE)
file.copy(from=f2, to=f2_doc, overwrite=TRUE)
file.copy(from=f3, to=f3_doc, overwrite=TRUE)

## =================
## Sample size
## =================

## Create table to store results
samp_size_tab <- data.frame(matrix(NA, nrow=nctry, ncol=6))
names(samp_size_tab) <- c("area_cont", "area_ctry", "area_name", "area_code",
                          "nfor", "ndef")

## Loop on countries
for (i in 1:nctry) {
    ## File path
    iso <- iso3[i]
    continent <- as.character(ctry_df$cont_run[ctry_df$iso3==iso])
    dir <- file.path(dir_fdb, dataset, continent)
    ## Area info
    area_cont <- as.character(ctry_df$area_cont[ctry_df$iso3==iso])
    area_ctry <- as.character(ctry_df$area_ctry[ctry_df$iso3==iso])
    area_name <- as.character(ctry_df$area_name[ctry_df$iso3==iso])
    area_code <- as.character(ctry_df$area_code[ctry_df$iso3==iso])
    ## Sample size
    f_name <- file.path(dir, iso, "output/sample_size.csv")
    samp_size_df <- read.table(f_name, header=TRUE, sep=",", stringsAsFactors=FALSE)
    nfor <- samp_size_df$n[samp_size_df$var=="nfor"]
    ndef <- samp_size_df$n[samp_size_df$var=="ndefor"]
    ## Fill in the table
    samp_size_tab[i, 1:4] <- cbind(area_cont, area_ctry, area_name, area_code)
    samp_size_tab[i, 5:6] <- cbind(nfor, ndef)
}

## Equivalence in ha
samp_size_tab$nforHa <- round(samp_size_tab$nfor*30*30/10000)
samp_size_tab$ndefHa <- round(samp_size_tab$ndef*30*30/10000)

## Sort continents, select col, and add total
samp_size_tab2 <- samp_size_tab %>%
    # Id
    mutate(id=ifelse(area_cont=="America", 1, ifelse(area_cont=="Africa", 2, 3))) %>%
    arrange(id, area_name) %>%
    # Select columns
    dplyr::select(area_cont, area_name, area_code, nfor, ndef, nforHa, ndefHa) %>%
    # Add total
    dplyr::add_row(area_cont="All continents", area_name="TOTAL", area_code="", 
            nfor=sum(.$nfor), ndef=sum(.$ndef), nforHa=sum(.$nforHa), ndefHa=sum(.$ndefHa))

## Save results
f <- here("Analysis", dataset, "samp_size.csv")
write.table(samp_size_tab2, file=f, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "samp_size.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## ===================================================
## Balanced data-set given for2010 for each study-area
## ===================================================

## Sample size
f <- here("Analysis", dataset, "samp_size.csv")
samp_size <- read_delim(f, delim=",") %>%
  dplyr::filter(area_cont!="All continents")

## Forest cover
f <- here("Analysis", dataset, "forest_cover_change_mean.csv")
fcc_tab <- read_delim(f, delim=",")

## nsamp_est function
nsamp_est_f <- function(nsamp, weight) {
  weight_max <- max(weight)
  nsamp_max <- unique(nsamp[weight==weight_max])
  nsamp_est <- nsamp_max*(weight/weight_max)
  return(as.integer(round(nsamp_est)))
}

## Balanced data-set
samp_df <- samp_size %>%
  mutate(nsamp=nfor+ndef) %>%
  select(-nfor, -ndef, -nforHa, -ndefHa) %>%
  mutate(for2010=fcc_tab$for2010) %>%
  mutate(weight=for2010/sum(for2010)) %>%
  mutate(nsamp_est=nsamp_est_f(nsamp, weight)) %>%
  mutate(nsamp_prop=pmin(nsamp, nsamp_est)) %>%
  arrange(area_code)

## Sample from full data-set
## see here: https://jennybc.github.io/purrr-tutorial/ls12_different-sized-samples.html
library(purrr)  # for map() function
library(tidyr)  # for nest() function
set.seed(4321)
f <- here("Analysis", dataset, "data_allctry.csv")
data_allctry <- read_delim(f, delim=",")
data_allctry_prop <- data_allctry %>%
  group_by(area_code) %>%
  nest() %>%
  ungroup() %>%
  arrange(area_code) %>%  # Check order
  mutate(nsamp_prop=samp_df$nsamp_prop) %>%
  mutate(samp=map2(data, nsamp_prop, sample_n)) %>%
  select(-data) %>%
  unnest(samp)

## Save balanced data-set
f <- here("Analysis", dataset, "data_allctry_prop.csv")
write_delim(data_allctry_prop, f, delim=",")

## ===============================================
## Projecting percentage of forest loss per region
## ===============================================

## Simulations with uncertainty
sim <- c("mean", "min", "max")
nsim <- length(sim)

## Loop on simulations
for (j in 1:nsim) {
  
  ## Simulation id
  s <- sim[j]
  
  ## Load previous fcc table
  f <- here("Analysis", dataset, glue("forest_cover_change_{s}.csv"))
  fcc_df <- read.table(f, header=TRUE, sep=",")
  
  ## All study-areas
  fcc_hist <- fcc_df %>%
    dplyr::select(1:3, for2000, for2010, for2020, andef)
  
  ## Brazil
  fcc_hist_bra <- fcc_hist %>%
    dplyr::filter(area_ctry=="Brazil") %>%
    dplyr::group_by(area_cont, area_ctry) %>%
    dplyr::summarise_if(is.numeric, sum) %>%
    dplyr::mutate(area_name="Brazil") %>%
    dplyr::relocate(area_name, .after=area_ctry)
  
  ## Replace Brazilian states data with data from whole Brazil.
  ## This is necessary as the annual deforestation is constant
  ## at the country scale but not at the federal state scale.
  fcc_hist <- fcc_hist %>%
    dplyr::filter(area_ctry!="Brazil") %>%
    dplyr::bind_rows(fcc_hist_bra)
  
  ## Historical deforestation per continent
  fcc_hist_cont <- fcc_hist %>%
    dplyr::group_by(area_cont) %>%
    dplyr::summarise_if(is.numeric, sum) %>%
    dplyr::select(-andef) %>%
    tidyr::pivot_longer(c(for2000, for2010, for2020),
                        names_to="year", values_to="fc") %>%
    dplyr::rename(cont=area_cont) %>%
    dplyr::mutate(year=as.integer(substr(year, 4, 7)))
  
  ## Project forest cover until year 2400
  fc_proj_cont_long <- data.frame(cont=rep(c("Africa","America","Asia"), 381),
                                  year=rep(c(2020:2400), each=3),
                                  fc=NA)
  
  ## Loop on year (starting from 0 for year 2020)
  for (i in 0:380) {
    fc_proj <- pmax(0, fcc_hist$for2020-(i*fcc_hist$andef))
    fc_yr_cont <- fcc_hist %>%
      dplyr::mutate(fc_proj=fc_proj) %>%
      dplyr::select(area_cont, fc_proj) %>%
      dplyr::group_by(area_cont) %>%
      dplyr::summarise_all(sum) %>%
      dplyr::select(fc_proj) %>%
      dplyr::pull()
    fc_proj_cont_long$fc[fc_proj_cont_long$year==2020+i] <- fc_yr_cont
  }
  
  ## Bind with historical deforestation
  fc_cont <- fc_proj_cont_long %>%
    ## Remove 2020 to avoid repetition
    dplyr::filter(year!=2020) %>%
    dplyr::bind_rows(fcc_hist_cont) %>%
    dplyr::arrange(year, cont)
  
  ## Percentage of loss compared to fc2000
  fc_perc_cont <- fc_cont %>% 
    dplyr::group_by(cont) %>% 
    dplyr::mutate(perc=100*(fc[year==2000]-fc)/fc[year==2000]) %>%
    dplyr::ungroup()
  ## Save results
  f_name <- glue("perc_loss_cont_{s}.csv")
  f <- here("Analysis", dataset, f_name)
  write.table(fc_perc_cont, file=f, sep=",", row.names=FALSE)
  ## Copy for manuscript
  f_doc <- here("Manuscript", "Supplementary_Materials", "tables", f_name)
  file.copy(from=f, to=f_doc, overwrite=TRUE)

}

## ===================================
## Plot change in percentage with time
## ====================================

## Load data
f <- here("Analysis", dataset, "perc_loss_cont_mean.csv")
fc_perc_cont <- read.table(f, header=TRUE, sep=",")
df_hist <- fc_perc_cont %>%
    dplyr::filter(year %in% c(2000, 2010, 2020))
df_proj <- fc_perc_cont %>%
    dplyr::filter(!(year %in% c(2000, 2010, 2020)))

## Simulations with uncertainty
# min
f <- here("Analysis", dataset, "perc_loss_cont_min.csv")
fc_perc_cont_min <- read.table(f, header=TRUE, sep=",")
df_proj_min <- fc_perc_cont_min %>%
  dplyr::filter(!(year %in% c(2000, 2010, 2020)))
# max
f <- here("Analysis", dataset, "perc_loss_cont_max.csv")
fc_perc_cont_max <- read.table(f, header=TRUE, sep=",")
df_proj_max <- fc_perc_cont_max %>%
  dplyr::filter(!(year %in% c(2000, 2010, 2020)))
# ci
df_proj_ci <- df_proj %>%
  mutate(perc_min=df_proj_min[["perc"]],
         perc_max=df_proj_max[["perc"]])

## mytheme
mytheme <- theme(
    axis.title=element_text(size=12),
    axis.text=element_text(size=10),
    legend.title=element_text(size=12),
    legend.text=element_text(size=10),
    legend.position=c(0.95, 0.05),
    legend.justification=c(1, 0),
    legend.background=element_rect(fill="transparent"))

## Plot
p <- ggplot(aes(x=year, y=perc, group=cont, col=cont), data=df_hist) +
    geom_point(size=1) +
    geom_ribbon(aes(ymin=perc_min, ymax=perc_max, group=cont, fill=cont),
                alpha=0.2, data=df_proj_ci, linetype=0) + 
    geom_line(data=df_proj_ci, size=0.8) +
    xlab("Year") + ylab("Percentage of forest cover loss\n(in comparison with year 2000)") +
    scale_color_manual(values=wes_palette("Moonrise2")[c(3, 2, 1)],
                       name="Continents",
                       breaks=c("America", "Africa", "Asia"),
                       labels=c("America", "Africa", "Asia")) +
    scale_fill_manual(values=wes_palette("Moonrise2")[c(3, 2, 1)],
                      name="Continents",
                      breaks=c("America", "Africa", "Asia"),
                      labels=c("America", "Africa", "Asia")) +
    ylim(0,100) +
    geom_hline(yintercept=75) +
    theme_bw() + mytheme

## Save results
f <- here("Analysis", dataset, "perc_loss_cont.png")
ggsave(f, p, width=16.6, height=10, units="cm", dpi=300)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "figures",
              "perc_loss_cont.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## =========================================================================
## Forest cover projections including yr75dis per region
## yrdis75: year during which 75% of for2000 wil have disappeared per region
## =========================================================================

## Simulations with uncertainty
sim <- c("mean", "min", "max")
nsim <- length(sim)

## Loop on simulations
for (j in 1:nsim) {
  
  ## Simulation id
  s <- sim[j]

  ## Load data
  f_name <- glue("perc_loss_cont_{s}.csv")
  f <- here("Analysis", dataset, f_name)
  fc_perc_cont <- read.table(f, header=TRUE, sep=",")
  
  ## yr75dis
  yr75dis_cont <- fc_perc_cont %>%
    dplyr::filter(perc>=75) %>%
    dplyr::group_by(cont) %>%
    dplyr::filter(year==min(year)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(yrdis=year-1) %>%
    dplyr::arrange(cont)
  
  ## loss21
  loss21_cont <- fc_perc_cont %>%
    dplyr::filter(year==2100)%>%
    dplyr::arrange(cont)
  
  ## Results in one table for continents
  fc_proj_cont <- fc_proj_cont_long %>%
    dplyr::filter(year %in% c(2040, 2050, 2060, 2080, 2100)) %>%
    tidyr::pivot_wider(names_from=year,
                       values_from=fc, names_prefix="for") %>%
    dplyr::mutate(loss21=loss21_cont$perc) %>%
    dplyr::mutate(yr75dis=yr75dis_cont$yrdis) %>%
    dplyr::mutate(id=c(2,1,3)) %>%
    dplyr::arrange(id) %>%
    dplyr::select(-id)
  
  ## Results for Brazil
  fcc_bra <- fcc_df %>%
    dplyr::filter(area_ctry=="Brazil") %>%
    dplyr::select(-pdef, -yrdis) %>%
    dplyr::group_by(area_ctry) %>%
    dplyr::summarise_if(is.numeric, sum) %>%
    dplyr::rename(cont=area_ctry) %>%
    dplyr::mutate(loss21=100*(for2000-for2100)/for2000) %>%
    dplyr::mutate(yr75dis=floor(2000+(for2000*0.75)/andef)) %>%
    dplyr::select(cont, for2040, for2050, for2060, for2080, for2100,
                  loss21, yr75dis)
  
  ## Results for India
  fcc_ind <- fcc_df %>%
    dplyr::filter(area_ctry=="India") %>%
    dplyr::select(-pdef, -yrdis) %>%
    dplyr::group_by(area_ctry) %>%
    dplyr::summarise_if(is.numeric, sum) %>%
    dplyr::rename(cont=area_ctry) %>%
    dplyr::mutate(loss21=100*(for2000-for2100)/for2000)
  ## yr75dis for India
  df <- fcc_df %>%
    dplyr::filter(area_ctry=="India") %>%
    dplyr::select(area_name, for2000, for2010, for2020, andef)
  for2000_ind <- sum(df$for2000)
  perc_ind <- vector()
  ## Loop from i=1 for 2021
  for (i in 1:80) {
    ## Deforestation from 2020
    fc_proj_ind <- sum(pmax(0, df$for2020-i*df$andef))
    ## Percentage of forest loss compared with 2000
    perc_ind[i] <- 100*(for2000_ind-fc_proj_ind)/for2000_ind
  }
  yr75dis_ind <- 2000+(which(perc_ind>=75)[1])-1
  ## Add yr75dis_ind to table for India
  fcc_ind <- fcc_ind %>%
    dplyr::mutate(yr75dis=yr75dis_ind) %>%
    dplyr::select(cont, for2040, for2050, for2060, for2080, for2100,
                  loss21, yr75dis)
  
  ## DRC and Indonesia
  fcc_DRC_IDN <- fcc_df %>% 
    dplyr::filter(area_ctry %in% c("DRC", "Indonesia")) %>%
    dplyr::rename(cont=area_ctry) %>%
    dplyr::mutate(loss21=100*(for2000-for2100)/for2000) %>%
    dplyr::mutate(yr75dis=floor(2000+(for2000*0.75)/andef)) %>%
    dplyr::select(cont, for2040, for2050, for2060, for2080, for2100,
                  loss21, yr75dis)
  
  ## All continents
  fc_allcont <- fc_cont %>% 
    dplyr::group_by(year) %>%
    dplyr::summarise_if(is.numeric, sum) %>%
    dplyr::mutate(perc=100*(fc[year==2000]-fc)/fc[year==2000])
  for2000_allcont <- fc_allcont$fc[fc_allcont$year==2000]
  
  yr75dis_allcont <- fc_allcont %>%
    dplyr::filter(perc>=75) %>%
    dplyr::filter(year==min(year)) %>%
    dplyr::mutate(yr75dis=year-1) %>%
    dplyr::select(yr75dis) %>%
    pull()
  
  fc_proj_allcont <- fc_proj_cont %>%
    dplyr::summarise_if(is.numeric, sum) %>%
    dplyr::mutate(cont="All continents") %>%
    dplyr::mutate(loss21=100*(for2000_allcont-for2100)/for2000_allcont) %>%
    dplyr::mutate(yr75dis=yr75dis_allcont) %>%
    dplyr::relocate(cont, .before=for2040)
  
  ## Combine regions
  fc_proj_regions <- fcc_ind %>%
    dplyr::bind_rows(fcc_DRC_IDN) %>%
    dplyr::bind_rows(fcc_bra) %>%
    dplyr::bind_rows(fc_proj_cont) %>%
    dplyr::bind_rows(fc_proj_allcont)
  
  ## Save results
  f_name <- glue("fcc_proj_region_{s}.csv")
  f <- here("Analysis", dataset, f_name)
  write.table(fc_proj_regions, file=f, sep=",", row.names=FALSE)
  ## Copy for manuscript
  f_doc <- here("Manuscript", "Supplementary_Materials", "tables", f_name)
  file.copy(from=f, to=f_doc, overwrite=TRUE)
  f_doc <- here("Manuscript", "Article", "tables", f_name)
  file.copy(from=f, to=f_doc, overwrite=TRUE)
  
}

## ===================
## Parameter estimates
## ===================

## Create table to store results
par_tab <- data.frame(matrix(NA, nrow=nctry, ncol=14))
names(par_tab) <- c("area_cont", "area_ctry", "area_name", "area_code",
                    "int", "pa", "alt",
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
    ## Area info
    area_cont <- as.character(ctry_df$area_cont[ctry_df$iso3==iso])
    area_ctry <- as.character(ctry_df$area_ctry[ctry_df$iso3==iso])
    area_name <- as.character(ctry_df$area_name[ctry_df$iso3==iso])
    area_code <- as.character(ctry_df$area_code[ctry_df$iso3==iso])
    ## Parameter estimates
    f_name <- file.path(dir, iso, "output/summary_hSDM.txt")
    par <- read.table(f_name, skip=4)
    names(par) <- c("Var", "Mean", "Sd", "CI_low", "CI_high")
    ## Fill in the table
    par_tab[i, 1:4] <- cbind(area_cont, area_ctry, area_name, area_code)
    for (r in 1:(nrow(par)-1)) {
        j <- which(long_var_names==par$Var[r])
        par_tab[i, j+4] <- par$Mean[r]
    }
}

## Rearrange study areas per continent
par_tab <- par_tab %>%
    dplyr::mutate(id=ifelse(area_cont=="America", 1, ifelse(area_cont=="Africa", 2, 3))) %>%
    dplyr::arrange(id, area_name) %>%
    dplyr::select(-id)

## Save results
f <- here("Analysis", dataset, "parameter_estimates.csv")
write.table(par_tab, file=f, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "parameter_estimates.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## =============================================
## Parameter estimate weighted mean by continent
## =============================================

## Load parameter estimates
f <- here("Analysis", dataset, "parameter_estimates.csv")
par_tab <- read.table(f, header=TRUE, sep=",")

## Add forest cover in 2010
f <- here("Analysis", dataset, "forest_cover_change_mean.csv")
fcc_tab <- read.table(f, header=TRUE, sep=",")
par_tab <- par_tab %>%
    dplyr::mutate(for2010=fcc_tab$for2010) %>%
    dplyr::relocate(for2010, .after=area_code)

## Replace NA with 0 as NA mean no effect.
par_tab <- par_tab %>%
    replace(., is.na(.), 0)

## Brazil
par_bra <- par_tab %>%
    dplyr::filter(area_ctry=="Brazil") %>%
    dplyr::select(-area_cont, -area_code, -area_name) %>%
    dplyr::mutate(across(.cols=!c(area_ctry, for2010),
                         .fns=~.x*for2010/sum(for2010))) %>%
    dplyr::select(-for2010) %>%
    dplyr::group_by(area_ctry) %>%
    dplyr::summarize(across(.cols=everything(),
                            .fns=function(x){sum(x, na.rm=TRUE)}))
## India
par_ind <- par_tab %>%
    dplyr::filter(area_ctry=="India") %>%
    dplyr::select(-area_cont, -area_code, -area_name) %>%
    dplyr::mutate(across(.cols=!c(area_ctry, for2010),
                         .fns=~.x*for2010/sum(for2010))) %>%
    dplyr::select(-for2010) %>%
    dplyr::group_by(area_ctry) %>%
    dplyr::summarize(across(.cols=everything(),
                            .fns=function(x){sum(x, na.rm=TRUE)}))

## All study areas
par_all <- par_tab %>%
    dplyr::select(-area_cont, -area_code, -area_name) %>%
    dplyr::mutate(across(.cols=!c(area_ctry, for2010),
                         .fns=~.x*for2010/sum(for2010))) %>%
    dplyr::select(-for2010) %>%
    dplyr::mutate(area_ctry="All continents") %>%
    dplyr::group_by(area_ctry) %>%
    dplyr::summarize(across(.cols=everything(),
                            .fns=function(x){sum(x, na.rm=TRUE)}))

## By continent
par_cont <- par_tab %>%
    dplyr::select(-area_ctry, -area_code, -area_name) %>%
    dplyr::group_by(area_cont) %>%
    dplyr::mutate(across(.cols=!c(for2010),
                         .fns=~.x*for2010/sum(for2010))) %>%
    dplyr::select(-for2010) %>%
    dplyr::summarize(across(.cols=everything(),
                            .fns=function(x){sum(x, na.rm=TRUE)})) %>%
    dplyr::mutate(id=ifelse(area_cont=="America", 1, ifelse(area_cont=="Africa", 2, 3))) %>%
    dplyr::arrange(id) %>%
    dplyr::select(-id) %>%
    dplyr::rename(area_ctry=area_cont)

## Combine regions
par_regions <- par_ind %>%
    dplyr::bind_rows(par_bra) %>%
    dplyr::bind_rows(par_cont) %>%
    dplyr::bind_rows(par_all)

## Save results
f <- here("Analysis", dataset, "weighted_param_region.csv")
write.table(par_regions, file=f, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "weighted_param_region.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## =============================
## Mean and SD for each variable
## =============================

## Create table to store results
mean_sd_tab <- data.frame(matrix(NA, nrow=nctry, ncol=18))
names(mean_sd_tab) <- c("area_cont", "area_ctry", "area_name", "area_code",
                    "alt_m", "alt_sd", "slope_m", "slope_sd", 
                    "ddefor_m", "ddefor_sd", "dedge_m", "dedge_sd",
                    "driver_m", "driver_sd", "droad_m", "droad_sd",
                    "dtown_m",  "dtown_sd")

## Loop on countries
for (i in 1:nctry) {
    iso <- iso3[i]
    continent <- as.character(ctry_df$cont_run[ctry_df$iso3==iso])
    dir <- file.path(dir_fdb, dataset, continent)
    ## Area info
    area_cont <- as.character(ctry_df$area_cont[ctry_df$iso3==iso])
    area_ctry <- as.character(ctry_df$area_ctry[ctry_df$iso3==iso])
    area_name <- as.character(ctry_df$area_name[ctry_df$iso3==iso])
    area_code <- as.character(ctry_df$area_code[ctry_df$iso3==iso])
    ## Sample
    f_name <- file.path(dir, iso, "output/sample.txt")
    sample <- read.table(f_name, header=TRUE, sep=",")
    ## Corrections
    if (iso=="BRA-RN") {sample <- sample %>% dplyr::filter(dist_defor <= 10000)}
    if (iso=="VIR") {sample <- sample %>% dplyr::mutate(dist_river=NA)}
    # Mean and SD
    Mean <- apply(sample, 2, mean, na.rm=TRUE)
    SD <- apply(sample, 2, sd, na.rm=TRUE)
    mat <- rbind(Mean, SD)[,c(1,9,2:6)]
    ## Fill in the table
    mean_sd_tab[i, 1:4] <- cbind(area_cont, area_ctry, area_name, area_code)
    mean_sd_tab[i, 5:18] <- c(mat)
}

## Rearrange study areas per continent
mean_sd_tab <- mean_sd_tab %>%
    dplyr::mutate(id=ifelse(area_cont=="America", 1, ifelse(area_cont=="Africa", 2, 3))) %>%
    dplyr::arrange(id, area_name) %>%
    dplyr::select(-id)

## Save results
f <- here("Analysis", dataset, "mean_sd_var.csv")
write.table(mean_sd_tab, file=f, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "mean_sd_var.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## ==========================
## Back-transformed parameter
## ==========================

## Mean and sd
f <- here("Analysis", dataset, "mean_sd_var.csv")
df_mu_sd <- read.table(f, header=TRUE, sep=",")
head(df_mu_sd)

## Parameter estimates
f <- here("Analysis", dataset, "parameter_estimates.csv")
df_par <- read.table(f, header=TRUE, sep=",")
head(df_par)
  
## Back-transformed parameters
df_bt_par <- df_par
intercept <- df_par$int
for (i in 1:7) {
  mu <- df_mu_sd[, 4+(i*2)-1]
  sd <- df_mu_sd[, 4+(i*2)]
  # Slope parameters (x1000, and x100 for slope and dist_edge)
  coeff <- ifelse(i==2, 100, 1000)
  df_bt_par[, 6+i] <- coeff*df_par[, 6+i] / sd
  # Intercept
  par <- ifelse(is.na(df_par[, 6+i]), 0, df_par[, 6+i])
  ratio_mu_sd <- ifelse(is.na(mu/sd), 0, mu/sd)
  intercept <- intercept - par*ratio_mu_sd
}
df_bt_par$int <- intercept

## Save results
f <- here("Analysis", dataset, "backtransformed_parameters.csv")
write.table(df_bt_par, file=f, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "backtransformed_parameters.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## ==============================================================
## Back-transformed parameter estimate weighted mean by continent
## ==============================================================

## Load parameter estimates
f <- here("Analysis", dataset, "backtransformed_parameters.csv")
par_tab <- read.table(f, header=TRUE, sep=",")

## Add forest cover in 2010
f <- here("Analysis", dataset, "forest_cover_change_mean.csv")
fcc_tab <- read.table(f, header=TRUE, sep=",")
par_tab <- par_tab %>%
    dplyr::mutate(for2010=fcc_tab$for2010) %>%
    dplyr::relocate(for2010, .after=area_code)

## Replace NA with 0 as NA mean no effect.
par_tab <- par_tab %>%
    replace(., is.na(.), 0)

## Brazil
par_bra <- par_tab %>%
    dplyr::filter(area_ctry=="Brazil") %>%
    dplyr::select(-area_cont, -area_code, -area_name) %>%
    dplyr::mutate(across(.cols=!c(area_ctry, for2010),
                         .fns=~.x*for2010/sum(for2010))) %>%
    dplyr::select(-for2010) %>%
    dplyr::group_by(area_ctry) %>%
    dplyr::summarize(across(.cols=everything(),
                            .fns=function(x){sum(x, na.rm=TRUE)}))
## India
par_ind <- par_tab %>%
    dplyr::filter(area_ctry=="India") %>%
    dplyr::select(-area_cont, -area_code, -area_name) %>%
    dplyr::mutate(across(.cols=!c(area_ctry, for2010),
                         .fns=~.x*for2010/sum(for2010))) %>%
    dplyr::select(-for2010) %>%
    dplyr::group_by(area_ctry) %>%
    dplyr::summarize(across(.cols=everything(),
                            .fns=function(x){sum(x, na.rm=TRUE)}))

## All study areas
par_all <- par_tab %>%
    dplyr::select(-area_cont, -area_code, -area_name) %>%
    dplyr::mutate(across(.cols=!c(area_ctry, for2010),
                         .fns=~.x*for2010/sum(for2010))) %>%
    dplyr::select(-for2010) %>%
    dplyr::mutate(area_ctry="All continents") %>%
    dplyr::group_by(area_ctry) %>%
    dplyr::summarize(across(.cols=everything(),
                            .fns=function(x){sum(x, na.rm=TRUE)}))

## By continent
par_cont <- par_tab %>%
    dplyr::select(-area_ctry, -area_code, -area_name) %>%
    dplyr::group_by(area_cont) %>%
    dplyr::mutate(across(.cols=!c(for2010),
                         .fns=~.x*for2010/sum(for2010))) %>%
    dplyr::select(-for2010) %>%
    dplyr::summarize(across(.cols=everything(),
                            .fns=function(x){sum(x, na.rm=TRUE)})) %>%
    dplyr::mutate(id=ifelse(area_cont=="America", 1, ifelse(area_cont=="Africa", 2, 3))) %>%
    dplyr::arrange(id) %>%
    dplyr::select(-id) %>%
    dplyr::rename(area_ctry=area_cont)

## Combine regions
par_regions <- par_ind %>%
    dplyr::bind_rows(par_bra) %>%
    dplyr::bind_rows(par_cont) %>%
    dplyr::bind_rows(par_all)

## Save results
f <- here("Analysis", dataset, "backtransformed_weighted_param_region.csv")
write.table(par_regions, file=f, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables",
              "backtransformed_weighted_param_region.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## ==================================================
## Deforestation probability -- variable relationship
## ==================================================

##-- Distance to road and forest edge --##

## Load data
f <- here("Analysis", dataset, "data_allctry_prop.csv")
data <- read_delim(f, delim=",")

## Percentiles
perc <- seq(0, 100, by=10)
nperc <- length(perc)

## Result table with local means
theta_lmean <- list()

## Compute theta and se by bins
y <- 1-data$fcc23  # Transform: defor=1, forest=0
data$dist_road_km <- data$dist_road/1000
data$dist_edge_km <- data$dist_edge/1000
varname <- c("dist_road_km", "dist_edge_km")
for (i in 1:length(varname)) {
  v <- varname[i]
  theta <- rep(0, nperc-1)
  se <- rep(0, nperc-1)
  x <- rep(0, nperc-1)
  quantiles <- quantile(data[[v]], probs=perc/100, na.rm=TRUE)
  ## Model icar
  theta_icar <- data$icar
  theta_icar_mean <- rep(0, nperc-1)
  ## Loop on percentiles
  for (j in 1:(nperc-1)) {
    inf <- quantiles[j]
    sup <- quantiles[j + 1]
    x[j] <- inf + (sup - inf) / 2
    # Observations in bin
    w <- (data[v] >= inf) & (data[v] < sup)
    if (j==(nperc-1)) {w <- (data[v] >= inf) & (data[v] <= sup)}
    y_bin <- y[w]
    # Local mean and se
    s <- sum(y_bin)  # success
    n <- length(y_bin)  # trials
    theta[j] <- ifelse(n != 0, s/n, NaN)
    ph <- (s + 1 / 2) / (n + 1)
    se[j] <- sqrt(ph * (1 - ph) / (n + 1))
    # icar
    t_bin <- theta_icar[w]
    theta_icar_mean[j] <- mean(t_bin)
  }
  ## Fill the list
  theta_lmean[[i]] <- data.frame(x=x, theta_obs=theta, theta_icar=theta_icar_mean)
}

## Combine data-set
theta_road <- theta_lmean[[1]]
theta_edge <- theta_lmean[[2]]
theta_dist <- bind_rows(data.frame(var="road", theta_road),
                       data.frame(var="forest edge", theta_edge))
## mytheme
mytheme <- theme(
  axis.title=element_text(size=12),
  axis.text=element_text(size=9),
  legend.title=element_text(size=12),
  legend.text=element_text(size=10),
  legend.position=c(0.9, 0.9),
  legend.justification=c(1, 1))

## Plot
p_dist <- ggplot(aes(x=x, y=theta_obs, group=var, col=var), data=theta_dist) +
  geom_point(size=1.2) +
  geom_line(aes(x=x, y=theta_icar, group=var, col=var), size=0.7) +
  xlab("Distance (km)") +
  ylab("Spatial probability of deforestation") +
  scale_color_manual(values=wes_palette("Moonrise2")[c(2, 1)],
                     name="Distance to:",
                     breaks=c("road", "forest edge"),
                     labels=c("road", "forest edge")) +
  coord_cartesian(xlim=c(0,22), ylim=c(0, 1)) +
  scale_y_continuous(breaks=seq(0, 1, 0.25)) +
  scale_x_continuous(breaks=c(0, 1, 5, 10, 15, 20, 22)) +
  theme_bw() + mytheme + theme(axis.title.y=element_blank())

##-- Protected areas --##

## Table of results
theta_pa <- data.frame(pa=c("Unprotected","Protected"), obs=NA, icar=NA)

## No PA
w0 <- which(data["pa"]==0)
nw0 <- length(w0)
theta_pa$obs[1] <- sum(y[w0]==1)/nw0
theta_pa$icar[1] <- mean(data$icar[w0])

## PA
w1 <- which(data["pa"]==1)
nw1 <- length(w1)
theta_pa$obs[2] <- sum(y[w1]==1)/nw1
theta_pa$icar[2] <- mean(data$icar[w1])

## PLot
p_pa <- ggplot(aes(x=pa, y=icar), data=theta_pa) +
  geom_col(colour="#bcc6c3", fill="#bcc6c3", size=1, width=0.25) +
  geom_point(aes(x=pa, y=obs), data=theta_pa, colour=wes_palette("Moonrise2")[1]) +
  xlab("Protection status") +
  ylab("Spatial probability of deforestation") +
  coord_cartesian(ylim=c(0, 1)) +
  scale_x_discrete(limits=c("Unprotected", "Protected")) +
  scale_y_continuous(breaks=seq(0, 1, 0.25)) +
  theme_bw() + mytheme
  
##-- Arrange and save plots --##
library(gridExtra)
textwidth <- 16.6
f <- here("Analysis", dataset, "proba-var.png")
png(filename=f, width=textwidth, height=textwidth*0.6, units="cm", res=300)
grid.arrange(p_pa, p_dist, ncol=2, nrow=1, widths=c(1.5, 3), heights=c(1))
dev.off()
## Copy for manuscript
f_doc <- here("Manuscript", "Article", "figures", "proba-var.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## ==================================
## Correlation plot between variables
## ==================================

## Logistic regression for correlation between continuous variables and protected areas.
var_names <- c("altitude", "slope", "dist_defor", "dist_edge",
               "dist_river", "dist_road", "dist_town")
corr_pa <- vector()
for (i in 1:7) {
  formula <- paste0("pa~scale(",var_names[i],")")
  mod_corr_pa <- glm(formula, data=data)
  corr_pa[i] <- round(coefficients(mod_corr_pa)[2], 2)
}

require("ggcorrplot")
## Correlation matrix
corr <- round(cor(data %>% dplyr::select(pa, altitude, slope, dist_defor:dist_town), use="pairwise.complete.obs"), 2)
corr[1, ] <- c(1, corr_pa)
corr[, 1] <- c(1, corr_pa)
## Variable names
colnames(corr) <- c("pa", "elev", "slope", "ddefor", "dedge", "driver", "droad", "dtown")
row.names(corr) <- colnames(corr)
## Palette
pal <- wes_palette("Moonrise2")
## Plot
f <- here("Analysis", dataset, "corr-var.png")
png(filename=f, width=textwidth, height=0.6*textwidth, units="cm", res=300)
ggcorrplot(corr, hc.order = TRUE, type = "lower",
   outline.col = "white", lab=TRUE, legend.title="Correlation",
   ggtheme = ggplot2::theme_bw, digits=2,
   colors = c(pal[1], "white", pal[2]))
dev.off()
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "figures", "corr-var.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## ===================
## PA effect
## ===================

## Create table to store results
parea_tab <- data.frame(matrix(NA, nrow=nctry, ncol=8))
names(parea_tab) <- c("area_cont", "area_ctry", "area_name", "area_code",
                      "Mean", "Sd", "CI_low", "CI_high")

## Loop on countries
for (i in 1:nctry) {
    iso <- iso3[i]
    continent <- as.character(ctry_df$cont_run[ctry_df$iso3==iso])
    dir <- file.path(dir_fdb, dataset, continent)
    ## Area info
    area_cont <- as.character(ctry_df$area_cont[ctry_df$iso3==iso])
    area_ctry <- as.character(ctry_df$area_ctry[ctry_df$iso3==iso])
    area_name <- as.character(ctry_df$area_name[ctry_df$iso3==iso])
    area_code <- as.character(ctry_df$area_code[ctry_df$iso3==iso])
    ## Parameter estimates
    f_name <- file.path(dir, iso, "/output/summary_hSDM.txt")
    par <- read.table(f_name, skip=4)
    names(par) <- c("Var", "Mean", "Sd", "CI_low", "CI_high")
    ## Fill in the table
    parea_tab[i, 1:4] <- cbind(area_cont, area_ctry, area_name, area_code)
    j <- which(par$Var=="C(pa)[T.1.0]")
    if (length(j)>0) {
        parea_tab[i, 5:8] <- par[j, 2:5]
    }
}

## Rearrange study areas per continent
parea_tab <- parea_tab %>%
    dplyr::mutate(id=ifelse(area_cont=="America", 1, ifelse(area_cont=="Africa", 2, 3))) %>%
    dplyr::arrange(id, area_name) %>%
    dplyr::select(-id)

## Add forest cover in 2100
f <- here("Analysis", dataset, "forest_cover_change_mean.csv")
fcc_tab <- read.table(f, header=TRUE, sep=",")
parea_tab <- parea_tab %>%
    dplyr::mutate(for2010=fcc_tab$for2010) %>%
    dplyr::relocate(for2010, .after=area_code)

## Save results
f <- here("Analysis", dataset, "parea_estimates.csv")
write.table(parea_tab, file=f, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "parea_estimates.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## ===================
## Road effect
## ===================

## Create table to store results
road_tab <- data.frame(matrix(NA, nrow=nctry, ncol=8))
names(road_tab) <- c("area_cont", "area_ctry", "area_name", "area_code",
                     "Mean", "Sd", "CI_low", "CI_high")

## Loop on countries
for (i in 1:nctry) {
    iso <- iso3[i]
    continent <- as.character(ctry_df$cont_run[ctry_df$iso3==iso])
    dir <- file.path(dir_fdb, dataset, continent)
    ## Area info
    area_cont <- as.character(ctry_df$area_cont[ctry_df$iso3==iso])
    area_ctry <- as.character(ctry_df$area_ctry[ctry_df$iso3==iso])
    area_name <- as.character(ctry_df$area_name[ctry_df$iso3==iso])
    area_code <- as.character(ctry_df$area_code[ctry_df$iso3==iso])
    ## Parameter estimates
    f_name <- file.path(dir, iso, "/output/summary_hSDM.txt")
    par <- read.table(f_name, skip=4)
    names(par) <- c("Var", "Mean", "Sd", "CI_low", "CI_high")
    ## Fill in the table
    road_tab[i, 1:4] <- cbind(area_cont, area_ctry, area_name, area_code)
    j <- which(par$Var=="scale(dist_road)")
    if (length(j)>0) {
        road_tab[i, 5:8] <- par[j, 2:5]
    }
}

## Rearrange study areas per continent
road_tab <- road_tab %>%
    dplyr::mutate(id=ifelse(area_cont=="America", 1, ifelse(area_cont=="Africa", 2, 3))) %>%
    dplyr::arrange(id, area_name) %>%
    dplyr::select(-id)

## Add forest cover in 2100
f <- here("Analysis", dataset, "forest_cover_change_mean.csv")
fcc_tab <- read.table(f, header=TRUE, sep=",")
road_tab <- road_tab %>%
    dplyr::mutate(for2010=fcc_tab$for2010) %>%
    dplyr::relocate(for2010, .after=area_code)

## Save results
f <- here("Analysis", dataset, "road_estimates.csv")
write.table(road_tab, file=f, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "road_estimates.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## =========================
## Significance PA and roads
## =========================

## Load data-sets
if (!exists("parea_tab")) {
  f <- here("Analysis", dataset,"results", "parea_estimates.csv")
  parea_tab <- read.table(f, header=TRUE, sep=",")
}
if (!exists("road_tab")) {
  f <- here("Analysis", dataset,"results", "road_estimates.csv")
  road_tab <- read.table(f, header=TRUE, sep=",")
}

## Significance PA
parea_tab <- parea_tab %>%
    mutate(sign=ifelse(CI_low * CI_high > 0 & !is.na(Mean), 1, 0))
## Number of countries for which the effet of protected areas is significant
sum_sign_PA <- sum(parea_tab$sign==1)
## Percentage of countries for which the effet of protected areas is significant
perc_sign_PA <- 100*sum_sign_PA/nrow(parea_tab)
## Weighted percentage with forest size in 2010
f <- here("Analysis", dataset, "forest_cover_change_mean.csv")
fcc_tab <- read.table(f, header=TRUE, sep=",")
weights <- fcc_tab$for2010
perc_sign_w_PA <- 100*sum((parea_tab$sign==1)*weights)/sum(weights)

## Significance road
road_tab <- road_tab %>%
    mutate(sign=ifelse(CI_low * CI_high > 0 & !is.na(Mean), 1, 0))
## Number of countries for which the effet of roads is significant
sum_sign_road <- sum(road_tab$sign==1)
## Percentage of countries for which the effet of roads is significant
perc_sign_road <- 100*sum_sign_road/nrow(road_tab)
## Weighted percentage with forest size in 2010
f <- here("Analysis", dataset, "forest_cover_change_mean.csv")
fcc_tab <- read.table(f, header=TRUE, sep=",")
weights <- fcc_tab$for2010
perc_sign_w_road <- 100*sum((road_tab$sign==1)*weights)/sum(weights)

## Save results
nctry <- length(fcc_tab$area_code)
nctry_sign <- c(sum_sign_PA, sum_sign_road)
perc <- c(perc_sign_PA, perc_sign_road)
perc_w <- c(perc_sign_w_PA, perc_sign_w_road) 
sign_PA_road <- data.frame(var=c("PA","road"), nctry=nctry, nctry_sign=nctry_sign,
                           perc=round(perc), perc_w=round(perc_w))

## Save results
f <- here("Analysis", dataset, "sign_PA_road.csv")
write.table(sign_PA_road, file=f, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "sign_PA_road.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## =====================================
## Carbon emissions (in tonnes=10e6 g C)
## =====================================

## Simulations with uncertainty
sim <- c("mean", "min", "max")
nsim <- length(sim)

## Loop on simulations
for (j in 1:nsim) {
  
  ## Simulation id
  s <- sim[j]

  ## Create table to store results
  Cem_tab <- data.frame(matrix(NA, nrow=nctry, ncol=16), stringsAsFactors=FALSE)
  C_var <- paste0("C", 2020 + c(0, 10, 15, 20, 30, 35, 40, 50, 60, 65, 70, 80))
  names(Cem_tab) <- c("area_cont", "area_ctry", "area_name", "area_code", C_var)
  
  ## Loop on countries
  for (i in 1:nctry) {
    iso <- iso3[i]
    continent <- as.character(ctry_df$cont_run[ctry_df$iso3==iso])
    dir <- file.path(dir_fdb, dataset, continent)
    ## Area info
    area_cont <- as.character(ctry_df$area_cont[ctry_df$iso3==iso])
    area_ctry <- as.character(ctry_df$area_ctry[ctry_df$iso3==iso])
    area_name <- as.character(ctry_df$area_name[ctry_df$iso3==iso])
    area_code <- as.character(ctry_df$area_code[ctry_df$iso3==iso])
    ## Carbon emissions
    f_name <- file.path(dir, iso, glue("output/{s}/C_emissions.csv"))
    Cem_df <- read.table(f_name, header=TRUE, sep=",", stringsAsFactors=FALSE)
    ## Fill in the table
    Cem_tab[i, 1:4] <- cbind(area_cont, area_ctry, area_name, area_code)
    Cem_tab[i, 5:16] <- Cem_df$C
  }
  
  ## Rearrange study areas per continent
  Cem_tab <- Cem_tab %>%
    dplyr::mutate(id=ifelse(area_cont=="America", 1, ifelse(area_cont=="Africa", 2, 3))) %>%
    dplyr::arrange(id, area_name) %>%
    dplyr::select(-id)
  
  ## Save results
  f <- here("Analysis", dataset, glue("C_emissions_{s}.csv"))
  write.table(Cem_tab, file=f, sep=",", row.names=FALSE)
  ## Copy for manuscript
  f_doc <- here("Manuscript", "Supplementary_Materials", "tables", glue("C_emissions_{s}.csv"))
  file.copy(from=f, to=f_doc, overwrite=TRUE)
  ## Copy for data
  f_doc <- here("Manuscript", "Supplementary_Data", "tables", glue("carbon_emissions_{s}.csv"))
  file.copy(from=f, to=f_doc, overwrite=TRUE)

}

## ========================
## Summary carbon emissions
## ========================

## Simulations with uncertainty
sim <- c("mean", "min", "max")
nsim <- length(sim)

## Loop on simulations
for (j in 1:nsim) {
  
  ## Simulation id
  s <- sim[j]

  ## Load data
  f <- here("Analysis", dataset, glue("C_emissions_{s}.csv"))
  Cem_tab <- read.table(f, header=TRUE, sep=",")
  
  ## Summarize results
  Cem_tab2 <- Cem_tab %>%
    group_by(area_cont) %>%
    summarize(across(starts_with("C"), sum)) %>%
    bind_rows(data.frame(area_cont="All continents",
                         summarize(Cem_tab, across(starts_with("C"), sum)),
                         stringsAsFactors=FALSE)) %>%
    mutate(across(C2020:C2100, function(x){x*1e-9}))  # Results in PgC
  
  ## Save results
  f_name <- glue("C_emissions_summary_{s}.csv")
  f <- here("Analysis", dataset, f_name)
  write.table(Cem_tab2, file=f, sep=",", row.names=FALSE)
  ## Copy for manuscript
  f_doc <- here("Manuscript", "Supplementary_Materials", "tables", f_name)
  file.copy(from=f, to=f_doc, overwrite=TRUE)

}

## In 2020-2050, 18.4 billions of tonnes of C emmitted = 18.1 PgC (0.60 PgC/year on 2020-2050)
##
## References for comparison:
## Baccini et al. (2017) only 0.8 PgC/year for 2003–2014. 
## Land use and land-cover change (LULCC) are believed to release between 0.81 and 1.14 PgC/yr.
## Baccini 2012 et al. reported 1 PgC/year on the period 2000-2010.
## See Van Der Werf, G. R. et al. CO2 emissions from forest loss. Nature Geosci. 2, 737–738 (2009).
## Friedlingstein, P. et al. Update on CO2 emissions. Nature Geosci. 3, 811–812 (2010).
## Le Quéré, C. et al. Trends in the sources and sinks of carbon dioxide. Nature Geosci. 2, 831–836 (2009).

## ======================
## Carbon emission trends
## ======================

## Simulations with uncertainty
sim <- c("mean", "min", "max")
nsim <- length(sim)

## Loop on simulations
for (j in 1:nsim) {
  
  ## Simulation id
  s <- sim[j]
  
  ## Load data
  f <- here("Analysis", dataset, glue("C_emissions_summary_{s}.csv"))
  Cem_tab2 <- read.table(f, header=TRUE, sep=",")
  
  ## Simulation id
  s <- sim[j]

  ## Compute carbon emission trends in the future
  C_trend <- Cem_tab2 %>%
    dplyr::mutate(T10_20=C2020/10,
                  T20_30=C2030/10,
                  T30_40=(C2040-C2030)/10,
                  T40_50=(C2050-C2040)/10,
                  T50_60=(C2060-C2050)/10,
                  T60_70=(C2070-C2060)/10,
                  T70_80=(C2080-C2070)/10,
                  T80_90=(C2090-C2080)/10,
                  T90_100=(C2100-C2090)/10) %>%
    dplyr::select(area_cont, T10_20:T90_100)
  
  ## Deforestation => increase in C source with time (deforestation of forest with higher carbon stocks).
  ## Climate change => decrease in C sink with time (higher mortality), see Hubau2020.
  ## The result is that forests will likely become a major C source in the future. 
  
  ## Save results
  f_name <- glue("C_trend_{s}.csv")
  f <- here("Analysis", dataset, f_name)
  write.table(C_trend, file=f, sep=",", row.names=FALSE)
  ## Copy for manuscript
  f_doc <- here("Manuscript", "Supplementary_Materials", "tables", f_name)
  file.copy(from=f, to=f_doc, overwrite=TRUE)

}

## ==============================
## Figure: carbon emission trends
## ==============================

## Load data
# Mean
f <- here("Analysis", dataset, glue("C_trend_mean.csv"))
C_trend_mean <- read.table(f, header=TRUE, sep=",")
## Transform dataset in long format
C_long_mean <- C_trend_mean %>%
  tidyr::pivot_longer(cols=T10_20:T90_100, names_to="T_int",
                      values_to="C_em") %>%
  dplyr::mutate(year=rep(seq(2015, 2095, by=10), 4))
# Min
f <- here("Analysis", dataset, glue("C_trend_min.csv"))
C_trend_min <- read.table(f, header=TRUE, sep=",")
## Transform dataset in long format
C_long_min <- C_trend_min %>%
  tidyr::pivot_longer(cols=T10_20:T90_100, names_to="T_int",
                      values_to="C_em") %>%
  dplyr::mutate(year=rep(seq(2015, 2095, by=10), 4))
# Max
f <- here("Analysis", dataset, glue("C_trend_max.csv"))
C_trend_max <- read.table(f, header=TRUE, sep=",")
## Transform dataset in long format
C_long_max <- C_trend_max %>%
  tidyr::pivot_longer(cols=T10_20:T90_100, names_to="T_int",
                      values_to="C_em") %>%
  dplyr::mutate(year=rep(seq(2015, 2095, by=10), 4))

## Historical data
C_hist <- C_long_mean %>% dplyr::filter(year==2015)

## Projected data
C_proj <- C_long_mean %>%
  mutate(C_em_min=C_long_min[["C_em"]],
         C_em_max=C_long_max[["C_em"]])

## mytheme
mytheme <- theme(
  axis.title=element_text(size=12),
  axis.text=element_text(size=10),
  legend.title=element_text(size=12),
  legend.text=element_text(size=10),
  legend.position=c(0.01, 0.99),
  legend.justification=c(0, 1),
  legend.background=element_rect(fill="transparent"))

## Plot
p <- ggplot(aes(x=year, y=C_em, group=area_cont, col=area_cont, fill=area_cont), data=C_proj) +
  geom_ribbon(aes(ymin=C_em_min, ymax=C_em_max, group=area_cont, fill=area_cont),
              alpha=0.2, data=C_proj, linetype=0) + 
  geom_line(aes(group=area_cont, col=area_cont), data=C_proj, size=0.8) +
  geom_point(aes(group=area_cont, col=area_cont), data=C_hist, size=1) +
  xlab("Year") + ylab("Annual carbon emissions (Pg/yr) associated to \ndeforestation of moist tropical forests") +
  scale_color_manual(values=wes_palette("Moonrise2")[c(4, 3, 2, 1)],
                     name="Continents",
                     breaks=c("All continents", "America", "Africa", "Asia"),
                     labels=c("All continents", "America", "Africa", "Asia")) +
  scale_fill_manual(values=wes_palette("Moonrise2")[c(4, 3, 2, 1)],
                    name="Continents",
                    breaks=c("All continents", "America", "Africa", "Asia"),
                    labels=c("All continents", "America", "Africa", "Asia")) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.25)) +
  scale_x_continuous(limits=c(2010, 2100), breaks=seq(2010,2100,by=10)) +
  theme_bw() + mytheme

## Save results
f <- here("Analysis", dataset, "C_trend.png")
ggsave(f, p, width=16.6, height=10, units="cm", dpi=300)
## Copy for manuscript
f_doc <- here("Manuscript", "Article", "figures",
              "C_trend.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## =======================================
## Table: species in biodiversity hotspots
## =======================================

f <- here("Analysis", "data", "species_biodiversity_hotspots.csv")
df <- read_delim(f, delim=",")
df2 <- df %>%
  dplyr::filter(Hotspot %in% c("Mesoamerica", "Guinean Forests of West Africa",
                              "Horn of Africa", "Madagascar and Indian Ocean Islands", "Indo-Burma",
                              "Western Ghats and Sri Lanka")) %>%
  dplyr::summarize(Plants=sum(Plants_E), Vertebrates=sum(Birds_E, Reptiles_E, Amphibians_E, Freshwater_fishes_E, Mammals_E))

## Save results
f <- here("Analysis", dataset, "species_loss.csv")
write_delim(df2, f, delim=",")
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "species_loss.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)
## Copy also the dataset
f <- here("Analysis", "data", "species_biodiversity_hotspots.csv")
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "species_biodiversity_hotspots.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## =======================================
## Uncertainty
## =======================================

## ---------------------------------------
## Forest cover change
## ---------------------------------------

# Data
f_mean <- here("Analysis", dataset, "forest_cover_change_mean.csv")
f_min <- here("Analysis", dataset, "forest_cover_change_min.csv")
f_max <- here("Analysis", dataset, "forest_cover_change_max.csv")
df_mean <- read_delim(f_mean, delim=",") %>% mutate(proj="mean")
df_min <- read_delim(f_min, delim=",") %>% mutate(proj="low")
df_max <- read_delim(f_max, delim=",") %>% mutate(proj="high")

# Combine datasets
df <- df_min %>%
  bind_rows(df_mean) %>%
  bind_rows(df_max) %>%
  # Id
  mutate(id_cont=ifelse(area_cont=="America", 1, ifelse(area_cont=="Africa", 2, 3))) %>%
  mutate(id_d=ifelse(proj=="low", 1, ifelse(proj=="mean", 2, 3))) %>%
  # Arrange
  arrange(id_cont, area_name, id_d) %>%
  relocate(proj, .after=area_code)

## Save results
f <- here("Analysis", dataset, "forest_cover_change_ci.csv")
write_delim(df, f, delim=",")
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Data", "tables", "forest_cover_change_ci.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## ---------------------------------------
## Carbon emissions
## ---------------------------------------

# Data
f_mean <- here("Analysis", dataset, "C_emissions_mean.csv")
f_min <- here("Analysis", dataset, "C_emissions_min.csv")
f_max <- here("Analysis", dataset, "C_emissions_max.csv")
df_mean <- read_delim(f_mean, delim=",") %>% mutate(proj="mean")
df_min <- read_delim(f_min, delim=",") %>% mutate(proj="low")
df_max <- read_delim(f_max, delim=",") %>% mutate(proj="high")

# Combine datasets
df <- df_min %>%
  bind_rows(df_mean) %>%
  bind_rows(df_max) %>%
  # Id
  mutate(id_cont=ifelse(area_cont=="America", 1, ifelse(area_cont=="Africa", 2, 3))) %>%
  mutate(id_d=ifelse(proj=="low", 1, ifelse(proj=="mean", 2, 3))) %>%
  # Arrange
  arrange(id_cont, area_name, id_d) %>%
  relocate(proj, .after=area_code)

## Save results
f <- here("Analysis", dataset, "C_emissions_ci.csv")
write_delim(df, f, delim=",")
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Data", "tables", "C_emissions_ci.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## ==================================
## End Of File
## ==================================