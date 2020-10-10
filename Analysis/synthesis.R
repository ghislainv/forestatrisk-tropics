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

## Set working directory
setwd(here())

## Dataset
##dataset <- "gfc2020_70" 
dataset <- "jrc2020"
dir.create(here("Analysis", dataset, "results"), recursive=TRUE)

## Result directory
dir_fdb <- "/home/forestatrisk-tropics"

## =================
## Countries info
## =================

## Load country info (encoding pb for ctry names)
ctry_df <- read.csv2(here("Analysis","ctry_run.csv"), header=TRUE, sep=";", encoding="UTF-8")

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
## write.table(sample_allctry_tab, file=here("Analysis", dataset, "results/sample_allctry.csv"), sep=",", row.names=FALSE)

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
fcc_tab <- read.table(here("Analysis", dataset, "results", "forest_cover_change.csv"),
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
f1 <- here("Analysis", dataset, "results", "performance_index.csv")
f2 <- here("Analysis", dataset, "results", "perf_mod.csv")
f3 <- here("Analysis", dataset, "results", "perf_cont_mod.csv")
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
f <- here("Analysis", dataset, "results", "samp_size.csv")
write.table(samp_size_tab2, file=f, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "samp_size.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## ===================
## Forest cover change
## ===================

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
    # Year during which forest will have disappeared
    mutate(yrdis=floor(2020 + for2020/andef))

## Corrections for Brazil with deforestation diffusion
if (dataset=="gfc2020_70") {
    fname_BRA <- "Brazil/fcc_BRA_gfc.csv"
} else if (dataset=="jrc2020") {
    fname_BRA <- "Brazil/fcc_BRA_jrc.csv"
}
fcc_BRA <- read.table(file.path(dir_fdb, dataset, fname_BRA), sep=",",
                      header=TRUE, stringsAsFactors=FALSE)
## Check order and replace with correct values for Brazil
codes_fcc <- paste0("BRA-",fcc_tab2$area_code[fcc_tab2$area_ctry=="Brazil"])
codes_BRA <- fcc_BRA$iso3
if (all(codes_fcc==codes_BRA)) {
    fcc_tab2[fcc_tab2$area_ctry=="Brazil", c(12:ncol(fcc_tab2))] <- round(fcc_BRA[, c(seq(7, 27, by=2), 30)])
}

## Sort continents, and select col
fcc_tab3 <- fcc_tab2 %>%
    # Id
    mutate(id=ifelse(area_cont=="America", 1, ifelse(area_cont=="Africa", 2, 3))) %>%
    arrange(id, area_name) %>%
    # Select columns
    dplyr::select(area_cont, area_ctry, area_name, area_code, for2000:yrdis)

## Save results
f <- here("Analysis", dataset, "results", "forest_cover_change.csv")
write.table(fcc_tab3, file=f, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "forest_cover_change.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## ====================================================
## Forest cover change summarized per region
## ====================================================

## Load previous fcc table
f <- here("Analysis", dataset, "results", "forest_cover_change.csv")
fcc_df <- read.table(f, header=TRUE, sep=",")

## For each continent
fcc_cont <- fcc_df %>%
    dplyr::group_by(area_cont) %>%
    dplyr::summarise_if(is.numeric, list(sum=sum, max=max)) %>%
    dplyr::select(area_cont, for2000_sum:andef_sum, for2030_sum:for2100_sum, yrdis_max) %>%
    dplyr::mutate(id_cont=c(2, 1, 3)) %>%
    dplyr::arrange(id_cont) %>%
    dplyr::select(-id_cont)

## For all continents
fcc_all <- fcc_df %>%
    dplyr::summarise_if(is.numeric, list(sum=sum, max=max)) %>%
    dplyr::select(for2000_sum:andef_sum, for2030_sum:for2100_sum, yrdis_max) %>%
    dplyr::mutate(area_cont="All continents") %>%
    dplyr::relocate(area_cont, .before=for2000_sum)

## For Brazil
fcc_bra <- fcc_df %>%
    dplyr::filter(area_ctry=="Brazil") %>%
    dplyr::group_by(area_ctry) %>%
    dplyr::summarise_if(is.numeric, list(sum=sum, max=max)) %>%
    dplyr::select(area_ctry, for2000_sum:andef_sum, for2030_sum:for2100_sum, yrdis_max) %>%
    dplyr::rename(area_cont=area_ctry)

## For India
fcc_ind <- fcc_df %>%
    dplyr::filter(area_ctry=="India") %>%
    dplyr::group_by(area_ctry) %>%
    dplyr::summarise_if(is.numeric, list(sum=sum, max=max)) %>%
    dplyr::select(area_ctry, for2000_sum:andef_sum, for2030_sum:for2100_sum, yrdis_max) %>%
    dplyr::rename(area_cont=area_ctry)

## Combine
TI <- 2020-2010  ## Time-interval
fcc_comb <- fcc_ind %>%
    # Add Brazil
    rbind(fcc_bra) %>%
    # Add continents
    rbind(fcc_cont) %>%
    rbind(fcc_all) %>%
    # Rename
    dplyr::rename_at(.vars=vars(starts_with("for")), .funs=substr, start=1, stop=7) %>%
    dplyr::rename(andef=andef_sum, yrdis=yrdis_max) %>%
    # Compute pdef
    dplyr::mutate(pdef=round(100*(1-(1-(for2010-for2020)/for2010)^(1/TI)), 1)) %>%
    # Compute loss21
    dplyr::mutate(loss21=100*(for2000-for2100)/for2000) %>%
    # Arrange columns
    dplyr::relocate(pdef, .after=andef) %>%
    dplyr::relocate(loss21, .before=yrdis)

## Save results
f <- here("Analysis", dataset, "results", "fcc_hist_region.csv")
write.table(fcc_comb, file=f, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "fcc_hist_region.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## ===============================================
## Projecting percentage of forest loss per region
## ===============================================

## Load previous fcc table
f <- here("Analysis", dataset, "results", "forest_cover_change.csv")
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
f <- here("Analysis", dataset, "results", "perc_loss_cont.csv")
write.table(fc_perc_cont, file=f, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "perc_loss_cont.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## ===================================
## Plot change in percentage with time
## ====================================

## Load data
fc_perc_cont <- read.table(f, header=TRUE, sep=",")
df_hist <- fc_perc_cont %>%
    dplyr::filter(year %in% c(2000, 2010, 2020))
df_proj <- fc_perc_cont %>%
    dplyr::filter(!(year %in% c(2000, 2010, 2020)))

## mytheme
mytheme <- theme(
    axis.text=element_text(size=20),
    axis.title=element_text(size=20),
    legend.title=element_text(size=20),
    legend.text=element_text(size=20),
    legend.position=c(0.95, 0.05),
    legend.justification=c(1, 0))

## Plot
p <- ggplot(aes(x=year, y=perc, group=cont, col=cont), data=df_hist) +
    geom_point(size=2) +
    geom_line(data=df_proj, size=1.2) +
    xlab("Year") + ylab("Percentage of forest cover loss\n(in comparison with year 2000)") +
    scale_color_discrete(name="Continents",
                         breaks=c("America", "Africa", "Asia"),
                         labels=c("America", "Africa", "Asia")) +
    ylim(0,100) +
    geom_hline(yintercept=75) +
    theme_bw() + mytheme

## Save results
f <- here("Analysis", dataset, "results", "perc_loss_cont.png")
ggsave(f, p)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "figures",
              "perc_loss_cont.png")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## =========================================================================
## Forest cover projections including yr75dis per region
## yrdis75: year during which 75% of for2000 wil have disappeared per region
## =========================================================================

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
    dplyr::filter(year %in% c(2040, 2060, 2080, 2100)) %>%
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
    dplyr::select(cont, for2040, for2060, for2080, for2100,
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
    dplyr::select(cont, for2040, for2060, for2080, for2100,
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
    rbind(fcc_bra) %>%
    rbind(fc_proj_cont) %>%
    rbind(fc_proj_allcont)

## Save results
f <- here("Analysis", dataset, "results", "fcc_proj_region.csv")
write.table(fc_proj_regions, file=f, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "fcc_proj_region.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## ===================
## Model parameters
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
    f_name <- file.path(dir, iso, "/output/summary_hSDM.txt")
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
f <- here("Analysis", dataset, "results", "parameter_estimates.csv")
write.table(par_tab, file=f, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "parameter_estimates.csv")
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

## Save results
f <- here("Analysis", dataset,"results", "parea_estimates.csv")
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

## Save results
f <- here("Analysis", dataset, "results", "road_estimates.csv")
write.table(road_tab, file=f, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "road_estimates.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## =========================
## Significance PA and roads
## =========================

## Significance PA
parea_tab <- parea_tab %>%
    mutate(sign=ifelse(CI_low * CI_high > 0 & !is.na(Mean), 1, 0))
## Percentage of country for which the effet of protected areas is significant
perc_sign_PA <- 100*sum(parea_tab$sign==1)/nrow(parea_tab)
## Weighted percentage with forest size in 2010
f <- here("Analysis", dataset, "results", "forest_cover_change.csv")
fcc_tab <- read.table(f, header=TRUE, sep=",")
weights <- fcc_tab$for2010
perc_sign_w_PA <- 100*sum((parea_tab$sign==1)*weights)/sum(weights)

## Significance road
road_tab <- road_tab %>%
    mutate(sign=ifelse(CI_low * CI_high > 0 & !is.na(Mean), 1, 0))
## Percentage of country for which the effet of protected areas is significant
perc_sign_road <- 100*sum(road_tab$sign==1)/nrow(road_tab)
## Weighted percentage with forest size in 2010
f <- here("Analysis", dataset, "results", "forest_cover_change.csv")
fcc_tab <- read.table(f, header=TRUE, sep=",")
weights <- fcc_tab$for2010
perc_sign_w_road <- 100*sum((road_tab$sign==1)*weights)/sum(weights)

## Save results
nctry <- length(fcc_tab$iso3)
perc <- c(perc_sign_PA, perc_sign_road)
perc_w <- c(perc_sign_w_PA, perc_sign_w_road) 
sign_PA_road <- data.frame(var=c("PA","road"), nctry=nctry, perc=round(perc), perc_w=round(perc_w))

## Save results
f <- here("Analysis", dataset, "results", "sign_PA_road.csv")
write.table(sign_PA_road, file=, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "sign_PA_road.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

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
f <- here("Analysis", dataset, "results", "C_emissions.csv")
write.table(Cem_tab, file=f, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "C_emissions.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

## ========================
## Summary carbon emissions
## ========================

## Summarize results
Cem_tab2 <- Cem_tab %>%
    mutate(cont=ifelse(cont=="Brazil", "America", cont)) %>%
    group_by(cont) %>%
    summarize_at(all_of(C_var), sum) %>%
    bind_rows(data.frame(cont="TOTAL",
                         summarize_at(Cem_tab, vars(C_var), sum),
                         stringsAsFactors=FALSE)) %>%
    mutate_at(vars(C2020:C2100), function(x){x*1e-9})  # Results in PgC

## Save results
f <- here("Analysis", dataset, "results", "C_emissions_summary.csv")
write.table(Cem_tab2, file=f, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "C_emissions_summary.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

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
    dplyr::mutate(T10_20=C2020/10,
           T20_30=C2030/10,
           T30_40=(C2040-C2030)/10,
           T40_50=(C2050-C2040)/10,
           T50_60=(C2060-C2050)/10,
           T60_70=(C2070-C2060)/10,
           T70_80=(C2080-C2070)/10,
           T80_90=(C2090-C2080)/10,
           T90_100=(C2100-C2090)/10) %>%
    dplyr::select(cont, T10_20:T90_100)

## Carbon emissions should continue to increase: from 0.66 PgC/yr on 2019-2035 to 0.814 PgC/yr on 2050-2085.
## Deforestation of forest areas with higher carbon stocks in the future.
## Higher carbon stocks because of environmental/elevation gradient (see Asner article) and because of remote, less degraded forests.
## Will decrease in Asia because many countries with no more forests after 2085.

## Deforestation => increase in C source with time (deforestation of forest with higher carbon stocks).
## Climate change => decrease in C sink with time (higher mortality), see Hubau2020.
## The result is that forests will likely become a major C source in the future. 

## Save results
f <- here("Analysis", dataset, "results", "C_trend.csv")
write.table(C_trend, file=f, sep=",", row.names=FALSE)
## Copy for manuscript
f_doc <- here("Manuscript", "Supplementary_Materials", "tables", "C_emissions.csv")
file.copy(from=f, to=f_doc, overwrite=TRUE)

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
