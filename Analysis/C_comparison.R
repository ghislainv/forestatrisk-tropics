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
require(ggplot2)
require(readr)

## Set working directory
setwd(here())

## Dataset
dataset <- "jrc2020"
dir.create(here("Analysis", dataset), recursive=TRUE)

## Load country info (encoding pb for ctry names)
ctry_df <- read_delim(here("Analysis", "data", "ctry_run.csv"), delim=";")

## Load C emissions from WHRC (in kg C)
C_whrc <- read_csv(here("Analysis", dataset, "C_emissions_mean_whrc.csv"))

## Import results by Clement
jrc_zarin <- read_csv(here("Analysis", "data", "STAT_Annual_M1_Zarin_CarbonLoss_Country.csv"))
jrc_zarin <- jrc_zarin %>%
    dplyr::filter(Year %in% 2010:2019) %>%
    group_by(ZoneName) %>%
    summarize(C=sum(S_Defor)) %>%
    left_join(ctry_df, by=c("ZoneName"="ctry_run")) %>%
    select("ZoneName", "C", "cont_run", "iso3") %>%
    arrange(cont_run, ZoneName)

## Export to manually fill NA
write_csv(jrc_zarin, file=here("Analysis", dataset, "jrc_zarin.csv"))

## Reimport
jrc_zarin <- read_csv(here("Analysis", dataset, "jrc_zarin.csv"))

## Join and set C unit in PgC
jrc_zarin <- jrc_zarin %>%
    left_join(C_whrc, by=c("iso3"="area_code")) %>%
    select(-(C2030:C2110)) %>%
    mutate(C=1e3*C, C2020=C2020/1e3)

## Plot
ggplot(data=jrc_zarin, aes(C, C2020)) +
    geom_point() +
    geom_abline(intercept=0, slope=1)
ggsave(here("Analysis", dataset, "jrc_zarin.png"))

## EOF
