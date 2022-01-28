#!/usr/bin/R

## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
## web             :https://ghislainv.github.io
## license         :GPLv3
## ==============================================================================

# ================================================
# Make a vector data-set combining country borders 
# and 1 degree grid.
# ================================================

# GRASS GIS 7.x.x is needed to run this script
# https://grass.osgeo.org/

# Libraries
library(rgrass7)
library(readr)
library(dplyr)
library(here)
library(glue)
library(ggplot2)

# Variables
dataset <- "jrc2020"
dir.create(here("Intensity", "output"))

Sys.setenv(PROJ_LIB="/usr/share/proj")

# =====================================
# Reproject in Lat-Long
# =====================================

continents <- c("America", "Africa", "Asia")

for (i in 1:length(continents)) {
  cont <- continents[i]
  input_file <- here("Maps", dataset, cont, glue("borders_{cont}.gpkg"))
  output_file <- here("Intensity", "output", glue("borders_{cont}_latlong.gpkg"))
  system(glue("ogr2ogr -s_srs EPSG:3395 -t_srs EPSG:4326 -f GPKG {output_file} {input_file}"))
}

# =====================================
# Create new grass location in Lat-Long
# =====================================

# dir.create("Intensity/grassdata")
# system("grass78 -c EPSG:4326 -e Intensity/grassdata/intensity")  # Ignore errors

# Connect R to grass location
# Make sure that /usr/lib/grass72/lib is in your PATH in RStudio
Sys.setenv(LD_LIBRARY_PATH=paste("/usr/lib/grass78/lib", Sys.getenv("LD_LIBRARY_PATH"), sep=":"))
initGRASS(gisBase="/usr/lib/grass78", home=tempdir(), 
          gisDbase="grassdata",
          location="intensity", mapset="PERMANENT",
          override=TRUE)

#======================================
# Import and compute overlays

# Import ctry country data into GRASS
# Snapping is necessary
# See https://en.wikipedia.org/wiki/Decimal_degrees 
# for correspondance between dd and meters (0.0001 dd ~ 4.3 m)

# Create 1 degree size grid:
system("g.region n=90 s=-90 w=-180 e=180 res=1 -p")
system("v.mkgrid map=grid grid=180,360")

for (i in 1:length(continents)) {
  
  # Import continent borders
  cont <- continents[i]
  input_file <- here("Intensity", "output", glue("borders_{cont}_latlong.gpkg"))
  snap <- ifelse(cont=="America", 0.0001, 1e-5)
  system(glue("v.in.ogr --o -o input={input_file} output=ctry_{cont} snap={snap}"))

  # Overlay country boundaries and grid
  system(glue("v.overlay --o ainput=ctry_{cont} binput=grid operator=and output=ctry_{cont}_grid"))
  
  # Rename columns
  system(glue("v.info -c ctry_{cont}_grid"))
  system(glue("v.db.renamecolumn ctry_{cont}_grid column=a_GID_0,area_code"))
  system(glue("v.db.renamecolumn ctry_{cont}_grid column=a_NAME_0,area_name"))
  
  # Export as shp
  output_file <- here("Intensity", "output", glue("borders_{cont}_grid.shp"))
  system(glue("v.out.ogr --o format=ESRI_Shapefile input=ctry_{cont}_grid output={output_file}"))
}

## 1. This vector data-sets have been imported as tables in Google Earth Engine
## with IDs "users/ghislainv/forestatrisk-tropics-intensity/borders_continents[i]_grid".

## 2. Then gee_fcc.py has been executed.
## system("python gee_fcc.py")

## 3. Finally, resulting datasets ForestArea_grid_{cont}.csv have been
## imported in folder "output"

# ==============================
# Import table obtained with GEE
# ==============================

# Add area info
ctry_run <- read_delim(here("Analysis", "data", "ctry_run.csv"), delim=";")
area_info <- ctry_run %>%
  select("iso3", "area_cont", "area_ctry", "area_name", "area_code")

# Area per cell grid
continents <- c("Afr", "Ame", "Asi")
ncont <- length(continents)
fc_df <- tibble()

# Loop on continent
for (i in 1:ncont) {
  cont <- continents[i]
  fc <- read_csv(here("Intensity", "output", glue("ForestArea_grid_{cont}.csv")))
  fc_df <- bind_rows(fc_df, fc)
}

# Compute total forest area per study area and year
fc_df2 <- fc_df %>%
  group_by(area_code) %>%
  summarise(across(starts_with("fc"), sum), .groups="keep") %>%
  mutate(across(starts_with("fc"), function(.x){round(.x/10000)})) %>%
  left_join(area_info, by="area_code") %>%
  select(iso3, area_cont:area_name, area_code, fc2000:fc2021) %>%
  mutate(id=ifelse(area_cont=="America", 1, ifelse(area_cont=="Africa", 2, 3))) %>%
  arrange(id, area_ctry) %>%
  select(-id)

# Save data
write_delim(fc_df2, here("Intensity", "output", "fc_gee_2000_2021.csv"), delim=",")

# Compute annual deforestation
d_df <- fc_df2 %>%
  select(-fc2021)
d_df[, 6:26] <- fc_df2[, 6:26]-fc_df2[, 7:27]
names(d_df)[6:26] <- paste0("d", 2000:2020)

# Save data
write_delim(d_df, here("Intensity", "output", "defor_gee_2000_2020.csv"), delim=",")

# ==============================
# Compute mean and 90% quantiles
# Keep d in ha/yr here
# ==============================

# Standard-error function
# Must take into account NAs in x for number of observations
se <- function(x) {
    x <- na.omit(x)
    return(sqrt(var(x)/length(x)))
}

# Load data
d_df <- read_delim(here("Intensity", "output", "defor_gee_2000_2020.csv"), delim=",")

# Transform data in long format
d_long <- d_df %>%
  tidyr::pivot_longer(cols=c(d2000:d2020),
                      names_to="year",
                      names_prefix="d",
                      values_to="defor")

# Remove outliers for Comoros and Gambia (unrealistic low or high annual deforestation estimates)
d_long <- d_long %>%
    mutate(defor=ifelse(area_code=="COM" & (defor<=10 | defor>=2000), NA, defor)) %>%
    mutate(defor=ifelse(area_code=="GMB" & (defor<=10 | defor>=2000), NA, defor))

# Compute mean and 95% confidence intervals of deforestation rates per country
# !! Ten years period on 2010-2019:
# !! We do not consider 2020 for which defor is underestimated (no difference between defor/degrad).
d_uncertainty <- d_long %>%
	filter(year %in% c(2010:2019)) %>%
	group_by(area_name) %>%
	summarize(iso3=unique(iso3),
	          area_cont=unique(area_cont),
	          area_ctry=unique(area_ctry),
	          area_code=unique(area_code),
	          d_se=round(se(defor)),
	          d_mean=round(mean(defor, na.rm=TRUE)),
                  d_min=round(mean(defor, na.rm=TRUE)-1.96*se(defor)),
                  d_max=round(mean(defor, na.rm=TRUE)+1.96*se(defor)),
                  .groups="keep") %>%
  relocate(area_name, .before=area_code) %>%
  mutate(id=ifelse(area_cont=="America", 1, ifelse(area_cont=="Africa", 2, 3))) %>%
  arrange(id, area_ctry) %>%
  select(-id)

# Save data
f <- here("Intensity", "output", "d_uncertainty.csv")
write_delim(d_uncertainty, f, delim=",")

# ==============================
# Plot
# d in Kha/yr for plots
# ==============================

# Dataframe for histogram
df_hist <- d_long %>%
  mutate(defor=defor/1000) %>% # !! defor in Kha/yr !!
  mutate(area_plot=ifelse(area_ctry=="India", paste0("IND-", area_code), area_code)) %>%
  mutate(area_plot=ifelse(area_ctry=="Brazil", paste0("BRA-", area_code), area_plot)) %>%
  mutate(cont_plot=ifelse(area_ctry=="Brazil", "Brazil", area_cont)) %>%
  dplyr::filter(area_code!="STP")
  
# Dataframe for confidence interval
df_ci <- d_uncertainty %>%
  mutate(across(starts_with("d_"), function(.x){round(.x/1000, 3)})) %>%
  mutate(area_plot=ifelse(area_ctry=="India", paste0("IND-", area_code), area_code)) %>%
  mutate(area_plot=ifelse(area_ctry=="Brazil", paste0("BRA-", area_code), area_plot)) %>%
  mutate(cont_plot=ifelse(area_ctry=="Brazil", "Brazil", area_cont)) %>%
  dplyr::filter(area_code!="STP")

# Some variables
continents <- c("America", "Africa", "Asia", "Brazil")
ncont <- length(continents)
text_width <- 16.6
fig_width <- text_width

# Loop on continents
for (i in 1:ncont) {
  
  # Continent
  cont <- continents[i]

  # Dataframe for histogram
  df_hist_cont <- df_hist %>% dplyr::filter(cont_plot==cont)
  # Dataframe for confidence interval
  df_ci_cont <- df_ci %>% dplyr::filter(cont_plot==cont)

  # nctry and npages
  study_areas <- unique(df_hist_cont$area_plot)
  nctry <- length(study_areas)
  npan_by_row <- 4  # Number of panels per row
  npan_by_col <- npan_by_row + 1
  npan_by_page <- npan_by_row*npan_by_col
  npages <- ceiling(nctry/(npan_by_page))

  # List of countries per pages
  area_list <- list()
  for (p in 1:npages) {
    inf <- (p-1)*npan_by_page+1
    sup <- min(nctry, p*npan_by_page)
    area_list[[p]] <- study_areas[inf:sup]
  }
  
  # Loop on pages
  for (p in 1:npages) {

    # Select countries for page
    df_hist_page <- df_hist_cont %>% dplyr::filter(area_plot %in% area_list[[p]])
    df_ci_page <- df_ci_cont %>% dplyr::filter(area_plot %in% area_list[[p]])
    
    # npan_in_page
    npan_in_page <- length(unique(df_hist_page$area_plot))
    nrow <- ceiling(npan_in_page/npan_by_row)
    
    # fig_height
    fig_height <- (fig_width/npan_by_row)*nrow
    
    # Plot histograms of disturbance
    p_hist <- ggplot(df_hist_page, aes(x=defor, group=area_plot)) + 
      geom_histogram(aes(y=..density..), bins=10, color=grey(0.5), fill="white") +
      facet_wrap(~area_plot, ncol=npan_by_row, nrow=nrow, scales="free") +
      geom_vline(data=df_ci_page, aes(xintercept=d_mean)) +
      geom_vline(data=df_ci_page, aes(xintercept=d_min), linetype="dashed") +
      geom_vline(data=df_ci_page, aes(xintercept=d_max), linetype="dashed") +
      xlab("Deforestation (Kha/yr)")
    
    # Save
    ggsave(here("Intensity", "output", glue("hist_defor_{cont}_{p}.png")),
           p_hist,
           width=fig_width, height=fig_height, units="cm")
  }
}

## # Comparing version of TMF
## df1 <- read_csv(here("Intensity", "output", "d_uncertainty.bak.csv"))
## df2 <- read_csv(here("Intensity", "output", "d_uncertainty.csv"))

## png(file=here("Intensity", "output", "comp_d_mean.png"))
## plot(df1$d_mean, df2$d_mean,
##      xlab="d_mean per country -- v0_2019",
##      ylab="dmean per country -- v1_2020")
## abline(a=0, b=1, col="red")
## dev.off()

## # Comparing version of TMF
## df1 <- read_csv(here("Intensity", "output", "fc_gee_2000_2020.csv"))
## df2 <- read_csv(here("Intensity", "output", "fc_gee_2000_2021.csv"))

## png(file=here("Intensity", "output", "comp_fc2020.png"))
## plot(df1$fc2020, df2$fc2020,
##      xlab="fc2020 per country -- v0_2019",
##      ylab="fc2020 per country -- v1_2020")
## abline(a=0, b=1, col="red")
## dev.off()

# EOF
