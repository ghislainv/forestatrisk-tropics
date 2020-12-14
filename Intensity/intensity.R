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
initGRASS(gisBase="/usr/lib/grass78",home=tempdir(), 
          gisDbase="grassdata",
          location="intensity",mapset="PERMANENT",
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

## This vector data-set as been imported as a table in Google Earth Engine
## with id "users/ghislainv/deforintensity/ctry_grid".

# ==============================
# Import table obtained with GEE
# ==============================

# Add area info
ctry_run <- read_delim(here("Analysis", "data", "ctry_run.csv"), delim=";")
area_info <- ctry_run %>%
  select("area_cont", "area_ctry", "area_name", "area_code")

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

# Compute total forest area per study area nd year
fc_df2 <- fc_df %>%
  group_by(area_code) %>%
  summarise(across(starts_with("fc"), sum), .groups="keep") %>%
  mutate(across(starts_with("fc"), function(.x){round(.x/10000)})) %>%
  left_join(area_info, by="area_code") %>%
  select(area_cont:area_name, area_code, fc2000:fc2020) %>%
  mutate(id=ifelse(area_cont=="America", 1, ifelse(area_cont=="Africa", 2, 3))) %>%
  arrange(id, area_ctry) %>%
  select(-id)

# Save data
write_delim(fc_df2, here("Intensity", "output", "fc_gee_2000_2020.csv"), delim=",")

# Compute annual deforestation
d_df <- fc_df2 %>%
  select(-fc2020)
d_df[, 5:24] <- fc_df2[,5:24]-fc_df2[,6:25]
names(d_df)[5:24] <- paste0("d", 2000:2019)

# Save data
write_delim(d_df, here("Intensity", "output", "defor_gee_2000_2020.csv"), delim=",")

# ==============================
# Compute mean and 95% quantiles
# Keep d in ha/yr here
# ==============================

# We assume defor+1 follows a logNormal distribution
# !! log(x) not define in 0 => defor+1 

# Functions
# Mean and variance of log-normal variable
# see: http://ww2.amstat.org/publications/jse/v13n1/olsson.html
log_theta <- function (x) {
	X <- log(x)
	mu_X <- mean(X, na.rm=TRUE)
	S2_X <- var(X, na.rm=TRUE)
	return(mu_X+S2_X/2)
}

var_log_theta <- function (x) {
	X <- log(x)
	n <- length(x)
	S2_X <- var(X, na.rm=TRUE)
	return(S2_X/n + S2_X^2/(2*(n-1)))
}

# Transform data in long format
d_long <- d_df %>%
  tidyr::pivot_longer(cols=c(d2000:d2019),
                      names_to="year",
                      names_prefix="d",
                      values_to="defor")

# Compute mean and 90% confidence intervals of deforestation rates per country
# t-distribution with n-1 = 9 degrees of freedom
d_uncertainty <- d_long %>%
	filter(year %in% c(2010:2019)) %>%
	group_by(area_name) %>%
	summarize(area_cont=unique(area_cont),
	          area_ctry=unique(area_ctry),
	          area_code=unique(area_code),
	          d_mean=exp(log_theta(defor+1))-1,
						d_min=exp(log_theta(defor+1)+qt(0.05,9)*sqrt(var_log_theta(defor+1)))-1,
						d_max=exp(log_theta(defor+1)+qt(0.95,9)*sqrt(var_log_theta(defor+1)))-1,
						.groups="keep") %>%
  relocate(area_name, .before=area_code) %>%
  mutate(id=ifelse(area_cont=="America", 1, ifelse(area_cont=="Africa", 2, 3))) %>%
  arrange(id, area_ctry) %>%
  select(-id)

# Save data
write_delim(d_uncertainty, here("Intensity", "output", "d_uncertainty.csv"), delim=",")

# ==============================
# Plot
# d in Kha/yr for plots
# ==============================

continents <- c("America", "Africa", "Asia")
ncont <- length(continents)

# Loop on continents
for (i in 1:ncont) {
  
  # Continent
  cont <- continents[i]
  # No country
  noctry <- ifelse(cont=="America", "Brazil", "Sao Tome and P.")
  
  # Dataframe for histogram
  df_hist <- d_long %>%
    dplyr::filter(area_cont==cont & area_ctry!=noctry) %>%
    mutate(defor=defor/1000)  # !! defor in Kha/yr !!
  
  # Dataframe for densities
  df_dens <- d_long %>%
    dplyr::filter(area_cont==cont & area_ctry!=noctry) %>%
    mutate(defor=defor/1000) %>% # !! defor in Kha/yr !!
    mutate(ldefor=log(defor+1/1000)) %>%
    group_by(area_code) %>%
    summarise(meanlog=mean(ldefor, na.rm=TRUE),sdlog=sd(ldefor, na.rm=TRUE),
              xmin=min(defor),xmax=max(defor), .groups="keep") %>%
    do(data.frame(defor=seq(.$xmin,.$xmax,length.out=100),
                  logDens=dlnorm(seq(.$xmin,.$xmax,length.out=100),.$meanlog,.$sdlog)))
  
  # Dataframe for confidence interval
  df_ci <- d_uncertainty %>%
    dplyr::filter(area_cont==cont & area_ctry!=noctry) %>%
    mutate(across(starts_with("d_"), function(.x){round(.x/1000, 3)}))
  
  # Plot histograms of disturbance
  p_hist <- ggplot(df_hist, aes(x=defor, group=area_code)) + 
    geom_histogram(aes(y=..density..), bins=10, color=grey(0.5), fill="white") +
    facet_wrap(~area_code, ncol=5, scales="free") +
    geom_line(data=df_dens, aes(y=logDens)) +
    geom_vline(data=df_ci, aes(xintercept=d_mean)) +
    geom_vline(data=df_ci, aes(xintercept=d_min), linetype="dashed") +
    geom_vline(data=df_ci, aes(xintercept=d_max), linetype="dashed") +
    xlab("Deforestation (Kha/yr)")
  ggsave(here("Intensity", "output", glue("hist_defor_{cont}.pdf")), p_hist, width=12, height=12)
  
}

# For Brazil
cont <- "Brazil"
# Dataframe for histogram
df_hist <- d_long %>%
  dplyr::filter(area_ctry=="Brazil") %>%
  mutate(defor=defor/1000)  # !! defor in Kha/yr !!

# Dataframe for densities
df_dens <- d_long %>%
  dplyr::filter(area_ctry=="Brazil") %>%
  mutate(defor=defor/1000) %>% # !! defor in Kha/yr !!
  mutate(ldefor=log(defor+1/1000)) %>%
  group_by(area_code) %>%
  summarise(meanlog=mean(ldefor, na.rm=TRUE),sdlog=sd(ldefor, na.rm=TRUE),
            xmin=min(defor),xmax=max(defor), .groups="keep") %>%
  do(data.frame(defor=seq(.$xmin,.$xmax,length.out=100),
                logDens=dlnorm(seq(.$xmin,.$xmax,length.out=100),.$meanlog,.$sdlog)))

# Dataframe for confidence interval
df_ci <- d_uncertainty %>%
  dplyr::filter(area_ctry=="Brazil") %>%
  mutate(across(starts_with("d_"), function(.x){round(.x/1000, 3)}))

# Plot histograms of disturbance
p_hist <- ggplot(df_hist, aes(x=defor, group=area_code)) + 
  geom_histogram(aes(y=..density..), bins=10, color=grey(0.5), fill="white") +
  facet_wrap(~area_code, ncol=5, scales="free") +
  geom_line(data=df_dens, aes(y=logDens)) +
  geom_vline(data=df_ci, aes(xintercept=d_mean)) +
  geom_vline(data=df_ci, aes(xintercept=d_min), linetype="dashed") +
  geom_vline(data=df_ci, aes(xintercept=d_max), linetype="dashed") +
  xlab("Deforestation (Kha/yr)")
ggsave(here("Intensity", "output", glue("hist_defor_{cont}.pdf")), p_hist, width=12, height=12)

# EOF