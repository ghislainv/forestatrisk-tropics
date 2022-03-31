#!/usr/bin/python
# -*- coding: utf-8 -*-

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# python_version  :2.7
# license         :GPLv3
# ==============================================================================

# PRIOR TO EXECUTING THE FOLLOWING SCRIPT, YOU MUST HAVE CREDENTIALS FOR
# 1. Google Earth Engine
# 2. rclone with Google Drive: https://rclone.org/drive/

# Standard library imports
import os
import pkg_resources
import shutil  # for rmtree
# import sys

# Third party imports
import ee
import forestatrisk as far
import pandas as pd

# index_ctry = int(sys.argv[1])-1

# ==================
# Settings
# Earth engine
ee.Initialize()
# GDAL
os.environ["GDAL_CACHEMAX"] = "1024"
# ==================

# Country isocode
file_ctry_run = pkg_resources.resource_filename("forestatrisk",
                                                "data/ctry_run.csv")
data_ctry_run = pd.read_csv(file_ctry_run, sep=";", header=0)
iso3 = list(data_ctry_run.iso3)
nctry = len(iso3)  # 120

# Temporary
iso3 = "COD"
cont = "Africa"
os.chdir("/home/ghislain/Bureau/fartest/COD")
# End temporary


# Function for multiprocessing
def run_country(iso3):

    # GDAL temp directory
    far.make_dir("/share/nas2-amap/gvieilledent/tmp/tmp_" + iso3)
    os.environ["CPL_TMPDIR"] = "/share/nas2-amap/gvieilledent/tmp/tmp_" + iso3

    # Set original working directory
    cont = data_ctry_run.cont_run[data_ctry_run["iso3"] == iso3].iloc[0]
    owd = "/share/nas2-amap/gvieilledent/jrc2020/" + cont
    os.chdir(owd)
    far.make_dir(iso3)
    os.chdir(os.path.join(owd, iso3))

    # Albers Equal Area projections
    if cont == "Africa":
        proj = ("+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 "
                "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    elif cont == "Asia":
        proj = ("+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 "
                "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    else:
        proj = ("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 "
                "+x_0=0 +y_0=0 +ellps=aust_SA +units=m no_defs")

    # Extract data for country on GEE
    far.data.country_biomass_run(
        iso3=iso3,
        proj="EPSG:4326",
        output_dir="data_raw",
        keep_dir=True,
        gdrive_remote_rclone="gdrive_gv",
        gdrive_folder="GEE_biomass_whrc")

    # Download data locally from Google Drive
    far.data.country_biomass_download(
        iso3=iso3,
        gdrive_remote_rclone="gdrive_gv",
        gdrive_folder="GEE_biomass_whrc",
        output_dir="data_raw")

    # Mosaic and resample
    far.data.country_biomass_compute(
        iso3=iso3,
        input_dir="data_raw",
        output_dir="data/emissions",
        proj=proj)

    # Carbon emissions
    carbon = far.emissions(input_stocks="data/emissions/biomass_whrc.tif",
                           input_forest="data/fcc23.tif")  # 745223502
    carbon_2050 = far.emissions(input_stocks="data/emissions/biomass_whrc.tif",
                                input_forest="data/fcc_2050.tif")  # 2521508823
    carbon_2100 = far.emissions(input_stocks="data/emissions/biomass_whrc.tif",
                                input_forest="data/fcc_2100.tif")  # 7413733535

    # Remove GDAL temp directory
    shutil.rmtree("/share/nas2-amap/gvieilledent/tmp/tmp_" + iso3)

    # Return country iso code
    return iso3


# Run country
for i in range(nctry):
    run_country(iso3[i])

# run_country(iso3[index_ctry])

# End
