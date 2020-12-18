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
# 3. WDPA: https://www.protectedplanet.net/

import sys
import os
import shutil  # for rmtree
import re  # regular expressions
import pkg_resources
import pandas as pd
import forestatrisk as far
from run_modelling_steps import run_modelling_steps

index_ctry = int(sys.argv[1])-1

# ==================
# Settings
# Earth engine
import ee
ee.Initialize()
# WDPA API
from dotenv import load_dotenv
load_dotenv("/home/gvieilledent/Code/forestatrisk-tropics/.env")
from pywdpa import get_token
get_token()
# GDAL
os.environ["GDAL_CACHEMAX"] = "1024"
# ==================

# Country isocode
file_ctry_run = pkg_resources.resource_filename("forestatrisk",
                                                "data/ctry_run.csv")
data_ctry_run = pd.read_csv(file_ctry_run, sep=";", header=0)
iso3 = list(data_ctry_run.iso3)
nctry = len(iso3)  # 120


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

    f = "data/forest/forest_t1.tif"
    if os.path.exists(f) is not True:

        # # Download data
        # far.data.country_download(
        #     iso3,
        #     gdrive_remote_rclone="gdrive_gv",
        #     gdrive_folder="GEE-forestatrisk-tropics-jrc-2020",
        #     output_dir="data_raw")

        # Download forest data
        far.data.country_forest_download(
            iso3,
            gdrive_remote_rclone="gdrive_gv",
            gdrive_folder="GEE-forestatrisk-tropics-jrc-2020",
            output_dir="data_raw")

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

        # Compute variables
        far.data.country_compute(
            iso3,
            temp_dir="data_raw",
            output_dir="data",
            proj=proj,
            data_country=False,
            data_forest=True,
            keep_temp_dir=True)

    # # If not Brazil
    # p = re.compile("BRA-.*")
    # m = p.match(iso3)
    # if m is None:
    #     # Model and Forecast
    #     run_modelling_steps(iso3, fcc_source="jrc")

    # Remove GDAL temp directory
    shutil.rmtree("/share/nas2-amap/gvieilledent/tmp/tmp_" + iso3)

    # Return country iso code
    return(iso3)

# Run country
run_country(iso3[index_ctry])

# End
