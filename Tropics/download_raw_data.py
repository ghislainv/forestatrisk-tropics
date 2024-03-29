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

import os
import sys
# import shutil  # for rmtree
import pkg_resources

from dotenv import load_dotenv
import ee
import pandas as pd
from pywdpa import get_token
import forestatrisk as far

# from run_modelling_steps import run_modelling_steps

index_ctry = int(sys.argv[1])-1

# ==================
# Settings
# Earth engine
ee.Initialize()
# WDPA API
load_dotenv("/home/gvieilledent/Code/forestatrisk-tropics/.env")
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

    # Set original working directory
    cont = data_ctry_run.cont_run[data_ctry_run["iso3"] == iso3].iloc[0]
    owd = "/share/nas2-amap/gvieilledent/jrc2020/" + cont
    os.chdir(owd)
    far.make_dir(iso3)
    os.chdir(os.path.join(owd, iso3))

    # Download data
    far.data.country_download(
        iso3,
        gdrive_remote_rclone="gdrive_gv",
        gdrive_folder="GEE-forestatrisk-tropics-jrc-2020",
        output_dir="data_raw")

    # # Exceptions WDPA
    # if cont == "Brazil":
    #     if iso3 == "BRA-AC":
    #         # Download data
    #         far.data.country_wdpa(iso3="BRA", output_dir="data_raw")
    #         # Rename
    #         for ext in [".dbf", ".prj", ".shp", ".shx"]:
    #             os.rename("data_raw/pa_BRA" + ext,
    #                       "data_raw/pa_" + iso3 + ext)
    #     else:
    #         # Copy data
    #         for ext in [".dbf", ".prj", ".shp", ".shx"]:
    #             shutil.copy(os.path.join(owd, "BRA-AC", "data_raw",
    #                                      "pa_BRA-AC" + ext),
    #                         "data_raw/pa_" + iso3 + ext)
    # if cont == "Asia":
    #     if iso3 == "AUS-QLD":
    #         # Download data
    #         far.data.country_wdpa(iso3="AUS", output_dir="data_raw")
    #         # Rename
    #         for ext in [".dbf", ".prj", ".shp", ".shx"]:
    #             os.rename("data_raw/pa_AUS" + ext,
    #                       "data_raw/pa_" + iso3 + ext)
    #     if iso3 == "IND-AND":
    #         # Download data
    #         far.data.country_wdpa(iso3="IND", output_dir="data_raw")
    #         # Rename
    #         for ext in [".dbf", ".prj", ".shp", ".shx"]:
    #             os.rename("data_raw/pa_IND" + ext,
    #                       "data_raw/pa_" + iso3 + ext)
    #     if iso3 in ["IND-EAST", "IND-WEST"]:
    #         # Copy data
    #         for ext in [".dbf", ".prj", ".shp", ".shx"]:
    #             shutil.copy(os.path.join(owd, "IND-AND", "data_raw",
    #                                      "pa_IND-AND" + ext),
    #                         "data_raw/pa_" + iso3 + ext)

    # Return country iso code
    return iso3


# Run country
# for i in range(nctry):
#     run_country(iso3[i])

run_country(iso3[index_ctry])

# End
