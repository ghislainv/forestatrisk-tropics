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
import pkg_resources
# import re  # regular expressions
import shutil  # for rmtree
import sys

from dotenv import load_dotenv
import ee
import forestatrisk as far
import pandas as pd
from pywdpa import get_token

from run_forecast import run_forecast

ctry_id = [2, 95, 103, 112, 40, 46, 18, 20, 87, 25, 92, 67, 31]

# index_ctry = int(sys.argv[1]) - 1
index_ctry = ctry_id[int(sys.argv[1]) - 1] - 1

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

    # GDAL temp directory
    far.make_dir("/share/nas2-amap/gvieilledent/tmp/tmp_" + iso3)
    os.environ["CPL_TMPDIR"] = "/share/nas2-amap/gvieilledent/tmp/tmp_" + iso3

    # Set original working directory
    cont = data_ctry_run.cont_run[data_ctry_run["iso3"] == iso3].iloc[0]
    owd = "/share/nas2-amap/gvieilledent/jrc2020/" + cont
    os.chdir(owd)
    far.make_dir(iso3)
    os.chdir(os.path.join(owd, iso3))

    # Forecast
    run_forecast(iso3, fcc_source="jrc")

    # Remove GDAL temp directory
    shutil.rmtree("/share/nas2-amap/gvieilledent/tmp/tmp_" + iso3)

    # Return country iso code
    return iso3


# Run country
# for i in range(nctry):
#     run_country(iso3[i])

run_country(iso3[index_ctry])

# End
