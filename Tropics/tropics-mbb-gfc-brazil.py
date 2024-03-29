#!/usr/bin/env python
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

# Country isocode for Brazil
file_ctry_run = pkg_resources.resource_filename("forestatrisk",
                                                "data/ctry_run.csv")
data_ctry_run = pd.read_csv(file_ctry_run, sep=";", header=0)
iso3 = data_ctry_run.loc[data_ctry_run["cont_run"]=="Brazil", "iso3"].tolist()
nctry = len(iso3)  # 26


# Function for multiprocessing
def run_country(iso3):

    # GDAL temp directory
    far.make_dir("/share/nas2-amap/gvieilledent/tmp/tmp_" + iso3)
    os.environ["CPL_TMPDIR"] = "/share/nas2-amap/gvieilledent/tmp/tmp_" + iso3

    # Set original working directory
    cont = data_ctry_run.cont_run[data_ctry_run["iso3"] == iso3].iloc[0]
    owd = "/share/nas2-amap/gvieilledent/gfc2020_70/" + cont
    os.chdir(owd)
    far.make_dir(iso3)
    os.chdir(os.path.join(owd, iso3))

    # Model and Forecast
    run_modelling_steps(iso3, fcc_source="gfc")

    # Remove GDAL tmp directory
    shutil.rmtree("/share/nas2-amap/gvieilledent/tmp/tmp_" + iso3)

    # Return country iso code
    return(iso3)

# Run country
run_country(iso3[index_ctry])

# End
