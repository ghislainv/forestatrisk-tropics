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
import subprocess
import pkg_resources
import pandas as pd
import forestatrisk as far

#index_ctry = int(sys.argv[1])-1

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
    owd = "/share/nas2-amap/gvieilledent/gfc2019_70/" + cont
    os.chdir(owd)
    far.make_dir(iso3)
    os.chdir(os.path.join(owd, iso3))

    # Copy borders
    far.make_dir("data_raw")
    in_dir = os.path.join("/share/nas2-amap/gvieilledent/gfc2019_50",
                          cont, iso3, "data_raw/")
    out_dir = os.path.join("/share/nas2-amap/gvieilledent/gfc2019_70",
                           cont, iso3, "data_raw/")
    in_f = in_dir + "gadm36_" + iso3 + "_0.*"
    cmd = " ".join(["cp", in_f, out_dir])
    subprocess.call(cmd, shell=True)

    # Compute gee forest data
    far.data.country_forest_run(
        iso3, proj="EPSG:3395",
        output_dir="data_raw",
        keep_dir=True,
        fcc_source="gfc", perc=70,
        gdrive_remote_rclone="gdrive_gv",
        gdrive_folder="GEE-forestatrisk-tropics-gfc-70")

    # Return country iso code
    return(iso3)

# Run country
for i in range(nctry):
    run_country(iso3[i])

# End
