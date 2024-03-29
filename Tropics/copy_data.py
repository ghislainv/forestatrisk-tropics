#!/usr/bin/python
# -*- coding: utf-8 -*-

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# python_version  :>=3
# license         :GPLv3
# ==============================================================================

# PRIOR TO EXECUTING THE FOLLOWING SCRIPT, YOU MUST HAVE CREDENTIALS FOR
# 1. Google Earth Engine
# 2. rclone with Google Drive: https://rclone.org/drive/
# 3. WDPA: https://www.protectedplanet.net/

import os
# import sys
import subprocess
import pkg_resources
import pandas as pd
import forestatrisk as far

# index_ctry = int(sys.argv[1])-1

# # ==================
# # Settings
# # Earth engine
# import ee
# ee.Initialize()
# # WDPA API
# from dotenv import load_dotenv
# load_dotenv("/home/gvieilledent/Code/forestatrisk-tropics/.env")
# from pywdpa import get_token
# get_token()
# # ==================

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
    far.make_dir("data_raw")

    # Input/output directories
    in_dir = os.path.join("/share/nas2-amap/gvieilledent/jrc2020_bak",
                          cont, iso3, "data_raw/")
    out_dir = os.path.join("/share/nas2-amap/gvieilledent/jrc2020",
                           cont, iso3, "data_raw/")

    # # Copy GADM borders
    # in_f = in_dir + "gadm36_" + iso3 + "_0.*"
    # cmd = " ".join(["cp", in_f, out_dir])
    # subprocess.call(cmd, shell=True)

    # Copy SRTM
    in_f = in_dir + "SRTM_V41_*.zip"
    cmd = " ".join(["cp", in_f, out_dir])
    subprocess.call(cmd, shell=True)

    # Copy WDPA
    in_f = in_dir + "pa_" + iso3 + ".*"
    cmd = " ".join(["cp", in_f, out_dir])
    subprocess.call(cmd, shell=True)

    # Copy OSM
    in_f = in_dir + "country.osm.pbf"
    cmd = " ".join(["cp", in_f, out_dir])
    subprocess.call(cmd, shell=True)

    # # Copy country data
    # in_dir = os.path.join("/share/nas2-amap/gvieilledent/jrc2020_bak",
    #                       cont, iso3, "data")
    # out_dir = os.path.join("/share/nas2-amap/gvieilledent/jrc2020",
    #                        cont, iso3, "data")
    # cmd = " ".join(["rclone sync", in_dir, out_dir])
    # subprocess.call(cmd, shell=True)

    # Return country iso code
    return(iso3)


# Run country
for i in range(nctry):
    run_country(iso3[i])

# run_country(iso3[index_ctry])

# End
