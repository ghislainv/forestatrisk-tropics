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

# Standard library imports
import os
import pkg_resources
import shutil  # for rmtree
import sys

# Third party imports
import forestatrisk as far
import pandas as pd

# Local imports
from run_forecast import run_forecast

# Arguments
index_ctry = int(sys.argv[1]) - 1

# ==================
# Settings
# GDAL
os.environ["GDAL_CACHEMAX"] = "1024"
# ==================

# Directories
cluster = "meso"  # (can be fdb, mbb, meso)
if cluster == "meso":
    work_dir = "/storage/replicated/cirad/projects/AMAP/vieilledentg/jrc2020/"
    temp_dir = "/lustre/vieilledentg/tmp/"
if cluster == "mbb":
    work_dir = "/share/nas2-amap/gvieilledent/jrc2020/"
    temp_dir = "/share/nas2-amap/gvieilledent/tmp/"
if cluster == "fdb":
    work_dir = "/home/forestatrisk-tropics/jrc2020/"
    temp_dir = "/home/forestatrisk-tropics/tmp/"

# Country isocode
file_ctry_run = pkg_resources.resource_filename("forestatrisk",
                                                "data/ctry_run.csv")
data_ctry_run = pd.read_csv(file_ctry_run, sep=";", header=0)
iso3 = list(data_ctry_run.iso3)
nctry = len(iso3)  # 120


# Function for multiprocessing
def run_country(iso3):

    # GDAL temp directory
    far.make_dir(temp_dir + "tmp_" + iso3)
    os.environ["CPL_TMPDIR"] = temp_dir + "tmp_" + iso3

    # Set original working directory
    cont = data_ctry_run.cont_run[data_ctry_run["iso3"] == iso3].iloc[0]
    owd = work_dir + cont
    os.chdir(owd)
    far.make_dir(iso3)
    os.chdir(os.path.join(owd, iso3))

    # Forecast
    run_forecast(iso3, fcc_source="jrc")

    # Remove GDAL temp directory
    shutil.rmtree(temp_dir + "tmp_" + iso3)

    # Print country iso code
    print(iso3)


# Run country
# for i in range(nctry):
#     run_country(iso3[i])

run_country(iso3[index_ctry])

# End
