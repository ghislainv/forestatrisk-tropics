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
# ==================

# Set working directory to nas
owd = "/share/nas2-amap/gvieilledent/Africa"
os.chdir(owd)

# Country isocode
file_ctry_run = pkg_resources.resource_filename("forestatrisk",
                                                "data/ctry_run.csv")
data_ctry_run = pd.read_csv(file_ctry_run, sep=";", header=0)
iso3 = list(data_ctry_run.iso3[data_ctry_run.cont_run == "Africa"])
iso3.sort()
# print(iso3)
# iso3 = ['AGO', 'BDI', 'BEN', 'CAF', 'CIV', 'CMR', 'COD', 'COG', 'COM',
#         'ETH', 'GAB', 'GHA', 'GIN', 'GMB', 'GNB', 'GNQ', 'KEN', 'LBR',
#         'MDG', 'MOZ', 'MUS', 'MWI', 'MYT', 'NGA', 'REU', 'RWA', 'SEN',
#         'SLE', 'SSD', 'STP', 'TGO', 'TZA', 'UGA', 'ZMB']
#iso3 = ["AGO", "COD"]


# Function for multiprocessing
def run_country(iso3):
    
    # Make new directory for country
    os.chdir(owd)
    far.make_dir(iso3)
    os.chdir(os.path.join(owd, iso3))
    
    # # Data
    # far.data.country(iso3=iso3,
    #                  proj="EPSG:3395",
    #                  data_country=True,
    #                  data_forest=True,
    #                  keep_data_raw=True,
    #                  fcc_source="jrc",
    #                  gdrive_remote_rclone="gdrive_gv",
    #                  gdrive_folder="GEE-forestatrisk-tropics")
    
    # Model and Forecast
    run_modelling_steps(fcc_source="jrc")
    
    # Return country iso code
    return(iso3)

# Run country
run_country(iso3[index_ctry])

# End
