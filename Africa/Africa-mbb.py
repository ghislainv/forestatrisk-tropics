#!/usr/bin/python
# -*- coding: utf-8 -*-

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# python_version  :2.7
# license         :GPLv3
# ==============================================================================

# PRIOR TO EXECUTING THE FOLLOWING SCRIPT, AUTHENTICATE TO
# 1. Google Earth Engine: earthengine authenticate
# 2. Google Drive with rclone: https://rclone.org/drive/

import sys
import os
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

# Change working directory to nas
os.chdir("/share/nas2-amap/gvieilledent/")

# Print iso3
iso3 = ['AGO', 'BDI', 'BEN', 'CAF', 'CIV', 'CMR', 'COD', 'COG', 'COM',
        'ETH', 'GAB', 'GHA', 'GIN', 'GMB', 'GNB', 'GNQ', 'KEN', 'LBR',
        'MDG', 'MOZ', 'MUS', 'MWI', 'MYT', 'NGA', 'REU', 'RWA', 'SEN',
        'SLE', 'SSD', 'STP', 'TGO', 'TZA', 'UGA', 'ZMB']

# Original working directory
owd = os.getcwd()

# Function for multiprocessing
def run_country(iso3):
    
    # Make new directory for country
    os.chdir(owd)
    far.make_dir(iso3)
    os.chdir(os.path.join(owd, iso3))
    
    # Data
    far.data.country(iso3=iso3,
                     proj="EPSG:3395",
                     data_country=True,
                     data_forest=True,
                     keep_data_raw=True,
                     fcc_source="jrc",
                     gdrive_remote_rclone="gdrive_gv",
                     gdrive_folder="GEE-forestatrisk-tropics")
    
    # Model and Forecast
    run_modelling_steps(fcc_source="jrc")
    
    # Return country iso code
    return(iso3)


# Run country
run_country(iso3[index_ctry])

# End
