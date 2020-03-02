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

import os
import forestatrisk as far
import pandas as pd
import pkg_resources
import psutil
import multiprocessing as mp

# Set PROJ_LIB
# os.environ["PROJ_LIB"] = "/home/ghislain/miniconda3/envs/forestatrisk/share/proj"

# ==================
# Settings
# Earth engine
import ee
ee.Initialize()
# WDPA API
from dotenv import load_dotenv
load_dotenv(os.path.expanduser("~/Code/forestatrisk-tropics/.env"))
from pywdpa import get_token
get_token()
# ==================

# Set wd
owd = "/home/ghislain/Code/forestatrisk-tropics/Asia"
os.chdir(owd)
from run_modelling_steps import run_modelling_steps

# Country isocode
file_ctry_run = pkg_resources.resource_filename("forestatrisk",
                                                "data/ctry_run.csv")
data_ctry_run = pd.read_csv(file_ctry_run, sep=";", header=0)
iso3 = list(data_ctry_run.iso3[data_ctry_run.cont_run == "Asia"])
iso3.sort()
# print(iso3)
# iso3 = ['BGD', 'BRN', 'BTN', 'FJI', 'IDN', 'IND', 'KHM', 'LAO', 'LKA', 'MMR',
#         'MYS', 'NCL', 'PHL', 'PNG', 'SGP', 'SLB', 'THA', 'TLS', 'VNM', 'VUT']
iso3 = ['IDN', 'IND', 'KHM', 'LAO', 'LKA', 'MMR',
        'MYS', 'NCL', 'PHL', 'PNG', 'SGP', 'SLB', 'THA', 'TLS', 'VNM', 'VUT']

# Number of countries
nctry = len(iso3)

# Number of cpu
total_cpu = psutil.cpu_count()
num_cpu = int(total_cpu * 0.75) if total_cpu > 2 else 1
# num_cpu = 3

#far.data.country_wdpa(iso3[0], output_dir=iso3[0] + "/data_raw")

# Function for multiprocessing
def run_country(iso3):

    # Make new directory for country
    os.chdir(owd)
    far.make_dir(iso3)
    os.chdir(os.path.join(owd, iso3))

    # Download data
    far.data.country_download(
        iso3,
        gdrive_remote_rclone="gdrive_gv",
        gdrive_folder="GEE-forestatrisk-tropics",
        output_dir="data_raw")
    
    # Compute variables
    far.data.country_compute(
        iso3,
        temp_dir="data_raw",
        output_dir="data",
        proj="EPSG:3395",
        data_country=True,
        data_forest=True,
        keep_temp_dir=True)
    
    # Model and Forecast
    run_modelling_steps(fcc_source="jrc")

    # Return country iso code
    return(iso3)

# For loop
for i in iso3:
    print("\nGEE for country: " + i)
    far.make_dir(i + "/data_raw")
    far.data.country_forest_run(
        iso3=i, proj="EPSG:3395",
        output_dir=i + "/data_raw",
        keep_dir=True,
        fcc_source="gfc", perc=50,
        gdrive_remote_rclone="gdrive_gv",
        gdrive_folder="GEE-forestatrisk-tropics-gfc-50"
    )
   #run_country(i)

# # Parallel computation
# pool = mp.Pool(processes=num_cpu)
# results = [pool.apply_async(run_country, args=(x,)) for x in iso3]
# pool.close()
# pool.join()
# output = [p.get() for p in results]
# print(output)

# Combine results

# Combine country borders
# os.system("find -type f -name *ctry_PROJ.shp \
# -exec ogr2ogr -update -append borders.shp {} \;")

# Spatial probability
# os.system("find -type f -name *prob.tif > list_prob.txt")
# os.system("gdalbuildvrt -input_file_list list_prob.txt prob.vrt")
# os.system("gdal_translate -co 'COMPRESS=LZW' -co 'PREDICTOR=2' -co 'BIGTIFF=YES' \
# prob.vrt prob.tif")
# Build overview
# os.system("gdaladdo -ro -r nearest prob.tif 16")
# Plot
# dfp.plot.prob("prob.tif", output_file="prob.png",
#               borders="borders.shp", zoom=None, dpi=300,
#               lw=0.5, c="grey")

# Forest cover in 2050
# os.system("find -type f -name *fcc_40yr.tif > list_fcc_40yr.txt")
# os.system("gdalbuildvrt -input_file_list list_fcc_40yr.txt fcc_40yr.vrt")
# os.system("gdal_translate -co 'COMPRESS=LZW' -co 'PREDICTOR=2' -co 'BIGTIFF=YES' \
# fcc_40yr.vrt fcc_40yr.tif")
# Build overview
# os.system("gdaladdo -ro -r nearest --config COMPRESS_OVERVIEW LZW \
# --config PREDICTOR_OVERVIEW 2 \
# --config BIGTIFF_OVERVIEW YES \
# fcc_40yr.tif 16")
# Plot
# dfp.plot.fcc("fcc_40yr.tif", output_file="fcc_40yr.png",
#              borders="borders.shp", overview=False, zoom=None, dpi=300,
#              lw=0.5, c="grey")

# Upload results to Google Cloud Storage with gsutil
# os.system("gsutil -o GSUtil:parallel_composite_upload_threshold=150M \
# cp fcc_40yr.tif gs://deforestprob/output_roadless/fcc_40yr.tif")
# os.system("gsutil cp fcc_40yr.png \
# gs://deforestprob/output_roadless/fcc_40yr.png")
# os.system("gsutil -o GSUtil:parallel_composite_upload_threshold=150M \
# cp prob.tif gs://deforestprob/output_roadless/prob.tif")
# os.system("gsutil cp prob.png gs://deforestprob/output_roadless/prob.png")

# End
